# $Id: GOR4.pm,v 1.0 2003/07/ 11
#
# BioPerl module for Bio::Tools::Analysis::Protein::GOR4
#
# Copyright Richard Adams
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1  NAME

Bio::Tools::Analysis::Protein::GOR4 - a wrapper around GOR4 protein
secondary structure prediction server

=head1  SYNOPSIS

  use Bio::Tools::Analysis::Protein::GOR4;
  #get a Bio::Seq or Bio::PrimarySeq
  use Bio::PrimarySeq;
  $seq = Bio::PrimarySeq->new
    (-seq=>'IKLCVHHJHJHJHJHJHJHNLAILAKAHLIELALAL',
     -primary_id=>'test'); # a Bio::PrimarySeqI object

  my $gor4 = Bio::Tools::Analysis::Protein::GOR4->new (-seq=>$seq);
  $gor4->run;
  print $gor4->result;# #raw text to standard error

=head1  DESCRIPTION

A module to remotely retrieve predictions of protein secondary
structure.  Each residue in the protein receives a score representing
the likelihood of existing in each of three different states (helix,
coil or sheet), e.g.,

  my $analysis_object = Bio::Tools::SimpleAnalysis::Protein::GOR4->
      new(-seq => $seq);

creates a new object

  $analysis_object->run;

submits the query to the server and obtains raw text output

Given an amino acid sequence the results can be obtained in 4 formats,
determined by the argument to the result method

=over 4

=item 1

The raw text of the program output

  my $rawdata = $analysis_object->result;

=item 2

An reference to an array of hashes of scores for each state and the
assigned state.

  my $data_ref = $analysis_object->result('parsed');
  print "score for helix at residue 2 is $data_ref->[1]{'helix'}\n";
  print "predicted struc  at residue 2 is $data_ref->[1]{'struc}\n";

=item 3

An array of Bio::SeqFeature::Generic objects where each feature is a
predicted unit of secondary structure. Only stretches of helix/sheet
predictions for longer than 4 residues are defined as helices. See 
Bio::Tools::Analysis::Domcut.pm for examples of how to add sequence
features.

  my @fts = $analysis_object->result(Bio::SeqFeatureI);
  for my $ft (@fts) {
      print " From ",  $ft->start, " to  ",$ft->end, " struc: " ,
             ($ft->each_tag_value('type'))[0]  ,"\n";
  }

=item 4

A Bio::Seq::Meta::Array implementing sequence.

This is a Bio::Seq object that can also hold data about each residue
in the sequence In this case, the sequence can be associated with a
single array of GOR4 prediction scores.  e.g.,

  my $meta_sequence = $analysis_object->result('all');
  print "helix scores from residues 10-20 are ",
      $meta_sequence->named_submeta_text("GOR4_helix",10,20), "\n";

Meta sequence names are : GOR4_helix, GOR4_sheet, GOR4_coil,
GOR4_struc, representing the scores for each residue.

Many methods common to all analyses are inherited from
Bio::Tools::Analysis::SimpleAnalysisBase.

=back

=head1 SEE ALSO

L<Bio::SimpleAnalysisI>, 
L<Bio::Tools::Analysis::SimpleAnalysisBase>, 
L<Bio::Seq::Meta::Array>, 
L<Bio::WebAgent>

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHORS

Richard Adams, Richard.Adams@ed.ac.uk,

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

use strict;

package Bio::Tools::Analysis::Protein::GOR4;

use IO::String;
use Bio::SeqIO;
use HTTP::Request::Common qw(POST);
use Bio::SeqFeature::Generic;
use Bio::Seq::Meta::Array;


use base qw(Bio::Tools::Analysis::SimpleAnalysisBase);

use constant MIN_STRUC_LEN => 3;
my $URL = 'http://npsa-pbil.ibcp.fr/cgi-bin/secpred_gor4.pl';
my $ANALYSIS_NAME = 'GOR4';
my $ANALYSIS_SPEC = {name => 'Gor4', type => 'Protein'};
my $INPUT_SPEC    = [
                     {mandatory =>'true',
                      type      => 'Bio::PrimarySeqI',
                      'name'    => 'seq',
                  },
                 ];
my  $RESULT_SPEC =
    {
     ''                 => 'bulk',              # same as undef
     'Bio::SeqFeatureI' => 'ARRAY of Bio::SeqFeature::Generic',
     raw                => '[ {struc =>, helix=> ,sheet=>, coil=>}]',
     meta                => 'Bio::Seq::Meta::Array object',
    };

=head2 result

 Name    : result
 Usage   : $job->result (...)
 Returns : a result created by running an analysis
 Args    : see keys of $RESULT_SPEC

The method returns a result of an executed job. If the job was
terminated by an error the result may contain an error message instead
of the real data.

This implementation returns differently processed data depending on
argument:

=over 3

=item undef

Returns the raw ASCII data stream but without HTML tags

=item 'Bio::SeqFeatureI'

The argument string defines the type of bioperl objects returned in an
array.  The objects are L<Bio::SeqFeature::Generic>.  Feature primary
tag is "2ary".  Feature tags are "type" (which can be helix, sheet or
coil) "method" (GOR4).

=item 'parsed'

Array of hash references of { helix =E<gt>, sheet =E<gt> , coil =E<gt> , struc=E<gt>}.

=item 'meta'

A Bio::Seq::Meta::Array object. Scores can be accessed using methods
from this class. Meta sequence names are GOR4_helix, GOR4_sheet,
GOR4_coil, GOR4_struc.


=back


=cut


sub result {
    my ($self,$value) = @_;

    my @scores;
    my @fts;

    if ($value ) {
        #parse into basic raw form, store this as well as '_result'
        if (!exists($self->{'_parsed'}) ) {
            my $result = IO::String->new($self->{'_result'});
            while (my $line = <$result>) {
                next unless $line =~ /^\w\s/; # or for sopma/hnn  /^[A-Z]\s/
                $line =~/(\w)\s+(\d+)\s+(\d+)\s+(\d+)/; # or for so
                push @scores, { struc => $1,
                                helix => $2,
                                sheet => $3,
                                coil  => $4,
                              };
            }
            $self->{'_parsed'} = \@scores;
        }
        if ($value eq 'Bio::SeqFeatureI') {
            $self->_get_2ary_coords();
            for my $type (keys %{$self->{'_parsed_coords'}} ) {
                next if $type =~  /\w{2,}/; #if not H,C,E or T
                for my $loc (@{$self->{'_parsed_coords'}{$type}} ) {
                    push @fts, Bio::SeqFeature::Generic->new
                        (-start   => $loc->{'start'},
                         -end     => $loc->{'end'},
                         -source  => 'GOR4',
                         -primary => 'Region',
                         -tag => {
                                  type => $type,
                                  method => $self->analysis_name,
                                 });
                }               #end of array of strucs of type
            }                   # end of all 2nd struc elements
            delete $self->{'_parsed_coords'}; #remove temp data
            return @fts;
        }                       #endif BioSeqFeature

        elsif ($value eq 'meta') {
            #1st of all make 3 or 4 arrays of scores for each type from column data
            my %type_scores;
            for my $aa (@{$self->{'_parsed'}}) {
                push @{$type_scores{'struc'}}, $aa->{'struc'};
                push @{$type_scores{'helix'}}, $aa->{'helix'};
                push @{$type_scores{'sheet'}}, $aa->{'sheet'};
                push @{$type_scores{'coil'}}, $aa->{'coil'};
            }
			
			## bless if necessary ##
			if (!$self->seq->isa("Bio::Seq::Meta::Array")){
            	bless ($self->seq, "Bio::Seq::Meta::Array");
				}
            $self->seq->isa("Bio::Seq::MetaI")
                || $self->throw("$self is not a Bio::Seq::MetaI");
            $Bio::Seq::Meta::Array::DEFAULT_NAME = 'GOR4_struc';

            ## now make meta_Sequence
            for my $struc_type (keys %type_scores) {
                my $meta_name = "GOR4". "_" . "$struc_type";
                my @meta = map{$_->{$struc_type}} @{$self->{'_parsed'}};
                if (grep{$_ eq $meta_name}$self->seq->meta_names ) {
                    $self->warn ("$meta_name already exists , not overwriting!");
                    next;
                }
                $self->seq->named_meta($meta_name,\@meta );
            }
            # return  seq array object implementing meta sequence #
            return $self->seq;

        } 
		else  {
           return $self->{'_parsed'};
         }
    }  #endif ($value)

    #return raw result if no return fomrt stated
    return $self->{'_result'};
}

sub _get_2ary_coords {
    #helper sub for result;
    ##extracts runs of structure > MIN_STRUC_LENresidues or less if Turn:
    #i.e., helical prediction for 1 residue isn't very meaningful...
    ## and poulates array of hashes with start/end values.
    ##keys of $Result are 'H' 'T' 'C' 'E'.
    #could be put into a secondary base class if need be
    my ($self) = @_;
    my @prot = @{$self->{'_parsed'}};
    my %Result;
    for (my $index = 0; $index <= $#prot; $index++) {

        my $type        = $prot[$index]{'struc'};
        next unless $type =~ /[HTCE]/;
        my $length = 1;
        for (my $j = $index + 1; $j <= $#prot; $j++) {
            my $test = $prot[$j];
            if ($test->{'struc'} eq $type) {
                $length++;
            } elsif (  $length > MIN_STRUC_LEN  ||
                       ($length <= MIN_STRUC_LEN && $type eq 'T') ) {
                push @{$Result{$type}}, {start => $index + 1 ,  end => $j};
                $index += $length -1;
                last;
            } else {
                $index += $length - 1;
                last;
            }
        }
    }
    $self->{'_parsed_coords'} = \%Result; #temp assignment
}

sub _init {
    my $self = shift;
    $self->url($URL);
    $self->{'_ANALYSIS_SPEC'} =$ANALYSIS_SPEC;
    $self->{'_INPUT_SPEC'} =$INPUT_SPEC;
    $self->{'_RESULT_SPEC'} =$RESULT_SPEC;
    $self->{'_ANALYSIS_NAME'} =$ANALYSIS_NAME;
    return $self;
}


sub  _run {
    my $self  = shift;
    $self->delay(1);
    # delay repeated calls by default by 3 sec, set delay() to change
    $self->sleep;
    $self->status('TERMINATED_BY_ERROR');
    my $request = POST $self->url,
        Content_Type => 'form-data',
            Content  => [title => "",
                         notice => $self->seq->seq,
                         ali_width => 70,
                        ];

    my $content = $self->request($request);
    my $text = $content->content;
    return unless $text;
    my ($next) = $text =~ /Prediction.*?=(.*?)>/;
    return unless $next;
    my $out = 'http://npsa-pbil.ibcp.fr/'.$next;
    my $req2 = HTTP::Request->new(GET=>$out);
    my $resp2 = $self->request($req2);
	$self->status('COMPLETED') if $resp2 ne '';
    $self->{'_result'} = $resp2->content;
}




1;
