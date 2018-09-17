# $Id: Domcut.pm,v 1.0 2003/07/ 11
#
# BioPerl module for Bio::Tools::Analysis::Protein::Domcut
#
# Copyright Richard Adams
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Analysis::Protein::Domcut -  a wrapper around Domcut server

=head1  SYNOPSIS

  use   Bio::Tools::Analysis::Protein::Domcut;
  #get a  Bio::PrimarySeq
  use Bio::PrimarySeq;
  my $seq = Bio::PrimarySeq->new
     (-seq=>'IKLCVNLAILAKAHLIELALAL',
     -primary_id=>'test'); # a Bio::PrimarySeqI object

  my $domcut = Bio::Tools::Analysis::Protein::Domcut->new (-seq=>$seq);
  $domcut->run;
  print $domcut->result;# #raw text to standard out

=head1  DESCRIPTION

A module to remotely retrieve predictions of protein domain
boundaries.  Each residue in the protein receives a score, those
better than the significance threshold and at a local minimum receive
a rank - i.e., the best minimum is rank 1, the second best minimum is
rank2 etc. These correspond to domain boundaries.  e.g.,

  my $analysis_object = Bio::Tools::Analysis::Protein::Domcut->new
     (-seq => $seq);

creates a new object. The sequence supplied must be a Bio::PrimarySeq and not
a Bio::Seq object. 

  $analysis_object->run;

submits the query to the server and obtains raw text output

Given an amino acid sequence the results can be obtained in 4 formats,
determined by the argument to the result method

=over 4

=item 1

The raw text of the program output

  my $rawdata = $analysis_object->result;

=item 2

A reference to an array of hashes of scores for each state and the
assigned state. Each element in the array is a residue (indexed from 0).

  my $data_ref = $analysis_object->result('parsed');
  print "score for helix at residue 2 is $data_ref->[1]{'helix'}\n";
  print "predicted struc  at residue 2 is $data_ref->[1]{'struc}\n";

=item 3

An array of Bio::SeqFeature::Generic objects where each feature is a
predicted unit of secondary structure. Only stretches of helix/sheet
predictions for longer than 4 residues are defined as helices.
So, in order to add features to an existing Bio::Seq object;

  # get a Bio::Seq object
  my $seqobj;
  my $tool = Bio::Tools::Analysis::Protein::Domcut->new
      ( -seq => $seqobj->primary_seq);
  $tool->run;

  my @fts = $tool->result(Bio::SeqFeatureI);

  $seqobj->add_SeqFeature(@fts);

  # if you want  meta sequences as well :
  my $meta = $tool->result('meta');
  $seqobj->primary_seq($meta);

  # can access meta data in a Bio::Seq object via a 
  # call to primary_seq:

  print $seq4->primary_seq->named_submeta_text('Domcut', 1,2), "\n";

=item 4

A Bio::Seq::Meta::Array implementing sequence.

This is a Bio::Seq object that can also hold data about each residue
in the sequence. In this case, the sequence can be associated with a
single array of Domcut prediction scores.  e.g.,

  my $meta_sequence = $analysis_object->result('meta');
  print "scores from residues 10 -20 are ",
      $meta_sequence->submeta_text(10,20), "\n";

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
package Bio::Tools::Analysis::Protein::Domcut;
use IO::String;
use Bio::SeqIO;
use HTTP::Request::Common qw(GET);
use Bio::SeqFeature::Generic;
use Bio::Seq::Meta::Array;

use base qw(Bio::Tools::Analysis::SimpleAnalysisBase);

my $URL = 'http://www.Bork.EMBL-Heidelberg.DE/Docu/mikita/domplot.cgi?';
my $ANALYSIS_NAME = 'Domcut';
my $ANALYSIS_SPEC =
    {
     'name'        => 'Domcut',
     'type'        => 'protein', #compulsory entry as is used for seq checking
     'version'     => 'n/a',
     'supplier'    => 'Ohara lab, Laboratory of DNA technology, 
                       Kazusa DNA Research Institute, 1532-3 Yana,
                       Kisarazu, Japan',
     'description' => 'to predict domain boundaries in proteins',
     'reference'   => 'Bioinformatics 19, 673-674 (2003)',
    };


my $INPUT_SPEC =
    [
     {
      'mandatory' => 'true',
      'type'      => 'Bio::PrimarySeqI',
      'name'      => 'seq',
     },
    ];

my  $RESULT_SPEC =
    {
     ''                 => 'bulk',              # same as undef
     'Bio::SeqFeatureI' => 'ARRAY of Bio::SeqFeature::Generic',
     'parsed'           => "Array of {'score' =>, 'rank'=> ]",
     'meta'             => 'Bio::Seq::Meta::Array object'
    };

=head2 result

 Name    : result
 Purpose : To retrieve results of analysis in one of several formats.
 Usage   : $job->result (...)
 Returns : a result created by running an analysis
 Args    : various - see keysin $RESULT_SPEC. 

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
array.  The objects are L<Bio::SeqFeature::Generic>. Tagnames are 'score' 
and 'rank'.

=item 'parsed'

Array of array references of [score, rank].

=item 'all'

A Bio::Seq::Meta::Array object. Scores can be accessed using methods
from this class. Meta sequence name is Domcut.

=back

=cut


sub result {
    my ($self,$value) = @_;
    my @scores;
    my @fts;

    if ($value ) {
        # parse raw text if not already done so
        if (!exists($self->{'_parsed'})) {
            my $result = IO::String->new($self->{'_result'});
            while (my $line = <$result>) {
                next if $line =~/#/;
                $line =~/(\-?\d\.\d+)\s+(\d+)?/;
                push @scores, {score => $1,
                               rank  => ($2)?$2:'' ,
                              };
            }
            #hold parsed results in object, saves having to reparse each time
            $self->{'_parsed'} = \@scores;
        }
        #make aarray of Bio::SeqFeature::Generic objects
        if ($value eq 'Bio::SeqFeatureI') {
            my $i = 0;          #array index (= aa num -1)
			my $in_trough = 0;
			my ($st, $end, $rank, $min_score, $min_locus) = (0,0,0,0,0);
			my $seqlen = $self->seq->length();
            for my $score (@{$self->{'_parsed'}}) {

				##start a potential trough
				if ($in_trough == 0 && $score->{'score'} < -0.09) {
					$in_trough = 1;
					$st        = $i+1;
					}

				## in a trough, is it ranked?
				elsif ( $in_trough == 1 && $score->{'score'} < -0.09 && $i +1 < $seqlen){
					if ($score->{'rank'} ) {
						$rank      = $score->{'rank'};
						$min_score = $score->{'score'};
					    $min_locus = $i + 1;
						}
				}
							
				## end of trough or end of sequence, make into feature
                ## if possible
				elsif ($in_trough == 1 && ($score->{'score'} > -0.09 ||
						 $i +1  == $seqlen) ){
					if ($rank != 0) {
                    	push @fts, Bio::SeqFeature::Generic->new (
                         	-start   => $st,
                            -end     => $i +1, #current position
                         	-primary => 'Linker',
							-source  => 'Domcut',
                         	-tag => {
                                  score   => $min_score,
                                  rank    => $rank,
								  residue => $min_locus,
                                 },
                        );
					}
					##and reset parameters ##
					($st, $in_trough, $min_locus, $min_score, $rank) = (0,0,0,0,0);
                }
                $i++;
            }
            return @fts;
        }
        ## convert parsed data into a meta array format
        elsif ($value eq 'meta') {

			## only need to bless  once
			if (! $self->seq->isa("Bio::Seq::MetaI")){
				bless ($self->seq, "Bio::Seq::Meta::Array");
				}
            $self->seq->isa("Bio::Seq::MetaI")
                || $self->throw("$self is not a Bio::Seq::MetaI");
            my $meta_name = "Domcut";

            #test that sequence does not have already a meta seq with same name
            if (grep{$_ eq $meta_name}$self->seq->meta_names ) {
                $self->warn ("$meta_name already exists , not overwriting!");
                next;
            }

            ### or should be an instance variable?? ##
            $Bio::Seq::Meta::Array::DEFAULT_NAME = 'Domcut';
            my @meta = map{$_->{'score'}} @{$self->{'_parsed'}};
            $self->seq->named_meta($meta_name,\@meta );

            # return  seq array object implementing meta sequence #
            return $self->seq;

        }
        #       return ref to array of predictions;
        elsif ($value eq 'parsed') {
            return $self->{'_parsed'};
        }
    }
    #else if no arguments return raw text
    return $self->{'_result'};
}

sub _init {
    my $self = shift;
    $self->url($URL);
    $self->{'_ANALYSIS_SPEC'} = $ANALYSIS_SPEC;
    $self->{'_INPUT_SPEC'}    = $INPUT_SPEC;
    $self->{'_RESULT_SPEC'}   = $RESULT_SPEC;
    $self->{'_ANALYSIS_NAME'} = $ANALYSIS_NAME;
    return $self;
}

sub _run {
    my $self  = shift;
    my $seq_fasta = $self->seq->seq;
    $self->delay(1);
    # delay repeated calls by default by 3 sec, set delay() to change
    $self->sleep;
    $self->status('TERMINATED_BY_ERROR');
    my $rqst = GET $self->url . "&seqnam=". "&sequence=".
                      $seq_fasta. "&outform=dat";

    my $content = $self->request($rqst);
    my $text = $content->content; #1st reponse
    $self->{'_result'} = $text;
    $self->status('COMPLETED') if $text ne '';
}

1;
