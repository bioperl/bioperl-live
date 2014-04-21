# $Id: Sopma.pm,v 1.0 2003/07/ 11
#
# BioPerl module for Bio::Tools::Analysis::Protein::Sopma
#
# Copyright Richard Adams
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code


=head1 NAME

Bio::Tools::Analysis::Protein::Sopma - a wrapper around the
Sopma protein secondary structure prediction server

=head1  SYNOPSIS

  use Bio::Tools::Analysis::Protein::Sopma;
  #get a Bio::Seq or Bio::PrimarySeq
  my $seq;

  my $sopma = Bio::Tools::Analysis::Protein::Sopma->new
      (-seq=>$seq, states=>4);
  $sopma->run;
  print $sopma->result;# #raw text to standard error

=head1  DESCRIPTION

A module to remotely retrieve predictions of protein secondary
structure.  Each residue in the protein receives a score representing
the likelihood of existing in each of four different states (helix,
coil, turn or sheet), e.g.,

  my $analysis_object = Bio::Tools::SimpleAnalysis::Protein::Sopma->new
      ( -seq          => $seq,
        -states       => 4,
        -window_width => 15,
      );

creates a new object.  Compulsory argument -seq.  Optional arguments
-states, -window_width,-similarity_threshold. These arguments can also be
set by direct methods , e.g.,

  $analysis_object->states(4);
  $analysis_object->run;

submits the query to the server and obtains raw text output. Given an
amino acid sequence the results can be obtained in 4 formats,
determined by the argument to the result method:

=over 4

=item 1

The raw text of the program output.

  my $rawdata = $analysis_object->result;

=item 2

A reference to an array of hashes of scores for each state and the
assigned state.

  my $data_ref = $analysis_object->result('parsed');
  print "score for helix at residue 2 is $data_ref->[1]{'helix'}\n";
  print "predicted struc  at residue 2 is $data_ref->[1]{'struc}\n";

Hash keys are 'helix', 'struc', 'sheet', 'coil', 'turn'.

=item 3

An array of Bio::SeqFeature::Generic objects where each feature is a
predicted unit of secondary structure. Only stretches of helix/sheet
predictions for longer than 4 residues are defined as helices/sheets.

  my @fts = $analysis_object->result(Bio::SeqFeatureI);
  for my $ft (@fts) {
      print " From ",  $ft->start, " to  ",$ft->end, " struc: " ,
             ($ft->each_tag_value('type'))[0]  ,"\n";
  }

=item 4

A Bio::Seq::Meta::Array implementing sequence.

This is a Bio::Seq object that can also hold data about each residue
in the sequence.  In this case, the sequence can be associated with a
arrays of Sopma prediction scores.  e.g.,

  my $meta_sequence = $analysis_object->result('meta');
  print "scores from residues 10 -20 are ",
      $meta_sequence->named_submeta_text("Sopma_helix",10,20), "\n";

Meta sequence names are : Sopma_helix, Sopma_sheet, Sopma_turn,
Sopma_coil, Sopma_struc, representing the scores for each residue.

Many methods common to all analyses are inherited from
Bio::Tools::Analysis::SimpleAnalysisBase.

=back

=head1 SEE ALSO

L<Bio::SimpleAnalysisI>, 
L<Bio::Tools::Analysis::SimpleAnalysisBase>
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

=cut

use strict;

package Bio::Tools::Analysis::Protein::Sopma;

use IO::String;
use Bio::SeqIO;
use HTTP::Request::Common qw (POST);
use Bio::SeqFeature::Generic;
use Bio::Seq::Meta::Array;


use base qw(Bio::Tools::Analysis::SimpleAnalysisBase);

#extends array for 2struc.
my $URL = 'http://npsa-pbil.ibcp.fr/cgi-bin/secpred_sopma.pl';
my $ANALYSIS_NAME= 'Sopma';
my $ANALYSIS_SPEC= {name => 'Sopma', type => 'Protein'};
my $INPUT_SPEC = [
                  {mandatory=>'true',
                   type     => 'Bio::PrimarySeqI',
                   'name'   => 'seq',
                  },
                  {mandatory =>'false',
                   type      => 'integer',
                   name      => 'similarity_threshold',
                   default   => 8,
                  },
                  {mandatory  =>'false',
                   type       => 'integer',
                   name       => 'window_width',
                   default    => 17,
                  },
                  {mandatory  =>'false',
                   type       => 'integer',
                   name       => 'states',
                   default    => 4,
                  },
                 ];
my  $RESULT_SPEC =
    {
     ''   => 'bulk',              # same as undef
     raw  => '[{struc=>, helix=>, turn=>, coil=>, sheet=>}]',
     meta => 'Bio::Seq::Meta::Array object',
     'Bio::SeqFeatureI' => 'ARRAY of Bio::SeqFeature::Generic',
    };
use constant MIN_STRUC_LEN => 3; 


=head2  similarity_threshold

  Useage  : $job->similarity_threshold(...)
  Returns : The  similarity threshold used in the analysis
  Args    : None (retrieves value) or  an integer (default = 8) 
            that sets the similarity threshold .

This method gets/sets the  similarity threshold for the prediction.

=cut

sub similarity_threshold {
    my ($self, $value) = @_;
    if ($value) {
        $self->throw ("similarity_threshold must be integer")
            unless $value =~ /^\d+$/;
        $self->{'_similarity_threshold'} = $value;
    }
    $self->{'_similarity_threshold'} ||= $self->input_spec->[1]{'default'};
    return $self->{'_similarity_threshold'};
}

=head2  window_width

  Usage    : $job->window_width(...)
  Returns  : The window width used in the analysis
  Args     : None (retrieves value) or  an integer (default = 17)
             that sets the window width.

This method gets/sets the window width for the prediction, .  If
attempted to set longer than the sequence, warns of error.

=cut

sub window_width {
    my ($self, $value) = @_;
    if ($value) {
        $self->throw ("window_width must be integer")
            unless $value =~ /^\d+$/;
        $self->{'_window_width'} = $value;
    }
    $self->{'_window_width'} ||= $self->input_spec->[2]{'default'};
    $self->warn ("window width longer than sequence!")
        unless $self->{'_window_width'} < $self->seq->length;
    return $self->{'_window_width'};
}

=head2  states

  Usage    : $job->states(...)
  Returns  : The number of secondary structure prediction states
  Args     : None (retrieves value) or either '3' or '4' to set
             prior to running analysis.

This method gets/sets the number of states for the prediction, either
3 or 4 (includes turns).

=cut

sub states {
    my ($self, $value) = @_;
    if ($value) {
        $self->throw ("number of states must be 3 or 4")
            unless $value == 3 or $value ==4;
        $self->{'_states'} = $value;
    }
    $self->{'_states'} ||= $self->input_spec->[3]{'default'};
    return $self->{'_states'};
}

=head2 result

  Usage   : $job->result (...)
  Returns : a result created by running an analysis
  Args    : various

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
tag is "2ary".  Feature tags are "type" (which can be helix, sheet
coil, or turn if 4 state prediction requested) "method" (Sopma)

=item 'parsed'

Array of hash references of scores/structure assignations 
{ helix =E<gt> , sheet =E<gt> , coil =E<gt> , struc=E<gt>}.

=item 'all'

A Bio::Seq::Meta::Array object. Scores can be accessed using methods
from this class. Meta sequence names are Sopma_helix, Sopma_sheet,
Sopma_coil, Sopma_turn (if defined), and Sopma_struc.


=back


=cut

sub result {
    my ($self,$value, $run_id) = @_;

    my @score;
    my @fts;

    if ($value ) {
        if (!exists($self->{'_parsed'} )) {
            my $result = IO::String->new($self->{'_result'});
            while (my $line = <$result>) {
                next unless $line =~ /^[HCET]\s/; # or for sopma/hnn  /^[A-Z]\s/
                $line =~/^([A-Z])\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/; # or for so
                push @score, { struc => $1,
                               helix => $2,
                               sheet => $3,
                               coil => $5,
                             };
                #include turn if 4states are requested
                $score[$#score]{'turn'} = $4 if $self->states == 4;
                #can optimize by duplicating code here
            }
            $self->{'_parsed'} = \@score;
        }
        if ($value eq 'Bio::SeqFeatureI') {
            $self->_get_2ary_coords();
            for my $type (keys %{$self->{'_parsed_coords'}} ) {
                next if $type =~  /\w{2,}/; #if not H,C,E or T

				## these 2 are added to distinguish features on same
               ## sequence run with different params
				my $tag_hash = {
								type   => $type,
                                method => $self->analysis_name,
								};
				$self->_add_params_to_result($tag_hash);

				## now make feature object
                for my $loc (@{$self->{'_parsed_coords'}{$type}} ) {
                    push  @fts,   Bio::SeqFeature::Generic->new
                        (-start   => $loc->{'start'},
                         -end     => $loc->{'end'},
                         -source  => 'Sopma',
                         -primary => 'Domain',
                         -tag => $tag_hash,
                                 );
                }               #end of array of strucs of type
            }                   # end of all 2nd struc elements
            delete $self->{'_parsed_coords'}; #remove temp data
            return @fts;
        }                       #endif BioSeqFeature

        elsif ($value eq 'meta') {
            #1st of all make 3 or 4 arrays of scores for each type from column data
            my %type_scores;
            for my $aa (@{$self->{'_parsed'}}) {
                for my $type (qw(struc helix sheet  coil)) {
                    push @{$type_scores{$type}}, $aa->{$type};
                }
                push @{$type_scores{'turn'}}, $aa->{'turn'} if  exists $aa->{'turn'};
            }
			
			## convert to meta sequence array ##
			if (!$self->seq->isa("Bio::Seq::Meta::Array")) {
           		 bless ($self->seq, "Bio::Seq::Meta::Array");
				}
            $self->seq->isa("Bio::Seq::MetaI")
                || $self->throw("$self is not a Bio::Seq::MetaI");


            $Bio::Seq::Meta::Array::DEFAULT_NAME = 'Sopma_struc';
            for my $struc_type (keys %type_scores) {
                my $meta_name = "Sopma". "_" . "$struc_type";
				if ($run_id) {
					$meta_name .= "|$run_id";
				}
                my @meta = map{$_->{$struc_type}} @{$self->{'_parsed'}};
                if (grep{$_ eq $meta_name}$self->seq->meta_names >0) {
                    $self->warn ("$meta_name already exists , not overwriting!");
                    next;
                }
                $self->seq->named_meta($meta_name,\@meta );
            }
            # return  seq array object implementing meta sequence #
            return $self->seq;

        }
		## else return parsed data if $value is defined
		 else {
            return $self->{'_parsed'};
        }

    }                           #endif ($value)
    #return raw result if no return format stated
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

sub _get_2ary_coords {
    #helper sub for result;
    ##extracts runs of structure > MIN_STRUC_LENresidues or less if Turn:
    #i.e., helical prediction for 1 residue isn't very meaningful...
    ## and poulates array of hashes with start/end values.
    ##keys of $Result are 'H' 'T' 'C' 'E'. 
    my ($self) = @_;
    my @prot = @{$self->{'_parsed'}};
    my %Result;
    for (my $index = 0; $index <= $#prot; $index++) {

        my $type        = $prot[$index]{'struc'};
        next unless $type && $type =~ /[HTCE]/;
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

sub  _run {
    my $self  = shift;
    $self->delay(1);
    # delay repeated calls by default by 3 sec, set delay() to change
    $self->sleep;
    $self->status('TERMINATED_BY_ERROR');
    my $request = POST 'http://npsa-pbil.ibcp.fr/cgi-bin/secpred_sopma.pl',
        Content_Type => 'form-data',
            Content  => [title     => "",
                         notice    => $self->seq->seq,
                         ali_width => 70,
                         states    => $self->states,
                         threshold => $self->similarity_threshold ,
                         width     => $self->window_width,
                        ];

    my $text = $self->request($request)->content;
    return $self unless $text;

    #### get text only version of results ## 
    my ($next) = $text =~ /Prediction.*?=(.*?)>/;
    my $out    = "http://npsa-pbil.ibcp.fr/". "$next";
    my $req2   = HTTP::Request->new(GET=>$out);
    my $resp2  = $self->request ($req2);
    $self->{'_result'} = $resp2->content;
    $self->status('COMPLETED') if $resp2 ne '';
    return $self;
}

sub _add_params_to_result{
	## called when making Seqfeature objects
	my ($self, $tag_hash) = @_;
	my $hash;
	## adds input parameter values to SeqFeatureI results where multiple
    ##  parameter values are possible. Only adds value if not default. 
	map{$hash->{$_->{'name'}} = $_}@{$self->input_spec()};

	for my $p (keys %$hash) {
		if (!ref($self->$p) && $self->$p ne $hash->{$p}{'default'}) {
			$tag_hash->{$p} = $self->$p;
		}
	}
				 
}


1;
