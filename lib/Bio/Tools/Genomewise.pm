#
# BioPerl module for Bio::Tools::Genomewise
#
# Copyright Jason Stajich <jason-at-bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Genomewise - Results of one Genomewise run

=head1 SYNOPSIS

  use Bio::Tools::Genomewise;
  my $gw = Bio::Tools::Genomewise(-file=>"genomewise.out");

  while (my $gene = $gw->next_prediction){
      my @transcripts = $gene->transcripts;
      foreach my $t(@transcripts){
        my @exons =  $t->exons;
        foreach my $e(@exons){
            print $e->start." ".$e->end."\n";
        }
      }
  }

=head1 DESCRIPTION

This is the parser for the output of Genewise. It takes either a file
handle or a file name and returns a
Bio::SeqFeature::Gene::GeneStructure object.  You will need to specify
the proper target sequence id on the object with the
$feature-E<gt>seq_id($seqid).

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

=head1 AUTHOR - Fugu Team, Jason Stajich 

 Email: fugui-at-worf.fugu-sg.org
        jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Genomewise;
use vars qw($Srctag);
use strict;

use Bio::Tools::AnalysisResult;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Gene::Exon;
use Bio::SeqFeature::FeaturePair;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::GeneStructure;

use base qw(Bio::Tools::Genewise);

$Srctag = 'genomewise';

=head2 new

 Title   : new
 Usage   : $obj->new(-file=>"genewise.out");
           $obj->new(-fh=>\*GW);
 Function: Constructor for genomewise wrapper. Takes either a file or filehandle
 Example :
 Returns : L<Bio::Tools::Genomewise>

=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  return $self;
}

=head2 _get_strand

 Title   : _get_strand
 Usage   : $obj->_get_strand
 Function: takes start and end values, swap them if start>end and returns end
 Example :
 Returns :$start,$end,$strand

=cut

=head2 score

 Title   : score
 Usage   : $obj->score
 Function: get/set for score info
 Example :
 Returns : a score value

=cut

=head2 _prot_id

 Title   : _prot_id
 Usage   : $obj->_prot_id
 Function: get/set for protein id 
 Example :
 Returns :a protein id

=cut

=head2 _target_id

 Title   : _target_id
 Usage   : $obj->_target_id
 Function: get/set for genomic sequence id
 Example :
 Returns :a target id

=cut


=head2 next_prediction

 Title   : next_prediction
 Usage   : while($gene = $genewise->next_prediction()) {
                  # do something
           }
 Function: Returns the gene structure prediction of the Genomewise result
           file. Call this method repeatedly until FALSE is returned.

 Example :
 Returns : a Bio::SeqFeature::Gene::GeneStructure object
 Args    :

=cut


sub next_prediction {
    my ($self) = @_;

    my $genes;
    while ($_ = $self->_readline) {
	$self->debug( $_ );
	last if m{^//};

	if( /^Gene\s+\d+\s*$/ ) {
	    $genes = Bio::SeqFeature::Gene::GeneStructure->new
		(-source => $Srctag,
		 -seq_id => $self->_target_id, # if this had been specified
		 );
	    $_ = $self->_readline;
	    $self->debug( $_ );

	    unless ( /^Gene\s+(\d+)\s+(\d+)\s*$/ ) {
		$self->warn("Unparseable genomewise output");
		last;
	    }
	    my $transcript = Bio::SeqFeature::Gene::Transcript->new
		(-source => $Srctag,
		 -seq_id => $self->_target_id, # if this had been specified
		 -start  => $1,
		 -end    => $2,
		 );
	    my $nbr = 1;
	    while( $_ = $self->_readline ) {    
		$self->debug( $_ );

		unless( m/^\s+Exon\s+(\d+)\s+(\d+)\s+phase\s+(\d+)/ ){
		    $self->_pushback($_);
		    last;
		}
		my ($e_start,$e_end,$phase,$e_strand) = ($1,$2,$3);
		
		($e_start,$e_end,$e_strand) = $self->_get_strand($e_start,
								 $e_end);
		$transcript->strand($e_strand) unless $transcript->strand != 0;
		
		my $exon = Bio::SeqFeature::Gene::Exon->new 
		    (-seq_id=>$self->_target_id,
		     -source => $Srctag,
		     -start=>$e_start, 
		     -end=>$e_end, 
		     -frame => $phase,
		     -strand=>$e_strand);
		$exon->add_tag_value("Exon",$nbr++);
		$exon->add_tag_value('phase',$phase);
		$transcript->add_exon($exon);
	    }
	    $genes->add_transcript($transcript);
	    last; # only process a single gene at a time
	}
    }
    return $genes;
}

1;
