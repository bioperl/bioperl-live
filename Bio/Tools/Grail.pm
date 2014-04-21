#
# BioPerl module for Bio::Tools::Grail
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Grail - Results of one Grail run

=head1 SYNOPSIS

   $grail = Bio::Tools::Grail->new(-file => 'result.grail');
   # filehandle:
   $grail = Bio::Tools::Grail->new( -fh  => \*INPUT );

   # parse the results
   while($gene = $grail->next_prediction()) {
       # $gene is an instance of Bio::Tools::Prediction::Gene

       # $gene->exons() returns an array of 
       # Bio::Tools::Prediction::Exon objects
       # all exons:
       @exon_arr = $gene->exons();

       # initial exons only
       @init_exons = $gene->exons('Initial');
       # internal exons only
       @intrl_exons = $gene->exons('Internal');
       # terminal exons only
       @term_exons = $gene->exons('Terminal');
       # singleton exons only -- should be same as $gene->exons() because
       # there are no other exons supposed to exist in this structure
       @single_exons = $gene->exons('Single');
   }

   # essential if you gave a filename at initialization (otherwise the file
   # will stay open)
   $genscan->close();

=head1 DESCRIPTION

The Grail module provides a parser for Grail gene structure prediction
output.


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

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::Grail;
use strict;

use Bio::Tools::Prediction::Gene;
use Bio::Tools::Prediction::Exon;
use Symbol;

use base qw(Bio::Root::IO Bio::Root::Root);

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->_initialize_io(@args);

  return $self;
}

=head2 next_prediction

 Title   : next_prediction
 Usage   : while($gene = $grail->next_prediction()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the Grail result
           file. Call this method repeatedly until FALSE is returned.

 Example :
 Returns : A Bio::Tools::Prediction::Gene object.
 Args    :

=cut

sub next_prediction {
    my ($self) = @_;
    
    # get next gene structure
    my $gene = $self->_prediction();

    if($gene) {
	# fill in predicted protein, and if available the predicted CDS
	#
	my ($id, $seq);
	# use the seq stack if there's a seq on it
	my $seqobj = pop(@{$self->{'_seqstack'}});
	if(! $seqobj) {
	    # otherwise read from input stream
	    ($id, $seq) = $self->_read_fasta_seq();
	    $seqobj = Bio::PrimarySeq->new('-seq' => $seq,
					   '-display_id' => $id,
					   '-alphabet' => "protein");
	}
	# check that prediction number matches the prediction number
	# indicated in the sequence id (there may be incomplete gene
	# predictions that contain only signals with no associated protein
	# and CDS, like promoters, poly-A sites etc)
	$gene->primary_tag() =~ /[^0-9]([0-9]+)$/;
	my $prednr = $1;
	if($seqobj->display_id() !~ /_predicted_\w+_$prednr\|/) {
	    # this is not our sequence, so push back for the next prediction
	    push(@{$self->{'_seqstack'}}, $seqobj);
	} else {
	    $gene->predicted_protein($seqobj);
	    # CDS prediction, too?
	    if($self->_has_cds()) {
		($id, $seq) = $self->_read_fasta_seq();
		$seqobj = Bio::PrimarySeq->new('-seq' => $seq,
					       '-display_id' => $id,
					       '-alphabet' => "dna");
		$gene->predicted_cds($seqobj);
	    }
	}
    }
    return $gene;
}

=head2 _parse_predictions

 Title   : _parse_predictions()
 Usage   : $obj->_parse_predictions()
 Function: Parses the prediction section. Automatically called by
           next_prediction() if not yet done.
 Example :
 Returns : 

=cut

sub _parse_predictions {
    my ($self) = @_;

    # code needs to go here
    
    $self->_predictions_parsed(1);
}

=head2 _prediction

 Title   : _prediction()
 Usage   : $gene = $obj->_prediction()
 Function: internal
 Example :
 Returns : 

=cut

sub _prediction {
    my ($self) = @_;

    return unless(exists($self->{'_preds'}) && @{$self->{'_preds'}});
    return shift(@{$self->{'_preds'}});
}

=head2 _add_prediction

 Title   : _add_prediction()
 Usage   : $obj->_add_prediction($gene)
 Function: internal
 Example :
 Returns : 

=cut

sub _add_prediction {
    my ($self, $gene) = @_;

    if(! exists($self->{'_preds'})) {
	$self->{'_preds'} = [];
    }
    push(@{$self->{'_preds'}}, $gene);
}

=head2 _predictions_parsed

 Title   : _predictions_parsed
 Usage   : $obj->_predictions_parsed
 Function: internal
 Example :
 Returns : TRUE or FALSE

=cut

sub _predictions_parsed {
    my ($self, $val) = @_;

    $self->{'_preds_parsed'} = $val if $val;
    if(! exists($self->{'_preds_parsed'})) {
	$self->{'_preds_parsed'} = 0;
    }
    return $self->{'_preds_parsed'};
}

=head2 _has_cds

 Title   : _has_cds()
 Usage   : $obj->_has_cds()
 Function: Whether or not the result contains the predicted CDSs, too.
 Example :
 Returns : TRUE or FALSE

=cut

sub _has_cds {
    my ($self, $val) = @_;

    $self->{'_has_cds'} = $val if $val;
    if(! exists($self->{'_has_cds'})) {
	$self->{'_has_cds'} = 0;
    }
    return $self->{'_has_cds'};
}

1;
