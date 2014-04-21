#
# BioPerl module for Bio::Tools::MZEF
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp-at-gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::MZEF - Results of one MZEF run

=head1 SYNOPSIS

   $mzef = Bio::Tools::MZEF->new(-file => 'result.mzef');
   # filehandle:
   $mzef = Bio::Tools::MZEF->new( -fh  => \*INPUT );
   # to indicate that the sequence was reversed prior to feeding it to MZEF
   # and that you want to have this reflected in the strand() attribute of 
   # the exons, as well have the coordinates translated to the non-reversed
   # sequence
   $mzef = Bio::Tools::MZEF->new( -file   => 'result.mzef',
                                  -strand => -1 );

   # parse the results
   # note: this class is-a Bio::Tools::AnalysisResult which implements
   # Bio::SeqAnalysisParserI, i.e., $genscan->next_feature() is the same
   while($gene = $mzef->next_prediction()) {
       # $gene is an instance of Bio::Tools::Prediction::Gene

       # $gene->exons() returns an array of 
       # Bio::Tools::Prediction::Exon objects
       # all exons:
       @exon_arr = $gene->exons();

       # internal exons only
       @intrl_exons = $gene->exons('Internal');
       # note that presently MZEF predicts only internal exons!
   }

   # essential if you gave a filename at initialization (otherwise the file
   # will stay open)
   $mzef->close();

=head1 DESCRIPTION

The MZEF module provides a parser for MZEF gene structure prediction
output.

This module inherits off L<Bio::Tools::AnalysisResult> and therefore
implements L<Bio::SeqAnalysisParserI>.

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

=head1 AUTHOR - Hilmar Lapp

Email hlapp-at-gmx.net (or hilmar.lapp-at-pharma.novartis.com)

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::MZEF;
use strict;

use Bio::Tools::Prediction::Gene;
use Bio::Tools::Prediction::Exon;

use base qw(Bio::Tools::AnalysisResult);

sub _initialize_state {
    my($self,@args) = @_;

    # first call the inherited method!
    my $make = $self->SUPER::_initialize_state(@args);

    # handle our own parameters
    my ($strand, $params) =
	$self->_rearrange([qw(STRAND
			      )],
			  @args);

    # our private state variables
    $strand = 1 unless defined($strand);
    $self->{'_strand'} = $strand;
    $self->{'_preds_parsed'} = 0;
    $self->{'_has_cds'} = 0;
    # array of pre-parsed predictions
    $self->{'_preds'} = [];
}

=head2 analysis_method

 Usage     : $mzef->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
             /mzef/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method { 
#-------------
    my ($self, $method) = @_;  
    if($method && ($method !~ /mzef/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}

=head2 next_feature

 Title   : next_feature
 Usage   : while($gene = $mzef->next_feature()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the MZEF result
           file. Call this method repeatedly until FALSE is returned.

           The returned object is actually a SeqFeatureI implementing object.
           This method is required for classes implementing the
           SeqAnalysisParserI interface, and is merely an alias for 
           next_prediction() at present.

           Note that with the present version of MZEF there will only be one
           object returned, because MZEF does not predict individual genes
           but just potential internal exons.
 Example :
 Returns : A Bio::Tools::Prediction::Gene object.
 Args    :

=cut

sub next_feature {
    my ($self,@args) = @_;
    # even though next_prediction doesn't expect any args (and this method
    # does neither), we pass on args in order to be prepared if this changes
    # ever
    return $self->next_prediction(@args);
}

=head2 next_prediction

 Title   : next_prediction
 Usage   : while($gene = $mzef->next_prediction()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the MZEF result
           file. Call this method repeatedly until FALSE is returned.

           Note that with the present version of MZEF there will only be one
           object returned, because MZEF does not predict individual genes
           but just potential internal exons.
 Example :
 Returns : A Bio::Tools::Prediction::Gene object.
 Args    :

=cut

sub next_prediction {
    my ($self) = @_;
    my $gene;

    # if the prediction section hasn't been parsed yet, we do this now
    $self->_parse_predictions() unless $self->_predictions_parsed();

    # return the next gene structure (transcript)
    return $self->_prediction();
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
    my ($method); # set but not used presently
    my $exon_tag = "InternalExon";
    my $gene;
    # my $seqname; # name given in output is poorly formatted
    my $seqlen;
    my $prednr = 1;

    while(defined($_ = $self->_readline())) {
	if(/^\s*(\d+)\s*-\s*(\d+)\s+/) {
	    # exon or signal
	    if(! defined($gene)) {
		$gene = Bio::Tools::Prediction::Gene->new(
                                       '-primary' => "GenePrediction$prednr",
				       '-source' => 'MZEF');
	    }
	    # we handle start-end first because may not be space delimited
	    # for large numbers
	    my ($start,$end) = ($1,$2);
	    s/^\s*(\d+)\s*-\s*(\d+)\s+//;
	    # split the rest into fields
	    chomp();
	    # format: Coordinates P Fr1 Fr2 Fr3 Orf 3ss Cds 5ss
	    # index:              0   1   2   3   4   5   6   7
	    my @flds = split(' ', $_);
	    # create the feature object depending on the type of signal --
	    # which is always an (internal) exon for MZEF
	    my $predobj = Bio::Tools::Prediction::Exon->new();
	    # set common fields
	    $predobj->source_tag('MZEF');
	    $predobj->significance($flds[0]);
	    $predobj->score($flds[0]); # what shall we set as overall score?
	    $predobj->strand($self->{'_strand'}); # MZEF searches only one
	    if($predobj->strand() == 1) {
		$predobj->start($start);
		$predobj->end($end);
	    } else {
		$predobj->start($seqlen-$end+1);
		$predobj->end($seqlen-$start+1);
	    }
	    # set scores
	    $predobj->start_signal_score($flds[5]);
	    $predobj->end_signal_score($flds[7]);
	    $predobj->coding_signal_score($flds[6]);
	    # frame -- we simply extract the one with highest score from the
	    # orf field, and store the individual scores for now
	    my $frm = index($flds[4], "1");
	    $predobj->frame(($frm < 0) ? undef : $frm);
	    $predobj->primary_tag($exon_tag);
	    $predobj->is_coding(1);
	    # add to gene structure (should be done only when start and end
	    # are set, in order to allow for proper expansion of the range)
	    $gene->add_exon($predobj);		
	    next;
	}
	if(/^\s*Internal .*(MZEF)/) {
	    $self->analysis_method($1);
	    next;
	}
	if(/^\s*File_Name:\s+(\S+)\s+Sequence_length:\s+(\d+)/) {
	    # $seqname = $1; # this is too poor currently (file name truncated
                             # to 10 chars) in order to be sensible enough
	    $seqlen = $2;
	    next;
	}
    }
    # $gene->seq_id($seqname);
    $self->_add_prediction($gene) if defined($gene);
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


1;
