#
# BioPerl module for Bio::Tools::Genscan
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Genscan - Results of one Genscan run

=head1 SYNOPSIS

   use Bio::Tools::Genscan;

   $genscan = Bio::Tools::Genscan->new(-file => 'result.genscan');
   # filehandle:
   $genscan = Bio::Tools::Genscan->new( -fh  => \*INPUT );

   # parse the results
   # note: this class is-a Bio::Tools::AnalysisResult which implements
   # Bio::SeqAnalysisParserI, i.e., $genscan->next_feature() is the same
   while($gene = $genscan->next_prediction()) {
       # $gene is an instance of Bio::Tools::Prediction::Gene, which inherits
       # off Bio::SeqFeature::Gene::Transcript.
       #
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
       # singleton exons: 
       ($single_exon) = $gene->exons();
   }

   # essential if you gave a filename at initialization (otherwise the file
   # will stay open)
   $genscan->close();

=head1 DESCRIPTION

The Genscan module provides a parser for Genscan gene structure prediction
output. It parses one gene prediction into a Bio::SeqFeature::Gene::Transcript-
derived object.

This module also implements the Bio::SeqAnalysisParserI interface, and thus
can be used wherever such an object fits. See L<Bio::SeqAnalysisParserI>.

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

Email hlapp@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Genscan;
use strict;
use Symbol;

use Bio::Root::Root;
use Bio::Tools::Prediction::Gene;
use Bio::Tools::Prediction::Exon;

use base qw(Bio::Tools::AnalysisResult);

my %ExonTags = ('Init' => 'Initial',
		'Intr' => 'Internal',
		'Term' => 'Terminal',
		'Sngl' => '');
    
sub _initialize_state {
    my ($self,@args) = @_;
    
    # first call the inherited method!
    $self->SUPER::_initialize_state(@args);

    # our private state variables
    $self->{'_preds_parsed'} = 0;
    $self->{'_has_cds'} = 0;
    # array of pre-parsed predictions
    $self->{'_preds'} = [];
    # seq stack
    $self->{'_seqstack'} = [];
}

=head2 analysis_method

 Usage     : $genscan->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
             /genscan/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method { 
#-------------
    my ($self, $method) = @_;  
    if($method && ($method !~ /genscan/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}

=head2 next_feature

 Title   : next_feature
 Usage   : while($gene = $genscan->next_feature()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the Genscan result
           file. Call this method repeatedly until FALSE is returned.

           The returned object is actually a SeqFeatureI implementing object.
           This method is required for classes implementing the
           SeqAnalysisParserI interface, and is merely an alias for 
           next_prediction() at present.

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
 Usage   : while($gene = $genscan->next_prediction()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the Genscan result
           file. Call this method repeatedly until FALSE is returned.

 Example :
 Returns : A Bio::Tools::Prediction::Gene object.
 Args    :

=cut

sub next_prediction {
    my ($self) = @_;
    my $gene;

    # if the prediction section hasn't been parsed yet, we do this now
    $self->_parse_predictions() unless $self->_predictions_parsed();

    # get next gene structure
    $gene = $self->_prediction();

    if($gene) {
	# fill in predicted protein, and if available the predicted CDS
	#
	my ($id, $seq);
	# use the seq stack if there's a seq on it
	my $seqobj = pop(@{$self->{'_seqstack'}});
	if(! $seqobj) {
	    # otherwise read from input stream
	    ($id, $seq) = $self->_read_fasta_seq();
	    # there may be no sequence at all, or none any more
	    if($id && $seq) {
		$seqobj = Bio::PrimarySeq->new('-seq' => $seq,
					       '-display_id' => $id,
					       '-alphabet' => "protein");
	    }
	}
	if($seqobj) {
	    # check that prediction number matches the prediction number
	    # indicated in the sequence id (there may be incomplete gene
	    # predictions that contain only signals with no associated protein
	    # and CDS, like promoters, poly-A sites etc)
	    $gene->primary_tag() =~ /[^0-9]([0-9]+)$/;
	    my $prednr = $1;
	    if($seqobj->display_id() !~ /_predicted_\w+_$prednr\|/) {
		# this is not our sequence, so push back for next prediction
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
    my $gene;
    my $seqname;

    while(defined($_ = $self->_readline())) {
	if(/^\s*(\d+)\.(\d+)/) {
	    # exon or signal
	    my $prednr = $1;
	    my $signalnr = $2; # not used presently
	    if(! defined($gene)) {
		$gene = Bio::Tools::Prediction::Gene->new(
                                       '-primary' => "GenePrediction$prednr",
				       '-source' => 'Genscan');
	    }
	    # split into fields
	    chomp();
	    my @flds = split(' ', $_);
	    # create the feature object depending on the type of signal
	    my $predobj;
	    my $is_exon = grep {$_ eq $flds[1];} (keys(%ExonTags));
	    if($is_exon) {
		$predobj = Bio::Tools::Prediction::Exon->new();
	    } else {
		# PolyA site, or Promoter
		$predobj = Bio::SeqFeature::Generic->new();
	    }
	    # set common fields
	    $predobj->source_tag('Genscan');
	    $predobj->score($flds[$#flds]);
	    $predobj->strand((($flds[2] eq '+') ? 1 : -1));
	    my ($start, $end) = @flds[(3,4)];
	    if($predobj->strand() == 1) {
		$predobj->start($start);
		$predobj->end($end);
	    } else {
		$predobj->end($start);
		$predobj->start($end);
	    }
	    # add to gene structure (should be done only when start and end
	    # are set, in order to allow for proper expansion of the range)
	    if($is_exon) {
		# first, set fields unique to exons
		$predobj->start_signal_score($flds[8]);
		$predobj->end_signal_score($flds[9]);
		$predobj->coding_signal_score($flds[10]);
		$predobj->significance($flds[11]);
		$predobj->primary_tag($ExonTags{$flds[1]} . 'Exon');
		$predobj->is_coding(1);
		# Figure out the frame of this exon. This is NOT the frame
		# given by Genscan, which is the absolute frame of the base
		# starting the first predicted complete codon. By comparing
		# to the absolute frame of the first base we can compute the
		# offset of the first complete codon to the first base of the
		# exon, which determines the frame of the exon.
		my $cod_offset;
		if($predobj->strand() == 1) {
		    $cod_offset = $flds[6] - (($predobj->start()-1) % 3);
		    # Possible values are -2, -1, 0, 1, 2. -1 and -2 correspond
		    # to offsets 2 and 1, resp. Offset 3 is the same as 0.
		    $cod_offset += 3 if($cod_offset < 1);		    
		} else {
		    # On the reverse strand the Genscan frame also refers to
		    # the first base of the first complete codon, but viewed
		    # from forward, which is the third base viewed from
		    # reverse.
		    $cod_offset = $flds[6] - (($predobj->end()-3) % 3);
		    # Possible values are -2, -1, 0, 1, 2. Due to the reverse
		    # situation, {2,-1} and {1,-2} correspond to offsets
		    # 1 and 2, resp. Offset 3 is the same as 0.
		    $cod_offset -= 3 if($cod_offset >= 0);
		    $cod_offset = -$cod_offset;
		}
		# Offsets 2 and 1 correspond to frame 1 and 2 (frame of exon
		# is the frame of the first base relative to the exon, or the
		# number of bases the first codon is missing).
		$predobj->frame(3 - $cod_offset);
		# then add to gene structure object
		$gene->add_exon($predobj, $ExonTags{$flds[1]});		
	    } elsif($flds[1] eq 'PlyA') {
		$predobj->primary_tag("PolyAsite");
		$gene->poly_A_site($predobj);
	    } elsif($flds[1] eq 'Prom') {
		$predobj->primary_tag("Promoter");
		$gene->add_promoter($predobj);
	    }
	    next;
	}
	if(/^\s*$/ && defined($gene)) {
	    # current gene is completed
	    $gene->seq_id($seqname);
	    $self->_add_prediction($gene);
	    $gene = undef;
	    next;
	}
	if(/^(GENSCAN)\s+(\S+)/) {
	    $self->analysis_method($1);
	    $self->analysis_method_version($2);
	    next;
	}
	if(/^Sequence\s+(\S+)\s*:/) {
	    $seqname = $1;
	    next;
	}
        
	if(/^Parameter matrix:\s+(\S+)/i) {
	    $self->analysis_subject($1);
	   next;
	}
	
	if(/^Predicted coding/) {
	    $self->_has_cds(1);
	    next;
	}
	/^>/ && do {
	    # section of predicted sequences
	    $self->_pushback($_);
	    last;
	};
    }
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

=head2 _read_fasta_seq

 Title   : _read_fasta_seq()
 Usage   : ($id,$seqstr) = $obj->_read_fasta_seq();
 Function: Simple but specialised FASTA format sequence reader. Uses
           $self->_readline() to retrieve input, and is able to strip off
           the traling description lines.
 Example :
 Returns : An array of two elements.

=cut

sub _read_fasta_seq {
    my ($self) = @_;
    my ($id, $seq);
    local $/ = ">";
    
    my $entry = $self->_readline();
    if($entry) {
	$entry =~ s/^>//;
	# complete the entry if the first line came from a pushback buffer
	while($entry !~ />$/) {
	    last unless $_ = $self->_readline();
	    $entry .= $_;
	}
	# delete everything onwards from an intervening empty line (at the
	# end there might be statistics stuff)
	$entry =~ s/\n\n.*$//s;
	# id and sequence
	if($entry =~ /^(\S+)\n([^>]+)/) {
	    $id = $1;
	    $seq = $2;
	} else {
	    $self->throw("Can't parse Genscan predicted sequence entry");
	}
	$seq =~ s/\s//g; # Remove whitespace
    }
    return ($id, $seq);
}

1;
