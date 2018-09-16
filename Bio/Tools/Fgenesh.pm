#
# BioPerl module for Bio::Tools::Fgenesh
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Christopher Dwan (chris@dwan.org)
#
# Copied, lock stock & barrel from Genscan.pm
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Fgenesh - parse results of one Fgenesh run

=head1 SYNOPSIS

   use Bio::Tools::Fgenesh;

   $fgenesh = Bio::Tools::Fgenesh->new(-file => 'result.fgenesh');
   # filehandle:
   $fgenesh = Bio::Tools::Fgenesh->new( -fh  => \*INPUT );

   # parse the results
   # note: this class is-a Bio::Tools::AnalysisResult which implements
   # Bio::SeqAnalysisParserI, i.e., $fgensh->next_feature() is the same
   while($gene = $fgenesh->next_prediction()) {
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
   $fgenesh->close();

=head1 DESCRIPTION

The Fgenesh module provides a parser for Fgenesh (version 2) gene structure 
prediction output. It parses one gene prediction into a 
Bio::SeqFeature::Gene::Transcript- derived object.

This module also implements the L<Bio::SeqAnalysisParserI> interface, and thus
can be used wherever such an object fits. 

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

=head1 AUTHOR - Chris Dwan

Email chris-at-dwan.org

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Fgenesh;
use strict;
use Symbol;

use Bio::Root::Root;
use Bio::Tools::Prediction::Gene;
use Bio::Tools::Prediction::Exon;

use base qw(Bio::Tools::AnalysisResult);

my %ExonTags = ('CDSf' => 'Initial',
		'CDSi' => 'Internal',
		'CDSl' => 'Terminal',
		'CDSo' => 'Singleton');
    
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
    if($method && ($method !~ /fgenesh/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}

=head2 next_feature

 Title   : next_feature
 Usage   : while($gene = $fgenesh->next_feature()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the Fgenesh result
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
 Usage   : while($gene = $fgenesh->next_prediction()) { ... }
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

	# use the seq stack if there's a seq on it
	my $seqobj = pop(@{$self->{'_seqstack'}});
        my ($id, $seq);
	unless ($seqobj) {
	   ($id, $seq) = $self->_read_fasta_seq();
           my $alphabet;
           if (($id =~ /mrna/) || ($id =~ /cds/)) { $alphabet = 'dna'; }
           else { $alphabet = 'protein'; }
           $seqobj = Bio::PrimarySeq->new('-seq'        => $seq,
                                          '-display_id' => $id,
                                          '-alphabet'   => $alphabet); 
        }
	if ($seqobj) {

	    # check that prediction number matches the prediction number
	    # indicated in the sequence id (there may be incomplete gene
	    # predictions that contain only signals with no associated protein
            # prediction.

	    $gene->primary_tag() =~ /[^0-9]([0-9]+)$/;
	    my $prednr = $1;
	    if ($id !~ /_predicted_(\w+)_$prednr/) {
		# this is not our sequence, so push back for next prediction
		push(@{$self->{'_seqstack'}}, $seqobj);
	    } else {
                if ($1 eq "protein") {
		  $gene->predicted_protein($seqobj);
                } elsif (($1 eq "mrna") || ($1 eq "cds")) {
                  $self->_has_cds(1);
                  $gene->predicted_cds($seqobj);
                  
                  # Have to go back in and get the protein...
                  ($id, $seq) = $self->_read_fasta_seq();
                  if ($id =~ /_cds_/) { 
                    ($id, $seq) = $self->_read_fasta_seq(); 
                  }
 
		  $seqobj = Bio::PrimarySeq->new('-seq' => $seq,
			    		         '-display_id' => $id,
						 '-alphabet' => "protein");
		  $gene->predicted_protein($seqobj);
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

	if(/^\s*(\d+)\s+([+\-])/) {
            my $line = $_;

	    # exon or signal
	    my $prednr = $1;
            my $strand = ($2 eq '+') ? 1 : -1;

	    if(! defined($gene)) {
		$gene = Bio::Tools::Prediction::Gene->new(
                                       '-primary' => "GenePrediction$prednr",
				       '-source' => 'Fgenesh');
	    }
	    # split into fields
	    chomp();
	    my @flds = split(/\s+/, ' ' . $line);
	    ## NB - the above adds leading whitespace before the gene
	    ## number in case there was none (as quick patch to code
	    ## below which expects it but it is not present after 999
	    ## predictions!) This allows >999 predictions to be parsed.

	    # create the feature object depending on the type of signal
	    my $predobj;
	    my $is_exon = grep {$line =~ $_} keys(%ExonTags);
            my ($start, $end);
	    if($is_exon) {
		$predobj = Bio::Tools::Prediction::Exon->new();
                $predobj->score($flds[8]);
                $start   = $flds[5];
                $end     = $flds[7];
	    } else {
		# PolyA site, or TSS 
		$predobj = Bio::SeqFeature::Generic->new();
                $predobj->score($flds[5]);
                $start   = $flds[4];
                $end     = $flds[4];
	    }
	    # set common fields
	    $predobj->source_tag('Fgenesh');
	    $predobj->strand($strand);

# Following tactical commenting-out made by
# malcolm.cook@stowers-institute.org since coordinate reversal is
# apparently vestigial copy/paste detritus from Genscan.pm origins of
# this module and this is NOT needed for fgenesh (at least in v
# 2.1.4).

#	    if($predobj->strand() == 1) {
		$predobj->start($start);
		$predobj->end($end);
#	    } else {
#		$predobj->end($start);
#		$predobj->start($end);
#	    }

            # print STDERR "start $start end $end\n";
	    # add to gene structure (should be done only when start and end
	    # are set, in order to allow for proper expansion of the range)
	    if($is_exon) {
		# first, set fields unique to exons
		$predobj->primary_tag($ExonTags{$flds[4]} . 'Exon');
		$predobj->is_coding(1);
		my $cod_offset;
		if($predobj->strand() == 1) {
		    $cod_offset = ($flds[9] - $predobj->start()) % 3;
		    # Possible values are -2, -1, 0, 1, 2. -1 and -2 correspond
		    # to offsets 2 and 1, resp. Offset 3 is the same as 0.
		    $cod_offset += 3 if($cod_offset < 1);		    
		} else {
		    # On the reverse strand the Genscan frame also refers to
		    # the first base of the first complete codon, but viewed
		    # from forward, which is the third base viewed from
		    # reverse.
		    $cod_offset = ($flds[11] - $predobj->end()) % 3;
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
                # print STDERR "  frame is " . $predobj->frame() . "\n";
		# then add to gene structure object
		$gene->add_exon($predobj, $ExonTags{$flds[1]});		
	    } elsif($flds[3] eq 'PolA') {
		$predobj->primary_tag("PolyAsite");
		$gene->poly_A_site($predobj);
	    } elsif($flds[3] eq 'TSS') {
	        $predobj->primary_tag("Promoter"); # (hey! a TSS is NOT a promoter... what's going on here?...
		$gene->add_promoter($predobj);
                #I'd like to do this (for now):
		#$predobj->primary_tag("TSS"); #this is not the right model, but, it IS a feature at least.
                #but the followg errs out
		#$gene->add_SeqFeature($predobj); #err: MSG: start is undefined when bounds at Bio::SeqFeature::Generic::add_SeqFeature 671 check since gene has no start yet
	    }
	    else {
	      $self->throw("unrecognized prediction line: " . $line);
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

	if(/^(FGENESH)\s+([\d\.]+)/) {
	    $self->analysis_method($1);
	    $self->analysis_method_version($2);
            if (/\s(\S+)\sgenomic DNA/) {
              $self->analysis_subject($1);
            }
	    next;
	}

	if(/^\s*Seq name:\s+(\S+)/) {
	    $seqname = $1;
	    next;
	}
        
	/^Predicted protein/ && do {
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
 Returns : An array of two elements: fasta_id & sequence

=cut

sub _read_fasta_seq {
    my ($self) = @_;
    my ($id, $seq);
    #local $/ = ">";
    
    my $entry = $self->_readline();
    # print " ^^ $entry\n";
    return unless ($entry);
    $entry = $self->_readline() if ($entry =~ /^Predicted protein/);
    # print " ^^ $entry\n";

    # Pick up the header / id.
    if ($entry =~ /^>FGENESH:/) {
      if ($entry =~ /^>FGENESH:\s+(\d+)/) {
         # print STDERR "  this is a predicted gene\n";
         $id  = "_predicted_protein_" . $1;
      } elsif ($entry =~ /^>FGENESH:\[mRNA\]\s*(\d+)/) {
	# print STDERR "  this is an mRNA\n";
         $id  = "_predicted_mrna_" . $1;
      } elsif ($entry =~ /^>FGENESH:\[exon\]\s+Gene:\s*(\d+)/) {
         $id  = "_predicted_cds_"  . $1;
      }
      $seq = "";
      $entry = $self->_readline();
    }

    my $done = 0;
    while (!$done) {
       # print "*** $entry\n";
       if (($entry =~ /^>FGENESH:\[exon\]/) && ($id =~ /^_predicted_cds_/)) {
         # print STDERR "  -- informed about an exon header...\n";
         $entry = $self->_readline();
       } else {
         $seq .= $entry;
         # print STDERR "  Added $entry\n";
       }

       last unless $entry  = $self->_readline();
       if (($entry =~ /^>/) && 
           (!(($entry =~ /^>FGENESH:\[exon\]/) && ($id =~ /^_predicted_cds_/)))) {
          $self->_pushback($entry); last;    
       }
    }

    # id and sequence
    $seq =~ s/\s//g; # Remove whitespace
    return ($id, $seq);
}

1;
