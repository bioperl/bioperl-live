# $ Id $
#
# bioperl module for Bio::LiveSeq::Mutator
#
# Cared for by Joseph Insana <insana@ebi.ac.uk> <jinsane@alaric.atnet.it>
#
# Copyright Joseph Insana
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME


  Bio::LiveSeq::Mutator - Package mutating LiveSequences

=head1 SYNOPSIS


=head1 DESCRIPTION

  This package holds a series of methods that can be used to mutate LiveSequence
  objects. The methods will generally require a mutmatrix (mutation matrix), i.e.
  a series of change() parameters describing subsequent mutations to be
  performed on a specific Gene or on other LiveSequence entities.
  Example of a MutMatrix:

     0 32 4
    ac 26 2
  gatt 15 0


  meaning:
     delete 4 nucleotides in the Gene's 1st Transcript (e.g.) departing from
       the one at positions 32 (according to Transcript's numbering)
     double point mutation of nucleotides 26 and 27 into "a" and "c"
     insertion of "gatt" before nucleotide at position 15

  The positions are checked all at the beginning, so they will not point to
  different nucleotides after indels has been made.

  I.e. the positions will be changed to "labels" internally and the labels
  will be used to issue the mutations in the LiveSequence objects.


=head1 FEEDBACK

=head2 Mailing Lists

  User feedback is an integral part of the evolution of this and other
  Bioperl modules. Send your comments and suggestions preferably to one
  of the Bioperl mailing lists.  Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

  report bugs to the Bioperl bug tracking system to help us keep track
  the bugs and their resolution.  Bug reports can be submitted via
  email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Joseph A.L. Insana

  Email:  Insana@ebi.ac.uk, jinsane@alaric.atnet.it

  Address: 

     EMBL Outstation, European Bioinformatics Institute
     Wellcome Trust Genome Campus, Hinxton
     Cambs. CB10 1SD, United Kingdom 

=head1 APPENDIX

  The rest of the documentation details each of the object
  methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LiveSeq::Mutator;
$VERSION=2.74;

# Version history:
# Thu Apr 13 18:09:55 BST 2000 v 1.0 begun
# Thu Apr 13 19:15:28 BST 2000 v 1.1 first tentative version, to be tested
# Fri Apr 14 03:13:18 BST 2000 v 1.2 debugged; generalized change_transcript to change_object
# Fri Apr 14 14:51:21 BST 2000 v 2.0 change_gene created
# Wed May  3 10:32:19 BST 2000 v 2.1 _pos2label_matrix now returns also a modified posmutmatrix with incorrect lines from input posmutmatrix deleted
# Wed May  3 10:46:06 BST 2000 v 2.2 change_* will now call create_mut_objs every time a change is issued -> this will allow complete mutation analysis/tracking
# Thu May  4 14:39:36 BST 2000 v 2.22 added before_mutation variable, i.e. possibility of passing information from create_mut*_before to create_mut*_after
# Tue May  9 15:21:12 BST 2000 v 2.3 put some code in create_* and in change_gene for creating SeqDiff, DNAMutation and RNAChange objects
# Thu May 11 17:11:10 BST 2000 v 2.31 started creating AAChange objects
# Tue May 16 17:58:09 BST 2000 v 2.4 praelabel and postlabel calculated (limits of the change) / _is_rna_affected coded (workaround in it not coded yet)
# Tue May 16 18:19:10 BST 2000 v 2.42 up/down streamseq code for RNAChange
# Wed May 17 16:59:34 BST 2000 v 2.5 setting codon_pos and exons_modified informations
# Mon May 22 14:52:19 BST 2000 v 2.52 aachange->allele_ori added
# Wed Jun  7 18:05:43 BST 2000 v 2.6 Heikki's update with code from Analyzer
# Fri Jun 23 15:45:00 BST 2000 v 2.61 Heikki's update
# Thu Jun 29 16:45:56 BST 2000 v 2.64 dnamut_start now decremented if position is before ATG
# LOT of changes in between by Heikki
# Wed Sep 20 16:26:20 BST 2000 v 2.7 to fix bug in factor7 (point mut g 65), at the beginning of an exon (the second), introduced the new variable RNApraelabel. Probably a RNApostlabel is also required. Maybe the "affected exon" does not work correctly when mutating at edges of an exon!!!!
# Mon Sep 25 14:56:15 BST 2000 v 2.74a tentatively added RNApostlabel but now this is alphaversion, has to be fully tested for possible inconsistency
# Mon Sep 25 15:06:43 BST 2000 v 2.74b decided to take away RNApostlabel fix because it creates other problems and is not the real fix. AAChange now uses RNApraelabel and works.

use strict;
use Carp qw(cluck croak carp);
use vars qw($VERSION @ISA);
use Bio::Variation::SeqDiff;
use Bio::Variation::DNAMutation;
use Bio::Variation::RNAChange;
use Bio::Variation::AAChange;
use Bio::Variation::Allele;

=head2 change_object

 Title   : change_object
 Usage   : $results = Bio::LiveSeq::Mutator::change_object(-object =>
                                                               $transcript,
                                                -mutmatrix => \@mutmatrix)
 Function: Returns a reference to a hash containing the results of the changes
           performed according to the instructions present in the mutmatrix,
           on the object (e.g. a transcript) specified.
           The numbering the positions in the MutMatrix must refer to is
           the numbering according to the object (i.e. entry numbering for
           DNA objects, transcript/codingregion numbering for Transcript
           objects).
 Args    : a reference to a LiveSeq object
           a reference to a mutmatrix array of arrays
 Returns : a reference to a hash
 Errorcode: 0

=cut

sub change_object {
    my %args = @_;
    my ($object,$mutmatrix)=($args{-object},$args{-mutmatrix});
    _change_object($object,$mutmatrix);
}

sub _change_object {
    my ($obj,$mutmatrix)=@_;
    my %results;

    unless (($obj)&&(ref($obj->{seq}) eq "Bio::LiveSeq::DNA")) {
	carp "Input object $obj, does not have a LiveSeq::DNA attached";
	return (0);
    }
    unless (($mutmatrix)&&(ref($mutmatrix) eq "ARRAY")) {
	carp "No mutmatrix array given";
	return (0);
    }

    my ($labelmutmatrix,$newposmutmatrix)=_pos2label_matrix($obj,$mutmatrix);
    my @labelmutmatrix=@{$labelmutmatrix};
    my @posmutmatrix=@{$newposmutmatrix};
    
    unless (@labelmutmatrix) { # if empty matrix
	carp "No changes to be performed. Errors in the mutmatrix!";
	return (0);
    }

    my $translation=$obj->get_Translation;
    if ($translation) { # if the object has a translation
	$results{ori_transl_seq}=$translation->seq;
    }
    $results{ori_dna_seq}=$obj->{seq}->seq;
    $results{ori_obj_seq}=$obj->seq;

    # issue changes
    my $k=0;
    my $changearray; my $before_mutation;
    foreach $changearray (@labelmutmatrix) {
	$before_mutation=create_mut_objs_before($obj,$changearray,$posmutmatrix[$k]);
	$obj->labelchange(@{$changearray});
	create_mut_objs_after($obj,$changearray,$posmutmatrix[$k],$before_mutation);
	$k++;
    }
    if ($translation) { # if the object has a translation
	$results{mut_transl_seq}=$translation->seq;
    }
    $results{mut_dna_seq}=$obj->{seq}->seq;
    $results{mut_obj_seq}=$obj->seq;

    return \%results;
}


# arguments: reference to LiveSeq object, reference to a position mutmatrix
# returns: reference to label mutmatrix, reference to new position mutmatrix
#          i.e. to a posmutmatrix deleted of incorrect lines
sub _pos2label_matrix {
    my $obj=$_[0];
    my @posmutmatrix=@{$_[1]};
    my @labelmutmatrix; my @newposmutmatrix;
    my @labels; my @changearray;
    my ($changearray,$label,$newseq,$position,$length);
    foreach $changearray (@posmutmatrix) {
	@changearray=@{$changearray};
	if((scalar(@changearray) < 2)||(scalar(@changearray) > 3)) {
	    carp "Error in mutmatrix syntax, ignoring the line -> @changearray\n";
	} else { # change the following if new columns added in mutmatrix
	    ($newseq,$position,$length)=@changearray;
	    $label=$obj->label($position); # get the label for the position
	    unless ($label > 0) { # label not found or error
		carp "No valid label found at the position defined in the mutmatrix!\nIgnoring the line -> @changearray\n";
	    } else {
		push(@labelmutmatrix,[$newseq,$label,$length]);
		push(@newposmutmatrix,$changearray);
	    }
	}
    }
    return (\@labelmutmatrix,\@newposmutmatrix);
}

=head2 change_gene

 Title   : change_gene
 Usage   : $results = Bio::LiveSeq::Mutator::change_gene(-gene => $gene,
                                                -mutmatrix => \@mutmatrix,
                                                -numbering => "entry")
 Function: Returns a reference to a hash containing the results of the changes
           performed according to the instructions present in the mutmatrix.
           The -numbering argument decides what is to be changed and what
           numbering scheme is to be used.
            -numbering => "entry" will change at the DNA level, using the
                          numbering of the database entry.
            -numbering => "coding" will change at the cDNA level, using the
                          numbering of the 1st Transcript
           Alternative transcripts can be used by specifying "coding 2" or
           "coding 3"....
 Args    : a reference to a LiveSeq object
           a reference to a mutmatrix array of arrays
           a string specifying a numbering scheme
 Returns : Bio::Variation::SeqDiff object
 Errorcode: 0

=cut

sub change_gene {
    my %args = @_;
    my ($gene,$mutmatrix,$numbering)=($args{-gene},$args{-mutmatrix},$args{-numbering});
    my $tmp;
    #print  `date`, " start\n";
    unless (($gene)&&(ref($gene) eq "Bio::LiveSeq::Gene")) {
	carp "Input object $gene is not a LiveSeq::Gene";
	return (0);
    }
    unless (($mutmatrix)&&(ref($mutmatrix) eq "ARRAY")) {
	carp "No mutmatrix array given";
	return (0);
    }

    unless ($numbering) {
	carp "Numbering scheme not specified, defaulting to 1st Transcript";
	$numbering="coding";
    }

    # SeqDiff object creation
    my $seqDiff  =  Bio::Variation::SeqDiff->new();
    $seqDiff->moltype($gene->get_DNA->moltype);
    my @transcripts=@{$gene->get_Transcripts};
    my $obj;
    if (($numbering eq "coding")||($numbering eq "coding 1")) {
	$obj=$transcripts[0];
    } elsif ($numbering eq "entry") {
	$obj=$transcripts[0]->{seq};
    } elsif (index($numbering,"coding ") == 0) { # string begins with "coding "
	my $transnumber=$numbering;
	$transnumber =~ s/coding //;
	$transnumber -= 1; # 1 -> 0, 2 -> 1
	if (($transnumber >= 0)&&($transnumber <= $#transcripts)) {
	    $obj=$transcripts[$transnumber];
	} else {
	    carp "The alternative Transcript - n.",$transnumber+1,
	    "- requested for the gene does not exist. Reverting to the 1st Transcript";
	    $obj=$transcripts[0];
	}
    } else {
	carp "Numbering scheme $numbering not known....";
	return (0);
    }
    #print  `date`, " seqDiff done\n";
    my ($labelmutmatrix,$newposmutmatrix)=_pos2label_matrix($obj,$mutmatrix);
    my @labelmutmatrix=@{$labelmutmatrix};
    my @posmutmatrix=@{$newposmutmatrix};
    unless (@labelmutmatrix) { # if empty matrix
	carp "No changes to be performed. Errors in the mutmatrix!";
	return (0);
    }
    #print  `date`, " seqDiff mutmatrix converted\n";
    my ($DNAobj, $RNAobj,$coord_system,$atg_label,$atg_offset);
    
    if ($obj->isa("Bio::LiveSeq::Transcript")) {
	$coord_system="coding";
	$RNAobj=$obj;
	$DNAobj=$obj->{seq};
	$seqDiff->rna_ori( $obj->seq );
	$tmp = $obj->get_Translation->seq;
	$seqDiff->aa_ori($obj->get_Translation->seq);
    } else {
	$coord_system="entry";
	$DNAobj=$obj;
	$RNAobj=$transcripts[0];
	$seqDiff->rna_ori($RNAobj->seq);
	$seqDiff->aa_ori($RNAobj->get_Translation->seq);
    }
    # set the numbering scheme used
    $seqDiff->numbering($coord_system);

    # set the accession number to the SeqDiff object
    $seqDiff->id($DNAobj->accession_number);

    $atg_label=$RNAobj->start;
    # the atg_offset takes in account that DNA object could be a subset of the
    # whole entry (via the light_weight loader)
    $atg_offset=$DNAobj->position($atg_label)+($DNAobj->start)-1;
    $seqDiff->offset($atg_offset-1);

    my $cds_end = $gene->downbound - $atg_offset + 1;
    $seqDiff->cds_end($cds_end);
    #print "offset: $atg_offset label: $atg_label \n"; # debug

    $seqDiff->dna_ori($DNAobj->seq);
    $seqDiff->rna_ori($RNAobj->seq);
    $DNAobj->coordinate_start($atg_label);
    
    # need to add more than one rna & aa
    #foreach $transcript (@transcripts) {
    #  $seqDiff{"ori_transcript_${i}_seq"}=$transcript->seq;
    #  $seqDiff{"ori_translation_${i}_seq"}=$transcript->get_Translation->seq;
    #}
    #print  `date`, " starting to issue changes\n";
    # issue changes
    my $k=0;
    my ($changearray,$before_mutation,$transpos);
    foreach $changearray (@labelmutmatrix) {
	if ($coord_system eq "coding") {
	    $transpos=$posmutmatrix[$k]->[1]; # transpos given by user
	} else {
	    $transpos=$RNAobj->position($changearray->[1],$atg_label); # transpos of label / It will be 0 if mutating an intron, negative if upstream of ATG
	}
	#print  `date`, " create_mutobj_before\n";

	$before_mutation=create_mut_objs_before
	    ($k+1,$seqDiff,$changearray,$transpos,
	     $DNAobj,$RNAobj,$atg_label,$coord_system);

	#print  `date`, " create_mutobj_before done, only labelchange left\n";

	$obj->labelchange(@{$changearray});

	#print  `date`, " labelchange done\n";

	create_mut_objs_after
	    ($k+1,$seqDiff,$changearray,$transpos,$before_mutation,
	     $DNAobj,$RNAobj,$atg_label,$coord_system);
	$k++;
    }
    #print  `date`, " getting mutated sequences\n";
    $seqDiff->dna_mut($DNAobj->seq);
    $seqDiff->rna_mut($RNAobj->seq);
    $seqDiff->aa_mut($obj->get_Translation->seq);
    #print  `date`, " getting mutated sequences is done\n";
    #$seqDiff{mut_dna_seq}=$gene->get_DNA->seq;
    #my $i=1;
    #foreach $transcript (@transcripts) {
    #  $seqDiff{"mut_transcript_${i}_seq"}=$transcript->seq;
    #  $seqDiff{"mut_translation_${i}_seq"}=$transcript->get_Translation->seq;
    #}
    #print `date`, "done\n";
  return $seqDiff;
}



# args: reference to label changearray, reference to position changearray
# Function: take care of the creation of mutation objects, with (if needed)
# information BEFORE the change takes place
sub create_mut_objs_before {

### Variable parsing
    my ($mut_number,$seqDiff,$labelchangearray,$transpos, $DNAobj,
	$RNAobj,$atg_label,$coord_system)=@_;
    my @labelchange=@{$labelchangearray};
    my %before_mutation;
    my ($seq, $label, $len ) = @labelchange;
    my ($praelabel,$postlabel,$lastmutlabel); #lastmutlabel used for aa_allele_ori

### Interpration of variables and conversion of the mutation array, creating
### new useful label

    $praelabel=$DNAobj->label(-1,$label); # 1 before label
    if (defined($len)) {
	if ($len == 0) {
	    $postlabel=$label;
	    $lastmutlabel=$label;
	} elsif ($len == 1) {
	    $lastmutlabel=$label; # last nucleotide affected
	    $postlabel=$DNAobj->label(2,$lastmutlabel); # $len after label
	} else {
	    $lastmutlabel=$DNAobj->label($len,$label); # last nucleotide affected
	    $postlabel=$DNAobj->label(2,$lastmutlabel); # $len after label
	}
    } else {
	$len = length($seq);
	$lastmutlabel=$label;
	$postlabel=$DNAobj->label(2,$label); # 1 after label
    }
    #print "\tDEBUG mutation n. $mut_number label $label lastmutlabel $lastmutlabel postlabel $postlabel\n";
    #unless(defined($len)) { 
    #$len = length($seq);
    #}
    #print  `date`, " initializations done\n";

### DNAmutation and alleles object creation

    my $entry_offset=$seqDiff->offset+1;


    my $dnamut_start=$label-$entry_offset+1;
    # if negative DNA positions (before ATG)
    if ($dnamut_start <= 0) { $dnamut_start--; }

    my $dnamut_end;
    $len == 0 ? $dnamut_end = $dnamut_start : $dnamut_end = $dnamut_start+$len -1;

    my $dnamut = Bio::Variation::DNAMutation->new(-start => $dnamut_start, 
						  -end => $dnamut_end,

						  #-length => $len,
						  #-allele_mut => $seq 
						  );
    $dnamut->isMutation(1);
    my $da_m = Bio::Variation::Allele->new;
    $da_m->seq($seq) if $seq;
    $dnamut->allele_mut($da_m);
    $dnamut->add_Allele($da_m);
    # allele_ori
    my $allele_ori=$DNAobj->labelsubseq($praelabel,undef,$postlabel); # get seq
    chop $allele_ori; # away the postlabel nucleotide
    $allele_ori=substr($allele_ori,1); # away the praelabel nucleotide
    my $da_o = Bio::Variation::Allele->new;
    $da_o->seq($allele_ori) if $allele_ori;
    $dnamut->allele_ori($da_o);
    $len == 0 ?  $dnamut->length($len) :  $dnamut->length(CORE::length $allele_ori); 

    #print " --------------- $dnamut_start -$len-  $dnamut_end -\n";
    $seqDiff->add_Variant($dnamut);
    $dnamut->mut_number($mut_number);
    $before_mutation{dnamut}=$dnamut; # passing the relevant information to "after"
    
    # setting proof
    if ($coord_system eq "entry") {
	$dnamut->proof('experimental');
    } else {
	$dnamut->proof('computed');
    }

    # how many nucleotides to store upstream and downstream of the change
    my $flanklen = 25;

    #print  `date`, " flanking sequnces start\n";
    my $upLabel = $DNAobj->label(1-$flanklen,$praelabel); # this could be unavailable!
    my $upstreamseq;
    if ($upLabel > 0) {
	$upstreamseq = $DNAobj->labelsubseq($upLabel, undef, $praelabel);
    } else { # from start (less than $flanklen nucleotides)
	$upstreamseq = $DNAobj->labelsubseq($DNAobj->start, undef, $praelabel);
    }
    $dnamut->upStreamSeq($upstreamseq);

    my $dnstreamseq = $DNAobj->labelsubseq($postlabel, $flanklen);
    $dnamut->dnStreamSeq($dnstreamseq); # $flanklen or less nucleotides


    #  print  `date`, " _is_rna_affected\n";

    #
    # Exon/Intron structure


### Check if mutation propagates to RNA (and AA) level
### Workaround has to be written for problems with inserted labels

    my @exons=$RNAobj->all_Exons;
    my $RNAstart=$RNAobj->start;
    my $RNAend=$RNAobj->end;
    my ($firstexon,$before,$after,$i);
    my ($rnaAffected) = 0;

    # check for inserted labels (that require follows instead of <,>)
    my $DNAend=$RNAobj->{seq}->end;
    if ($praelabel > $DNAend or $postlabel > $DNAend) { 
	#this means one of the two labels is an inserted one 
	#(coming from a previous mutation. This would falsify all <,> 
	#checks, so the follow() has to be used
	croak "Workaround not implemented yet!! prae $praelabel post $postlabel dnaend $DNAend\n\n";
    } else {
	my $strand=$exons[0]->strand;
	if (($strand == 1 and $postlabel <= $RNAstart) or 
	    ($strand != 1 and $postlabel >= $RNAstart)) {
	    carp "RNA not affected because change occurs before RNAstart";
	}
	elsif (($strand == 1 and $praelabel >= $RNAend) or 
	       ($strand != 1 and $praelabel <= $RNAend)) {
	    carp "RNA not affected because change occurs after RNAend";
	}
	elsif (scalar @exons == 1) {
	    #if @exons now empty -> no introns, just one exon
	    $rnaAffected = 1; # then RNA is affected!
	} else {	    
	    # otherwise check for change occurring inside an intron
	    $firstexon=shift(@exons);
	    $before=$firstexon->end;
	    
	    foreach $i (0..$#exons) {
		$after=$exons[$i]->start;
		if ( ($strand == 1 and ($praelabel < $before || $praelabel >= $after ||
		      $postlabel > $after)) or 
		     ($strand != 1 and ($praelabel > $before || $praelabel <= $after ||
		      $postlabel < $after)) ){ # check this if really correct
		    $rnaAffected = 1; 
		    # $i is number of exon and can be used for proximity check
		}
		$before=$exons[$i]->end;
	    }
	    
	}
    }
	   #carp "RNA not affected because change occurs inside an intron";
	    #return(0); # if still not returned, then not affected, return 0

    #
    # RNA         checking if start of mutation is inside 
    #             the Transcript boundaries (i.e. in coding region)
    #
    #if (_is_rna_affected($RNAobj,$praelabel,$postlabel)) {
    if ($rnaAffected) {

### Creation of RNA variation objects
	
	my $rnapos_end;
	$len == 0 ?  $rnapos_end = $transpos : $rnapos_end = $transpos + $len -1;
	my $rnachange = Bio::Variation::RNAChange->new(-start => $transpos, 
						       -end =>  $rnapos_end
						       );
	$rnachange->isMutation(1);

	# setting proof
	if ($coord_system eq "coding") {
	    $rnachange->proof('experimental');
	} else {
	    $rnachange->proof('computed');
	}

	$seqDiff->add_Variant($rnachange);
	$rnachange->DNAMutation($dnamut);
	$dnamut->RNAChange($rnachange);
	$rnachange->mut_number($mut_number);
	$before_mutation{rnachange}=$rnachange;

#    print  `date`, " before RNAObj query\n";
	# setting the codon_position of the "start" nucleotide of the change
	$rnachange->codon_pos(($RNAobj->frame($label))+1); # codon_pos=frame+1

	my @exons=$RNAobj->all_Exons;
	$before_mutation{exons}=\@exons;
	#print  `date`, " before flank, after exons. RNAObj query\n";
	# if cannot retrieve from Transcript, Transcript::upstream_seq will be used
	# before "fac7 g 65" bug discovered
	# $upLabel=$RNAobj->label(1-$flanklen,$praelabel);
	my $RNApraelabel=$RNAobj->label(-1,$label); # to fix fac7g65 bug
	# for the fix, all praelabel used in the next block have been changed to RNApraelabel
	$upLabel=$RNAobj->label(1-$flanklen,$RNApraelabel);
	if ($RNAobj->valid($upLabel)) {
	    $upstreamseq=$RNAobj->labelsubseq($upLabel, undef, $RNApraelabel);
	} else {
	    $upstreamseq=$RNAobj->labelsubseq($RNAobj->start, undef, $RNApraelabel);
	    my $lacking=$flanklen-length($upstreamseq); # how many missing
	    my $upstream_atg=$exons[0]->subseq(-$lacking,-1);
	    $upstreamseq=$upstream_atg . $upstreamseq;
	}

	#
	# Exon/Intron 
	#
	#first only for point mutations, elaborate later
	if ($transpos == $rnapos_end) {
	    my ($affected_exon,$exon_number)=$RNAobj->in_which_Exon($RNApraelabel); # changed from praelabel to RNApraelabel (does it work like this??)
	    $rnachange->exons_modified($exon_number); 
	    $dnamut->region('exon');
	    $dnamut->region_value($exon_number);
	}
	#
	# Flanking sequences
	#
	$rnachange->upStreamSeq($upstreamseq);
	#print  `date`, " before RNAObj dnstream\n";
	# won't work OK if postlabel NOT in Transcript
	# now added RNApostlabel but this has to be fully tested
	# for the fix, all postlabel used in the next block have been changed to RNApostlabel
	# removed for now.... (bad fix)
	#my $RNApostlabel=$RNAobj->label(2,$label); # to fix fac7g65 bug
	$dnstreamseq=$RNAobj->labelsubseq($postlabel, $flanklen);
	if ($dnstreamseq == -1) { # if out of transcript was requested
	    my $lastexon=$exons[-1];
	    my $lastexonlength=$lastexon->length;
	    $dnstreamseq=$RNAobj->labelsubseq($postlabel); # retrieves till RNAend
	    my $lacking=$flanklen-length($dnstreamseq); # how many missing
	    my $downstream_stop=$lastexon->subseq($lastexonlength+1,undef,$lacking);
	    $dnstreamseq .= $downstream_stop;
	} else {
	    $rnachange->dnStreamSeq($dnstreamseq);
	}
	#print  `date`, " before AAobject creation\n";
	# AAChange creation
	my $AAobj=$RNAobj->get_Translation;
	# storage of praelabel here, to be used in create_mut_objs_after
	my $aachange = Bio::Variation::AAChange->new(-start => $RNApraelabel
						     );
	$aachange->isMutation(1);
	$aachange->proof('computed');

	$seqDiff->add_Variant($aachange);
	$rnachange->AAChange($aachange);
	$aachange->RNAChange($rnachange);

	$aachange->mut_number($mut_number);
	$before_mutation{aachange}=$aachange;

	my $ra_o = Bio::Variation::Allele->new;
	$ra_o->seq($allele_ori) if $allele_ori;
	$rnachange->allele_ori($ra_o);
	
	$rnachange->length(CORE::length $allele_ori);

	my $ra_m = Bio::Variation::Allele->new;
	$ra_m->seq($seq) if $seq;
	$rnachange->allele_mut($ra_m);
	$rnachange->add_Allele($ra_m);
	
	#$rnachange->allele_mut($seq);
	$rnachange->end($rnachange->start) if $rnachange->length == 0;

	# this holds the aminoacid sequence that will be affected by the mutation
	my $aa_allele_ori=$AAobj->labelsubseq($label,undef,$lastmutlabel);

	my $aa_o = Bio::Variation::Allele->new;
	$aa_o->seq($aa_allele_ori) if $aa_allele_ori;
	$aachange->allele_ori($aa_o);
	#$aachange->allele_ori($aa_allele_ori);
	
	my $aa_length_ori = length($aa_allele_ori);
	$aachange->length($aa_length_ori); #print "==========$aa_length_ori\n";
	$aachange->end($aachange->start + $aa_length_ori - 1 );
    }
    #
    # Create RNAChange object for 5' and 3' UTR mutations
    #
    elsif ($seqDiff->offset != 0 ) {
	#print  `date`, " creating RNAChange\n";
	my $rnachange = Bio::Variation::RNAChange->new; 

	my $ra_o = Bio::Variation::Allele->new;
	$ra_o->seq($allele_ori) if $allele_ori;
	$rnachange->allele_ori($ra_o);
	#$rnachange->allele_ori($dnamut->allele_ori);

	my $ra_m = Bio::Variation::Allele->new;
	$ra_m->seq($dnamut->allele_mut->seq) if $dnamut->allele_mut->seq;
	$rnachange->allele_mut($ra_m);
	$rnachange->add_Allele($ra_m);
	#$rnachange->allele_mut($dnamut->allele_mut);

	$rnachange->upStreamSeq($dnamut->upStreamSeq);
	$rnachange->dnStreamSeq($dnamut->dnStreamSeq);
	$rnachange->start($dnamut->start);
	$rnachange->end($dnamut->end);
	$rnachange->length($dnamut->length);
	$rnachange->mut_number($dnamut->mut_number);
	#print  `date`, " creating RNAChange\n";
	# setting proof
	if ($coord_system eq "coding") {
	    $rnachange->proof('experimental');
	} else {
	    $rnachange->proof('computed');
	}
	# setting region
	if ($rnachange->end < 0) {
	    $rnachange->region('5\'UTR');
	} else {
	    $rnachange->region('3\'UTR');
	}
	#print  `date`, " RNAChange done\n";
	$seqDiff->add_Variant($rnachange);
	$rnachange->DNAMutation($dnamut);
	$dnamut->RNAChange($rnachange);
	$before_mutation{rnachange}=$rnachange;
    }
    else {
	carp "Mutation starts outside coding region, RNAChange object not created";
    }
    #  print "labelchange: @labelchange\ntranspos: $transpos\n"; #debug
    #print  `date`, " end of before mutation\n";

  return (\%before_mutation);
}



# args: reference to label changearray, reference to position changearray
# Function: take care of the creation of mutation objects, with
# information AFTER the change takes place
sub create_mut_objs_after {
    my ($mut_number,$seqDiff,$labelchangearray,$transpos, 
	$before_mutation,$DNAobj,$RNAobj,$atg_label,$coord_system)=@_;

    # the $before_mutation can be a reference to a hash containing 
    #any information coming from the create_mut_objs_before function

    my @labelchange=@{$labelchangearray};
    my $dnamut=$before_mutation->{dnamut};
    my $rnachange=$before_mutation->{rnachange};
    
    if ($rnachange->region eq 'coding') {
	my $aachange=$before_mutation->{aachange};
	my ($AAobj,$aa_start_praelabel,$aa_start,$mut_translation);
	$AAobj=$RNAobj->get_Translation;
	$aa_start_praelabel=$aachange->start;
	$aa_start=$AAobj->position($RNAobj->label(2,$aa_start_praelabel));
	$aachange->start($aa_start);
	$mut_translation=$AAobj->seq;

	# this now takes in account possible praeinsertions
	my $aa_m = Bio::Variation::Allele->new;
	$aa_m->seq(substr($mut_translation,$aa_start-1)) if substr($mut_translation,$aa_start-1);
	$aachange->allele_mut($aa_m);
	$aachange->add_Allele($aa_m);
	#$aachange->allele_mut(substr($mut_translation,$aa_start-1)); 
	#$aachange->allele_mut($mut_translation); 
	my ($rlenori, $rlenmut);
        $rlenori = CORE::length($aachange->RNAChange->allele_ori->seq);
	$rlenmut = CORE::length($aachange->RNAChange->allele_mut->seq);
	#point mutation
	if ($rlenori == 1 and $rlenmut == 1 and $aachange->allele_ori->seq ne '*') {
	    my $alleleseq;
	    if ($aachange->allele_mut->seq) {
		$alleleseq = substr($aachange->allele_mut->seq, 0, 1);
		$aachange->allele_mut->seq($alleleseq);
	    }
	    $aachange->end($aachange->start);
	    $aachange->length(1);
	}
	#inframe mutatation
	elsif (($rlenori+$rlenmut)%3 == 0 ) {

	    if ($aachange->RNAChange->codon_pos == 1){
		if ($aachange->RNAChange->allele_mut->seq eq  '') {
		    $aachange->allele_mut->seq('');
		}
		elsif ($aachange->RNAChange->allele_ori->seq eq '' ) {
		    $aachange->allele_mut->seq(substr $aachange->allele_mut->seq, 0, 
					  length ($aachange->RNAChange->allele_mut->seq) / 3);
		    $aachange->allele_ori->seq('');
		    $aachange->length(0);
		}
	    } else {
		#elsif  ($aachange->RNAChange->codon_pos == 2){
		if (not $aachange->RNAChange->allele_mut->seq ) {
		    $aachange->allele_mut-seq(substr $aachange->allele_mut->seq, 0, 1);
		}
		elsif (not $aachange->RNAChange->allele_ori->seq) {
		    $aachange->allele_mut->seq(substr $aachange->allele_mut->seq, 0, 
					  length ($aachange->RNAChange->allele_mut->seq) / 3 +1);
		}
	    }
	} else {
	    #frameshift
	    #my $pos = index $aachange->allele_mut
	    #$aachange->allele_mut(substr($aachange->allele_mut, 0, 1));
	    $aachange->length(CORE::length($aachange->allele_ori->seq));
	    my $aaend = $aachange->start + $aachange->length -1;
	    $aachange->end($aachange->start);
	    
	}

	# splicing site deletion check
	my @beforeexons=@{$before_mutation->{exons}};
	my @afterexons=$RNAobj->all_Exons;
	my $i;
	if (scalar(@beforeexons) ne scalar(@afterexons)) {
	    carp "Exons have been modified at mutation n.$mut_number!";
	    $rnachange->exons_modified(1);
	} else {
	  EXONCHECK:
	    foreach $i (0..$#beforeexons) {
		if ($beforeexons[$i] ne $afterexons[$i]) {
		    carp "Exons have been modified at mutation n.$mut_number!";
		    $rnachange->exons_modified(1);
		    last EXONCHECK;
		}
	    }
	}
    }
    return 1;
}

# this function can be modified to detect proximity of mutation to splicesite
# this function cheats: it uses labels
sub _is_rna_affected {
    my ($RNAobj,$praelabel,$postlabel)=@_;

    my @exons=$RNAobj->all_Exons;
    my $RNAstart=$RNAobj->start;
    my $RNAend=$RNAobj->end;
    my ($firstexon,$before,$after,$i);

    # check for inserted labels (that require follows instead of <,>)
    my $DNAend=$RNAobj->{seq}->end;
    if (($praelabel > $DNAend)||($postlabel > $DNAend)) { 
	#this means one of the two labels is an inserted one 
	#(coming from a previous mutation. This would falsify all <,> 
	#checks, so the follow() has to be used
	croak "Workaround not implemented yet!! prae $praelabel post $postlabel dnaend $DNAend\n\n";
    } else {
	my $strand=$exons[0]->strand;
	if ($strand == 1) {
	    if ($postlabel <= $RNAstart) {
		carp "RNA not affected because change occurs before RNAstart";
		return(0);
	    }
	    if ($praelabel >= $RNAend) {
		carp "RNA not affected because change occurs after RNAend";
		return(0);
	    }
	    $firstexon=shift(@exons);
	    $before=$firstexon->end;

	    unless (@exons) { 
		#if @exons now empty -> no introns, just one exon
		return(1); # then RNA is affected!
	    }
	    
	    # otherwise check for change occurring inside an intron
	    foreach $i (0..$#exons) {
		$after=$exons[$i]->start;
		if (($praelabel < $before)||($praelabel >= $after)||
		    ($postlabel > $after)) { # check this if really correct
		    return(1); 
                    # $i is number of exon and can be used for proximity check
		}
		$before=$exons[$i]->end;
	    }
	    carp "RNA not affected because change occurs inside an intron";
	    return(0); # if still not returned, then not affected, return 0
	} else { # reverse strand, all (<,>) swapped
	    if ($postlabel >= $RNAstart) {
		carp "RNA not affected because change occurs before RNAstart, revstrand";
		return(0);
	    }
	    if ($praelabel <= $RNAend) {
		carp "RNA not affected because change occurs after RNAend, revstrand";
		return(0);
	    }
	    $firstexon=shift(@exons);
	    $before=$firstexon->end;
	    
	    unless (@exons) { # if @exons now empty -> no introns, just one exon
		return(1); # then RNA is affected!
	    }
	    
	    # otherwise check for change occurring inside an intron
	    foreach $i (1..$#exons) {
		$after=$exons[$i]->start;
		if (($praelabel > $before)||($praelabel <= $after)||
		    ($postlabel < $after)) { # check this if really correct
		    return(1); # $i is number of exon and can be used for proximity check
		}
		$before=$exons[$i]->end;
	    }
	    carp "RNA not affected because change occurs inside an intron, revstrand";
	    return(0); # if still not returned, then not affected, return 0
	}
    }
}

sub _is_rna_affected2 {
    my ($RNAobj,$praelabel,$postlabel)=@_;

    my @exons=$RNAobj->all_Exons;
    my $RNAstart=$RNAobj->start;
    my $RNAend=$RNAobj->end;
    my ($firstexon,$before,$after,$i);

    # check for inserted labels (that require follows instead of <,>)
    my $DNAend=$RNAobj->{seq}->end;
    if (($praelabel > $DNAend)||($postlabel > $DNAend)) { 
	#this means one of the two labels is an inserted one 
	#(coming from a previous mutation. This would falsify all <,> 
	#checks, so the follow() has to be used
	croak "Workaround not implemented yet!! prae $praelabel post $postlabel dnaend $DNAend\n\n";
    } else {
	my $strand=$exons[0]->strand;
	if ($strand == 1) {
	    if (($strand == 1 and $postlabel <= $RNAstart) or $postlabel >= $RNAstart) {
		carp "RNA not affected because change occurs before RNAstart";
		return(0);
	    }
	    if (($strand == 1 and $praelabel >= $RNAend) or $praelabel <= $RNAend) {
		carp "RNA not affected because change occurs after RNAend";
		return(0);
	    }
	    $firstexon=shift(@exons);
	    $before=$firstexon->end;

	    unless (@exons) { 
		#if @exons now empty -> no introns, just one exon
		return(1); # then RNA is affected!
	    }
	    
	    # otherwise check for change occurring inside an intron
	    if ($strand == 1) {
		foreach $i (0..$#exons) {
		    $after=$exons[$i]->start;
		    if (($praelabel < $before)||($praelabel >= $after)||
			($postlabel > $after)) { # check this if really correct
			return(1); 
			# $i is number of exon and can be used for proximity check
		    }
		    $before=$exons[$i]->end;
		}
	    } else {
		foreach $i (1..$#exons) {
		    $after=$exons[$i]->start;
		    if (($praelabel > $before)||($praelabel <= $after)||
			($postlabel < $after)) { # check this if really correct
			return(1); # $i is number of exon and can be used for proximity check
		    }
		    $before=$exons[$i]->end;
		}
	    }
	    carp "RNA not affected because change occurs inside an intron";
	    return(0); # if still not returned, then not affected, return 0
	}
    }
}

1
