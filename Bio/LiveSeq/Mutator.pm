# $Id$
#
# bioperl module for Bio::LiveSeq::Mutator
#
# Cared for by Joseph Insana <insana@ebi.ac.uk> <jinsana@gmx.net>
#
# Copyright Joseph Insana
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

  Bio::LiveSeq::Mutator - Package mutating LiveSequences

=head1 SYNOPSIS

  # $gene is a Bio::LiveSeq::Gene object
  my $mutate = Bio::LiveSeq::Mutator->new('-gene' => $gene,
  					  '-numbering' => "coding"
  					   );
  # $mut is a Bio::LiveSeq::Mutation object
  $mutate->add_Mutation($mut);
  # $results is a Bio::Variation::SeqDiff object
  my $results=$mutate->change_gene();
  if ($results) {
      my $out = Bio::Variation::IO->new( '-format' => 'flat');
      $out->write($results);
  }

=head1 DESCRIPTION

  This class mutates L<Bio::LiveSeq::Gene> objects and returns
  L<Bio::Variation::SeqDiff> object. Mutations are described as
  L<Bio::LiveSeq::Mutation> objects.

=head1 FEEDBACK


User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                         - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Heikki Lehvaslaiho & Joseph A.L. Insana

  Email:  heikki@ebi.ac.uk
          insana@ebi.ac.uk, jinsana@gmx.net

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
use vars qw(@ISA);
use strict;

use Carp qw(cluck croak carp);
use vars qw($VERSION @ISA);
use Bio::Variation::SeqDiff;
use Bio::Variation::DNAMutation;
use Bio::Variation::RNAChange;
use Bio::Variation::AAChange;
use Bio::Variation::Allele;

#use integer;
# Object preamble - inheritance

use Bio::Root::RootI;

@ISA = qw( Bio::Root::RootI );

sub new {
    my($class,@args) = @_;
    my $self;
    $self = {};
    bless $self, $class;

    my ($gene, $numbering) =
	    $self->_rearrange([qw(GENE
				  NUMBERING
				  )],
			      @args);

    $self->{ 'mutations' } = [];

    $gene && $self->gene($gene);
    $numbering && $self->numbering($numbering);

    #class constant;
    $self->{'flanklen'} = 25;
    return $self; # success - we hope!
}

=head2 gene

 Title   : gene
 Usage   : $mutobj = $obj->gene;
         : $mutobj = $obj->gene($objref);
 Function: 

           Returns or sets the link-reference to a
           L<Bio::LiveSeq::Gene> object. If no value has ben set, it
           will return undef

 Returns : an object reference  or undef
 Args    : a L<Bio::LiveSeq::Gene> 

=cut

sub gene {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::LiveSeq::Gene') ) {
	  $self->throw("Is not a Bio::LiveSeq::Gene object but a [$value]");
	  return undef;
      }
      else {
	  $self->{'gene'} = $value;
      }
  }
  unless (exists $self->{'gene'}) {
      return (undef);
  } else {
      return $self->{'gene'};
  }
}


=head2 numbering

 Title   : numbering
 Usage   : $obj->numbering();
 Function:

            Sets and returns coordinate system used in positioning the
            mutations.

 Example :
 Returns : string
 Args    : string (coding [transcript number] | entry)

=cut


sub numbering {
    my ($self,$value) = @_;
    if( defined $value) {
	if ($value =~ /(coding)( )?(\d+)?/ or $value eq 'entry') {
	    $self->{'numbering'} = $value;
	} else { # defaulting to 'coding'
	    $self->{'numbering'} = 'coding';
	}
    }
    unless (exists $self->{'numbering'}) {
	return 'coding';
    } else {
	return $self->{'numbering'};
    }
}

=head2 add_Mutation

 Title   : add_Mutation
 Usage   : $self->add_Mutation($ref)
 Function: adds a L<Bio::Liveseq::Mutation> object
 Example :
 Returns :
 Args    : a L<Bio::Liveseq::Mutation>

=cut

sub add_Mutation{
    my ($self,$value) = @_;
    if( $value->isa('Bio::Liveseq::Mutation') ) {
	my $com = ref $value;
	$self->throw("Is not a Mutation object but a [$com]" );
	return undef;
    }
    if (! $value->pos) {
	carp "No value for mutation position in the sequence!";
	return undef;
    }
    if (! $value->seq && ! $value->len) {
	carp "Either mutated sequence or length of the deletion must be given!";
	return undef;
    }
    push(@{$self->{'mutations'}},$value);
}

=head2 each_Mutation

 Title   : each_Mutation
 Usage   : foreach $ref ( $a->each_Mutation )
 Function: gets an array of L<Bio::Liveseq::Mutation> objects
 Example :
 Returns : array of Mutations
 Args    :


=cut

sub each_Mutation{
   my ($self) = @_;
   return @{$self->{'mutations'}};
}


=head2 mutation

 Title   : mutation
 Usage   : $mutobj = $obj->mutation;
         : $mutobj = $obj->mutation($objref);
 Function: 

           Returns or sets the link-reference to the current mutation
           object.  If the value is not set, it will return undef.
           Internal method.

 Returns : an object reference  or undef

=cut


sub mutation {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::LiveSeq::Mutation') ) {
	  $self->throw("Is not a Bio::LiveSeq::Mutation object but a [$value]");
	  return undef;
      }
      else {
	  $self->{'mutation'} = $value;
      }
  }
  unless (exists $self->{'mutation'}) {
      return (undef);
  } else {
      return $self->{'mutation'};
  }
}

=head2 DNA

 Title   : DNA
 Usage   : $mutobj = $obj->DNA;
         : $mutobj = $obj->DNA($objref);
 Function:

           Returns or sets the reference to the LiveSeq object holding
           the reference sequence. If there is no link, it will return
           undef.
           Internal method.

 Returns : an object reference or undef

=cut

sub DNA {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::LiveSeq::DNA') and ! $value->isa('Bio::LiveSeq::Transcript') ) {
	  $self->throw("Is not a Bio::LiveSeq::DNA/Transcript object but a [$value]");
	  return undef;
      }
      else {
	  $self->{'DNA'} = $value;
      }
  }
  unless (exists $self->{'DNA'}) {
      return (undef);
  } else {
      return $self->{'DNA'};
  }
}


=head2 RNA

 Title   : RNA
 Usage   : $mutobj = $obj->RNA;
         : $mutobj = $obj->RNA($objref);
 Function: 

           Returns or sets the reference to the LiveSeq object holding
           the reference sequence. If the value is not set, it will return
           undef.
           Internal method.

 Returns : an object reference  or undef

=cut


sub RNA {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::LiveSeq::Transcript') ) {
	  $self->throw("Is not a Bio::LiveSeq::RNA/Transcript object but a [$value]");
	  return undef;
      }
      else {
	  $self->{'RNA'} = $value;
      }
  }
  unless (exists $self->{'RNA'}) {
      return (undef);
  } else {
      return $self->{'RNA'};
  }
}


=head2 dnamut

 Title   : dnamut
 Usage   : $mutobj = $obj->dnamut;
         : $mutobj = $obj->dnamut($objref);
 Function: 

           Returns or sets the reference to the current DNAMutation object.
           If the value is not set, it will return undef.
           Internal method.

 Returns : an L<Bio::Variation::DNAMutation> or undef

=cut


sub dnamut {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::Variation::DNAMutation') ) {
	  $self->throw("Is not a Bio::Variation::DNAMutation object but a [$value]");
	  return undef;
      }
      else {
	  $self->{'dnamut'} = $value;
      }
  }
  unless (exists $self->{'dnamut'}) {
      return (undef);
  } else {
      return $self->{'dnamut'};
  }
}


=head2 rnachange

 Title   : rnachange
 Usage   : $mutobj = $obj->rnachange;
         : $mutobj = $obj->rnachange($objref);
 Function: 

           Returns or sets the reference to the current RNAChange object.
           If the value is not set, it will return undef.
           Internal method.

 Returns : an L<Bio::Variation::RNAChange> or undef

=cut


sub rnachange {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::Variation::RNAChange') ) {
	  $self->throw("Is not a Bio::Variation::RNAChange object but a [$value]");
	  return undef;
      }
      else {
	  $self->{'rnachange'} = $value;
      }
  }
  unless (exists $self->{'rnachange'}) {
      return (undef);
  } else {
      return $self->{'rnachange'};
  }
}


=head2 aachange

 Title   : aachange
 Usage   : $mutobj = $obj->aachange;
         : $mutobj = $obj->aachange($objref);
 Function: 

           Returns or sets the reference to the current AAChange object.
           If the value is not set, it will return undef.
           Internal method.

 Returns : an L<Bio::Variation::AAChange> or undef

=cut


sub aachange {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::Variation::AAChange') ) {
	  $self->throw("Is not a Bio::Variation::AAChange object but a [$value]");
	  return undef;
      }
      else {
	  $self->{'aachange'} = $value;
      }
  }
  unless (exists $self->{'aachange'}) {
      return (undef);
  } else {
      return $self->{'aachange'};
  }
}


=head2 exons

 Title   : exons
 Usage   : $mutobj = $obj->exons;
         : $mutobj = $obj->exons($objref);
 Function: 

           Returns or sets the reference to a current array of Exons.
           If the value is not set, it will return undef.
           Internal method.

 Returns : an array of L<Bio::LiveSeq::Exon> objects or undef

=cut


sub exons {
  my ($self,$value) = @_;
  if (defined $value) {
      $self->{'exons'} = $value;
  }
  unless (exists $self->{'exons'}) {
      return (undef);
  } else {
      return $self->{'exons'};
  }
}



=head2 change_gene

 Title   : change_gene
 Usage   : my $mutate = Bio::LiveSeq::Mutator->new(-gene => $gene,
						   numbering => "coding"
						   );
           $mutate->add_Mutation($mut);
           my $results=$mutate->change_gene();

 Function:

           Returns a Bio::Variation::SeqDiff object containing the
           results of the changes performed according to the
           instructions present in the mutmatrix.  The -numbering
           argument decides what molecule is being changed and what
           numbering scheme being used:

            -numbering => "entry"

               determines the DNA level, using the numbering of the
               beginning of the sequence

            -numbering => "coding"

               determines the cDNA level, using the numbering from the
               beginning of the 1st transcript

           Alternative transcripts can be used by specifying "coding 2" or
           "coding 3" ...

 Args    : Bio::LiveSeq::Gene object
           reference to a mutmatrix array of arrays
           string specifying a numbering scheme (defaults to 'coding')
 Returns : Bio::Variation::SeqDiff object or 0 on error

=cut

sub change_gene {
    my ($self) = @_;

    #
    # Sanity check
    #
    unless ($self->gene) {
	carp "Input object Bio::LiveSeq::Gene is not given\n";
	return 0;
    }
    #
    # Setting the reference sequence based on -numbering
    #
    my @transcripts=@{$self->gene->get_Transcripts};
    my $refseq; # will hold Bio::LiveSeq:Transcript object or Bio::LiveSeq::DNA
    if ($self->numbering =~ /(coding)( )?(\d+)?/ ) {
	$self->numbering($1);
	my $transnumber = $3;
	$transnumber-- if $3; # 1 -> 0, 2 -> 1
	if ($transnumber && $transnumber >= 0 && $transnumber <= $#transcripts) {
	    $refseq=$transcripts[$transnumber];
	} else {
	    $transnumber && carp "The alternative transcript number", $transnumber+1,
	    "- does not exist. Reverting to the 1st transcript\n";
	    $refseq=$transcripts[0];
	}
    } else {
	$refseq=$transcripts[0]->{'seq'};
    }
    #
    # Converting mutation positions to labels
    #
    carp "no mutations", return 0 unless  $self->_mutationpos2label($refseq);
    #
    # Recording the state: SeqDiff object creation  ?? transcript no.??
    #
    my $seqDiff = Bio::Variation::SeqDiff->new();
    $seqDiff->moltype($self->gene->get_DNA->moltype);
    $seqDiff->numbering($self->numbering);
    my ($DNAobj, $RNAobj);
    if ($refseq->isa("Bio::LiveSeq::Transcript")) {
	$self->RNA($refseq);
	$self->DNA($refseq->{'seq'});
	$seqDiff->rna_ori($refseq->seq );
	$seqDiff->aa_ori($refseq->get_Translation->seq);
    } else {
	$self->DNA($refseq);
	$self->RNA($transcripts[0]);
	$seqDiff->rna_ori($self->RNA->seq);
	$seqDiff->aa_ori($self->RNA->get_Translation->seq);
    }
    $seqDiff->dna_ori($self->DNA->seq);
    # put the accession number into the SeqDiff object ID
    $seqDiff->id($self->DNA->accession_number);

    # the atg_offset takes in account that DNA object could be a subset of the
    # whole entry (via the light_weight loader)
    my $atg_label=$self->RNA->start;
    my $atg_offset=$self->DNA->position($atg_label)+($self->DNA->start)-1;
    $seqDiff->offset($atg_offset - 1);
    my $cds_end = $self->gene->downbound - $atg_offset + 1;
    $seqDiff->cds_end($cds_end);
    $self->DNA->coordinate_start($atg_label);

    # need to add more than one rna & aa
    #foreach $transcript (@transcripts) {
    #  $seqDiff{"ori_transcript_${i}_seq"}=$transcript->seq;
    #  $seqDiff{"ori_translation_${i}_seq"}=$transcript->get_Translation->seq;
    #}

    # issue changes
    my $k;
    foreach my $mutation ($self->each_Mutation) {
#	print STDERR "+2:", $mutation->label, "+\n";
	next unless $mutation->label > 0;
	$self->mutation($mutation);

	$mutation->issue(++$k);
	#
	# current position on the trascript
	#
	if ($self->numbering =~ /coding/) {
	    $mutation->transpos($mutation->pos); # transpos given by user
	} else {
	    #transpos of label / It will be 0 if mutating an intron, negative if upstream of ATG
	    $mutation->transpos($self->RNA->position($mutation->label,$atg_label));
	}
	#
	# Calculate adjacent labels based on the position on the current sequence
	#
	$mutation->prelabel($self->DNA->label(-1, $mutation->label)); # 1 before label
	if ($mutation->len == 0) {
	    $mutation->postlabel($mutation->label);
	    $mutation->lastlabel($mutation->label);
	} elsif ($mutation->len == 1) {
	    $mutation->lastlabel($mutation->label); # last nucleotide affected
	    $mutation->postlabel($self->DNA->label(2,$mutation->lastlabel)); # $len after label
	} else {
	    $mutation->lastlabel($self->DNA->label($mutation->len,$mutation->label));
	    $mutation->postlabel($self->DNA->label(2,$mutation->lastlabel));
	}
	my $dnamut = $self->_set_DNAMutation($seqDiff);
	#
	#
	#
	if ($self->_rnaAffected) {
	    $self->_set_effects($seqDiff, $dnamut);
	}
	elsif ($seqDiff->offset != 0 ) {
	    $self->_untranslated ($seqDiff, $dnamut);
	}
	else {
	    carp "Mutation starts outside coding region, RNAChange object not created";
	}

	#########################################################################
	# Mutations are done here!                                              #
	$refseq->labelchange($mutation->seq, $mutation->label, $mutation->len); #
	#########################################################################
#	create_mut_objs_after ( $seqDiff,$before_mutation);
	$self->_post_mutation ($seqDiff);

	$self->dnamut(undef);
	$self->rnachange(undef);
	$self->aachange(undef);
	$self->exons(undef);
    }
    # record the final state of all three sequences
    $seqDiff->dna_mut($self->DNA->seq);
    $seqDiff->rna_mut($self->RNA->seq);
    $seqDiff->aa_mut($refseq->get_Translation->seq);

    #$seqDiff{mut_dna_seq}=$gene->get_DNA->seq;
    #my $i=1;
    #foreach $transcript (@transcripts) {
    #  $seqDiff{"mut_transcript_${i}_seq"}=$transcript->seq;
    #  $seqDiff{"mut_translation_${i}_seq"}=$transcript->get_Translation->seq;
    #}
    return $seqDiff;
}

=head2 _mutationpos2label

 Title   : _mutationpos2label
 Usage   :
 Function:
 Example :
 Returns : number of valid mutations
 Args    : LiveSeq sequence object

=cut

sub _mutationpos2label {
    my ($self, $refseq) = @_;
    my $count;
    my @bb = @{$self->{'mutations'}};
    my $cc = scalar @bb;
    #print STDERR "x,------$cc-----\n";
    foreach my $mut (@{$self->{'mutations'}}) {
	my $label = $refseq->label($mut->pos); # get the label for the position
	$mut->label($label), $count++ if $label > 0 ;
	#print STDERR "x", $mut->pos,'|' ,$mut->label, "\n";
    }
    return $count;
}

#
# Calculate labels around mutated nucleotide
#

=head2 _set_DNAMutation

 Title   : _set_DNAMutation
 Usage   :
 Function:

           Stores DNA level mutation attributes before mutation into
           L<Bio::Variation::DNAMutation> object.  Links it to SeqDiff
           object.

 Example :
 Returns : L<Bio::Variation::DNAMutation> object
 Args    : L<Bio::Variation::SeqDiff> object

=cut

sub _set_DNAMutation {
    my ($self, $seqDiff) = @_;

    my $dnamut_start = $self->mutation->label - $seqDiff->offset;#+1+1;
    # if negative DNA positions (before ATG)
    if ($dnamut_start <= 0) { $dnamut_start--; }
    my $dnamut_end;
    $self->mutation->len == 0 ?
	$dnamut_end = $dnamut_start :
	$dnamut_end = $dnamut_start+$self->mutation->len;
    my $dnamut = Bio::Variation::DNAMutation->new(-start => $dnamut_start,
						  -end => $dnamut_end,
						  );
    $dnamut->mut_number($self->mutation->issue);
    $dnamut->isMutation(1);
    my $da_m = Bio::Variation::Allele->new;
    $da_m->seq($self->mutation->seq) if $self->mutation->seq;
    $dnamut->allele_mut($da_m);
    $dnamut->add_Allele($da_m);
    # allele_ori
    my $allele_ori = $self->DNA->labelsubseq($self->mutation->prelabel,
					     undef,
					     $self->mutation->postlabel); # get seq
    chop $allele_ori; # chop the postlabel nucleotide
    $allele_ori=substr($allele_ori,1); # away the praelabel nucleotide
    my $da_o = Bio::Variation::Allele->new;
    $da_o->seq($allele_ori) if $allele_ori;
    $dnamut->allele_ori($da_o);
    $self->mutation->len == 0 ?
	$dnamut->length($self->mutation->len) : $dnamut->length(CORE::length $allele_ori);
    #print " --------------- $dnamut_start -$len-  $dnamut_end -\n";
    $seqDiff->add_Variant($dnamut);
    $self->dnamut($dnamut);
    $dnamut->mut_number($self->mutation->issue);
    # setting proof
    if ($seqDiff->numbering eq "entry") {
	 $dnamut->proof('experimental');
    } else {
	 $dnamut->proof('computed');
    }
    # how many nucleotides to store upstream and downstream of the change
    my $flanklen = $self->{'flanklen'};
    #print  `date`, " flanking sequences start\n";
    my $uplabel = $self->DNA->label(1-$flanklen,$self->mutation->prelabel); # this could be unavailable!

    my $upstreamseq;
    if ($uplabel > 0) {
	 $upstreamseq =
	     $self->DNA->labelsubseq($uplabel, undef, $self->mutation->prelabel);
    } else { # from start (less than $flanklen nucleotides)
	 $upstreamseq =
	     $self->DNA->labelsubseq($self->DNA->start, undef, $self->mutation->prelabel);
    }
    $dnamut->upStreamSeq($upstreamseq);
    my $dnstreamseq = $self->DNA->labelsubseq($self->mutation->postlabel, $flanklen);
    $dnamut->dnStreamSeq($dnstreamseq); # $flanklen or less nucleotides
    return $dnamut;
}


#
### Check if mutation propagates to RNA (and AA) level
### Workaround has to be written for problems with inserted labels
#
# returns a boolean value
#

sub _rnaAffected {
    my ($self) = @_;
    my @exons=$self->RNA->all_Exons;
    my $RNAstart=$self->RNA->start;
    my $RNAend=$self->RNA->end;
    my ($firstexon,$before,$after,$i);
    my ($rnaAffected) = 0;

    # check for inserted labels (that require follows instead of <,>)
    my $DNAend=$self->RNA->{'seq'}->end;
    if ($self->mutation->prelabel > $DNAend or $self->mutation->postlabel > $DNAend) {
	 #this means one of the two labels is an inserted one
	 #(coming from a previous mutation. This would falsify all <,>
	 #checks, so the follow() has to be used
	 carp "Attention, workaround not fully tested yet! Expect unpredictable results.\n";
	 if (($self->mutation->postlabel==$RNAstart) or (follows($self->mutation->postlabel,$RNAstart))) {
	     carp "RNA not affected because change occurs before RNAstart";
	 }
	 elsif (($RNAend==$self->mutation->praelabel) or (follows($RNAend,$self->mutation->prelabel))) {
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
		 if (follows($self->mutation->prelabel,$before) or
			($after==$self->mutation->prelabel) or
			follows($after,$self->mutation->prelabel) or
			follows($after,$self->mutation->postlabel)) {

		     $rnaAffected = 1;
		     # $i is number of exon and can be used for proximity check
		 }
		 $before=$exons[$i]->end;
	     }
	
	 }
    } else {
	 my $strand=$exons[0]->strand;
	 if (($strand == 1 and $self->mutation->postlabel <= $RNAstart) or
	     ($strand != 1 and $self->mutation->postlabel >= $RNAstart)) {
	     carp "RNA not affected because change occurs before RNAstart";
	 }
	 elsif (($strand == 1 and $self->mutation->prelabel >= $RNAend) or
		($strand != 1 and $self->mutation->prelabel <= $RNAend)) {
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
		 if ( ($strand == 1 and ($self->mutation->prelabel < $before ||
					 $self->mutation->prelabel >= $after ||
					 $self->mutation->postlabel > $after)) or
		      ($strand != 1 and ($self->mutation->prelabel > $before ||
					 $self->mutation->prelabel <= $after ||
					 $self->mutation->postlabel < $after)) ) {
                     # check this if really correct
		     $rnaAffected = 1;
		     # $i is number of exon and can be used for proximity check
		 }
		 $before=$exons[$i]->end;
	     }
	
	 }
    }
    #carp "RNA not affected because change occurs inside an intron";
    #return(0); # if still not returned, then not affected, return 0
    return $rnaAffected;
}

#
# ### Creation of RNA and AA variation objects
#

=head2 _set_effects

 Title   : _set_effects
 Usage   :
 Function:

           Stores RNA and AA level mutation attributes before mutation
           into L<Bio::Variation::RNACange> and
           L<Bio::Variation::AACange> objects.  Links them to
           SeqDiff object.

 Example :
 Returns :
 Args    : L<Bio::Variation::SeqDiff> object
           L<Bio::Variation::DNAMutation> object
=cut

sub _set_effects {
    my ($self, $seqDiff, $dnamut) = @_;
    my ($rnapos_end, $upstreamseq, $dnstreamseq);
    my $flanklen = $self->{'flanklen'};

    $self->mutation->len == 0 ?
	$rnapos_end = $self->mutation->transpos :
	$rnapos_end = $self->mutation->transpos + $self->mutation->len -1;
    my $rnachange = Bio::Variation::RNAChange->new(-start => $self->mutation->transpos,
						    -end =>  $rnapos_end
						    );
    $rnachange->isMutation(1);

    # setting proof
    if ($seqDiff->numbering eq "coding") {
	 $rnachange->proof('experimental');
    } else {
	 $rnachange->proof('computed');
    }

    $seqDiff->add_Variant($rnachange);
    $self->rnachange($rnachange);
    $rnachange->DNAMutation($dnamut);
    $dnamut->RNAChange($rnachange);
    $rnachange->mut_number($self->mutation->issue);

    # setting the codon_position of the "start" nucleotide of the change
    $rnachange->codon_pos(($self->RNA->frame($self->mutation->label))+1); # codon_pos=frame+1

    my @exons=$self->RNA->all_Exons;
    $self->exons(\@exons);
    #print  `date`, " before flank, after exons. RNAObj query\n";
    # if cannot retrieve from Transcript, Transcript::upstream_seq will be used
    # before "fac7 g 65" bug discovered
    # $uplabel=$self->RNA->label(1-$flanklen,$praelabel);
    my $RNApraelabel=$self->RNA->label(-1,$self->mutation->label); # to fix fac7g65 bug
    # for the fix, all praelabel used in the next block have been changed to RNApraelabel
    my $uplabel=$self->RNA->label(1-$flanklen,$RNApraelabel);
    if ($self->RNA->valid($uplabel)) {
	 $upstreamseq=$self->RNA->labelsubseq($uplabel, undef, $RNApraelabel);
    } else {
	 $upstreamseq=$self->RNA->labelsubseq($self->RNA->start, undef, $RNApraelabel);
	 my $lacking=$flanklen-length($upstreamseq); # how many missing
	 my $upstream_atg=$exons[0]->subseq(-$lacking,-1);
	 $upstreamseq=$upstream_atg . $upstreamseq;
    }

    #
    # Exon/Intron
    #
    #first only for point mutations, elaborate later
    if ($self->mutation->transpos == $rnapos_end) {
	 my ($affected_exon,$exon_number)=$self->RNA->in_which_Exon($RNApraelabel);
            # changed from praelabel to RNApraelabel (does it work like this??)
	 $rnachange->exons_modified($exon_number);
	 $dnamut->region('exon');
	 $dnamut->region_value($exon_number);
    }
    #
    # Flanking sequences
    #
    $rnachange->upStreamSeq($upstreamseq);

    # won't work OK if postlabel NOT in Transcript
    # now added RNApostlabel but this has to be /fully tested/
    # for the fix, all postlabel used in the next block have been changed to RNApostlabel
    my $RNApostlabel; # to fix fac7g64 bug
    if ($self->mutation->len == 0) {
      $RNApostlabel=$self->mutation->label;
    } else {
      $RNApostlabel=$self->RNA->label(2,$self->mutation->label);
    }
    $dnstreamseq=$self->RNA->labelsubseq($RNApostlabel, $flanklen);
    if ($dnstreamseq eq '-1') { # if out of transcript was requested
	 my $lastexon=$exons[-1];
	 my $lastexonlength=$lastexon->length;
	 $dnstreamseq=$self->RNA->labelsubseq($RNApostlabel); # retrieves till RNAend
	 my $lacking=$flanklen-length($dnstreamseq); # how many missing
	 my $downstream_stop=$lastexon->subseq($lastexonlength+1,undef,$lacking);
	 $dnstreamseq .= $downstream_stop;
    } else {
	 $rnachange->dnStreamSeq($dnstreamseq);
    }
    # AAChange creation
    my $AAobj=$self->RNA->get_Translation;
    # storage of praelabel here, to be used in create_mut_objs_after
    my $aachange = Bio::Variation::AAChange->new(-start => $RNApraelabel
						  );
    $aachange->isMutation(1);
    $aachange->proof('computed');

    $seqDiff->add_Variant($aachange);
    $self->aachange($aachange);
    $rnachange->AAChange($aachange);
    $aachange->RNAChange($rnachange);

    $aachange->mut_number($self->mutation->issue);
#    $before_mutation{aachange}=$aachange;

    my $ra_o = Bio::Variation::Allele->new;
    $ra_o->seq($dnamut->allele_ori->seq) if $dnamut->allele_ori->seq;
    $rnachange->allele_ori($ra_o);

    $rnachange->length(CORE::length $rnachange->allele_ori->seq);

    my $ra_m = Bio::Variation::Allele->new;
    $ra_m->seq($self->mutation->seq) if $self->mutation->seq;
    $rnachange->allele_mut($ra_m);
    $rnachange->add_Allele($ra_m);

    #$rnachange->allele_mut($seq);
    $rnachange->end($rnachange->start) if $rnachange->length == 0;

    # this holds the aminoacid sequence that will be affected by the mutation
    my $aa_allele_ori=$AAobj->labelsubseq($self->mutation->label,undef, $self->mutation->lastlabel);

    my $aa_o = Bio::Variation::Allele->new;
    $aa_o->seq($aa_allele_ori) if $aa_allele_ori;
    $aachange->allele_ori($aa_o);
    #$aachange->allele_ori($aa_allele_ori);

    my $aa_length_ori = length($aa_allele_ori);
    $aachange->length($aa_length_ori); #print "==========$aa_length_ori\n";
    $aachange->end($aachange->start + $aa_length_ori - 1 );
}

=head2 _untranslated

 Title   : _untranslated
 Usage   :
 Function:

           Stores RNA change attributes before mutation
           into L<Bio::Variation::RNAChange object.  Links it to
           SeqDiff object.

 Example :
 Returns :
 Args    : L<Bio::Variation::SeqDiff> object
           L<Bio::Variation::DNAMutation> object

=cut

sub  _untranslated {
    my ($self, $seqDiff, $dnamut) = @_;
    my $rnachange = Bio::Variation::RNAChange->new;
    my $ra_o = Bio::Variation::Allele->new;
    $ra_o->seq($dnamut->allele_ori->seq) if $dnamut->allele_ori->seq;
    $rnachange->allele_ori($ra_o);
    my $ra_m = Bio::Variation::Allele->new;
    $ra_m->seq($dnamut->allele_mut->seq) if $dnamut->allele_mut->seq;
    $rnachange->allele_mut($ra_m);
    $rnachange->add_Allele($ra_m);
    $rnachange->upStreamSeq($dnamut->upStreamSeq);
    $rnachange->dnStreamSeq($dnamut->dnStreamSeq);
    $rnachange->start($dnamut->start);
    $rnachange->end($dnamut->end);
    $rnachange->length($dnamut->length);
    $rnachange->mut_number($dnamut->mut_number);
    # setting proof
    if ($seqDiff->numbering eq "coding") {
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
    $seqDiff->add_Variant($rnachange);
    $self->rnachange($rnachange);
    $rnachange->DNAMutation($dnamut);
    $dnamut->RNAChange($rnachange);
}


# args: reference to label changearray, reference to position changearray
# Function: take care of the creation of mutation objects, with
# information AFTER the change takes place
sub _post_mutation {
    my ($self) = @_;

    if ($self->rnachange->region eq 'coding') {
	 my $aachange=$self->aachange;
	 my ($AAobj,$aa_start_praelabel,$aa_start,$mut_translation);
	 $AAobj=$self->RNA->get_Translation;
	 $aa_start_praelabel=$aachange->start;
	 $aa_start=$AAobj->position($self->RNA->label(2,$aa_start_praelabel));
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
		     $aachange->end($aachange->start + $aachange->length - 1 );
		 }
		 elsif ($aachange->RNAChange->allele_ori->seq eq '' ) {
		     $aachange->allele_mut->seq(substr $aachange->allele_mut->seq, 0,
					   length ($aachange->RNAChange->allele_mut->seq) / 3);
		     $aachange->allele_ori->seq('');
		     $aachange->end($aachange->start + $aachange->length - 1 );
		     $aachange->length(0);
		 }
	     } else {
		 #elsif  ($aachange->RNAChange->codon_pos == 2){
		 if (not $aachange->RNAChange->allele_mut->seq ) {
		     $aachange->allele_mut->seq(substr $aachange->allele_mut->seq, 0, 1);
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
	 my @beforeexons=@{$self->exons};
	 my @afterexons=$self->RNA->all_Exons;
	 my $i;
	 if (scalar(@beforeexons) ne scalar(@afterexons)) {
	     my $mut_number = $self->mutation->issue;
	     carp "Exons have been modified at mutation n.$mut_number!";
	     $self->rnachange->exons_modified(1);
	 } else {
	   EXONCHECK:
	     foreach $i (0..$#beforeexons) {
		 if ($beforeexons[$i] ne $afterexons[$i]) {
	     my $mut_number = $self->mutation->issue;
		     carp "Exons have been modified at mutation n.$mut_number!";
		     $self->rnachange->exons_modified(1);
		     last EXONCHECK;
		 }
	     }
	 }
    }
    return 1;
}

1;
