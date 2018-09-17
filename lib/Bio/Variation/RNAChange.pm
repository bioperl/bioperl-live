#
# BioPerl module for Bio::Variation::RNAChange
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Variation::RNAChange - Sequence change class for RNA level

=head1 SYNOPSIS

   $rnachange = Bio::Variation::RNAChange->new
       ('-start'         => $start,
        '-end'           => $end,
        '-length'        => $len,
        '-codon_pos'     => $cp,
        '-upStreamSeq'   => $upflank,
        '-dnStreamSeq'   => $dnflank,
        '-proof'         => $proof,
   	'-isMutation'    => 1,
        '-mut_number'    => $mut_number
       );
   $a1 = Bio::Variation::Allele->new;
   $a1->seq('a');
   $rnachange->allele_ori($a1);
   my $a2 = Bio::Variation::Allele->new;
   $a2->seq('t');
   $rnachange->add_Allele($a2);
   $rnachange->allele_mut($a2);

   print "The codon change is ", $rnachange->codon_ori, 
       ">", $rnachange->codon_mut, "\n"; 

   # add it to a SeqDiff container object
   $seqdiff->add_Variant($rnachange);

   # and create links to and from DNA level mutation objects
   $rnachange->DNAMutation($dnamut);
   $dnamut->RNAChange($rnachange);

=head1 DESCRIPTION

The instantiable class Bio::Variation::DNAMutation describes basic
sequence changes at RNA molecule level. It uses methods defined in
superclass Bio::Variation::VariantI. See L<Bio::Variation::VariantI>
for details.

You are normally expected to create a corresponding
Bio::Variation::DNAMutation object even if mutation is defined at
RNA level. The numbering follows then cDNA numbering.  Link the
DNAMutation object to the RNAChange object using the method
DNAMutation(). If the variation described by a RNAChange object is
translated, link the corresponding Bio::Variation::AAChange object
to it using method AAChange(). See L<Bio::Variation::DNAMutation> and
L<Bio::Variation::AAChange> for more information.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the 
Bioperl mailing lists  Your participation is much appreciated.

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

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Variation::RNAChange;
use strict;

# Object preamble - inheritance

use Bio::Tools::CodonTable;

use base qw(Bio::Variation::VariantI);

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($start, $end, $length, $strand, $primary, $source,
        $frame, $score, $gff_string,
        $allele_ori,  $allele_mut,  $upstreamseq,  $dnstreamseq,
	$label,  $status,  $proof,  $region,  $region_value, $region_dist, $numbering,
	$mut_number,  $isMutation,
	$codon_ori, $codon_mut, $codon_pos, $codon_table, $cds_end) =
	    $self->_rearrange([qw(START
				  END
				  LENGTH
				  STRAND
				  PRIMARY
				  SOURCE
				  FRAME
				  SCORE
				  GFF_STRING
				  ALLELE_ORI
				  ALLELE_MUT
				  UPSTREAMSEQ
				  DNSTREAMSEQ
				  LABEL
				  STATUS
				  PROOF
				  REGION
				  REGION_VALUE
				  REGION_DIST
				  NUMBERING
				  MUT_NUMBER
				  ISMUTATION
				  CODON_ORI
				  CODON_MUT
				  CODON_POS
				  TRANSLATION_TABLE
				  CDS_END
				  )],@args);
    
    $self->primary_tag("Variation");
    
    $self->{ 'alleles' } = [];
    
    $start && $self->start($start);
    $end   && $self->end($end);
    $length && $self->length($length);
    $strand && $self->strand($strand);
    $primary && $self->primary_tag($primary);
    $source  && $self->source_tag($source);
    $frame   && $self->frame($frame);
    $score   && $self->score($score);
    $gff_string && $self->_from_gff_string($gff_string);
    
    $allele_ori && $self->allele_ori($allele_ori);
    $allele_mut  && $self->allele_mut($allele_mut);
    $upstreamseq  && $self->upStreamSeq($upstreamseq);
    $dnstreamseq  && $self->dnStreamSeq($dnstreamseq);
    
    $label  && $self->label($label);
    $status  && $self->status($status);
    $proof && $self->proof($proof);
    $region  && $self->region($region);
    $region_value  && $self->region_value($region_value);
    $region_dist  && $self->region_dist($region_dist);
    $numbering && $self->numbering($numbering);
    $mut_number && $self->mut_number($mut_number);
    $isMutation && $self->isMutation($isMutation);
    
    $codon_ori  && $self->codon_ori($codon_ori);
    $codon_mut  && $self->codon_mut($codon_mut);
    $codon_pos  && $self->codon_pos($codon_pos);
    $codon_table && $self->codon_table($codon_table);
    $cds_end  && $self->cds_end($cds_end);
    return $self; # success - we hope!
}


=head2 codon_ori

 Title   : codon_ori
 Usage   : $obj->codon_ori();
 Function: 

            Sets and returns codon_ori triplet.  If value is not set,
            creates the codon triplet from the codon position and
            flanking sequences.  The string has to be three characters
            long. The character content is not checked.

 Example : 
 Returns : string
 Args    : string

=cut

sub codon_ori {
    my ($self,$value) = @_;
    if (defined $value) {
	if (length $value != 3) {
	    $self->warn("Codon string \"$value\" is not three characters long");
	}
	$self->{'codon_ori'} = $value;
    }
    elsif (! $self->{'codon_ori'}) {
	my $codon_ori = '';

	if ($self->region eq 'coding' && $self->start && $self->start  >= 1) {
	    
	    $self->warn('Codon position is not defined') 
		if not defined $self->codon_pos;
	    $self->warn('Upstream flanking sequence  is not defined') 
		if not defined $self->upStreamSeq;
	    $self->warn('Downstream flanking sequence  is not defined') 
		if not defined $self->dnStreamSeq;

	    my $cpos = $self->codon_pos; 
	    $codon_ori = substr($self->upStreamSeq, -$cpos +1  , $cpos-1);
	    $codon_ori .= substr($self->allele_ori->seq, 0, 4-$cpos) 
		if $self->allele_ori and $self->allele_ori->seq;
	    $codon_ori .= substr($self->dnStreamSeq, 0, 3-length($codon_ori));
	}
	$self->{'codon_ori'} = lc $codon_ori;
    }
    return $self->{'codon_ori'};
}


=head2 codon_mut

 Title   : codon_mut
 Usage   : $obj->codon_mut();
 Function: 

            Sets and returns codon_mut triplet.  If value is not
            set, creates the codon triplet from the codon position and
            flanking sequences. Return undef for other than point mutations.

 Example : 
 Returns : string
 Args    : string

=cut


sub codon_mut {
    my ($self,$value) = @_;
    if (defined $value) {
	if (length $value != 3 ) {
	    $self->warn("Codon string \"$value\" is not three characters long");
	}
	$self->{'codon_mut'} = $value;
    }
    else {
	my $codon_mut = '';
	if ($self->allele_ori->seq and $self->allele_mut->seq and
	  CORE::length($self->allele_ori->seq) == 1 and 
	  CORE::length($self->allele_mut->seq) == 1 and
	    $self->region eq 'coding' and $self->start >= 1) {

	    $self->warn('Codon position is not defined') 
		if not defined $self->codon_pos;
	    $self->warn('Upstream flanking sequnce  is not defined') 
		if not defined $self->upStreamSeq;
	    $self->warn('Downstream flanking sequnce  is not defined') 
		if not defined $self->dnStreamSeq;
	    $self->throw('Mutated allele is not defined') 
		if not defined $self->allele_mut;
	    
	    my $cpos = $self->codon_pos;
	    $codon_mut = substr($self->upStreamSeq, -$cpos +1  , $cpos-1);
	    $codon_mut .= substr($self->allele_mut->seq, 0, 4-$cpos) 
		if $self->allele_mut and $self->allele_mut->seq; 
	    $codon_mut .= substr($self->dnStreamSeq, 0, 3-length($codon_mut));
	    
	    $self->{'codon_mut'} = lc $codon_mut;
	}
    }
    return $self->{'codon_mut'};
}


=head2 codon_pos

 Title   : codon_pos
 Usage   : $obj->codon_pos();
 Function: 

            Sets and returns the position of the mutation start in the
            codon. If value is not set, returns false.

 Example : 
 Returns : 1,2,3
 Args    : none if get, the new value if set

=cut


sub codon_pos {
    my ($self,$value) = @_;
    if( defined $value) {
	if ( $value !~ /[123]/ ) {
	    $self->throw("'$value' is not a valid codon position");
	}
	$self->{'codon_pos'} = $value;
    }
    return $self->{'codon_pos'};
}


=head2 codon_table

 Title   : codon_table
 Usage   : $obj->codon_table();
 Function: 

            Sets and returns the codon table id of the RNA
            If value is not set, returns 1, 'universal' code, as the default.

 Example : 
 Returns : integer
 Args    : none if get, the new value if set

=cut


sub codon_table {
    my ($self,$value) = @_;
    if( defined $value) {
	if (  not $value =~ /^\d$/ ) {
	    $self->throw("'$value' is not a valid codon table ID\n".
			"Has to be a positive integer. Defaulting to 1\n");
	} else {
	    $self->{'codon_table'} = $value;
	}
    }
    if( ! exists $self->{'codon_table'} ) {
	return 1;
    } else {
	return $self->{'codon_table'};
    }
}


=head2 DNAMutation

 Title   : DNAMutation
 Usage   : $mutobj = $obj->DNAMutation;
         : $mutobj = $obj->DNAMutation($objref);
 Function: Returns or sets the link-reference to a mutation/change object.
           If there is no link, it will return undef
 Returns : an obj_ref or undef

=cut


sub DNAMutation {
    my ($self,$value) = @_;
    if (defined $value) {
	if( ! $value->isa('Bio::Variation::DNAMutation') ) {
	    $self->throw("Is not a Bio::Variation::DNAMutation object but a [$self]");
	    return;
	}
	else {
	    $self->{'DNAMutation'} = $value;
	}
    }
    unless (exists $self->{'DNAMutation'}) {
	return;
    } else {
	return $self->{'DNAMutation'};
    }
}


=head2 AAChange

 Title   : AAChange
 Usage   : $mutobj = $obj->AAChange;
         : $mutobj = $obj->AAChange($objref);
 Function: Returns or sets the link-reference to a mutation/change object.
           If there is no link, it will return undef
 Returns : an obj_ref or undef

=cut

sub AAChange {
    my ($self,$value) = @_;
    if (defined $value) {
	if( ! $value->isa('Bio::Variation::AAChange') ) {
	    $self->throw("Is not a Bio::Variation::AAChange object but a [$self]");
	return;
	}
	else {
	    $self->{'AAChange'} = $value;
	}
    }
    unless (exists $self->{'AAChange'}) {
	return;
    } else {
	return $self->{'AAChange'};
    }
}    


=head2 exons_modified

 Title   : exons_modified
 Usage   : $modified = $obj->exons_modified;
         : $modified = $obj->exons_modified(1);
 Function: Returns or sets information (example: a simple boolean flag) about
           the modification of exons as a result of a mutation.

=cut

sub exons_modified {
  my ($self,$value)=@_;
  if (defined($value)) {
    $self->{'exons_modified'}=$value;
  }
  return ($self->{'exons_modified'});
}

=head2 region

 Title   : region
 Usage   : $obj->region();
 Function: 

            Sets and returns the name of the sequence region type or
            protein domain at this location.  If value is not set,
            returns false.

 Example : 
 Returns : string
 Args    : string

=cut



sub region {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'region'} = $value;
    } 
    elsif (not defined $self->{'region'}) {

	$self->warn('Mutation start position is not defined') 
	    if not defined $self->start and $self->verbose;
	$self->warn('Mutation end position is not defined') 
	    if not defined $self->end and $self->verbose;
	$self->warn('Length of the CDS is not defined, the mutation can be beyond coding region!')
	    if not defined $self->cds_end and $self->verbose;
	
	$self->region('coding');
	if ($self->end && $self->end < 0 ){
	    $self->region('5\'UTR');
	}
	elsif ($self->start && $self->cds_end && $self->start > $self->cds_end ) {
	    $self->region('3\'UTR');
	}
    }
    return $self->{'region'};
}

=head2 cds_end

 Title   : cds_end
 Usage   : $cds_end = $obj->get_cds_end();
 Function: 

           Sets or returns the cds_end from the beginning of the DNA sequence
           to the coordinate start used to describe variants.
           Should be the location of the last nucleotide of the
           terminator codon of the gene.

 Example : 
 Returns : value of cds_end, a scalar
 Args    : 

=cut



sub cds_end {
    my ($self, $value) = @_;
    if (defined $value) {
	$self->warn("[$value] is not a good value for sequence position") 
	    if not $value =~ /^\d+$/ ;
	$self->{'cds_end'} = $value;
    } else {
	$self->{'cds_end'} = $self->SeqDiff->cds_end if $self->SeqDiff;
    }
    return $self->{'cds_end'};
}


=head2 label

 Title   : label
 Usage   : $obj->label();
 Function: 

            Sets and returns mutation event label(s).  If value is not
            set, or no argument is given returns false.  Each
            instantiable subclass of L<Bio::Variation::VariantI> needs
            to implement this method. Valid values are listed in
            'Mutation event controlled vocabulary' in
            http://www.ebi.ac.uk/mutations/recommendations/mutevent.html.

 Example : 
 Returns : string
 Args    : string

=cut

sub label {
    my ($self) = @_;
    my ($o, $m, $type);
    $o = $self->allele_ori->seq if $self->allele_ori and $self->allele_ori->seq;
    $m = $self->allele_mut->seq if $self->allele_mut and $self->allele_mut->seq;

    my $ct  = Bio::Tools::CodonTable -> new ( -id => $self->codon_table );
    if ($o and $m and CORE::length($o) == 1 and CORE::length($m) == 1) { 
	if (defined $self->AAChange) {
	    if ($self->start > 0 and $self->start < 4 ) {
		$type = 'initiation codon';
	    }
	    elsif ($self->codon_ori && $ct->is_ter_codon($self->codon_ori) ) {
		#AAChange->allele_ori and $self->AAChange->allele_ori->seq eq '*' ) {
		$type = 'termination codon';
	    }
	    elsif ($self->codon_mut && $ct->is_ter_codon($self->codon_mut) ) {
		#elsif ($self->AAChange->allele_mut and $self->AAChange->allele_mut->seq eq "*") {
		$type = 'nonsense';
	    } 
	    elsif ($o and $m and ($o eq $m or 
				  $self->AAChange->allele_ori->seq eq 
				  $self->AAChange->allele_mut->seq)) {
		$type = 'silent';
	    } else {
		$type = 'missense';
	    }
	} else {
	    $type = 'unknown';
	}
    }  else {
	my $len = 0;
	$len = CORE::length($o) if $o;
	$len -= CORE::length($m) if $m;
	if ($len%3 == 0 ) {
	    $type = 'inframe';
	} else {
	    $type = 'frameshift';
	}
	if (not $m ) {
	    $type .= ', '. 'deletion';
	}
	elsif (not $o ) {
	    $type .= ', '. 'insertion';
	}
	else {
	    $type .= ', '. 'complex';
	}	
	if ($self->codon_ori && $ct->is_ter_codon($self->codon_ori) ) {
	    $type .= ', '. 'termination codon';
	}
    }

    $self->{'label'} = $type;
    return $self->{'label'};
}


=head2 _change_codon_pos

 Title   : _change_codon_pos
 Usage   : $newCodonPos = _change_codon_pos($myCodonPos, 5)
 Function: 

           Keeps track of the codon position in a changeing sequence

 Returns : codon_pos = integer 1, 2 or 3
 Args    : valid codon position 
           signed integer offset to a new location in sequence

=cut


sub _change_codon_pos ($$)  {
    my ($cpos, $i) = @_;

    $cpos = ($cpos + $i%3)%3;
    if ($cpos > 3 ) {
	$cpos = $cpos - 3;
    }
    elsif ($cpos < 1 ) {
	$cpos = $cpos + 3;
    }
    return $cpos;
}

1;
