# $Id$
#
# BioPerl module for Bio::Variation::DNAMutation
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Variation::DNAMutation - DNA level mutation class

=head1 SYNOPSIS

    $dnamut = Bio::Variation::DNAMutation->new
        ('-start'         => $start,
         '-end'           => $end,
         '-length'        => $len,
         '-upStreamSeq'   => $upflank,
         '-dnStreamSeq'   => $dnflank,
         '-proof'         => $proof,
	 '-isMutation'    => 1,
         '-mut_number'    => $mut_number
        );
    $a1 = Bio::Variation::Allele->new;
    $a1->seq('a');
    $dnamut->allele_ori($a1);
    my $a2 = Bio::Variation::Allele->new;
    $a2->seq('t');
    $dnamut->add_Allele($a2);

    print "Restriction changes are ", $dnamut->restriction_changes, "\n";

    # add it to a SeqDiff container object
    $seqdiff->add_Variant($dnamut);


=head1 DESCRIPTION

The instantiable class Bio::Variation::DNAMutation describes basic
sequence changes in genomic DNA level. It uses methods defined in
superclass L<Bio::Variation::VariantI>.

If the variation described by a DNAMutation object is transcibed, link
the corresponding L<Bio::Variation::RNAChange> object to it using
method ->RNAChange.

=head1 FEEDBACK

=head2 Mailing Lists

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

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki@ebi.ac.uk
Address: 

     EMBL Outstation, European Bioinformatics Institute
     Wellcome Trust Genome Campus, Hinxton
     Cambs. CB10 1SD, United Kingdom 

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Variation::DNAMutation;
my $VERSION=1.0;
use vars qw(@ISA);
use strict;

# Object preamble - inheritance
use Bio::Variation::VariantI;

@ISA = qw( Bio::Variation::VariantI );

sub new {
    my($class,@args) = @_;
    my $self;
    $self = {};
    bless $self, $class;
    
    my ($start, $end, $length, $strand, $primary, $source, 
	$frame, $score, $gff_string,
	$allele_ori,  $allele_mut,  $upstreamseq,  $dnstreamseq,  
	$label,  $status,  $proof,  $region, $region_value, $numbering, 
	$cpg, $mut_number, $ismutation) =
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
				  NUMBERING
				  CPG
				  MUT_NUMBER
				  ISMUTATION
				  )],
			      @args);

    $self->SUPER::primary_tag("Variation");

    $self->{ 'alleles' } = [];

    $start && $self->SUPER::start($start);
    $end   && $self->SUPER::end($end);
    $length && $self->SUPER::length($length);
    $strand && $self->SUPER::strand($strand);
    $primary && $self->SUPER::primary_tag($primary);
    $source  && $self->SUPER::source_tag($source);
    $frame   && $self->SUPER::frame($frame);
    $score   && $self->SUPER::score($score);
    $gff_string && $self->SUPER::_from_gff_string($gff_string);
    
    $allele_ori && $self->SUPER::allele_ori($allele_ori);
    $allele_mut  && $self->SUPER::allele_mut($allele_mut);
    $upstreamseq  && $self->SUPER::upStreamSeq($upstreamseq);
    $dnstreamseq  && $self->SUPER::dnStreamSeq($dnstreamseq);
    
    $label  && $self->label($label);
    $status  && $self->SUPER::status($status);
    $proof && $self->SUPER::proof($proof);
    $region  && $self->SUPER::region($region);
    $region_value  && $self->SUPER::region_value($region_value);
    $numbering && $self->SUPER::numbering($numbering);
    $mut_number && $self->SUPER::mut_number($mut_number);
    $ismutation && SUPER::isMutation($ismutation);
    
    $cpg && $self->CpG($cpg);
    
    return $self; # success - we hope!
}


=head2 CpG

 Title   : CpG
 Usage   : $obj->CpG()
 Function: sets and returns boolean values for variation 
           hitting a CpG site.  Unset value return -1.
 Example : $obj->CpG()
 Returns : boolean
 Args    : optional true of false value


=cut


sub CpG {
   my ($obj,$value) = @_;
   if( defined $value) {
       $value ? ($value = 1) : ($value = 0);
       $obj->{'cpg'} = $value;
   }
    elsif (not defined $obj->{'label'}) {
	$obj->{'cpg'} = $obj->_CpG_value;
    }
   else {
       return $obj->{'cpg'};
   }
}



sub _CpG_value {
    my ($self) = @_;
    if ($self->allele_ori eq $self->allele_mut and length ($self->allele_ori) == 1 ) {
    
	# valid only for point mutations
	# CpG methylation-mediated deamination:
	#   CG -> TG | CG -> CA substitutions
	# implementation here is  less strict: if CpG dinucleotide was hit
	
	if ( ( ($self->allele_ori eq 'c') && (substr($self->upStreamSeq, 0, 1) eq 'g') ) ||
	     ( ($self->allele_ori eq 'g') && (substr($self->dnStreamSeq, -1, 1) eq 'c') ) ) {
	    return 1;
	}
	else {
	    return 0;
	}
    } else {
	$self->warn('CpG makes sense only in the context of point mutation');
	return undef;
    }
}


=head2 RNAChange

 Title   : RNAChange
 Usage   : $mutobj = $obj->RNAChange;
         : $mutobj = $obj->RNAChange($objref);
 Function: Returns or sets the link-reference to a mutation/change object.
           If there is no link, it will return undef
 Returns : an obj_ref or undef

=cut


sub RNAChange {
  my ($self,$value) = @_;
  if (defined $value) {
      if( ! $value->isa('Bio::Variation::RNAChange') ) {
	  $self->throw("Is not a Bio::Variation::RNAChange object but a [$self]");
	  return (undef);
      }
      else {
	  $self->{'RNAChange'} = $value;
      }
  }
  unless (exists $self->{'RNAChange'}) {
      return (undef);
  } else {
      return $self->{'RNAChange'};
  }
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
    my ($self, $value) = @_;
    my ($o, $m, $type);
    $o = $self->allele_ori->seq if $self->allele_ori and $self->allele_ori->seq;
    $m = $self->allele_mut->seq if $self->allele_mut and $self->allele_mut->seq;
    
    if (not $o and not $m ) {
	$self->warn("[DNAMutation, label] Both alleles should not be empty!\n");
	$type = 'no change'; # is this enough?
    }
    elsif ($o && $m && length($o) == length($m) && length($o) == 1) {
	$type = 'point';
	$type .= ", ". _point_type_label($o, $m);
    }
    elsif (not $o ) {
	$type = 'insertion';
    }
    elsif (not $m  ) {
	$type = 'deletion';
    }
    else {
	$type = 'complex';
    }
    $self->{'label'} = $type;
    return $self->{'label'};
}


sub _point_type_label {
    my ($o, $m) = @_;
    my ($type);
    my %transition = ('a' => 'g',
		   'g' => 'a',
		   'c' => 't',
		   't' => 'c');
    $o = lc $o;
    $m = lc $m;
    if ($o eq $m) {
	$type = 'no change';
    }
    elsif ($transition{$o} eq $m ) {
	$type = 'transition';
    }
    else {
	$type = 'transversion';
    }
}


=head2 sysname

 Title   : sysname
 Usage   : $self->sysname
 Function: 

           This subroutine creates a string corresponding to the
           'systematic name' of the mutation. Systematic name is
           specified in Antonorakis & MDI Nomenclature Working Group:
           Human Mutation 11:1-3, 1998. 
           http://www.interscience.wiley.com/jpages/1059-7794/nomenclature.html
 Returns : string

=cut


sub sysname {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'sysname'} = $value;
    } else {
	$self->warn('Mutation start position is not defined') 
	    if not defined $self->start;
	my $sysname = '';
	# show the moltype only if $self->SeqDiff->moltype is set;
	my $mol = '';
	if ($self->SeqDiff->moltype && $self->SeqDiff->moltype eq 'dna') {
	    $mol = 'g.';
	}
	elsif ($self->SeqDiff->moltype && $self->SeqDiff->moltype eq 'rna') {
	    $mol = 'c.';
	}
	my $sep;
	if ($self->isMutation) {
	    $sep = '>';
	} else {
	    $sep = '|';
	}
	my $sign = '+'; 
	$sign = '' if $self->start < 1;
	$sysname .=  $mol ;#if $mol;
	$sysname .= $sign. $self->start;

	my @alleles = $self->each_Allele;
	$self->allele_mut($alleles[0]);

	$sysname .= 'del' if $self->label =~ /deletion/;
	$sysname .= 'ins' if $self->label =~ /insertion/;
	$sysname .=  uc $self->allele_ori->seq if $self->allele_ori->seq;



	#push @alleles, $self->allele_mut if $self->allele_mut;
	foreach my $allele (@alleles) {
	    $self->allele_mut($allele);
	    $sysname .= $sep if $self->label =~ /point/;
	    $sysname .=  uc $self->allele_mut->seq if $self->allele_mut->seq;
	}
	$self->{'sysname'} = $sysname;
	#$self->{'sysname'} = $sign. $self->start. 
	#    uc $self->allele_ori->seq. $sep. uc $self->allele_mut->seq;
    }
    return $self->{'sysname'};
}

1;
