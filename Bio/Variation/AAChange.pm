# $Id$
#
# BioPerl module for Bio::Variation::AAChange
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Variation::AAChange - Sequence change class for polypeptides

=head1 SYNOPSIS

   $aamut = Bio::Variation::AAChange->new 
       ('-start'         => $start,
 	'-end'           => $end,
 	'-length'        => $len,
 	'-proof'         => $proof,
 	'-isMutation'    => 1,
 	'-mut_number'    => $mut_number             
 	);

   my $a1 = Bio::Variation::Allele->new;
   $a1->seq($ori) if $ori;
   $aamut->allele_ori($a1);
   my $a2 = Bio::Variation::Allele->new;
   $a2->seq($mut) if $mut;
   $aachange->add_Allele($a2);
   $aachange->allele_mut($a2);

   print  "\n"; 

   # add it to a SeqDiff container object
   $seqdiff->add_Variant($rnachange);

   # and create links to and from RNA level variant objects
   $aamut->RNAChange($rnachange);
   $rnachange->AAChange($rnachange);

=head1 DESCRIPTION

The instantiable class Bio::Variation::RNAChange describes basic
sequence changes at polypeptide  level. It uses methods defined in
superclass L<Bio::Variation::VariantI>.

If the variation described by a AAChange object has a known
L<Bio::Variation::RNAAChange> object, create the link with method
AAChange().

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


package Bio::Variation::AAChange;
my $VERSION=1.0;
use vars qw(@ISA $MATRIX);
use strict;

# Object preamble - inheritance
use Bio::Variation::VariantI;

@ISA = qw( Bio::Variation::VariantI );

BEGIN {

my $matrix = << "__MATRIX__";
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
__MATRIX__

    my %blosum = ();
    $matrix =~ /^ +(.+)$/m;
    my @aas = split / +/, $1;
    foreach my $aa (@aas) {
	my $tmp = $aa;
	$tmp = "\\$aa" if $aa eq '*';
	$matrix =~ /^($tmp) +([-+]?\d.*)$/m;
	my @scores = split / +/, $2 if defined $2;
	my $count = 0;
	foreach my $ak (@aas) {
	    $blosum{$aa}->{$aas[$count]} = $scores[$count];
	    $count++;
	}
    }
    sub _matrix;
    $MATRIX = \%blosum;
}

sub new {
    my($class,@args) = @_;
    my $self;
    $self = {};
    bless $self, $class;

    my ($start, $end, $length, $strand, $primary, $source, 
	$frame, $score, $gff_string,
	$allele_ori,  $allele_mut,  $upstreamseq,  $dnstreamseq,  
	$label,  $status,  $proof,  $re_changes,  $region, $region_value, 
        $region_dist, 
	$numbering,  $mut_number,  $ismutation) =
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
				  RE_CHANGES
				  REGION
				  REGION_VALUE
				  REGION_DIST
				  NUMBERING
				  MUT_NUMBER
				  ISMUTATION
				  )],@args);
    
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
    $upstreamseq  && $self->SUPER::upstreamseq($upstreamseq);
    $dnstreamseq  && $self->SUPER::dnstreamseq($dnstreamseq);

    $label  && $self->label($label);
    $status  && $self->SUPER::status($status);
    $proof && $self->SUPER::proof($proof);
    $region  && $self->SUPER::region($region);
    $region_value  && $self->SUPER::region_value($region_value);
    $region_dist  && $self->SUPER::region_dist($region_dist);
    $numbering && $self->SUPER::numbering($numbering);
    $mut_number && $self->SUPER::mut_number($mut_number);
    $ismutation && $self->SUPER::isMutation($ismutation);

    return $self; # success - we hope!
}

=head2 RNAChange

 Title   : RNAChange
 Usage   : $mutobj = $self->RNAChange;
         : $mutobj = $self->RNAChange($objref);
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
    my ($self) = @_;
    my ($o, $m, $type);
    $o = $self->allele_ori->seq if $self->allele_ori and $self->allele_ori->seq;
    $m = $self->allele_mut->seq if $self->allele_mut and $self->allele_mut->seq;

    if ($self->start == 1 ) {
	if ($o and substr($o, 0, 1) ne substr($m, 0, 1)) {
	    $type = 'no translation';
	}
	elsif ($o and $m and $o eq $m ) {
	    $type = 'silent';
	}
	# more ...
    }
    elsif ($o and substr($o, 0, 1) eq '*' ) {
	if ($m and substr($o, 0, 1) ne substr($m, 0, 1)) {
	    $type = 'post-elongation';
	}
	elsif ($m and $o eq $m ) {
	    $type = 'silent, conservative';
	}
    }
    elsif ($o and $m and $o eq $m) {
	$type = 'silent, conservative';
    }
    elsif ($m and $m eq '*') {
	$type = 'truncation';
    }
    elsif ($o and $m and $o eq $m) {
	$type = 'silent, conservative';
    }
    elsif (not $m or 
	   ($o and $m and  length($o) > length($m) and 
	    substr($m, -1, 1) ne '*')) {
	$type = 'deletion';
	if ($o and $m and $o !~ $m and $o !~ $m) {
	    $type .= ', complex'; 
	}
    }
    elsif (not $o or 
	   ($o and $m and length($o) < length($m) and 
	    substr($m, -1, 1) ne '*' ) ) {
	$type = 'insertion';	
	if ($o and $m and $o !~ $m and $o !~ $m) {
	    $type .= ', complex'; 
	}
    }
    elsif  ($o and $m and $o ne $m and 
	    length $o == 1 and  length $m  == 1 ) {
	$type = 'substitution';
	my $value = $self->similarity_score;
	if (defined $value) {
	    my $cons = ($value < 0) ? 'nonconservative' : 'conservative';
	    $type .= ", ". $cons;
	}
    } else {
	$type = 'out-of-frame translation, truncation';
    }
    $self->{'label'} = $type;
    return $self->{'label'};
}


=head2 similarity_score

 Title   : similarity_score
 Usage   : $self->similarity_score
 Function: Measure for evolutionary conservativeness
           of single amino substitutions. Uses BLOSUM62.
           Negative numbers are noncoservative changes.
 Returns : integer, undef if not single amino acid change

=cut

sub similarity_score {
    my ($self) = @_;
    my ($o, $m, $type);
    $o = $self->allele_ori->seq if $self->allele_ori and $self->allele_ori->seq;
    $m = $self->allele_mut->seq if $self->allele_mut and $self->allele_mut->seq;
    return undef unless $o and $m and length $o == 1 and length $m == 1;
    return undef unless $o =~ /[ARNDCQEGHILKMFPSTWYVBZX*]/i and $m =~ /[ARNDCQEGHILKMFPSTWYVBZX*]/i;
    return $MATRIX->{"\U$o"}->{"\U$m"};
}

=head2 trivname

 Title   : trivname
 Usage   : $self->trivname
 Function: 

           Given a Bio::Variation::AAChange object with linked
           Bio::Variation::RNAChange and Bio::Variation::DNAMutation
           objects, this subroutine creates a string corresponding to
           the 'trivial name' of the mutation. Trivial name is
           specified in Antonorakis & MDI Nomenclature Working Group:
           Human Mutation 11:1-3, 1998.
           http://www.interscience.wiley.com/jpages/1059-7794/nomenclature.html

 Returns : string

=cut


sub trivname {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'trivname'} = $value;
    } else {
	my ( $aaori, $aamut,$aamutsymbol, $aatermnumber, $aamutterm) = 
	    ('', '', '', '', '');
	my $o = $self->allele_ori->seq if $self->allele_ori and $self->allele_ori->seq;
	#my $m = $self->allele_mut->seq if $self->allele_mut and $self->allele_mut->seq;

	$aaori = substr ($o, 0, 1) if $o;
	$aaori =~ tr/\*/X/;

	my $sep;
	if ($self->isMutation) {
	    $sep = '>';
	} else {
	    $sep = '|';
	}
	my $trivname = $aaori. $self->start;
	$trivname .= $sep if $sep eq '|';

	my @alleles = $self->each_Allele;
	foreach my $allele (@alleles) {
	    my $m = $allele->seq if $allele->seq;

	    $self->allele_mut($allele);
	    #$trivname .=  $sep. uc $m if $m;	    
	    
	    $aamutterm = substr ($m, -1, 1) if $m;
	    if ($self->RNAChange->label =~ /initiation codon/ and 
		( $o and $m and $o ne $m)) {
		$aamut = 'X';
	    } 
	    elsif (CORE::length($o) == 1 and CORE::length($m) == 1 ) {
		$aamutsymbol = '';
		$aamut = $aamutterm;
	    }
	    elsif ($self->RNAChange->label =~ /deletion/) {
		$aamutsymbol = 'del';	    
		if ($aamutterm eq '*') {
		    $aatermnumber = $self->start + length($m) -1;
		    $aamut = 'X'. $aatermnumber;
		}	 
		if ($self->RNAChange  && $self->RNAChange->label =~ /inframe/){
		    $aamut = '-'. length($self->RNAChange->allele_ori->seq)/3 ;
		}		
	    }
	    elsif ($self->RNAChange->label =~ /insertion/ or 
		   $self->RNAChange->label =~ /complex/) {
		$aamutsymbol = 'ins';
		if (($aamutterm eq '*') && (length($m)-1 != 0)) {
		    $aatermnumber = $self->start + length($m)-1;
		    $aamut =  $aatermnumber. 'X';
		}
		if ($self->RNAChange->label =~ /inframe/){
		    $aamut = '+'. length($self->RNAChange->allele_mut->seq)/3 ;
		}
	    }
	    elsif ($self->label =~ /truncation/) {
		$aamut = $m; 
	    } else {
		$aamutsymbol = '';
		$aamut = $aamutterm;
	    }	
	    $aamut =~ tr/\*/X/;
	    $trivname .= $aamutsymbol. $aamut. $sep;
	}
	chop $trivname;
	$self->{'trivname'} = $trivname; 
    }
    return $self->{'trivname'};
}

1;
