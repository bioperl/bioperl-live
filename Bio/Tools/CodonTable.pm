# $Id$
#
# bioperl module for Bio::Tools::CodonTable
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::CodonTable - Bioperl codon table object

=head1 SYNOPSIS

  This is a read-only class for all known codon tables.  The IDs are
  the ones used by nucleotide sequence databases.  All common IUPAC
  ambiguity codes for DNA, RNA and animo acids are recognized.

  # to use
  use Bio::Tools::CodonTable;

  # defaults to ID 1 "Standard"
  $myCodonTable   = Bio::Tools::CodonTable->new();
  $myCodonTable2  = Bio::Tools::CodonTable -> new ( -id => 3 );

  # change codon table
  $myCodonTable->id(5);

  # examine codon table
  print  join (' ', "The name of the codon table no.", $myCodonTable->id(4),
	       "is:", $myCodonTable->name(), "\n");

  # translate a codon
  $aa = $myCodonTable->translate('ACU');
  $aa = $myCodonTable->translate('act');
  $aa = $myCodonTable->translate('ytr');

  # reverse translate an amino acid
  @codons = $myCodonTable->revtranslate('A');
  @codons = $myCodonTable->revtranslate('Ser');
  @codons = $myCodonTable->revtranslate('Glx');
  @codons = $myCodonTable->revtranslate('cYS', 'rna');

  #boolean tests
   print "Is a start\n"       if $myCodonTable->is_start_codon('ATG');
   print "Is a termianator\n" if $myCodonTable->is_ter_codon('tar');
   print "Is a unknown\n"     if $myCodonTable->is_unknown_codon('JTG');

=head1 DESCRIPTION

Codon tables are also called translation tables or genetics codes
since that is what they try to represent. A bit more complete picture
of the full complexity of codon usage in various taxonomic groups
presented at the NCBI Genetic Codes Home page.


CodonTable is a BioPerl class that knows all current translation
tables that are used by primary nucleotide sequence databases
(GenBank, EMBL and DDBJ). It provides methods to output information
about tables and relationships between codons and amino acids.

This class and its methods recognized all common IUPAC ambiguity codes
for DNA, RNA and animo acids. The translation method follows the
conventions in EMBL and TREMBL databases.

It is a nuisance to separate RNA and cDNA representations of nucleic
acid transcripts. The CodonTable object accepts codons of both type as
input and allows the user to set the mode for output when reverse
translating. Its default for output is DNA.

Note: This class deals with individual codons and amino acids, only.
      Call it from your own objects to translate and reverse translate
      longer sequences.


The amino acid codes are IUPAC recommendations for common amino acids:

          A           Ala            Alanine
          R           Arg            Arginine
          N           Asn            Asparagine
          D           Asp            Aspartic acid
          C           Cys            Cysteine
          Q           Gln            Glutamine
          E           Glu            Glutamic acid
          G           Gly            Glycine
          H           His            Histidine
          I           Ile            Isoleucine
          L           Leu            Leucine
          K           Lys            Lysine
          M           Met            Methionine
          F           Phe            Phenylalanine
          P           Pro            Proline
          S           Ser            Serine
          T           Thr            Threonine
          W           Trp            Tryptophan
          Y           Tyr            Tyrosine
          V           Val            Valine
          B           Asx            Aspartic acid or Asparagine
          Z           Glx            Glutamine or Glutamic acid
          X           Xaa            Any or unknown amino acid


It is worth noting that, "Bacterial" codon table no. 11 produces an
polypeptide that is, confusingly, identical to the standard one. The
only differences are in available initiator codons.


NCBI Genetic Codes home page:
     http://www.ncbi.nlm.nih.gov/htbin-post/Taxonomy/wprintgc?mode=c

EBI Translation Table Viewer:
     http://www.ebi.ac.uk/cgi-bin/mutations/trtables.cgi

Amended ASN.1 version with ids 16 and 21 is at:
     ftp://ftp.ebi.ac.uk/pub/databases/geneticcode/

Thank your for Matteo diTomasso for the original Perl implementation
of these tables.

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

package Bio::Tools::CodonTable;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI
use Bio::Root::RootI;
@ISA = qw(Bio::Root::RootI);

# first set internal values for all translation tables

my @names =  #id
    (
     'Standard', #1
     'Vertebrate Mitochondrial',#2
     'Yeast Mitochondrial',# 3
     'Mold, Protozoan, and CoelenterateMitochondrial and Mycoplasma/Spiroplasma',#4
     'Invertebrate Mitochondrial',#5
     'Ciliate, Dasycladacean and Hexamita Nuclear',# 6
       '', '',
     'Echinoderm Mitochondrial',#9
     'Euplotid Nuclear',#10
     '"Bacterial"',# 11
     'Alternative Yeast Nuclear',# 12
     'Ascidian Mitochondrial',# 13
     'Flatworm Mitochondrial',# 14
     'Blepharisma Nuclear',# 15
     'Chlorophycean Mitochondrial',# 16
       '', '',  '', '',
     'Trematode Mitochondrial',# 21
     'Scenedesmus obliquus Mitochondrial', #22
     'Thraustochytrium Mitochondrial' #23
     );

my @tables =
    qw(
       FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
       FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG
       FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG
       FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
       FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG
       FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
       '' ''
       FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
       FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
       FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
       FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
       FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG
       FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
       FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
       FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
       '' '' '' ''
       FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG   
       FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
       FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
       );


my @starts =
    qw(
       ---M---------------M---------------M----------------------------
       --------------------------------MMMM---------------M------------
       ----------------------------------MM----------------------------
       --MM---------------M------------MMMM---------------M------------
       ---M----------------------------MMMM---------------M------------
       -----------------------------------M----------------------------
       '' ''
       -----------------------------------M----------------------------
       -----------------------------------M----------------------------
       ---M---------------M------------MMMM---------------M------------
       -------------------M---------------M----------------------------
       -----------------------------------M----------------------------
       -----------------------------------M----------------------------
       -----------------------------------M----------------------------
       -----------------------------------M----------------------------
       '' ''  '' ''
       -----------------------------------M---------------M------------  
       -----------------------------------M----------------------------
       --------------------------------M--M---------------M------------
       );

my @nucs = qw(t c a g);
my $x = 0;
my ($codons, $trCol);
for my $i (@nucs) {
    for my $j (@nucs) {
        for my $k (@nucs) {
            my $codon = "$i$j$k";
            $codons->{$codon} = $x;
            $trCol->{$x} = $codon;
            $x++;
        }
    }
}

my  %onecode =
    ('Ala' => 'A',     'Asx' => 'B',
     'Cys' => 'C',     'Asp' => 'D',
     'Glu' => 'E',     'Phe' => 'F',
     'Gly' => 'G',     'His' => 'H',
     'Ile' => 'I',     'Lys' => 'K',
     'Leu' => 'L',     'Met' => 'M',
     'Asn' => 'N',     'Pro' => 'P',
     'Gln' => 'Q',     'Arg' => 'R',
     'Ser' => 'S',     'Thr' => 'T',
     'Val' => 'V',     'Trp' => 'W',
     'Xaa' => 'X',     'Tyr' => 'Y',
     'Glx' => 'Z',     'Ter' => '*'
     );

my %iupac_dna =
    ( 'a' => [qw( a       )],
      'c' => [qw( c       )],
      'g' => [qw( g       )],
      't' => [qw( t       )],
      'u' => [qw( t       )],
      'm' => [qw( a c     )],
      'r' => [qw( a g     )],
      'w' => [qw( a t     )],
      'k' => [qw( g t     )],
      'y' => [qw( c t     )],
      's' => [qw( c g     )],
      'v' => [qw( a c g   )],
      'h' => [qw( a c t   )],
      'd' => [qw( a g t   )],
      'b' => [qw( c g t   )],
      'n' => [qw( a c g t )],
      'x' => [qw( a c g t )],
     );

my %iupac_aa =
    ( A => [qw(A)],      B => [qw(D N)],
      C => [qw(C)],      D => [qw(D)],
      E => [qw(E)],      F => [qw(F)],
      G => [qw(G)],      H => [qw(H)],
      I => [qw(I)],      K => [qw(K)],
      L => [qw(L)],      M => [qw(M)],
      N => [qw(N)],      P => [qw(P)],
      Q => [qw(Q)],      R => [qw(R)],
      S => [qw(S)],      T => [qw(T)],
      U => [qw(U)],      V => [qw(V)],
      W => [qw(W)],      X => [qw(X)],
      Y => [qw(Y)],      Z => [qw(E Q)],
      '*' => ['*']
      );

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my($id) =
	$self->_rearrange([qw(ID
			     )],
			 @args);

    $id = 1 if ( ! $id );
    $id  && $self->id($id);
    return $self; # success - we hope!
}

=head2 id

 Title   : id
 Usage   : $obj->id(3); $id_integer = $obj->id();
 Function:

           Sets or returns the id of the translation table.  IDs are
           integers from 1 to 15, excluding 7 and 8 which have been
           removed as redundant. If an invalid ID is given the method
           returns 0, false.


 Example :
 Returns : value of id, a scalar, 0 if not a valid
 Args    : newvalue (optional)

=cut

sub id{
   my ($self,$value) = @_;
   if( defined $value) {
       if (  !(defined $tables[$value-1]) or $tables[$value-1] eq '') {
	   $self->warn("Not a valid codon table ID [$value] ");
	   $value = 0;
       }
       $self->{'id'} = $value;
   }
   return $self->{'id'};
}

=head2 name

 Title   : name
 Usage   : $obj->name()
 Function: returns the descriptive name of the translation table
 Example :
 Returns : A string
 Args    : None


=cut

sub name{
   my ($self) = @_;

   my ($id) = $self->{'id'};
   return $names[$id-1];

}

=head2 translate

 Title   : translate
 Usage   : $obj->translate('YTR')
 Function: Returns a string of one letter amino acid codes from 
           nucleotide sequence input. The imput can be of any length.

           Returns 'X' for unknown codons and codons that code for
           more than one amino acid. Returns an empty string if input
           is not three characters long. Exceptions for these are:

             - IUPAC amino acid code B for Aspartic Acid and
               Asparagine, is used.
             - IUPAC amino acid code Z for Glutamic Acid, Glutamine is
               used.
             - if the codon is two nucleotides long and if by adding
               an a third character 'N', it codes for a single amino
               acid (with exceptions above), return that, otherwise
               return empty string.

           Returns empty string for other input strings that are not
           three characters long.

 Example :
 Returns : a string of one letter ambiguous IUPAC amino acid codes
 Args    : ambiguous IUPAC nucleotide string


=cut

sub translate {
    my ($self, $seq) = @_;
    my $id = $self->id;

    my ($partial) = 0;
    $partial = 2 if length($seq) % 3 == 2;
    
    $seq = lc $seq; 
    $seq =~ tr/u/t/;
    my $protein = "";
    if ($seq =~ /[^actg]/ ) { #ambiguous chars
        for (my $i = 0; $i < (length($seq) - 2 ); $i+=3) {
            my $triplet = substr($seq, $i, 3);
	    if (exists $codons->{$triplet}) {
		$protein .= substr($tables[$id-1], $codons->{$triplet},1);
	    } else {
		$protein .= $self->_translate_ambiguous_codon($triplet);
	    }
	}
    } else { # simple, strict translation
	for (my $i = 0; $i < (length($seq) - 2 ); $i+=3) {
            my $triplet = substr($seq, $i, 3); 
            if (exists $codons->{$triplet}) {
                $protein .= substr($tables[$id-1], $codons->{$triplet}, 1);
	    } else {
                $protein .= 'X';
            }
        }
    }
    if ($partial == 2) { # 2 overhanging nucleotides
	my $triplet = substr($seq, ($partial -4)). "n";
	if (exists $codons->{$triplet}) {
	    my $aa = substr($tables[$id-1], $codons->{$triplet},1);       
	    $protein .= $aa;
	} else {
	    $protein .= $self->_translate_ambiguous_codon($triplet, $partial);
	}
    }
    return $protein;
}

sub _translate_ambiguous_codon {
    my ($self, $triplet, $partial) = @_;
    $partial ||= 0;
    my $id = $self->id;
    my $aa;
    my @codons = _unambiquous_codons($triplet);
    my %aas =();
    foreach my $codon (@codons) {
	$aas{substr($tables[$id-1],$codons->{$codon},1)} = 1;
    }
    my $count = scalar keys %aas;
    if ( $count == 1 ) {
	$aa = (keys %aas)[0];
    }
    elsif ( $count == 2 ) {
	if ($aas{'D'} and $aas{'N'}) {
	    $aa = 'B';
	}
	elsif ($aas{'E'} and $aas{'Q'}) {
	    $aa = 'Z';
	} else {
	    $partial ? ($aa = '') : ($aa = 'X');
	}
    } else {
	$partial ? ($aa = '') :  ($aa = 'X');
    }
    return $aa;
}

sub translate_old{
   my ($self, $value) = @_;
   my ($id) = $self->{'id'};
   my ($partial) = 0;
   my $result;

   $value  = lc $value;
   $value  =~ tr/u/t/;

   if (length $value != 3 and length $value != 2) {
       return '';
   }
   elsif ($value =~ /[^atgc]/i or length $value == 2 ) {
       if (length $value == 2 ) {
	   $value = $value. 'n';
	   $partial = 1;
       }
       my @codons = _unambiquous_codons($value);

       my %aas =();
       foreach my $codon (@codons) {
	   $aas{substr($tables[$id-1],$codons->{$codon},1)} = 1;	
       }
       #foreach my $x (keys %aas) {print "$x\n";} ###

       my $count = scalar keys %aas;
       if ( $count == 1 ) {
	   return (keys %aas)[0];
       }
       elsif ( $count == 2 ) {
	   if ($aas{'D'} and $aas{'N'}) {
	       return 'B';
	   }
	   elsif ($aas{'E'} and $aas{'Q'}) {
	       return 'Z';
	   } else {
	       $partial ? return '' :  return 'X';
	   }
       } else {
	   $partial ? return '' :  return 'X';
       }
   } else {
       return translate_strict (@_);
   }
}

=head2 translate_strict

 Title   : translate_strict
 Usage   : $obj->translate_strict('ACT')
 Function: returns one letter amino acid code for a codon input

           Fast and simple translation. User is responsible to resolve
           ambiguous nucleotide codes before calling this
           method. Returns 'X' for unknown codons and an empty string
           for input strings that are not three characters long.

           It is not recommended to use this method in a production
           environment. Use method translate, instead.

 Example :
 Returns : A string
 Args    : a codon = a three nucleotide character string


=cut

sub translate_strict{
   my ($self, $value) = @_;
   my ($id) = $self->{'id'};

   $value  = lc $value;
   $value  =~ tr/u/t/;

   if (length $value != 3 ) {
       return '';
   }
   elsif (!(defined $codons->{$value}))  {
       return 'X';
   }
   else {
       return substr($tables[$id-1],$codons->{$value},1);
   }
}

=head2 revtranslate

 Title   : revtranslate
 Usage   : $obj->revtranslate('G')
 Function: returns codons for an amino acid

           Returns an empty string for unknown amino acid
           codes. Ambiquous IUPAC codes Asx,B, (Asp,D; Asn,N) and
           Glx,Z (Glu,E; Gln,Q) are resolved. Both single and three
           letter amino acid codes are accepted. '*' and 'Ter' are
           used for terminator.

           By default, the output codons are shown in DNA.  If the
           output is needed in RNA (tr/t/u/), add a second argument
           'RNA'.

 Example : $obj->revtranslate('Gly', 'RNA')
 Returns : An array of three lower case letter strings i.e. codons
 Args    : amino acid, 'RNA'

=cut

sub revtranslate {
    my ($self, $value, $coding) = @_;
    my ($id) = $self->{'id'};
    my (@aas,  $p);
    my (@codons) = ();

    if (length($value) == 3 ) {
	$value = lc $value;
	$value = ucfirst $value;
	$value = $onecode{$value};
    }
    if ( defined $value and $value =~ /[ARNDCQEGHILKMFPSTWYVBZX*]/ and length($value) == 1 ) {
	$value = uc $value;
	@aas = @ {$iupac_aa{$value}} ;	
	foreach my $aa (@aas) {
	    #print $aa, " -2\n";
	    $aa = '\*' if $aa eq '*';
	    while ($tables[$id-1] =~ m/$aa/g) {
		$p = pos $tables[$id-1];
		push (@codons, $trCol->{--$p});
	    }
	}
    }

    if ($coding and uc ($coding) eq 'RNA') {
	for my $i (0..$#codons)  {
	    $codons[$i] =~ tr/t/u/;
	}
    }

    return @codons;
}

=head2 is_start_codon

 Title   : is_start_codon
 Usage   : $obj->is_start_codon('ATG')
 Function: returns true (1) for all codons that can be used as a
           translation start, false (0) for others.
 Example : $myCodonTable->is_start_codon('ATG')
 Returns : boolean
 Args    : codon


=cut

sub is_start_codon{
   my ($self, $value) = @_;
   my ($id) = $self->{'id'};

   $value  = lc $value;
   $value  =~ tr/u/t/;

   if (length $value != 3  )  {
       return 0;
   }
   else {
       my $result = 1;
       my @ms = map { substr($starts[$id-1],$codons->{$_},1) } _unambiquous_codons($value);
       foreach my $c (@ms) {
	   $result = 0 if $c ne 'M';
       }
       return $result;
   }
}



=head2 is_ter_codon

 Title   : is_ter_codon
 Usage   : $obj->is_ter_codon('GAA')
 Function: returns true (1) for all codons that can be used as a
           translation tarminator, false (0) for others.
 Example : $myCodonTable->is_ter_codon('ATG')
 Returns : boolean
 Args    : codon


=cut

sub is_ter_codon{
   my ($self, $value) = @_;
   my ($id) = $self->{'id'};

   $value  = lc $value;
   $value  =~ tr/u/t/;

   if (length $value != 3  )  {
       return 0;
   }
   else {
       my $result = 1;
       my @ms = map { substr($tables[$id-1],$codons->{$_},1) } _unambiquous_codons($value);
       foreach my $c (@ms) {
	   $result = 0 if $c ne '*';
       }
       return $result;
   }
}

=head2 is_unknown_codon

 Title   : is_unknown_codon
 Usage   : $obj->is_unknown_codon('GAJ')
 Function: returns false (0) for all codons that are valid,
	    true (1) for others.
 Example : $myCodonTable->is_unknown_codon('NTG')
 Returns : boolean
 Args    : codon


=cut

sub is_unknown_codon{
   my ($self, $value) = @_;
   my ($id) = $self->{'id'};

   $value  = lc $value;
   $value  =~ tr/u/t/;

   if (length $value != 3  )  {
       return 1;
   }
   else {
       my $result = 0;
       my @cs = map { substr($tables[$id-1],$codons->{$_},1) } _unambiquous_codons($value);
       $result = 1 if scalar @cs == 0;
       return $result;
   }
}

=head2 _unambiquous_codons

 Title   : _unambiquous_codons
 Usage   : @codons = _unambiquous_codons('ACN')
 Function:
 Example :
 Returns : array of strings (one letter unambiguous amino acid codes)
 Args    : a codon = a three IUPAC nucleotide character string

=cut

sub _unambiquous_codons{
    my ($value) = @_;
    my @nts = ();
    my @codons = ();
    my ($i, $j, $k);
    @nts = map { $iupac_dna{$_} }  split(//, $value);
    for my $i (@{$nts[0]}) {
	for my $j (@{$nts[1]}) {
	    for my $k (@{$nts[2]}) {
		push @codons, "$i$j$k";
	    }
	}
    }
    return @codons;
}


1;
