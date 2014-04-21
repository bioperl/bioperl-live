#
# bioperl module for Bio::Tools::CodonTable
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

Bio::Tools::CodonTable - Codon table object

=head1 SYNOPSIS

  # This is a read-only class for all known codon tables.  The IDs are
  # the ones used by nucleotide sequence databases.  All common IUPAC
  # ambiguity codes for DNA, RNA and amino acids are recognized.

  use Bio::Tools::CodonTable;

  # defaults to ID 1 "Standard"
  $myCodonTable   = Bio::Tools::CodonTable->new();
  $myCodonTable2  = Bio::Tools::CodonTable->new( -id => 3 );

  # change codon table
  $myCodonTable->id(5);

  # examine codon table
  print  join (' ', "The name of the codon table no.", $myCodonTable->id(4),
           "is:", $myCodonTable->name(), "\n");

  # print possible codon tables
  $tables = Bio::Tools::CodonTable->tables;
  while ( ($id,$name) = each %{$tables} ) {
    print "$id = $name\n";
  }

  # translate a codon
  $aa = $myCodonTable->translate('ACU');
  $aa = $myCodonTable->translate('act');
  $aa = $myCodonTable->translate('ytr');

  # reverse translate an amino acid
  @codons = $myCodonTable->revtranslate('A');
  @codons = $myCodonTable->revtranslate('Ser');
  @codons = $myCodonTable->revtranslate('Glx');
  @codons = $myCodonTable->revtranslate('cYS', 'rna');

  # reverse translate an entire amino acid sequence into a IUPAC
  # nucleotide string

  my $seqobj    = Bio::PrimarySeq->new(-seq => 'FHGERHEL');
  my $iupac_str = $myCodonTable->reverse_translate_all($seqobj);

  # boolean tests
  print "Is a start\n"       if $myCodonTable->is_start_codon('ATG');
  print "Is a terminator\n" if $myCodonTable->is_ter_codon('tar');
  print "Is a unknown\n"     if $myCodonTable->is_unknown_codon('JTG');

=head1 DESCRIPTION

Codon tables are also called translation tables or genetic codes
since that is what they represent. A bit more complete picture
of the full complexity of codon usage in various taxonomic groups
is presented at the NCBI Genetic Codes Home page.

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

Note: 

This class deals primarily with individual codons and amino
acids. However in the interest of speed you can L<translate>
longer sequence, too. The full complexity of protein translation
is tackled by L<Bio::PrimarySeqI::translate>.


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
          O           Pyl            Pyrrolysine (22nd amino acid)
          U           Sec            Selenocysteine (21st amino acid)
          S           Ser            Serine
          T           Thr            Threonine
          W           Trp            Tryptophan
          Y           Tyr            Tyrosine
          V           Val            Valine
          B           Asx            Aspartic acid or Asparagine
          Z           Glx            Glutamine or Glutamic acid
          J           Xle            Isoleucine or Valine (mass spec ambiguity)
          X           Xaa            Any or unknown amino acid


It is worth noting that, "Bacterial" codon table no. 11 produces an
polypeptide that is, confusingly, identical to the standard one. The
only differences are in available initiator codons.


NCBI Genetic Codes home page:
     http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c

EBI Translation Table Viewer:
     http://www.ebi.ac.uk/cgi-bin/mutations/trtables.cgi

Amended ASN.1 version with ids 16 and 21 is at:
     ftp://ftp.ebi.ac.uk/pub/databases/geneticcode/

Thanks to Matteo diTomasso for the original Perl implementation
of these tables.

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

package Bio::Tools::CodonTable;
use vars qw(@NAMES @TABLES @STARTS $TRCOL $CODONS %IUPAC_DNA $CODONGAP $GAP
                %IUPAC_AA %THREELETTERSYMBOLS $VALID_PROTEIN $TERMINATOR);
use strict;

# Object preamble - inherits from Bio::Root::Root
use Bio::Tools::IUPAC;
use Bio::SeqUtils;

use base qw(Bio::Root::Root);


# first set internal values for all translation tables

BEGIN { 
    use constant CODONSIZE => 3;
    $GAP = '-';
    $CODONGAP = $GAP x CODONSIZE;

    @NAMES =            #id
    (
     'Standard',        #1
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
     'Thraustochytrium Mitochondrial', #23
     'Strict', #24, option for only ATG start
     );

    @TABLES =
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
       FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
       );

   #           (bases used for these tables, for reference)
   # 1 TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
   # 2 TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
   # 3 TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

    @STARTS =
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
       -----------------------------------M----------------------------
       );

    my @nucs = qw(t c a g);
    my $x = 0;
    ($CODONS, $TRCOL) = ({}, {});
    for my $i (@nucs) {
    for my $j (@nucs) {
        for my $k (@nucs) {
        my $codon = "$i$j$k";
        $CODONS->{$codon} = $x;
        $TRCOL->{$x} = $codon;
        $x++;
        }
    }
    }
    %IUPAC_DNA = Bio::Tools::IUPAC->iupac_iub();
    %IUPAC_AA = Bio::Tools::IUPAC->iupac_iup();
    %THREELETTERSYMBOLS = Bio::SeqUtils->valid_aa(2);
    $VALID_PROTEIN = '['.join('',Bio::SeqUtils->valid_aa(0)).']';
    $TERMINATOR = '*';
}

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
 Function: Sets or returns the id of the translation table.  IDs are
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
       if (  !(defined $TABLES[$value-1]) or $TABLES[$value-1] eq '') {
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
   return $NAMES[$id-1];
}

=head2 tables

 Title   : tables
 Usage   : $obj->tables()  or  Bio::Tools::CodonTable->tables()
 Function: returns a hash reference where each key is a valid codon
           table id() number, and each value is the corresponding
           codon table name() string
 Example :
 Returns : A hashref
 Args    : None


=cut

sub tables{
  my %tables;
  for my $id (1 .. @NAMES) {
    my $name = $NAMES[$id-1];
    $tables{$id} = $name if $name;
  }
  return \%tables;
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
    my ($self, $seq, $complete_codon) = @_;
    $self->throw("Calling translate without a seq argument!") unless defined $seq;
    return '' unless $seq;

    my $id = $self->id;
    my ($partial) = 0;
    $partial = 2 if length($seq) % CODONSIZE == 2;
    
    $seq = lc $seq;
    $seq =~ tr/u/t/;
    my $protein = "";
    if ($seq =~ /[^actg]/ ) { #ambiguous chars
        for (my $i = 0; $i < (length($seq) - (CODONSIZE-1)); $i+= CODONSIZE) {
            my $triplet = substr($seq, $i, CODONSIZE);
        if( $triplet eq $CODONGAP ) {
        $protein .= $GAP;
        } elsif (exists $CODONS->{$triplet}) {
        $protein .= substr($TABLES[$id-1], 
                   $CODONS->{$triplet},1);
        } else {
        $protein .= $self->_translate_ambiguous_codon($triplet);
        }
    }
    } else { # simple, strict translation
    for (my $i = 0; $i < (length($seq) - (CODONSIZE -1)); $i+=CODONSIZE) {
            my $triplet = substr($seq, $i, CODONSIZE); 
            if( $triplet eq $CODONGAP ) {
        $protein .= $GAP;
        } if (exists $CODONS->{$triplet}) {
                $protein .= substr($TABLES[$id-1], $CODONS->{$triplet}, 1);
        } else {
                $protein .= 'X';
            }
        }
    }
    if ($partial == 2 && $complete_codon) { # 2 overhanging nucleotides
    my $triplet = substr($seq, ($partial -4)). "n";
    if( $triplet eq $CODONGAP ) {
        $protein .= $GAP;
    } elsif (exists $CODONS->{$triplet}) {
        my $aa = substr($TABLES[$id-1], $CODONS->{$triplet},1);       
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
    my @codons = $self->unambiguous_codons($triplet);
    my %aas =();
    foreach my $codon (@codons) {
    $aas{substr($TABLES[$id-1],$CODONS->{$codon},1)} = 1;
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
   my $id = $self->{'id'};

   $value  = lc $value;
   $value  =~ tr/u/t/;

   return '' unless length $value == 3;

   return 'X' unless defined $CODONS->{$value};

   return substr( $TABLES[$id-1], $CODONS->{$value}, 1 );
}

=head2 revtranslate

 Title   : revtranslate
 Usage   : $obj->revtranslate('G')
 Function: returns codons for an amino acid

           Returns an empty string for unknown amino acid
           codes. Ambiguous IUPAC codes Asx,B, (Asp,D; Asn,N) and
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
    my @codons;

    if (length($value) == 3 ) {
        $value = lc $value;
        $value = ucfirst $value;
        $value = $THREELETTERSYMBOLS{$value};
    }
    if ( defined $value and $value =~ /$VALID_PROTEIN/
          and length($value) == 1 ) {
        my $id = $self->{'id'};

        $value = uc $value;
        my @aas = @{$IUPAC_AA{$value}};
        foreach my $aa (@aas) {
            #print $aa, " -2\n";
            $aa = '\*' if $aa eq '*';
          while ($TABLES[$id-1] =~ m/$aa/g) {
              my $p = pos $TABLES[$id-1];
              push (@codons, $TRCOL->{--$p});
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

=head2 reverse_translate_all

 Title   : reverse_translate_all
 Usage   : my $iup_str = $cttable->reverse_translate_all($seq_object)
           my $iup_str = $cttable->reverse_translate_all($seq_object,
                                                         $cutable,
                                                         15);
 Function: reverse translates a protein sequence into IUPAC nucleotide
           sequence. An 'X' in the protein sequence is converted to 'NNN'
           in the nucleotide sequence.
 Returns : a string
 Args    : a Bio::PrimarySeqI compatible object (mandatory)
           a Bio::CodonUsage::Table object and a threshold if only
             codons with a relative frequency above the threshold are
             to be considered.
=cut

sub reverse_translate_all {

    my ($self, $obj, $cut, $threshold) = @_;

    ## check args are OK

    if (!$obj || !$obj->isa('Bio::PrimarySeqI')){
        $self->throw(" I need a Bio::PrimarySeqI object, not a [".
                        ref($obj) . "]");
        }
    if($obj->alphabet ne 'protein') {
        $self->throw("Cannot reverse translate, need an amino acid sequence .".
                     "This sequence is of type [" . $obj->alphabet ."]");
        }
    my @data;
    my @seq = split '', $obj->seq;

    ## if we're not supplying a codon usage table...
    if( !$cut && !$threshold) {
        ## get lists of possible codons for each aa. 
        for my $aa (@seq) {
            if ($aa =~ /x/i) {
                push @data, (['NNN']);
            }else {
                my @cods = $self->revtranslate($aa);
                push @data, \@cods;
            }
        }
    }else{
    #else we are supplying a codon usage table, we just want common codons
    #check args first. 
        if(!$cut->isa('Bio::CodonUsage::Table'))    {
            $self->throw("I need a Bio::CodonUsage::Table object, not a [".
                     ref($cut). "].");
            }
        my $cod_ref = $cut->probable_codons($threshold);
        for my $aa (@seq) {
            if ($aa =~ /x/i) {
                push @data, (['NNN']);
                next;
                }
            push @data, $cod_ref->{$aa};
        }
    }

    return $self->_make_iupac_string(\@data);

}

=head2 reverse_translate_best

 Title   : reverse_translate_best
 Usage   : my $str = $cttable->reverse_translate_best($seq_object,$cutable);
 Function: Reverse translates a protein sequence into plain nucleotide
           sequence (GATC), uses the most common codon for each amino acid
 Returns : A string
 Args    : A Bio::PrimarySeqI compatible object and a Bio::CodonUsage::Table object

=cut

sub reverse_translate_best {

    my ($self, $obj, $cut) = @_;

    if (!$obj || !$obj->isa('Bio::PrimarySeqI')){
        $self->throw(" I need a Bio::PrimarySeqI object, not a [".
                         ref($obj) . "]");
    }
    if ($obj->alphabet ne 'protein')    {
        $self->throw("Cannot reverse translate, need an amino acid sequence .".
                         "This sequence is of type [" . $obj->alphabet ."]");
    }
    if ( !$cut | !$cut->isa('Bio::CodonUsage::Table'))  {
        $self->throw("I need a Bio::CodonUsage::Table object, not a [".
                         ref($cut). "].");
    }

    my $str = '';
    my @seq = split '', $obj->seq;

    my $cod_ref = $cut->most_common_codons();

    for my $aa ( @seq ) {
        if ($aa =~ /x/i) {
            $str .= 'NNN';
            next;
        }
        if ( defined $cod_ref->{$aa} ) {
            $str .= $cod_ref->{$aa};
        } else {
            $self->throw("Input sequence contains invalid character: $aa");
        }
    }
   $str;
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
   shift->_codon_is( shift, \@STARTS, 'M' );
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
    shift->_codon_is( shift, \@TABLES, $TERMINATOR );
}

# desc: compares the passed value with a single entry in the given
#       codon table
# args: a value (typically a three-char string like 'atg'),
#       a reference to the appropriate set of codon tables,
#       a single-character value to check for at the position in the
#       given codon table
# ret:  boolean, true if the given codon table contains the $key at the
#       position corresponding to $value
sub _codon_is {
   my ($self, $value, $table, $key ) = @_;

   return 0 unless length $value == 3;

   $value  = lc $value;
   $value  =~ tr/u/t/;

   my $id = $self->{'id'};
   for my $c ( $self->unambiguous_codons($value) ) {
       my $m = substr( $table->[$id-1], $CODONS->{$c}, 1 );
       return 0 unless $m eq $key;
   }
   return 1;
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
   $value  = lc $value;
   $value  =~ tr/u/t/;
   return 1 unless $self->unambiguous_codons($value);
   return 0;
}

=head2 unambiguous_codons

 Title   : unambiguous_codons
 Usage   : @codons = $self->unambiguous_codons('ACN')
 Returns : array of strings (one-letter unambiguous amino acid codes)
 Args    : a codon = a three IUPAC nucleotide character string

=cut

sub unambiguous_codons{
    my ($self,$value) = @_;
    my @nts = map { $IUPAC_DNA{uc $_} }  split(//, $value);

    my @codons;
    for my $i ( @{$nts[0]} ) {
    for my $j ( @{$nts[1]} ) {
    for my $k ( @{$nts[2]} ) {
        push @codons, lc "$i$j$k";
    }}}
    return @codons;
}

=head2 _unambiquous_codons

deprecated, now an alias for unambiguous_codons

=cut

sub _unambiquous_codons {
    unambiguous_codons( undef, @_ );
}

=head2 add_table

 Title   : add_table
 Usage   : $newid = $ct->add_table($name, $table, $starts)
 Function: Add a custom Codon Table into the object.
           Know what you are doing, only the length of
           the argument strings is checked!
 Returns : the id of the new codon table
 Args    : name, a string, optional (can be empty)
           table, a string of 64 characters
           startcodons, a string of 64 characters, defaults to standard

=cut

sub add_table {
    my ($self, $name, $table, $starts) = @_;

    $name ||= 'Custom'. scalar @NAMES + 1;
    $starts ||= $STARTS[0]; 
    $self->throw('Suspect input!')
        unless length($table) == 64 and length($starts) == 64;

    push @NAMES, $name;
    push @TABLES, $table;
    push @STARTS, $starts;

    return scalar @NAMES;

}

sub _make_iupac_string {

    my ($self, $cod_ref) = @_;
    if(ref($cod_ref) ne 'ARRAY') {
        $self->throw(" I need a reference to a list of references to codons, ".
                     " not a [". ref($cod_ref) . "].");
        }
    my %iupac_hash   = Bio::Tools::IUPAC->iupac_rev_iub();
    my $iupac_string = ''; ## the string to be returned
    for my $aa (@$cod_ref) {

        ## scan through codon positions, record the differing values,   
        # then look up in the iub hash
        for my $index(0..2) {
            my %h;
            map { my $k = substr($_,$index,1);
                $h{$k}  = undef;} @$aa;
            my $lookup_key = join '', sort{$a cmp $b}keys %h;

            ## extend string 
            $iupac_string .= $iupac_hash{uc$lookup_key};
        }
    }
    return $iupac_string;

}


1;
