package Bio::Tools::CodonTable;

use utf8;
use strict;
use warnings;

use Bio::Tools::IUPAC;
use Bio::SeqUtils;

use base qw(Bio::Root::Root);

# ABSTRACT: Codon table object
# AUTHOR: Heikki Lehvaslaiho <heikki@bioperl.org>
# OWNER: Heikki Lehvaslaiho <heikki@bioperl.org>
# LICENSE: Perl_5

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
  print "Is a terminator\n"  if $myCodonTable->is_ter_codon('tar');
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
     (Last update of the Genetic Codes: Apr. 25, 2024)
     https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c

The "value notation" / "print form" ASN.1 version is at:
     ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt

Thanks to Matteo diTomasso for the original Perl implementation
of these tables.

=cut


# set internal values for all translation tables
use constant CODONSIZE => 3;
our $GAP = '-';
our $CODONGAP = $GAP x CODONSIZE;
our %IUPAC_DNA = Bio::Tools::IUPAC->iupac_iub();
our %IUPAC_AA = Bio::Tools::IUPAC->iupac_iup();
our %THREELETTERSYMBOLS = Bio::SeqUtils->valid_aa(2);
our $VALID_PROTEIN = '['.join('',Bio::SeqUtils->valid_aa(0)).']';
our $TERMINATOR = '*';

our (@NAMES, @TABLES, @STARTS);
# Parse the ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt file which
# is below __DATA__ in this module (see the end of the file).  This
# fills the @NAMES, @TABLES, and @STARTS variables.  To update to a
# new release of gc.prt, replace the content below __DATA__.
{
    # Init tables has with special option (id=0) for ATG-only start
    my %tables = (
        0 => {
            name => "Strict",
            ncbieaa => "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            sncbieaa => "----------**--*--------------------M----------------------------",
        },
    );

    while (defined(my $line = <DATA>)) {
        next if $line =~ /^\s*--/;  # skip comment lines
        if ($line =~ /^\s*\{\s*$/) {  # start of a table description
            my $name = "";
            my $id = 0;
            my $ncbieaa = "";
            my $sncbieaa = "";
            do {
                if ($line =~ /^\s*(name|id|ncbieaa|sncbieaa)\s+(.+)/) {
                    my $key = $1;
                    my $rem = $2;
                    if ($key eq "id") {
                        $rem =~ /^(\d+)/;
                        $id = int $1;
                    } else {
                        # The remaining keys --- name, ncbieaa, and
                        # sncbieaa --- are strings which may be
                        # multi-line (e.g., name for table with id 4).
                        # We are assuming that there is no " character
                        # inside the value so we keep appending lines
                        # until we find an end ".
                        while ($rem !~ /^"(.*)"/ && ! eof DATA) {
                            $rem .= <DATA>;
                        }
                        $rem =~ s/\n//g;
                        $rem =~ /^"(.*)"/;
                        my $str = $1;
                        if ($key eq "name" && ! $name) {
                            # ignore alternative names, e.g. SGC0,
                            # only keep the first name listed.
                            $name = $str;
                        } elsif ($key eq "ncbieaa") {
                            $ncbieaa = $str;
                        } elsif ($key eq "sncbieaa") {
                            $sncbieaa = $str;
                        }
                    }
                }
            } until (($line = <DATA>) =~ /^\s*}\s*,?$/);  # we reached the end of table description
            $tables{$id} = {
                name => $name,
                ncbieaa => $ncbieaa,
                sncbieaa => $sncbieaa
            };
        }
    }
    close DATA;
    # use Data::Dumper;
    # print Dumper %tables;

    # After parsing gc.prt, fill in @NAMES, @TABLES, and @STARTS
    my $highest_id = (sort {$a <=> $b} keys %tables)[-1];
    for (my $i = 0; $i < $highest_id; $i++) {
        if (defined $tables{$i}) {
            push @NAMES, $tables{$i}->{name};
            push @TABLES, $tables{$i}->{ncbieaa};
            push @STARTS, $tables{$i}->{sncbieaa};
        } else {
            push @NAMES, '';
            push @TABLES, '';
            push @STARTS, '';
        }
    }
}

our ($TRCOL, $CODONS);
{
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
}

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my($id) =
        $self->_rearrange([qw(ID
                 )],
             @args);

    $id = 1 if ( ! defined ( $id ) );
    $self->id($id);
    return $self; # success - we hope!
}

=head2 id

 Title   : id
 Usage   : $obj->id(3); $id_integer = $obj->id();
 Function: Sets or returns the id of the translation table.  IDs are
           integers from 0 (special ATG-only start) to 25, excluding
           7-8 and 17-20 which have been removed. If an invalid ID is
           given the method returns 1, the standard table.
 Example :
 Returns : value of id, a scalar, warn and fall back to 1 (standard table)
           if specified id is not valid
 Args    : newvalue (optional)

=cut

sub id{
    my ($self,$value) = @_;
    if( defined $value) {
        if (! defined $TABLES[$value] || $TABLES[$value] eq '' || $value < 0) {
            $self->warn("Not a valid codon table ID [$value], using [1] instead ");
            $value = 1;
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
   return $NAMES[$id];
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
  for my $id (0 .. $#NAMES) {
    my $name = $NAMES[$id];
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
                $protein .= substr($TABLES[$id],
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
            }
            if (exists $CODONS->{$triplet}) {
                $protein .= substr($TABLES[$id], $CODONS->{$triplet}, 1);
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
            my $aa = substr($TABLES[$id], $CODONS->{$triplet},1);
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
        $aas{substr($TABLES[$id],$CODONS->{$codon},1)} = 1;
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

   return substr( $TABLES[$id], $CODONS->{$value}, 1 );
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
    if (    defined $value and $value =~ /$VALID_PROTEIN/
        and length($value) == 1
        ) {
        my $id = $self->{'id'};

        $value = uc $value;
        my @aas = @{$IUPAC_AA{$value}};
        foreach my $aa (@aas) {
            #print $aa, " -2\n";
            $aa = '\*' if $aa eq '*';
            while ($TABLES[$id] =~ m/$aa/g) {
                my $p = pos $TABLES[$id];
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
   return $str;
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
   my ($self, $value) = @_;
   my $id = $self->{'id'};

   # We need to ensure U is mapped to T (ie. UAG)
   $value = uc $value;
   $value =~ tr/U/T/;

   if (length $value != 3  )  {
       # Incomplete codons are not stop codons
       return 0;
   } else {
       my $result = 0;

       # For all the possible codons, if any are not a stop
       # codon, fail immediately
       for my $c ( $self->unambiguous_codons($value) ) {
	   my $m = substr( $TABLES[$id], $CODONS->{$c}, 1 );
	   if($m eq $TERMINATOR) {
	       $result = 1;
	   } else {
	       return 0;
	   }
       }
       return $result;
   }
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
       my $m = substr( $table->[$id], $CODONS->{$c}, 1 );
       if ($m eq $key) { return 1; }
   }
   return 0;
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

    $name   ||= 'Custom' . $#NAMES + 1;
    $starts ||= $STARTS[1];
    $self->throw('Suspect input!')
        unless length($table) == 64 and length($starts) == 64;

    push @NAMES,  $name;
    push @TABLES, $table;
    push @STARTS, $starts;

    return $#NAMES;
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

# Follows the content of
# ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt, which is the NCBI
# genetic codon table in ASN.1 value notation / print format.  We do
# not have a ASN.1 decoder for value notation but it's easy enough to
# parse.

__DATA__
--**************************************************************************
--  This is the NCBI genetic code table
--  Initial base data set from Andrzej Elzanowski while at PIR International
--  Addition of Eubacterial and Alternative Yeast by J.Ostell at NCBI
--  Base 1-3 of each codon have been added as comments to facilitate
--    readability at the suggestion of Peter Rice, EMBL
--  Later additions by Taxonomy Group staff at NCBI
--
--  Version 4.6
--     Renamed genetic code 24 to Rhabdopleuridae Mitochondrial
--
--  Version 4.5
--     Added Cephalodiscidae mitochondrial genetic code 33
--
--  Version 4.4
--     Added GTG as start codon for genetic code 3
--     Added Balanophoraceae plastid genetic code 32
--
--  Version 4.3
--     Change to CTG -> Leu in genetic codes 27, 28, 29, 30
--
--  Version 4.2
--     Added Karyorelict nuclear genetic code 27
--     Added Condylostoma nuclear genetic code 28
--     Added Mesodinium nuclear genetic code 29
--     Added Peritrich nuclear genetic code 30
--     Added Blastocrithidia nuclear genetic code 31
--
--  Version 4.1
--     Added Pachysolen tannophilus nuclear genetic code 26
--
--  Version 4.0
--     Updated version to reflect numerous undocumented changes:
--     Corrected start codons for genetic code 25
--     Name of new genetic code is Candidate Division SR1 and Gracilibacteria
--     Added candidate division SR1 nuclear genetic code 25
--     Added GTG as start codon for genetic code 24
--     Corrected Pterobranchia Mitochondrial genetic code (24)
--     Added genetic code 24, Pterobranchia Mitochondrial
--     Genetic code 11 is now Bacterial, Archaeal and Plant Plastid
--     Fixed capitalization of mitochondrial in codes 22 and 23
--     Added GTG, ATA, and TTG as alternative start codons to code 13
--
--  Version 3.9
--     Code 14 differs from code 9 only by translating UAA to Tyr rather than
--     STOP.  A recent study (Telford et al, 2000) has found no evidence that
--     the codon UAA codes for Tyr in the flatworms, but other opinions exist.
--     There are very few GenBank records that are translated with code 14,
--     but a test translation shows that retranslating these records with code
--     9 can cause premature terminations.  Therefore, GenBank will maintain
--     code 14 until further information becomes available.
--
--  Version 3.8
--     Added GTG start to Echinoderm mitochondrial code, code 9
--
--  Version 3.7
--     Added code 23 Thraustochytrium mitochondrial code
--        formerly OGMP code 93
--        submitted by Gertraude Berger, Ph.D.
--
--  Version 3.6
--     Added code 22 TAG-Leu, TCA-stop
--        found in mitochondrial DNA of Scenedesmus obliquus
--        submitted by Gertraude Berger, Ph.D.
--        Organelle Genome Megasequencing Program, Univ Montreal
--
--  Version 3.5
--     Added code 21, Trematode Mitochondrial
--       (as deduced from: Garey & Wolstenholme,1989; Ohama et al, 1990)
--     Added code 16, Chlorophycean Mitochondrial
--       (TAG can translated to Leucine instaed to STOP in chlorophyceans
--        and fungi)
--
--  Version 3.4
--     Added CTG,TTG as allowed alternate start codons in Standard code.
--        Prats et al. 1989, Hann et al. 1992
--
--  Version 3.3 - 10/13/95
--     Added alternate intiation codon ATC to code 5
--        based on complete mitochondrial genome of honeybee
--        Crozier and Crozier (1993)
--
--  Version 3.2 - 6/24/95
--  Code       Comments
--   10        Alternative Ciliate Macronuclear renamed to Euplotid Macro...
--   15        Blepharisma Macro.. code added
--    5        Invertebrate Mito.. GTG allowed as alternate initiator
--   11        Eubacterial renamed to Bacterial as most alternate starts
--               have been found in Archea
--
--
--  Version 3.1 - 1995
--  Updated as per Andrzej Elzanowski at NCBI
--     Complete documentation in NCBI toolkit documentation
--  Note: 2 genetic codes have been deleted
--
--   Old id   Use id     - Notes
--
--   id 7      id 4      - Kinetoplast code now merged in code id 4
--   id 8      id 1      - all plant chloroplast differences due to RNA edit
--
--
--*************************************************************************

Genetic-code-table ::= {
 {
  name "Standard" ,
  name "SGC0" ,
  id 1 ,
  ncbieaa  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "---M------**--*----M---------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Vertebrate Mitochondrial" ,
  name "SGC1" ,
  id 2 ,
  ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
  sncbieaa "----------**--------------------MMMM----------**---M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Yeast Mitochondrial" ,
  name "SGC2" ,
  id 3 ,
  ncbieaa  "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "----------**----------------------MM---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
    name "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate
 Mitochondrial; Mycoplasma; Spiroplasma" ,
  name "SGC3" ,
  id 4 ,
  ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "--MM------**-------M------------MMMM---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Invertebrate Mitochondrial" ,
  name "SGC4" ,
  id 5 ,
  ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
  sncbieaa "---M------**--------------------MMMM---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear" ,
  name "SGC5" ,
  id 6 ,
  ncbieaa  "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "--------------*--------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Echinoderm Mitochondrial; Flatworm Mitochondrial" ,
  name "SGC8" ,
  id 9 ,
  ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
  sncbieaa "----------**-----------------------M---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Euplotid Nuclear" ,
  name "SGC9" ,
  id 10 ,
  ncbieaa  "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "----------**-----------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Bacterial, Archaeal and Plant Plastid" ,
  id 11 ,
  ncbieaa  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "---M------**--*----M------------MMMM---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Alternative Yeast Nuclear" ,
  id 12 ,
  ncbieaa  "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "----------**--*----M---------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Ascidian Mitochondrial" ,
  id 13 ,
  ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
  sncbieaa "---M------**----------------------MM---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Alternative Flatworm Mitochondrial" ,
  id 14 ,
  ncbieaa  "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
  sncbieaa "-----------*-----------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Blepharisma Macronuclear" ,
  id 15 ,
  ncbieaa  "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "----------*---*--------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Chlorophycean Mitochondrial" ,
  id 16 ,
  ncbieaa  "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "----------*---*--------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Trematode Mitochondrial" ,
  id 21 ,
  ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
  sncbieaa "----------**-----------------------M---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Scenedesmus obliquus Mitochondrial" ,
  id 22 ,
  ncbieaa  "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "------*---*---*--------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Thraustochytrium Mitochondrial" ,
  id 23 ,
  ncbieaa  "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "--*-------**--*-----------------M--M---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Rhabdopleuridae Mitochondrial" ,
  id 24 ,
  ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
  sncbieaa "---M------**-------M---------------M---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Candidate Division SR1 and Gracilibacteria" ,
  id 25 ,
  ncbieaa  "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "---M------**-----------------------M---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Pachysolen tannophilus Nuclear" ,
  id 26 ,
  ncbieaa  "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "----------**--*----M---------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Karyorelict Nuclear" ,
  id 27 ,
  ncbieaa  "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "--------------*--------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Condylostoma Nuclear" ,
  id 28 ,
  ncbieaa  "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "----------**--*--------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Mesodinium Nuclear" ,
  id 29 ,
  ncbieaa  "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "--------------*--------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Peritrich Nuclear" ,
  id 30 ,
  ncbieaa  "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "--------------*--------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Blastocrithidia Nuclear" ,
  id 31 ,
  ncbieaa  "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "----------**-----------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Balanophoraceae Plastid" ,
  id 32 ,
  ncbieaa  "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "---M------*---*----M------------MMMM---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Cephalodiscidae Mitochondrial" ,
  id 33 ,
  ncbieaa  "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
  sncbieaa "---M-------*-------M---------------M---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 }
}
