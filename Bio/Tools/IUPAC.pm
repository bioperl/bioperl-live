
#
# BioPerl module for IUPAC
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

IUPAC - Generates unique Seq objects from an ambiguous Seq object

=head1 SYNOPSIS

use Bio::Seq;
use Bio::Tools::IUPAC;

my $ambiseq = new Bio::Seq (-seq => 'ARTCGUTGR', -type => 'Dna');
my $stream  = new Bio::Tools::IUPAC($ambiseq);

while ($uniqueseq = $stream->next_seq()) {
    # process the unique Seq object.
}

=head1 DESCRIPTION

IUPAC is a tool that produces a stream of unique, "strict"-satisfying Seq
objects from an ambiquous Seq object (containing non-standard characters given
the meaning shown below)

        Extended Dna / Rna alphabet :
        (includes symbols for nucleotide ambiguity)
        ------------------------------------------
        Symbol       Meaning      Nucleic Acid
        ------------------------------------------
         A            A           Adenine
         C            C           Cytosine
         G            G           Guanine
         T            T           Thymine
         U            U           Uracil
         M          A or C
         R          A or G
         W          A or T
         S          C or G
         Y          C or T
         K          G or T
         V        A or C or G
         H        A or C or T
         D        A or G or T
         B        C or G or T
         X      G or A or T or C
         N      G or A or T or C

        IUPAC-IUB SYMBOLS FOR NUCLEOTIDE NOMENCLATURE:
          Cornish-Bowden (1985) Nucl. Acids Res. 13: 3021-3030.

-----------------------------------

       Amino Acid alphabet:
        ------------------------------------------
        Symbol           Meaning
        ------------------------------------------
        A        Alanine
        B        Aspartic Acid, Asparagine
        C        Cystine
        D        Aspartic Acid
        E        Glutamic Acid
        F        Phenylalanine
        G        Glycine
        H        Histidine
        I        Isoleucine
        K        Lysine
        L        Leucine
        M        Methionine
        N        Asparagine
        P        Proline
        Q        Glutamine
        R        Arginine
        S        Serine
        T        Threonine
        V        Valine
        W        Tryptophan
        X        Unknown
        Y        Tyrosine
        Z        Glutamic Acid, Glutamine
        *        Terminator


        IUPAC-IUP AMINO ACID SYMBOLS:
          Biochem J. 1984 Apr 15; 219(2): 345-373
          Eur J Biochem. 1993 Apr 1; 213(1): 2

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Tools::IUPAC;

use strict;
use vars qw(@ISA @EXPORT_OK $AUTOLOAD);

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;
use Exporter;

@ISA = qw(Bio::Root::Object Exporter);
@EXPORT_OK = qw();

=head2 new

 Title   : new
 Usage   : new Bio::Tools::IUPAC $seq;
 Function: returns a new seq stream (akin to SeqIO)
 Returns : a Bio::Tools::IUPAC stream object that will produce unique
           Seq objects on demand.
 Args    : an ambiguously coded Seq.pm object that has a specified 'type'


=cut

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
    my($self,@args) = @_;

    my $make = $self->SUPER::_initialize;

    # set stuff in self from @args

    my $seq = shift @args;
    $seq->isa('Bio::Seq') or $self->throw("Must supply a Seq.pm object to IUPAC!");
    $self->{'_SeqObj'} = $seq;
    if ($self->{'_SeqObj'}->type() =~ m/^[dr]na$/i ) { # nucleotide seq object
	$self->{'_iupac_iub'} = { A => [qw(A)],
				  C => [qw(C)],
				  G => [qw(G)],
				  T => [qw(T)],
				  U => [qw(U)],
				  M => [qw(A C)],
				  R => [qw(A G)],
				  W => [qw(A T)],
				  S => [qw(C G)],
				  Y => [qw(C T)],
				  K => [qw(G T)],
				  V => [qw(A C G)],
				  H => [qw(A C T)],
				  D => [qw(A G T)],
				  B => [qw(C G T)],
				  X => [qw(G A T C)],
				  N => [qw(G A T C)]
				};
	$self->{'_alpha'} = [ map { $self->{'_iupac_iub'}{uc($_)} } split('', $self->{'_SeqObj'}->str()) ];
    } elsif ($self->{'_SeqObj'}->type() =~ m/^amino$/i ) { # amino acid seq object
	$self->{'_iupac_iup'} = { A => [qw(A)],
				  B => [qw(D N)],
				  C => [qw(C)],
				  D => [qw(D)],
				  E => [qw(E)],
				  F => [qw(F)],
				  G => [qw(G)],
				  H => [qw(H)],
				  I => [qw(I)],
				  K => [qw(K)],
				  L => [qw(L)],
				  M => [qw(M)],
				  N => [qw(N)],
				  P => [qw(P)],
				  Q => [qw(Q)],
				  R => [qw(R)],
				  S => [qw(S)],
				  T => [qw(T)],
				  U => [qw(U)],
				  V => [qw(V)],
				  W => [qw(W)],
				  X => [qw(X)],
				  Y => [qw(Y)],
				  Z => [qw(E Q)],
				  '*' => ['*']
				};
	$self->{'_alpha'} = [ map { $self->{'_iupac_iup'}{uc($_)} } split('', $self->{'_SeqObj'}->str()) ];
    } else { # unknown type: we could make a guess, but let's not.
	$self->throw("You must specify the 'type' of sequence provided to IUPAC");
    }
    $self->{'_string'} = [(0) x length($self->{'_SeqObj'}->str())];
    scalar @{$self->{'_string'}} or $self->throw("Sequence has zero-length!");
    $self->{'_string'}->[0] = -1;
    return $make; # success - we hope!
}

=head2 next_seq

 Title   : next_seq
 Usage   : $iupac->next_seq()
 Function: returns the next unique Seq object
 Returns : a Seq.pm object
 Args    : none.


=cut

sub next_seq{
    my ($self,@args) = @_;

    for my $i ( 0 .. $#{$self->{'_string'}} ) {
	next unless $self->{'_string'}->[$i] || @{$self->{'_alpha'}->[$i]} > 1;
	if ( $self->{'_string'}->[$i] == $#{$self->{'_alpha'}->[$i]} ) { # rollover
	    if ( $i == $#{$self->{'_string'}} ) { # end of possibilities
		return undef;
	    } else {
		$self->{'_string'}->[$i] = 0;
		next;
	    }
	} else {
	    $self->{'_string'}->[$i]++;
	    my $j = -1;
	    $self->{'_SeqObj'}->setseq(join('', map { $j++; $self->{'_alpha'}->[$j]->[$_]; } @{$self->{'_string'}}));
	    my $desc = $self->{'_SeqObj'}->desc();
	    $self->{'_num'}++;
	    1 while $self->{'_num'} =~ s/(\d)(\d\d\d)(?!\d)/$1,$2/;
	    $desc =~ s/( \[Bio::Tools::IUPAC-generated\sunique sequence # [^\]]*\])|$/ \[Bio::Tools::IUPAC-generated unique sequence # $self->{'_num'}\]/;
	    $self->{'_SeqObj'}->desc($desc);
	    $self->{'_num'} =~ s/,//g;
	    return $self->{'_SeqObj'};
	}
    }
}

sub AUTOLOAD {

    my $self = shift @_;
    my $method = $AUTOLOAD;
    $method =~ s/.*:://;
    return $self->{'_SeqObj'}->$method(@_);
}

1;

