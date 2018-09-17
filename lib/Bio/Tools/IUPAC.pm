#
# BioPerl module for IUPAC
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::IUPAC - Generates unique sequence objects or regular expressions from
an ambiguous IUPAC sequence

=head1 SYNOPSIS

 use Bio::PrimarySeq;
 use Bio::Tools::IUPAC;

 # Get the IUPAC code for proteins
 my %iupac_prot = Bio::Tools::IUPAC->new->iupac_iup;

 # Create a sequence with degenerate residues
 my $ambiseq = Bio::PrimarySeq->new(-seq => 'ARTCGUTGN', -alphabet => 'dna');

 # Create all possible non-degenerate sequences
 my $iupac = Bio::Tools::IUPAC->new(-seq => $ambiseq);
 while ($uniqueseq = $iupac->next_seq()) {
     # process the unique Bio::Seq object.
 }

 # Get a regular expression that matches all possible sequences
 my $regexp = $iupac->regexp();

=head1 DESCRIPTION

Bio::Tools::IUPAC is a tool that manipulates sequences with ambiguous residues
following the IUPAC conventions. Non-standard characters have the meaning 
described below:

    IUPAC-IUB SYMBOLS FOR NUCLEOTIDE (DNA OR RNA) NOMENCLATURE:
      Cornish-Bowden (1985) Nucl. Acids Res. 13: 3021-3030

    ---------------------------------------------------------------
    Symbol       Meaning      Nucleic Acid
    ---------------------------------------------------------------
     A            A           Adenine
     C            C           Cytosine
     G            G           Guanine
     T            T           Thymine
     U            U           Uracil
     M          A or C        aMino
     R          A or G        puRine
     W          A or T        Weak
     S          C or G        Strong
     Y          C or T        pYrimidine
     K          G or T        Keto
     V        A or C or G     not T (closest unused char after T)
     H        A or C or T     not G (closest unused char after G)
     D        A or G or T     not C (closest unused char after C)
     B        C or G or T     not A (closest unused char after A)
     X      G or A or T or C  Unknown (very rarely used)
     N      G or A or T or C  Unknown (commonly used)


    IUPAC-IUP AMINO ACID SYMBOLS:
      Biochem J. 1984 Apr 15; 219(2): 345-373
      Eur J Biochem. 1993 Apr 1; 213(1): 2

    ------------------------------------------
    Symbol           Meaning
    ------------------------------------------
    A        Alanine
    B        Aspartic Acid, Asparagine
    C        Cysteine
    D        Aspartic Acid
    E        Glutamic Acid
    F        Phenylalanine
    G        Glycine
    H        Histidine
    I        Isoleucine
    J        Isoleucine/Leucine
    K        Lysine
    L        Leucine
    M        Methionine
    N        Asparagine
    O        Pyrrolysine
    P        Proline
    Q        Glutamine
    R        Arginine
    S        Serine
    T        Threonine
    U        Selenocysteine
    V        Valine
    W        Tryptophan
    X        Unknown
    Y        Tyrosine
    Z        Glutamic Acid, Glutamine
    *        Terminator

There are a few things Bio::Tools::IUPAC can do for you:

=over

=item *

report the IUPAC mapping between ambiguous and non-ambiguous residues

=item *

produce a stream of all possible corresponding unambiguous Bio::Seq objects given
an ambiguous sequence object

=item *

convert an ambiguous sequence object to a corresponding regular expression

=back

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

=head1 AUTHOR - Aaron Mackey

Email amackey-at-virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


package Bio::Tools::IUPAC;

use strict;
use base qw(Bio::Root::Root);
use vars qw(%IUB %IUB_AMB %REV_IUB %IUP %IUP_AMB $AUTOLOAD);

BEGIN {
    # Ambiguous nucleic residues are matched to unambiguous residues
    %IUB = (
        A => [qw(A)],
        C => [qw(C)],
        G => [qw(G)],
        T => [qw(T)],
        U => [qw(U)],
        M => [qw(A C)],
        R => [qw(A G)],
        S => [qw(C G)],
        W => [qw(A T)],
        Y => [qw(C T)],
        K => [qw(G T)],
        V => [qw(A C G)],
        H => [qw(A C T)],
        D => [qw(A G T)],
        B => [qw(C G T)],
        N => [qw(A C G T)],
        X => [qw(A C G T)],
    );

    # Same as %IUB but ambiguous residues are matched to ambiguous residues only
    %IUB_AMB = (
        M => [qw(M)],
        R => [qw(R)],
        W => [qw(W)],
        S => [qw(S)],
        Y => [qw(Y)],
        K => [qw(K)],
        V => [qw(M R S V)],
        H => [qw(H M W Y)],
        D => [qw(D K R W)],
        B => [qw(B K S Y)],
        N => [qw(B D H K M N R S V W Y)],
    );

    # The inverse of %IUB
    %REV_IUB = (
        A    => 'A',
        T    => 'T',
        U    => 'U',
        C    => 'C',
        G    => 'G',
        AC   => 'M',
        AG   => 'R',
        AT   => 'W',
        CG   => 'S',
        CT   => 'Y',
        GT   => 'K',
        ACG  => 'V',
        ACT  => 'H',
        AGT  => 'D',
        CGT  => 'B',
        ACGT => 'N',
        N    => 'N'
    );

    # Same thing with proteins now
    %IUP = (
        A => [qw(A)],
        B => [qw(D N)],
        C => [qw(C)],
        D => [qw(D)],
        E => [qw(E)],
        F => [qw(F)],
        G => [qw(G)],
        H => [qw(H)],
        I => [qw(I)],
        J => [qw(I L)],
        K => [qw(K)],
        L => [qw(L)],
        M => [qw(M)],
        N => [qw(N)],
        O => [qw(O)],
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
        '*' => [qw(*)],
    );

    %IUP_AMB = (
        B => [qw(B)],
        J => [qw(J)],
        Z => [qw(Z)],
    );

}


=head2 new

 Title   : new
 Usage   : Bio::Tools::IUPAC->new($seq);
 Function: Create a new IUPAC object, which acts as a sequence stream (akin to
           SeqIO)
 Args    : an ambiguously coded sequence object that has a specified 'alphabet'
 Returns : a Bio::Tools::IUPAC object.

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($seq) = $self->_rearrange([qw(SEQ)],@args);

    if ( (not defined $seq) && @args && ref($args[0]) ) {
        # parameter not passed as named parameter?
        $seq = $args[0];
    }

    if (defined $seq) {
        if (not $seq->isa('Bio::PrimarySeqI')) {
            $self->throw('Must supply a sequence object');
        }
        if (length $seq->seq == 0) {
            $self->throw('Sequence had zero-length');
        }
        $self->{'_seq'} = $seq;
    }

    return $self;
}


sub _initialize {
    my ($self) = @_;
    my %iupac = $self->iupac;
    $self->{'_alpha'} = [ map { $iupac{uc $_} } split('', $self->{'_seq'}->seq) ];
    $self->{'_string'} = [(0) x length($self->{'_seq'}->seq())];
    $self->{'_string'}->[0] = -1;
}


=head2 next_seq

 Title   : next_seq
 Usage   : $iupac->next_seq();
 Function: returns the next unique sequence object
 Args    : none.
 Returns : a Bio::Seq object

=cut

sub next_seq {
    my ($self) = @_;

    if (not exists $self->{'_string'}) {
        $self->_initialize();
    }

    for my $i ( 0 .. $#{$self->{'_string'}} ) {
        next unless $self->{'_string'}->[$i] || @{$self->{'_alpha'}->[$i]} > 1;
        if ( $self->{'_string'}->[$i] == $#{$self->{'_alpha'}->[$i]} ) { # rollover
            if ( $i == $#{$self->{'_string'}} ) { # end of possibilities
                return;
            } else {
                $self->{'_string'}->[$i] = 0;
                next;
            }
        } else {
            $self->{'_string'}->[$i]++;
            my $j = -1;
            my $seqstr = join('', map { $j++; $self->{'_alpha'}->[$j]->[$_]; } @{$self->{'_string'}});
            my $desc   = $self->{'_seq'}->desc() || '';
            $self->{'_num'}++;
            1 while $self->{'_num'} =~ s/(\d)(\d\d\d)(?!\d)/$1,$2/;
            $desc =~ s/( \[Bio::Tools::IUPAC-generated\sunique sequence # [^\]]*\])|$/ \[Bio::Tools::IUPAC-generated unique sequence # $self->{'_num'}\]/;
            $self->{'_num'} =~ s/,//g;

            # Return a fresh sequence object
            return Bio::PrimarySeq->new(-seq  => $seqstr, -desc => $desc);
        }
    }
}


=head2 iupac

 Title   : iupac
 Usage   : my %symbols = $iupac->iupac;
 Function: Returns a hash of symbols -> symbol components of the right type
           for the given sequence, i.e. it is the same as iupac_iup() if
           Bio::Tools::IUPAC was given a proteic sequence, or iupac_iub() if the 
           sequence was nucleic. For example, the key 'M' has the value ['A', 'C'].
 Args    : none
 Returns : Hash

=cut

sub iupac {
    my ($self) = @_;
    my $alphabet = lc( $self->{'_seq'}->alphabet() );
    if ( ($alphabet eq 'dna') or ($alphabet eq 'rna') ) {
        return %IUB; # nucleic
    } elsif ( $alphabet eq 'protein' ) {
        return %IUP; # proteic
    } else {
        $self->throw("The input sequence had the unknown alphabet '$alphabet'\n");
    }
}



=head2 iupac_amb

 Title   : iupac_amb
 Usage   : my %symbols = $iupac->iupac_amb;
 Function: Same as iupac() but only contains a mapping between ambiguous residues
           and the ambiguous residues they map to. For example, the key 'N' has
           the value ['R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N'],
           i.e. it matches all other ambiguous residues.
 Args    : none
 Returns : Hash

=cut

sub iupac_amb {
    my ($self) = @_;
    my $alphabet = lc( $self->{'_seq'}->alphabet() );
    if ( ($alphabet eq 'dna') or ($alphabet eq 'rna') ) {
        return %IUB_AMB; # nucleic
    } elsif ( $alphabet eq 'protein' ) {
        return %IUP_AMB; # proteic
    } else {
        $self->throw("The input sequence had the unknown alphabet '$alphabet'\n");
    }
}


=head2 iupac_iup

 Title   : iupac_iup
 Usage   : my %aasymbols = $iupac->iupac_iup;
 Function: Returns a hash of PROTEIN symbols -> non-ambiguous symbol components
 Args    : none
 Returns : Hash

=cut

sub iupac_iup {
   return %IUP;
}


=head2 iupac_iup_amb

 Title   : iupac_iup_amb
 Usage   : my %aasymbols = $iupac->iupac_iup_amb;
 Function: Returns a hash of PROTEIN symbols -> ambiguous symbol components
 Args    : none
 Returns : Hash

=cut

sub iupac_iup_amb {
   return %IUP_AMB;
}


=head2 iupac_iub

 Title   : iupac_iub
 Usage   : my %dnasymbols = $iupac->iupac_iub;
 Function: Returns a hash of DNA symbols -> non-ambiguous symbol components
 Args    : none
 Returns : Hash

=cut

sub iupac_iub {
   return %IUB;
}


=head2 iupac_iub_amb

 Title   : iupac_iub_amb
 Usage   : my %dnasymbols = $iupac->iupac_iub;
 Function: Returns a hash of DNA symbols -> ambiguous symbol components
 Args    : none
 Returns : Hash

=cut

sub iupac_iub_amb {
   return %IUB_AMB;
}


=head2 iupac_rev_iub

 Title   : iupac_rev_iub
 Usage   : my %dnasymbols = $iupac->iupac_rev_iub;
 Function: Returns a hash of nucleotide combinations -> IUPAC code
           (a reverse of the iupac_iub hash).
 Args    : none
 Returns : Hash

=cut

sub iupac_rev_iub {
   return %REV_IUB;
}


=head2 count

 Title   : count
 Usage   : my $total = $iupac->count();
 Function: Calculates the number of unique, unambiguous sequences that
           this ambiguous sequence could generate
 Args    : none
 Return  : int

=cut

sub count {
    my ($self) = @_;
    if (not exists $self->{'_string'}) {
        $self->_initialize();
    }
    my $count = 1;
    $count *= scalar(@$_) for (@{$self->{'_alpha'}});
    return $count;
}


=head2 regexp

 Title   : regexp
 Usage   : my $re = $iupac->regexp();
 Function: Converts the ambiguous sequence into a regular expression that
           matches all of the corresponding ambiguous and non-ambiguous sequences.
           You can further manipulate the resulting regular expression with the
           Bio::Tools::SeqPattern module. After you are done building your
           regular expression, you might want to compile it and make it case-
           insensitive:
              $re = qr/$re/i;
 Args    : 1 to match RNA: T and U characters will match interchangeably
 Return  : regular expression

=cut

sub regexp {
    my ($self, $match_rna) = @_;
    my $re;
    my $seq = $self->{'_seq'}->seq;
    my %iupac = $self->iupac;
    my %iupac_amb = $self->iupac_amb;
    for my $pos (0 .. length($seq)-1) {
        my $res = substr $seq, $pos, 1;
        my $iupacs = $iupac{$res};
        my $iupacs_amb = $iupac_amb{$res} || [];
        if (not defined $iupacs) {
            $self->throw("Primer sequence '$seq' is not a valid IUPAC sequence.".
                         " Offending character was '$res'.\n");
        }
        my $part = join '', (@$iupacs, @$iupacs_amb);
        if ($match_rna) {
            $part =~ s/T/TU/i || $part =~ s/U/TU/i;
        }
        if (length $part > 1) {
           $part = '['.$part.']';
        }
        $re .= $part;
    }
    return $re;
}


sub AUTOLOAD {
    my $self = shift @_;
    my $method = $AUTOLOAD;
    $method =~ s/.*:://;
    return $self->{'_seq'}->$method(@_)
        unless $method eq 'DESTROY';
}

1;

