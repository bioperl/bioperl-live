# 
# BioPerl module for Bio::Tools::Gel
# Copyright Allen Day <allenday@ucla.edu>
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Gel - Calculates relative electrophoretic migration distances

=head1 SYNOPSIS

    use Bio::PrimarySeq;
    use Bio::Restriction::Analysis;
    use Bio::Tools::Gel;

    # get a sequence
    my $d = 'AAAAAAAAAGAATTCTTTTTTTTTTTTTTGAATTCGGGGGGGGGGGGGGGGGGGG';
    my $seq1 = Bio::Seq->new(-id=>'groundhog day',-seq=>$d);

    # cut it with an enzyme
    my $ra=Bio::Restriction::Analysis->new(-seq=>$seq1);
    @cuts = $ra->fragments('EcoRI'), 3;

    # analyse the fragments in a gel
    my $gel = Bio::Tools::Gel->new(-seq=>\@cuts,-dilate=>10);
    my %bands = $gel->bands;
    foreach my $band (sort {$b <=> $a} keys %bands){
      print $band,"\t", sprintf("%.1f", $bands{$band}),"\n";
    }

    #prints:
    #20   27.0
    #25   26.0
    #10   30.0

=head1 DESCRIPTION

This takes a set of sequences or Bio::Seq objects, and calculates their
respective migration distances using:
    distance = dilation * (4 - log10(length(dna));

Source: Molecular Cloning, a Laboratory Manual. Sambrook, Fritsch, Maniatis. 
CSHL Press, 1989.

Bio::Tools::Gel currently calculates migration distances based solely on
the length of the nucleotide sequence.  Secondary or tertiary structure, 
curvature, and other biophysical attributes of a sequence are currently 
not considered.  Polypeptide migration is currently not supported.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Allen Day

Email allenday@ucla.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::Tools::Gel;
use strict;

use Bio::PrimarySeq;

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $gel = Bio::Tools::Gel->new(-seq => $sequence,-dilate => 3);
 Function: Initializes a new Gel
 Returns : Bio::Tools::Gel
 Args    : -seq      => Bio::Seq(s), scalar(s) or list of either/both 
                        (default: none)
           -dilate   => Expand band migration distances (default: 1)

=cut

sub new {
    my ($class, @args) = @_;

    my $self = $class->SUPER::new(@args);
    my ($seqs, $dilate) = $self->_rearrange([qw(SEQ DILATE)], @args);
    if( ! ref($seqs)  ) {
        $self->add_band([$seqs]);
    } elsif( ref($seqs) =~ /array/i ||
        $seqs->isa('Bio::PrimarySeqI') ) {
        $self->add_band($seqs);
    }
    $self->dilate($dilate || 1);
  
    return $self;
}


=head2 add_band

 Title   : add_band
 Usage   : $gel->add_band($seq);
 Function: Calls _add_band with a (possibly created) Bio::Seq object.
 Returns : 
 Args    : Bio::Seq, scalar sequence, or list of either/both.

=cut

sub add_band {
    my ($self, $args) = @_;

    foreach my $arg (@$args){
        my $seq;
        if( ! ref $arg ) {
            if( $arg =~ /^\d+/ ) {
                # $arg is a number
                $seq = Bio::PrimarySeq->new(-seq=>'N'x$arg, -id => $arg);
            } else {
                # $arg is a sequence string
                $seq = Bio::PrimarySeq->new(-seq=>$arg, -id=>length $arg);
            }
        } elsif( $arg->isa('Bio::PrimarySeqI') ) {
            # $arg is a sequence object
            $seq = $arg;
        }

        $self->_add_band($seq);
    }
    return 1;
}


=head2 _add_band

 Title   : _add_band
 Usage   : $gel->_add_band($seq);
 Function: Adds a new band to the gel.
 Returns : 
 Args    : Bio::Seq object

=cut

sub _add_band {
    my ($self, $arg) = @_;  
    if ( defined $arg) {
        push (@{$self->{'bands'}},$arg);
    }
    return 1;
}


=head2 dilate

 Title   : dilate
 Usage   : $gel->dilate(1);
 Function: Sets/retrieves the dilation factor.
 Returns : dilation factor 
 Args    : Float or none

=cut

sub dilate {
    my ($self, $arg) = @_;
    return $self->{dilate} unless $arg;
    $self->throw("-dilate should be numeric") if defined $arg and $arg =~ /[^e\d\.]/;
    $self->{dilate} = $arg;
    return $self->{dilate};
}


sub migrate {
    my ($self, $arg) = @_;
    $arg = $self unless $arg;
    if ( $arg ) {
        return 4 - log10($arg);
    } else {
        return 0;
    }
}


=head2 bands

 Title   : bands
 Usage   : $gel->bands;
 Function: Calculates migration distances of sequences.
 Returns : hash of (seq_id => distance)
 Args    : 

=cut

sub bands {
    my $self = shift;
    $self->throw("bands() is read-only") if @_;

    my %bands = ();

    foreach my $band (@{$self->{bands}}){
        my $distance = $self->dilate * migrate($band->length);
        $bands{$band->id} = $distance;
    }

    return %bands;
}


=head2 log10

 Title   : log10
 Usage   : log10($n);
 Function: returns base 10 log of $n.
 Returns : float
 Args    : float

=cut

# from "Programming Perl"
sub log10 {
    my $n = shift;
    return log($n)/log(10);
}

1;
