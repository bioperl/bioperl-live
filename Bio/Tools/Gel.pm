# BioPerl module for Bio::Tools::Gel
# Copyright Allen Day <allenday@ucla.edu>
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Gel - Calculates relative electrophoretic migration distances


=head1 SYNOPSIS

    #An example of a virtual restriction digest and subsequent gel run
    use Bio::Seq;
    use Bio::Tools::RestrictionEnzyme;
    use Bio::Tools::Gel;

    my $d = 'AAAAAAAAAGAATTCTTTTTTTTTTTTTTGAATTCGGGGGGGGGGGGGGGGGGGG';
    my $seq1 = Bio::Seq->new(-id=>'groundhog day',-seq=>$d);
    my $EcoRI = Bio::Tools::RestrictionEnzyme->new(-NAME=>'EcoRI');
    my @cuts = $EcoRI->cut_seq($seq);

    my $gel = Bio::Tools::Gel->new(-seq=>\@cuts,-dilate=>10);
    my %bands = $gel->bands;
    foreach my $band (keys %bands){
      print $band,"\t",$bands{$band},"\n";
    }

    #prints:
    #25      26.0205999132796
    #10      30
    #20      26.9897000433602


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

  bioperl-l@bioperl.org            - General discussion
http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Allen Day

Email allenday@ucla.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Gel;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::Seq;

@ISA = qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $gel = new Bio::Tools::Gel(-seq => $sequence,-dilate => 3);
 Function: Initializes a new Gel
 Returns : Bio::Tools::Gel
 Args    : -seq      => Bio::Seq(s), scalar(s) or list of either/both 
                        (default: none)
           -dilate   => Expand band migration distances (default: 1)

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->_init(@args);
  
  return $self;
}

=head2 _init

 Title   : _init
 Usage   : _init(@args);
 Function: Handles arguments to new()
 Returns : 
 Args    : @args from new()

=cut

sub _init {
  my($self,%args) = @_;

  my @seqs = ref $args{-seq} eq 'ARRAY' ? @{$args{-seq}}
                                        :   $args{-seq}; 
  $self->add_band(\@seqs);
  $self->dilate($args{-dilate} || 1);
}

=head2 add_band

 Title   : add_band
 Usage   : $gel->add_band($seq);
 Function: Calls _add_band with a (possibly created) Bio::Seq object.
 Returns : 
 Args    : Bio::Seq, scalar sequence, or list of either/both.

=cut

sub add_band {
  my($self,$args) = @_;

  foreach my $arg (@$args){

    my $seq = ref $arg eq 'Bio::Seq' ? $arg
            : Bio::Seq->new(-seq=>$arg,-id=>length($arg));

    $seq->validate_seq or $seq->throw("invalid symbol in sequence");
    $self->_add_band($seq);
  }
}

=head2 _add_band

 Title   : _add_band
 Usage   : $gel->_add_band($seq);
 Function: Adds a new band to the gel.
 Returns : 
 Args    : Bio::Seq object

=cut

sub _add_band {
  my($self,$arg) = @_;
  push @{$self->{bands}},$arg;
}

=head2 dilate

 Title   : dilate
 Usage   : $gel->dilate(1);
 Function: Sets/retrieves the dilation factor.
 Returns : dilation factor 
 Args    : Float or none

=cut

sub dilate {
  my($self,$arg) = @_;
  return $self->{dilate} unless $arg;
  $self->throw("-dilate should be numeric") if defined $arg and $arg =~ /[^e\d\.]/;
  $self->{dilate} = $arg;
  return $self->{dilate};
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
    my $distance = $self->dilate * (4 - log10($band->length));
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

#from programming perl
sub log10 {
  my $n = shift;
  return log($n)/log(10);
}

1;
