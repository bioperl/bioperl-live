package Bio::DB::GFF::Segment;

use strict;
use Bio::Root::RootI;

use vars qw($VERSION @ISA);
@ISA = qw(Bio::Root::RootI);
$VERSION = '0.20';

use overload '""' => 'asString';

*abs_start = \&start;
*abs_stop  = \&stop;

sub new {
  my $class = shift;
  my ($factory,$segname,$start,$stop,$segclass) = @_;

  $factory or $class->throw("->new(): provide a factory argument");
  $class = ref $class if ref $class;
  return bless { factory   => $factory,
		 sourceseq => $segname,
		 class     => $segclass,
		 start     => $start,
		 stop      => $stop,
		 strand    => 0,
	       },$class;
}

# read-only accessors
sub factory { shift->{factory} }

# start, stop, length
sub start  { shift->{start} }
sub stop   { shift->{stop}  }
sub length { abs($_[0]->{start} - $_[0]->{stop})+1 }
sub strand { 
  my $self = shift; 
  return $self->{stop} >= $self->{start} ? +1 : -1; 
}
sub refseq    { shift->{sourceseq} }
sub sourceseq { shift->{sourceseq} }
sub class     { shift->{class}     }

# deep copy of the thing
sub clone {
  my $self = shift;
  my %h = %$self;
  return bless \%h,ref($self);
}

sub subseq {
  my $self = shift;
  my ($newstart,$newstop) = @_;
  my ($refseq,$start,$stop,$class) = ($self->{sourceseq},
				      $self->{start},$self->{stop},
				      $self->class);

  # We deliberately force subseq to return objects of type RelSegment
  # Otherwise, when we get a subsequence from a Feature object,
  # its method and source go along for the ride, which is incorrect.
  my $new = $self->new_from_segment($self);
  if ($start <= $stop) {
    @{$new}{qw(start stop)} = ($start + $newstart - 1, $start + $newstop  - 1);
  } else {
    @{$new}{qw(start stop)} = ($start - ($newstart - 1), $start - ($newstop  - 1)),

  }

  $new;
}

sub asString {
  my $self = shift;
  my $label = $self->refseq;
  my $start = $self->start;
  my $stop  = $self->stop;
  return "$label:$start,$stop";
}

sub id { shift->asString }

sub dna {
  my $self = shift;
  my ($ref,$class,$start,$stop,$strand) 
    = @{$self}{qw(sourceseq class start stop strand)};
  ($start,$stop) = ($stop,$start) if $strand eq '-';
  $self->factory->dna($ref,$class,$start,$stop);
}


1;
__END__
# Below is stub documentation for your module. You better edit it!

=head1 NAME

Ace::Sequence::Mysql - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Ace::Sequence::Mysql;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for Ace::Sequence::Mysql, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.


=head1 AUTHOR

A. U. Thor, a.u.thor@a.galaxy.far.far.away

=head1 SEE ALSO

perl(1).

=cut
