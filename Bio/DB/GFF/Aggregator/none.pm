=head1 NAME

Bio::DB::GFF::Aggregator::none -- No aggregation

=head1 SYNOPSIS

  use Bio::DB::GFF;

  # Open the sequence database
  my $db      = Bio::DB::GFF->new( -adaptor => 'dbi:mysql',
                                   -dsn     => 'dbi:mysql:elegans42',
				   -aggregator => 'none'
				 );


=head1 DESCRIPTION

Bio::DB::GFF::Aggregator::none can be used to indicate that you do not
want any aggregation performed.  It is equivalent to providing undef
to the B<-aggregator> argument.  It overrides disaggregate() and
aggregate() so that they do exactly nothing.

=cut

package Bio::DB::GFF::Aggregator::none;

use strict;

use base qw(Bio::DB::GFF::Aggregator);

sub disaggregate {
  my $self  = shift;
  my $types = shift;
  # no change
}

sub aggregate {
  my $self = shift;
  my $features = shift;
  return;  # no change
}

1;
