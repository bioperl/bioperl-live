package Bio::Coordinate::Result::Gap;
use utf8;
use strict;
use warnings;
use parent qw(Bio::Location::Simple Bio::Coordinate::ResultI);

# ABSTRACT: Another name for L<Bio::Location::Simple>.
# AUTHOR:   Heikki Lehvaslaiho <heikki@bioperl.org>
# OWNER:    Heikki Lehvaslaiho
# LICENSE:  Perl_5

=head1 SYNOPSIS

  $loc = Bio::Coordinate::Result::Gap->new(-start=>10,
                                          -end=>30,
                                          -strand=>1);

=head1 DESCRIPTION

This is a location object for coordinate mapping results.
=cut

1;
