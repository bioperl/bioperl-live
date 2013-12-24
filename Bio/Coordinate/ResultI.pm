package Bio::Coordinate::ResultI;
use utf8;
use strict;
use warnings;
use parent qw(Bio::LocationI);

# ABSTRACT: Interface to identify coordinate mapper results.
# AUTHOR:   Heikki Lehvaslaiho <heikki@bioperl.org>
# OWNER:    Heikki Lehvaslaiho
# LICENSE:  Perl_5

=head1 SYNOPSIS

  # not to be used directly

=head1 DESCRIPTION

ResultI identifies Bio::LocationIs returned by
Bio::Coordinate::MapperI implementing classes from other locations.

=cut

1;
