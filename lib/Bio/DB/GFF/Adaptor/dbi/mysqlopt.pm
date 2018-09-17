package Bio::DB::GFF::Adaptor::dbi::mysqlopt;

=head1 NAME

Bio::DB::GFF::Adaptor::dbi::mysqlopt -- Deprecated database adaptor

=head1 SYNOPSIS

This adaptor has been superseded by Bio::DB::GFF::Adaptor::dbi::mysql.

See L<Bio::DB::GFF> and L<Bio::DB::GFF::Adaptor::dbi::mysql>

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bioperl>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2002 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

use strict;
use base qw(Bio::DB::GFF::Adaptor::dbi::mysql);

1;
