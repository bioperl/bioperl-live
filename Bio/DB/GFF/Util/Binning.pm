=head1 NAME

Bio::DB::GFF::Util::Binning - binning utility for Bio::DB::GFF index

=head1 SYNOPSIS

 use Bio::DB::GFF::Util::Binning qw(bin bin_bot bin_top);
 my $tier = bin($start,$stop,$min);

=head1 DESCRIPTION

This is a utility module that exports the functions bin(), bin_bot()
and bin_top().  These functions translate a range on the genome into a
named bin that is used as an index in the Bio::DB::GFF schema.  The
index makes certain range retrieval queries much faster.

=head1 API

The remainder of the document describes the function calls.  No calls
are exported by default, but must be imported explicitly.

=over 4

=cut

package Bio::DB::GFF::Util::Binning;

use strict;
require Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
@ISA = 'Exporter';
@EXPORT_OK = qw(bin bin_bot bin_top);
@EXPORT = @EXPORT_OK;

=item $bin_name = bin($start,$stop,$bin_size)

Given a start, stop and bin size on the genome, translate this
location into a bin name.  In a list context, returns the bin tier
name and the position that the bin begins.

=cut

sub bin {
  my ($start,$stop,$min) = @_;
  my $tier = $min;
  my ($bin_start,$bin_end);
  while (1) {
    $bin_start = int $start/$tier;
    $bin_end   = int $stop/$tier;
    last if $bin_start == $bin_end;
    $tier *= 10;
  }
  return wantarray ? ($tier,$bin_start) : bin_name($tier,$bin_start);
}

=item $bottom = bin_bot($tier,$start)

Given a tier name and a range start position, returns the lower end of
the bin range.

=cut

sub bin_bot {
  my $tier = shift;
  my $pos  = shift;
  bin_name($tier,int($pos/$tier));
}

=item $top = bin_top($tier,$end)

Given a tier name and the end of a range, returns the upper end of the
bin range.

=cut

*bin_top = \&bin_bot;

sub bin_name { sprintf("%d.%06d",@_) }

sub log10 {
  my $i = shift;
  log($i)/log(10);
}

1;

=back

=head1 BUGS

None known yet.

=head1 SEE ALSO

L<Bio::DB::GFF>,

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
