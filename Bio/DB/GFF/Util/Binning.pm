package Bio::DB::GFF::Util::Binning;

use strict;
require Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
@ISA = 'Exporter';
@EXPORT_OK = qw(bin bin_bot bin_top);
@EXPORT = @EXPORT_OK;

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

sub bin_bot {
  my $tier = shift;
  my $pos  = shift;
  bin_name($tier,int($pos/$tier));
}
*bin_top = \&bin_bot;

sub bin_name { sprintf("%d.%06d",@_) }

sub log10 {
  my $i = shift;
  log($i)/log(10);
}

1;
