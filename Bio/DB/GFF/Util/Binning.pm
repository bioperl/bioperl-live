package Bio::DB::GFF::Util::Binning;

use strict;
require Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
@ISA = 'Exporter';
@EXPORT_OK = qw(bin bin_bot bin_top);
@EXPORT = @EXPORT_OK;

sub bin {
  my ($start,$stop,$min) = @_;
  my $magnitude = int(1+log10(($stop - $start) || 1));
  my $tier = 10 ** $magnitude;
  $tier = $min if $tier < $min;
  my $bin_start = int $start/$tier;
  my $bin_end   = int $stop/$tier;
  $tier *= 10 if $bin_start != $bin_end;
  my $index = int $start/$tier;
  return wantarray ? ($tier,$index) : bin_name($tier,$index);
}

sub bin_bot {
  my $tier = shift;
  my $pos  = shift;
  my $bin_no = int($pos/$tier);
  bin_name($tier,$bin_no > 0 ? $bin_no-1 : 0);
}

sub bin_top {
  my $tier = shift;
  my $pos  = shift;
  bin_name($tier,int($pos/$tier));
}

sub bin_name { sprintf("%d.%06d",@_) }

sub log10 {
  my $i = shift;
  log($i)/log(10);
}

1;
