package Bio::DB::GFF::Aggregator::transcript;

use strict;
use Bio::DB::GFF::Aggregator;

use vars qw($VERSION @ISA);
@ISA = qw(Bio::DB::GFF::Aggregator);

$VERSION = '0.10';

sub method { 'transcript' }

sub part_names {
  return qw(intron exon CDS);
}

sub main_name {
#  return 'Sequence:curated';
  return 'Sequence';
}

1;
