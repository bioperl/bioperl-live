#!/usr/bin/perl
use strict;

use Bio::DB::Registry;
use Bio::SeqIO;
use strict;

my $registry = new Bio::DB::Registry();

print "services are ", join(',', $registry->services), "\n";
my $db = $registry->get_database("embl");
my $seq = $db->get_Seq_by_id("J02231");
my $out = new Bio::SeqIO;
$out->write_seq($seq);
