#!/usr/bin/perl

use lib '.';

use strict;
use Bio::DB::GFF;

my $db = Bio::DB::GFF->new(-adaptor=>'dbi::mysqlopt',
			   -dsn=>'dbi:mysql:elegans'
			  ) or die;

$db->initialize(1);  # drop and reinitialize schema!
# $db->lock_on_load(1);
$db->load("$ENV{HOME}/projects/mysql_gff");
print "done\n";

