# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use Test;
use strict;
BEGIN { plan tests => 3}
use Bio::Tools::Blast::HTML;
use Bio::Tools::Blast;

ok(1);
my $blastobj = new Bio::Tools::Blast();
ok(ref($blastobj),qr/Bio::Tools::Blast/);
my @reportarray;
$blastobj->to_html(-file=>'t/blast.report', 
		   -out=>\@reportarray);
ok(@reportarray, 428, 'did not produce any htmlified blast report lines');

