#-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use lib '.','./blib/lib';

BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan test => 9;
}
use Bio::DB::Fasta;

my $db = Bio::DB::Fasta->new('./t/data/dbfa',-reindex=>1);
ok($db);
ok($db->length('CEESC13F') > 0);
ok(length $db->seq('CEESC13F:1,10') == 10);
ok(length $db->seq('AW057119',1,10) == 10);
undef $db;
my (%h,$dna1,$dna2);

ok(tie(%h,'Bio::DB::Fasta','./t/data/dbfa'));
ok($h{'AW057146'});
ok($dna1 = $h{'AW057146:1,10'});
ok($dna2 = $h{'AW057146:10,1'});

my $revcom = reverse $dna1;
$revcom =~ tr/gatcGATC/ctagCTAG/;
ok($dna2 eq $revcom);

END {
  unlink "./t/data/dbfa/directory.index";
}
