use strict;
use Bio::Root::IO; # needed for the BEGIN block

BEGIN {	
    eval { require Test; };
    if( $@ ) { 
	use lib 't';
    }
    use Test;

    plan tests => 29;
    use vars qw($loaded $testout $expectedout);
    $testout = "blast.t.out";  # output from this script.
    $expectedout = Bio::Root::IO->catfile("t","expected.blast.out");
    unlink $testout;
    $^W = 0; 
}
END {
    ok(0) unless $loaded;
    unlink $testout;  # commented out since you may want to check it...
}

use Bio::Tools::Blast;
$loaded = 1;

my($blast,@hits,@inds,$cfile,$ufile);

open (OUT,">$testout");

ok $blast = Bio::Tools::Blast->new(-file   =>Bio::Root::IO->catfile("t","blast.report"),
					-signif => 1e-5,
					-parse  => 1,
					-stats  => 1,
					-check_all_hits => 1,
					);
ok $blast->display();
ok $blast->is_signif;
ok $blast->signif, '1.0e-05', "Signif: ".$blast->signif;
ok $blast->num_hits, 4;
ok $blast->length, 504;
ok $blast->program, 'TBLASTN';
ok $blast->query,  'gi|1401126';
ok $blast->hit->name, 'gb|U49928|HSU49928';
ok $blast->hit->length, 3096;

@hits  = $blast->hits;

ok $hits[0]->expect, '0.0';
ok $hits[1]->expect, '4e-07';
ok $hits[2]->expect, '1e-05';
ok $hits[1]->frac_identical, '0.25';
ok $hits[1]->hsp->frac_conserved, '0.43';
ok $hits[1]->hsp->score, 137;
ok $hits[1]->hsp->bits, '57.8';

ok scalar (@inds = $hits[1]->hsp->seq_inds('query', 'iden', 1));
ok $inds[0], '66-68';

ok print OUT $blast->table_labels;
ok print OUT $blast->table;
print OUT "\n\n";
ok print OUT $blast->table_labels_tiled;
ok print OUT $blast->table_tiled;
print OUT "\n\n";
close OUT;

ok (-s $blast->file);
ok ($cfile = $blast->compress_file);
ok (-s $cfile and -B $cfile), 1,"Can't compress Blast file";
ok ($ufile = $blast->uncompress_file);
ok (-s $ufile and -T $ufile), 1,"Can't uncompress Blast file";

print "checking expected output...\n";

ok system('diff', $testout, $expectedout), 0, "diff $testout $expectedout";
