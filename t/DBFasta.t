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
    plan test => 11;
}
use Bio::DB::Fasta;
use Bio::Root::IO;
my $db = Bio::DB::Fasta->new(Bio::Root::IO->catfile(qw(. t data dbfa)),
	                     -reindex=>1);
ok($db);
ok($db->length('CEESC13F') > 0);
ok(length $db->seq('CEESC13F:1,10') == 10);
ok(length $db->seq('AW057119',1,10) == 10);
my $primary_seq = $db->get_Seq_by_id('AW057119');
ok($primary_seq);
ok(length($primary_seq->seq) > 0);
undef $db;
undef $primary_seq;

my (%h,$dna1,$dna2);
ok(tie(%h,'Bio::DB::Fasta',Bio::Root::IO->catfile(qw(. t data dbfa))));
ok($h{'AW057146'});
ok($dna1 = $h{'AW057146:1,10'});
ok($dna2 = $h{'AW057146:10,1'});

my $revcom = reverse $dna1;
$revcom =~ tr/gatcGATC/ctagCTAG/;
ok($dna2 eq $revcom);

END {
  unlink Bio::Root::IO->catfile(qw(t data dbfa directory.index));
  unlink Bio::Root::IO->catfile(qw(t data dbfa directory.index.dir));
  unlink Bio::Root::IO->catfile(qw(t data dbfa directory.index.pag));
}
