#-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use lib '.','./blib/lib';
use vars qw($DEBUG $NUMTESTS $exit);

$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    eval {
	require Bio::DB::Fasta;
    };
    if ( $@ ) {
	warn("A depedendancy for Bio::DB::Fasta is not installed - skipping tests. $@\n") if $DEBUG;
	$exit = 1;
    }
    plan test => ($NUMTESTS = 12);
}
require Bio::DB::Fasta;
use Bio::Root::IO;

# this obfuscation is to deal with lockfiles by GDBM_File which can
# only be created on local filesystems apparently so will cause test
# to block and then fail when the testdir is on an NFS mounted system

use File::Copy;
my $io = Bio::Root::IO->new(-verbose => $DEBUG);
my $tempdir = $io->tempdir('CLEANUP' => 1);
my $test_dbdir = $io->catfile($tempdir, 'dbfa');
mkdir($test_dbdir); # make the directory
my $indir = $io->catfile(qw(. t data dbfa)); 
opendir(INDIR,$indir) || die("cannot open dir $indir");
# effectively do a cp -r but only copy the files that are in there, no subdirs
for my $file ( map { $io->catfile($indir,$_) } readdir(INDIR) ) {
    next unless (-f $file );
    copy($file, $test_dbdir);
}
closedir(INDIR);

# now use this temporary dir for the db file
my $db = Bio::DB::Fasta->new($test_dbdir, -reindex => 1);
ok($db);
ok($db->length('CEESC13F') > 0);
ok(length $db->seq('CEESC13F:1,10') == 10);
ok(length $db->seq('AW057119',1,10) == 10);
my $primary_seq = $db->get_Seq_by_id('AW057119');
ok($primary_seq);
ok(length($primary_seq->seq) > 0);
ok(!defined $db->get_Seq_by_id('foobarbaz'));
undef $db;
undef $primary_seq;

my (%h,$dna1,$dna2);
ok(tie(%h,'Bio::DB::Fasta',$test_dbdir));
ok($h{'AW057146'});
ok($dna1 = $h{'AW057146:1,10'});
ok($dna2 = $h{'AW057146:10,1'});

my $revcom = reverse $dna1;
$revcom =~ tr/gatcGATC/ctagCTAG/;
ok($dna2 eq $revcom);


END {
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('Unable to complete DBFasta tests',1);
	}
	# test dir is cleaned up automagically by tempdir(CLEANUP => 1)
}
