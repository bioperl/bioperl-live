# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);


BEGIN { 
    eval { require Test::More; };
    if( $@ ) {
		use lib 't/lib';
    }
    use Test::More;

    $NUMTESTS = 19;
    eval { 
		require DB_File; 
		require Bio::DB::Flat; 
		require Bio::Root::IO; 
    };
    if( $@ ) {
		plan skip_all => "DB_File not loaded. This means flat.t test cannot be executed. Skipping";
	} else {
	    plan tests => $NUMTESTS;	
    }
	use_ok('Bio::Root::IO');
	use_ok('Bio::DB::Flat');
	use_ok('Cwd');
}

my $testnum;
my $verbose = 0;

## End of black magic.

#First of all we need to create an flat db

my $cd = cwd();
my $tmpdir = Bio::Root::IO->catfile($cd,qw(t tmp));
&maketmpdir();
my $db = Bio::DB::Flat->new(-directory  => $tmpdir,
                            -index      => 'bdb',
			    -dbname     => 'mydb',
			    -format     => 'fasta',
			    -verbose    => $verbose,
                 	    -write_flag => 1
                            );
ok($db);
my $dir = Bio::Root::IO->catfile($cd,qw(t data AAC12660.fa));
my $result = $db->build_index(glob($dir));
ok($result);

#Now let's get the sequence out again
my $seq = $db->get_Seq_by_id('AAC12660');
ok($seq);
is($seq->length,504);
undef $db;
&cleanup();
&maketmpdir();
$db = Bio::DB::Flat->new(-directory  => $tmpdir,
                         -index      => 'bdb',
                         -format     => 'embl',
			 -dbname     => 'myembl',
                         -verbose    => $verbose,
                         -write_flag => 1
			 );

$dir= Bio::Root::IO->catfile($cd,qw(t data factor7.embl));

$result = $db->build_index(glob($dir));
ok($result);
$seq = $db->get_Seq_by_id('HSCFVII');
ok($seq);
is($seq->length,12850);

# deal with wantarray conditions
$seq = $db->get_Seq_by_acc('J02933');
ok($seq && ref($seq));
is($seq->length,12850);


undef $db;

&cleanup();
&maketmpdir();


$db = Bio::DB::Flat->new(-directory  => $tmpdir,
			 -index      => 'binarysearch',
			 -format     => 'fasta',
			 -dbname     => 'mybinfa',
			 -verbose    => $verbose,
			 -write_flag => 1
			 );

$dir= Bio::Root::IO->catfile($cd,qw(t data dbfa 1.fa));
$result = $db->build_index($dir);
ok($result);
$seq = $db->get_Seq_by_id('AW057119');
ok($seq);
is($seq->length,808);
undef $db;

&cleanup();
&maketmpdir();
$db = Bio::DB::Flat->new(-directory  => $tmpdir,
			 -index      => 'binarysearch',
			 -format     => 'swiss',
			 -dbname     => 'mybinswiss',
			 -verbose    => $verbose,
			 -write_flag => 1
			 );
$dir= Bio::Root::IO->catfile($cd,qw(t data swiss.dat));
$result = $db->build_index($dir);

ok($result);
$seq = $db->get_Seq_by_id('ACON_CAEEL');
ok($seq);
is($seq->length,788);

$seq = $db->get_Seq_by_id('ACON_CAEEL');
ok($seq && ref($seq));

undef $db;


&cleanup();

sub maketmpdir {
    mkdir ($tmpdir,0777);
}

sub cleanup {
    eval { 
	 Bio::Root::IO->rmtree($tmpdir) if( defined $tmpdir && -d $tmpdir);
    };
} 

END {
    &cleanup();
}
