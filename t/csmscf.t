# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#

use ExtUtils::testlib;
use strict;
require 'dumpvar.pl';

BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 11;
}


print("Checking if the Bio::SeqIO::csmscf module could be used, even though it shouldn't be directly use'd...\n");
        # test 1
use Bio::SeqIO::csmscf;
ok(1);

print("Checking to see if SeqWithQuality objects can be created from an scf file...\n");
	# my $in_scf = Bio::SeqIO->new(-file => "<t/data/chad100.scf" , '-format' => 'csmscf');
my $in_scf = Bio::SeqIO->new(-file => "<t/data/chad100.scf" , '-format' => 'csmscf');
ok(1);

my $swq = $in_scf->next_scf();

print("Checking to see that SeqIO::scf returned the right kind of object (SeqWithQuality)...\n");
ok (ref($swq) eq "Bio::Seq::SeqWithQuality");

print("Checking to see if the SeqWithQuality object contains the right stuff.\n");
print("sequence :  ".$swq->seq()."\n");
ok (length($swq->seq())>10);
my $qualities = join(' ',@{$swq->qual()});
print("qualities: $qualities\n");
ok (length($qualities)>10);
print("id       : ".$swq->id()."\n");
ok ($swq->id());

print("Now checking to see that you can retrieve traces for the individual colour channels in the scf...\n");
print("First, trying to retrieve a base channel that doesn't exist...\n");
eval { $in_scf->get_trace("h"); };
ok ($@ =~ /that wasn't A,T,G, or C/);
print("Now trying to request a valid base channel...\n");
my $a_channel = $in_scf->get_trace("a");
ok (length($a_channel) > 10);
my $c_channel = $in_scf->get_trace("c");
ok (length($c_channel) > 10);
my $g_channel = $in_scf->get_trace("g");
ok (length($g_channel) > 10);
my $t_channel = $in_scf->get_trace("t");
ok (length($t_channel) > 10);

	# everything ok? <deep breath> ok, now we test the writing components
	# 1. try to create an empty file
print("Trying to create a new scf file from the existing object (from above)...\n");

my $out_scf = Bio::SeqIO->new(-file => ">t/data/csmscf.scf" , '-format' => 'csmscf');
$out_scf->write_scf(	-SeqWithQuality	=>	$swq,
			-MACH		=>	'CSM sequence-o-matic 5000',
			-TPSW		=>	'trace processing software',
			-BCSW		=>	'basecalling software',
			-DATF		=>	'AM_Version=2.00',
			-DATN		=>	'a22c.alf',
			-CONV		=>	'Bioperl-scf.pm');




# dumpValue($in_scf);


