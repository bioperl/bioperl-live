# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#

use strict;

BEGIN {
    use vars qw($DEBUG);
    $DEBUG = $ENV{'BIOPERLDEBUG'};

    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 16;
}

END {
    unlink qw(write_scf.scf write_scf_no_sequence.scf 
	      write_scf_no_qualities.scf);
}


require 'dumpvar.pl';

print("Checking if the Bio::SeqIO::scf module could be used, even though it shouldn't be directly use'd...\n") if( $DEBUG);
        # test 1
use Bio::SeqIO::scf;
ok(1);

print("Checking to see if SeqWithQuality objects can be created from an scf file...\n") if( $DEBUG );
	# my $in_scf = Bio::SeqIO->new(-file => "<t/data/chad100.scf" , '-format' => 'csmscf');
my $in_scf = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data",
							       "chad100.scf"),
			     '-format' => 'scf',
			     '-verbose' => $DEBUG || 0);
ok(1);

my $swq = $in_scf->next_seq();

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
ok ($@ =~ /that wasn\'t A,T,G, or C/);

print("Now trying to request a valid base channel...\n");
my $a_channel = $in_scf->get_trace("a");
ok (length($a_channel) > 10);
my $c_channel = $in_scf->get_trace("c");
ok (length($c_channel) > 10);
my $g_channel = $in_scf->get_trace("g");
ok (length($g_channel) > 10);
my $t_channel = $in_scf->get_trace("t");
ok (length($t_channel) > 10);

print("Now checking to see if peak indices can be pulled for a v2 scf...\n");
my @indices = @{$in_scf->get_peak_indices()};
     # print ("@indices\n");
ok (scalar(@indices) == 761);

print("Now checking to see if peak indices can be pulled for a version 3 scf...\n");
my $in_scf_v3 = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data",
							       "version3.scf"),
			     '-format' => 'scf',
			     '-verbose' => $DEBUG || 0);

my $v3 = $in_scf_v3->next_seq();
@indices = @{$in_scf_v3->get_peak_indices()};
     # print("The peak indices for a v3 scf were there: @indices\n");
ok (scalar(@indices) == 1106);

print("Now getting the header from a scf to determine how many sample points there were (an example application)...\n");
my %header = %{$in_scf_v3->get_header()};
ok (%header->{bases} == 1106);
ok (%header->{samples} == 14107);


	# everything ok? <deep breath> ok, now we test the writing components
	# 1. try to create an empty file
print("Trying to create a new scf file from the existing object (from above)...\n");

my $out_scf = Bio::SeqIO->new('-file' => ">write_scf.scf",
			      '-format' => 'scf');
$out_scf->write_seq(-SeqWithQuality	=>	$swq,
		    -MACH		=>	'CSM sequence-o-matic 5000',
		    -TPSW		=>	'trace processing software',
		    -BCSW		=>	'basecalling software',
		    -DATF		=>	'AM_Version=2.00',
		    -DATN		=>	'a22c.alf',
		    -CONV		=>	'Bioperl-scf.pm');
ok( -e "write_scf.scf");

print("Trying to create an scf using null qualities.\n");

$swq = Bio::Seq::SeqWithQuality->new(-seq=>'ATCGTACGTACGTC',
				     -qual=>"");

$out_scf = Bio::SeqIO->new('-file' => ">write_scf_no_qualities.scf",
			   '-format' => 'scf');
$out_scf->write_seq(	-SeqWithQuality	=>	$swq,
			-MACH		=>	'CSM sequence-o-matic 5000',
			-TPSW		=>	'trace processing software',
			-BCSW		=>	'basecalling software',
			-DATF		=>	'AM_Version=2.00',
			-DATN		=>	'a22c.alf',
			-CONV		=>	'Bioperl-scf.pm');

	print("Trying to create an scf using null sequence but with qualities.\n");

$out_scf = Bio::SeqIO->new('-verbose' => 1,
			   '-file' => ">write_scf_no_sequence.scf",
			   '-format' => 'scf');

$swq = Bio::Seq::SeqWithQuality->new(-seq=>'',
				     -qual=>"10 20 30 40 50 20 10 30 40 50",
				     -alphabet=>'dna');

$out_scf->write_seq(	-SeqWithQuality	=>	$swq,
			-MACH		=>	'CSM sequence-o-matic 5000',
			-TPSW		=>	'trace processing software',
			-BCSW		=>	'basecalling software',
			-DATF		=>	'AM_Version=2.00',
			-DATN		=>	'a22c.alf',
			-CONV		=>	'Bioperl-scf.pm');

