# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#

use strict;

BEGIN {
    use vars qw($DEBUG $verbose);
    $DEBUG = $ENV{'BIOPERLDEBUG'};
    $verbose = $DEBUG ? 0 : -1;
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

# test 1
use Bio::SeqIO::scf;


my $in_scf = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data",
							       "chad100.scf"),
			     '-format' => 'scf',
			     '-verbose' => $verbose);
ok($in_scf);

my $swq = $in_scf->next_seq();

ok (ref($swq) eq "Bio::Seq::SeqWithQuality");

ok (length($swq->seq())>10);
my $qualities = join(' ',@{$swq->qual()});


ok (length($qualities)>10);
ok ($swq->id());

eval { $in_scf->get_trace("h"); };
ok ($@, qr/that wasn\'t A,T,G, or C/);

my $a_channel = $in_scf->get_trace("a");
ok (length($a_channel) > 10);
my $c_channel = $in_scf->get_trace("c");
ok (length($c_channel) > 10);
my $g_channel = $in_scf->get_trace("g");
ok (length($g_channel) > 10);
my $t_channel = $in_scf->get_trace("t");
ok (length($t_channel) > 10);

my @indices = @{$in_scf->get_peak_indices()};
ok (scalar(@indices), 761);

my $in_scf_v3 = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile
				("t","data",
				 "version3.scf"),
			     '-format' => 'scf',
			     '-verbose' => $verbose);

my $v3 = $in_scf_v3->next_seq();
@indices = @{$in_scf_v3->get_peak_indices()};
ok (scalar(@indices) == 1106);

my %header = %{$in_scf_v3->get_header()};
ok (%header->{bases} == 1106);
ok (%header->{samples} == 14107);


my $out_scf = Bio::SeqIO->new('-file' => ">write_scf.scf",
			      '-format' => 'scf',
			      '-verbose' => $verbose);
$out_scf->write_seq(-SeqWithQuality	=>	$swq,
		    -MACH		=>	'CSM sequence-o-matic 5000',
		    -TPSW		=>	'trace processing software',
		    -BCSW		=>	'basecalling software',
		    -DATF		=>	'AM_Version=2.00',
		    -DATN		=>	'a22c.alf',
		    -CONV		=>	'Bioperl-scf.pm');
ok( -e "write_scf.scf" && ! -z "write_scf.scf" );

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
