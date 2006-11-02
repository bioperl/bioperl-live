# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#

use strict;

BEGIN {
    use vars qw($DEBUG $verbose);
    $DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
    $verbose = $DEBUG ? 0 : -1;
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 15;
}

END {
     unlink qw(
               write_scf.scf
               write_scf_synthetic_traces.scf 
	       write_scf_subtrace.scf
               write_scf_version2.scf
     );
}

use Dumpvalue();

my $dumper = new Dumpvalue();
$dumper->veryCompact(1) if $DEBUG;

use Bio::SeqIO::scf;
use Bio::Seq::SequenceTrace;

my $in_scf = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data",
							       "chad100.scf"),
			     '-format' => 'scf',
			     '-verbose' => $verbose);
ok($in_scf);

my $swq = $in_scf->next_seq();

ok (ref($swq) eq "Bio::Seq::SequenceTrace");

ok (length($swq->seq())>10);
my $qualities = join(' ',@{$swq->qual()});


ok (length($qualities)>10);
my $id = $swq->id();
ok ($swq->id() eq "ML4942R");

my $a_channel = $swq->trace("a");
ok (scalar(@$a_channel) > 10);
my $c_channel = $swq->trace("c");
ok (length($c_channel) > 10);
my $g_channel = $swq->trace("g");
ok (length($g_channel) > 10);
my $t_channel = $swq->trace("t");
ok (length($t_channel) > 10);

my $ref = $swq->peak_indices();
my @indices = @$ref;
ok (scalar(@indices), 761);

warn("Now checking version3...\n") if $DEBUG;
my $in_scf_v3 = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile
				("t","data",
				 "version3.scf"),
			     '-format' => 'scf',
			     '-verbose' => $verbose);

my $v3 = $in_scf_v3->next_seq();
my $ind = $v3->peak_indices();
my @ff = @$ind;
# print("ind is $ind which is @ff\n");
@indices = @{$v3->peak_indices()};
ok (scalar(@indices) == 1106);

# $dumper->dumpValue($v3);


my %header = %{$in_scf_v3->get_header()};
ok $header{bases}, 1106;
ok $header{samples},  14107;

my %comments = %{$in_scf_v3->get_comments()};
ok $comments{'NAME'}, 'IIABP1D4373';
ok $comments{'CONV'},  'phred version=0.990722.h';

warn("Now testing the _writing_ of scfs\n") if $DEBUG;

my $out_scf = Bio::SeqIO->new('-file' => ">write_scf.scf",
			      '-format' => 'scf',
			      '-verbose' => $verbose);
exit;	# the new way
$out_scf->write_seq(
	-target	=>	$v3,
	-MACH		=>	'CSM sequence-o-matic 5000',
	-TPSW		=>	'trace processing software',
	-BCSW		=>	'basecalling software',
	-DATF		=>	'AM_Version=2.00',
	-DATN		=>	'a22c.alf',
	-CONV		=>	'Bioperl-scf.pm');



ok( -e "write_scf.scf" && ! -z "write_scf.scf" );


$out_scf = Bio::SeqIO->new('-verbose' => 1,
			   '-file' => ">write_scf_synthetic_traces.scf",
			   '-format' => 'scf');

$swq = Bio::Seq::Quality->new(-seq=>'ATCGATCGAA',
				     -qual=>"10 20 30 40 50 20 10 30 40 50",
				     -alphabet=>'dna');

my $trace = Bio::Seq::SequenceTrace->new(
                         -swq =>   $swq);

$out_scf->write_seq(	
               -target	=>	$trace,
			-MACH		=>	'CSM sequence-o-matic 5000',
			-TPSW		=>	'trace processing software',
			-BCSW		=>	'basecalling software',
			-DATF		=>	'AM_Version=2.00',
			-DATN		=>	'a22c.alf',
			-CONV		=>	'Bioperl-scf.pm');

warn("Trying to write an scf with a subset of a real scf...\n") if $DEBUG;
$out_scf = Bio::SeqIO->new('-verbose' => 1,
			   '-file' => ">write_scf_subtrace.scf",
			   '-format' => 'scf');
$in_scf_v3 = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile
				("t","data",
				 "version3.scf"),
			     '-format' => 'scf',
			     '-verbose' => $verbose);
$v3 = $in_scf_v3->next_seq();
# print("The full trace object is as follows:\n");
my $sub_v3 = $v3->sub_trace_object(5,50);

warn("The subtrace object is this:\n") if $DEBUG;
$dumper->dumpValue($sub_v3) if $DEBUG;

$out_scf->write_seq(
          -target   =>   $sub_v3
);


my $in_scf_v2 = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile
				("t","data",
				 "version2.scf"),
			     '-format' => 'scf',
			     '-verbose' => $verbose);
$v3 = $in_scf_v2->next_seq();
ok($v3);

$out_scf = Bio::SeqIO->new(-file   =>   ">write_scf_version2.scf",
                           -format =>   "scf");
$out_scf->write_seq( -target  => $v3,
                     -version => 2 );

# now some version 2 things.
