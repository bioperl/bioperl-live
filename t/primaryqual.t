# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#
# modeled after the t/Allele.t test script

use strict;
use vars qw($DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'};
my $verbose = $DEBUG ? 0 : -1;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 26;
}

END { 
    unlink qw(batch_write_qual.qual write_qual.qual);
	
}
# redirect STDERR to STDOUT
#open (STDERR, ">&STDOUT");
use Bio::Root::IO;
use Bio::SeqIO;

use Bio::Seq::PrimaryQual;
ok(1);
if( $DEBUG ) {
    print("Checking Bio::Seq::PrimaryQual methods...\n") ; 
    print("Note: These are more for Chad's benefit then for yours because I run these every time I mangle the modules to make sure that things did not break. So there.\n");
    print("Ready? There are 7 tests.\n");
    print("1. Checking that new() works...\n");
    print("1a) New with quality values...\n");
}
my $string_quals = "10 20 30 40 50 40 30 20 10";
print("Quals are $string_quals\n") if($DEBUG); 
my $qualobj = Bio::Seq::PrimaryQual->new( '-qual' => $string_quals,
					  '-id'  => 'QualityFragment-12',
					  '-accession_number' => 'X78121',
					  );
ok($qualobj);
print("1b) New with a reference to an array of quality values...\n") if($DEBUG);
my @q2 = split/ /,$string_quals;
$qualobj = Bio::Seq::PrimaryQual->new
    ( '-qual' => \@q2,
      '-primary_id'	=>	'chads primary_id',			
      '-desc'		=>	'chads desc',
      '-accession_number' => 'chads accession_number',
      '-id'		=>	'chads id'
      );

ok($qualobj);

my $rqual = $qualobj->qual();
ok(ref($rqual), qr/ARRAY/);
my $newqualstring = "50 90 1000 20 12 0 0";

$qualobj->qual($newqualstring);
my $retrieved_quality = $qualobj->qual();
my $retrieved_quality_string = join(' ',@$retrieved_quality);

ok($newqualstring, $retrieved_quality_string);


my @newqualarray = split/ /,$newqualstring;
$qualobj->qual(\@newqualarray);
$retrieved_quality = $qualobj->qual();
$retrieved_quality_string = join(' ',@$retrieved_quality);
ok($newqualstring, $retrieved_quality_string);

eval {
    $qualobj->qual("chad");
};
ok($@, qr/not look healthy/);

eval { $qualobj->qual(""); };
ok(!$@);

eval { $qualobj->qual(" 4"); };
ok(!$@);

$qualobj->qual("10 20 30 40 50 40 30 20 10");
my @subquals = @{$qualobj->subqual(3,6);};

my @false_comparator = qw(30 40 70 40);
my @true_comparator = qw(30 40 50 40);
ok(!&compare_arrays(\@subquals,\@true_comparator));

eval { $qualobj->subqual(-1,6); };
ok($@, qr/EX/ );
eval { $qualobj->subqual(1,6); };
ok(!$@);
eval { $qualobj->subqual(1,9); };
ok(!$@);
eval { $qualobj->subqual(9,1); };
ok($@, qr/EX/ );

ok($qualobj->display_id(), "chads id");
$qualobj->display_id("chads new display_id");
ok($qualobj->display_id(), "chads new display_id");

ok($qualobj->accession_number(), "chads accession_number");

$qualobj->accession_number("chads new accession_number");
ok($qualobj->accession_number(), "chads new accession_number");

ok($qualobj->primary_id(), "chads primary_id");

$qualobj->primary_id("chads new primary_id");
ok($qualobj->primary_id(), "chads new primary_id");
ok($qualobj->desc(), "chads desc");
$qualobj->desc("chads new desc");
ok($qualobj->desc(), "chads new desc");
ok($qualobj->id(), "chads new display_id");
$qualobj->id("chads new id");
ok($qualobj->id(), "chads new id");

my $in_qual  = Bio::SeqIO->new(-file => "<" . Bio::Root::IO->catfile("t","data","qualfile.qual") , '-format' => 'qual');
ok($in_qual);
my $pq = $in_qual->next_seq();
ok($pq->qual()->[99], '39'); # spot check boundary
ok($pq->qual()->[100], '39'); # spot check boundary

my $out_qual = Bio::SeqIO->new('-file' => ">write_qual.qual",
			       '-format' => 'qual',
			       '-verbose' => $verbose);
$out_qual->write_seq(-source	=>	$pq);

my $swq545 = Bio::Seq::SeqWithQuality->new (	-seq	=>	"ATA",
						-qual	=>	$pq
					);
$out_qual->write_seq(-source	=>	$swq545);


$in_qual = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data","qualfile.qual") , 
			   '-format' => 'qual',
			   '-verbose' => $verbose);

my $out_qual2 = Bio::SeqIO->new('-file' => ">batch_write_qual.qual",
				'-format' => 'qual',
				'-verbose' => $verbose);

while ( my $batch_qual = $in_qual->next_seq() ) {
    $out_qual2->write_seq(-source	=>	$batch_qual);
}

sub display {
    if($DEBUG ) {
 	my @quals;
	print("I saw these in qualfile.qual:\n") ;
	while ( my $qual = $in_qual->next_seq() ) {
	    # ::dumpValue($qual);
	    print($qual->display_id()."\n");
	    @quals = @{$qual->qual()};
	    print("(".scalar(@quals).") quality values.\n");
	}
    }
}

# dumpValue($qualobj);

sub compare_arrays {
        my ($a1,$a2) = @_;
        return 1 if (scalar(@{$a1}) != scalar(@{$a2}));
        my ($v1,$v2,$diff,$curr);
        for ($curr=0;$curr<scalar(@{$a1});$curr++){
                return 1 if ($a1->[$curr] ne $a2->[$curr]);
        }
        return 0;
}
