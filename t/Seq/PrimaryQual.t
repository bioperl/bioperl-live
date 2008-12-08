# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 35);
	
    use_ok('Bio::SeqIO');
    use_ok('Bio::Seq::Quality');
    use_ok('Bio::Seq::PrimaryQual');
}

my $DEBUG = test_debug();
my $verbose = -1 unless $DEBUG;

# redirect STDERR to STDOUT
open (STDERR, ">&STDOUT");

my $string_quals = "10 20 30 40 50 40 30 20 10";
print("Quals are $string_quals\n") if($DEBUG);
my $qualobj = Bio::Seq::PrimaryQual->new(
					  '-qual' => $string_quals,
					  '-id'  => 'QualityFragment-12',
					  '-accession_number' => 'X78121',
					  );
ok($qualobj);
is($qualobj->display_id, 'QualityFragment-12');
is($qualobj->accession_number, 'X78121');

my @q2 = split/ /,$string_quals;
$qualobj = Bio::Seq::PrimaryQual->new
    ( '-qual'             => \@q2,
      '-primary_id'	     => 'chads primary_id',
      '-desc'		        => 'chads desc',
      '-accession_number' => 'chads accession_number',
      '-id'		           => 'chads id',
		'-header'           => 'chads header'
      );

is($qualobj->primary_id, 'chads primary_id');
my $rqual = $qualobj->qual();
is(ref($rqual),"ARRAY");

my $newqualstring = "50 90 1000 20 12 0 0";

$qualobj->qual($newqualstring);
my $retrieved_quality = $qualobj->qual();
my $retrieved_quality_string = join(' ', @$retrieved_quality);
is($retrieved_quality_string,$newqualstring);

my @newqualarray = split/ /,$newqualstring;
$qualobj->qual(\@newqualarray);
$retrieved_quality = $qualobj->qual();
$retrieved_quality_string = join(' ',@$retrieved_quality);
is($retrieved_quality_string,$newqualstring);

eval {
    $qualobj->qual("chad");
};
like($@, qr/not look healthy/);

eval { $qualobj->qual(""); };
ok(!$@);

eval { $qualobj->qual(" 4"); };
ok(!$@);

$qualobj->qual("4 10");

is($qualobj->length(),2 );

$qualobj->qual("10 20 30 40 50 40 30 20 10");
my @subquals = @{$qualobj->subqual(3,6);};
is(@subquals, 4);
     # chad, note to self, evaluate border conditions
is ("30 20 10", join(' ',@{$qualobj->subqual(7,9)}));


my @false_comparator = qw(30 40 70 40);
my @true_comparator = qw(30 40 50 40);
ok(!&compare_arrays(\@subquals,\@true_comparator));

eval { $qualobj->subqual(-1,6); };
like($@, qr/EX/ );
eval { $qualobj->subqual(1,6); };
ok(!$@);
eval { $qualobj->subqual(1,9); };
ok(!$@);
eval { $qualobj->subqual(9,1); };
like($@, qr/EX/ );


is($qualobj->display_id(), "chads id");
$qualobj->display_id("chads new display_id");
is($qualobj->display_id(), "chads new display_id");

is($qualobj->accession_number(), "chads accession_number");
$qualobj->accession_number("chads new accession_number");
is($qualobj->accession_number(), "chads new accession_number");
is($qualobj->primary_id(), "chads primary_id");
$qualobj->primary_id("chads new primary_id");
is($qualobj->primary_id(), "chads new primary_id");

is($qualobj->desc(), "chads desc");
$qualobj->desc("chads new desc");
is($qualobj->desc(), "chads new desc");
is($qualobj->display_id(), "chads new display_id");
$qualobj->display_id("chads new id");
is($qualobj->display_id(), "chads new id");

is($qualobj->header(), "chads header");

my $in_qual  = Bio::SeqIO->new(-file => test_input_file('qualfile.qual') ,
			       '-format' => 'qual',
			       '-verbose' => $verbose);
ok($in_qual);
my $pq = $in_qual->next_seq();
is($pq->qual()->[99], '39'); # spot check boundary
is($pq->qual()->[100], '39'); # spot check boundary

my $out_qual = Bio::SeqIO->new('-file'    => ">".test_output_file(),
                               '-format'  => 'qual',
                               '-verbose' => $verbose);
$out_qual->write_seq(-source	=>	$pq);

my $swq545 = Bio::Seq::Quality->new (	-seq	=>	"ATA",
                                        -qual	=>	$pq
                                    );
$out_qual->write_seq(-source	=>	$swq545);

$in_qual = Bio::SeqIO->new('-file' => test_input_file('qualfile.qual') , 
			   '-format' => 'qual',
			   '-verbose' => $verbose);

my $out_qual2 = Bio::SeqIO->new('-file' => ">".test_output_file(),
				'-format'  => 'qual',
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
