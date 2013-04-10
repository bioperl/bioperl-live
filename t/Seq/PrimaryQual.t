# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 70);

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
ok my $qualobj = Bio::Seq::PrimaryQual->new(
    -qual             => $string_quals,
    -id               => 'QualityFragment-12',
    -accession_number => 'X78121',
);
is $qualobj->display_id, 'QualityFragment-12';
is $qualobj->accession_number, 'X78121';

my @q2 = split / /, $string_quals;
$qualobj = Bio::Seq::PrimaryQual->new(
    -qual             => \@q2,
    -primary_id       => 'chads primary_id',
    -desc             => 'chads desc',
    -accession_number => 'chads accession_number',
    -id               => 'chads id',
    -header           => 'chads header'
);

is $qualobj->primary_id, 'chads primary_id';
isa_ok $qualobj->qual(), 'ARRAY';

my $newqualstring = "50 90 1000 20 12 0 0";

ok $qualobj->qual($newqualstring);
is join(' ', @{$qualobj->qual()}), $newqualstring;

my @newqualarray = split / /,$newqualstring;
ok $qualobj->qual(\@newqualarray);
is join(' ', @{$qualobj->qual()}), $newqualstring;

is $qualobj->validate_qual($string_quals ), 1;
is $qualobj->validate_qual(""            ), 1;
is $qualobj->validate_qual("0"           ), 1;
is $qualobj->validate_qual(undef         ), 1;
is $qualobj->validate_qual("   "         ), 1;
is $qualobj->validate_qual("10 20 30 30" ), 1;
is $qualobj->validate_qual(" 20  9 5   " ), 1;
is $qualobj->validate_qual("+1 9.3 50e-1"), 1;
is $qualobj->validate_qual(" 4"          ), 1;
is $qualobj->validate_qual("chad"        ), 0;
is $qualobj->validate_qual("10 one"      ), 0;

ok $qualobj->qual("10 20 30 30");
ok $qualobj->qual("+1 9.3 50e-1");
throws_ok { $qualobj->qual("chad"); } qr/.+/;
throws_ok { $qualobj->validate_qual("chad", 1) } qr/.+/;

ok $qualobj->qual("4 10");
is $qualobj->length(), 2;

$qualobj->qual("10 20 30 40 50 40 30 20 10");
ok my @subquals = @{$qualobj->subqual(3,6);};
is @subquals, 4;
is "30 20 10", join(' ',@{$qualobj->subqual(7,9)});

throws_ok { $qualobj->subqual(-1,6); } qr/EX/;
ok $qualobj->subqual(1,6);
ok $qualobj->subqual(1,9);
throws_ok { $qualobj->subqual(9,1); } qr/EX/;


is $qualobj->display_id(), "chads id";
is $qualobj->display_id("chads new display_id"), "chads new display_id";
is $qualobj->display_id(), "chads new display_id";

is $qualobj->accession_number(), "chads accession_number";
is $qualobj->accession_number("chads new accession_number"), "chads new accession_number";
is $qualobj->accession_number(), "chads new accession_number";
is $qualobj->primary_id(), "chads primary_id";
is $qualobj->primary_id("chads new primary_id"), "chads new primary_id";
is $qualobj->primary_id(), "chads new primary_id";

is $qualobj->desc(), "chads desc";
is $qualobj->desc("chads new desc"), "chads new desc";
is $qualobj->desc(), "chads new desc";
is $qualobj->display_id(), "chads new display_id";
is $qualobj->display_id("chads new id"), ("chads new id");
is $qualobj->display_id(), "chads new id";

is $qualobj->header(), "chads header";

ok my $in_qual  = Bio::SeqIO->new(
    -file    => test_input_file('qualfile.qual'),
    -format  => 'qual',
    -verbose => $verbose,
);
ok my $pq = $in_qual->next_seq();
is $pq->qual()->[99] , '39'; # spot check boundary
is $pq->qual()->[100], '39'; # spot check boundary

ok my $out_qual = Bio::SeqIO->new(
    -file    => ">".test_output_file(),
    -format  => 'qual',
    -verbose => $verbose,
);
ok $out_qual->write_seq(-source => $pq);

ok my $swq545 = Bio::Seq::Quality->new (
    -seq  => "ATA",
    -qual => $pq
);
ok $out_qual->write_seq(-source => $swq545);

ok $in_qual = Bio::SeqIO->new(
    -file    => test_input_file('qualfile.qual'),
    -format  => 'qual',
    -verbose => $verbose,
);

ok my $out_qual2 = Bio::SeqIO->new(
    -file    => ">".test_output_file(),
    -format  => 'qual',
    -verbose => $verbose,
);

while ( my $batch_qual = $in_qual->next_seq() ) {
    ok $out_qual2->write_seq(-source => $batch_qual);
}

