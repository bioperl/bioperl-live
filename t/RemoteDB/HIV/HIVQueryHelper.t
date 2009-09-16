#-*-perl-*-
# testing HIVQueryHelper.pm and lanl-schema.xml indirectly
# $Id: HIVQueryHelper.t 231 2008-12-11 14:32:00Z maj $

use strict;
use warnings;

BEGIN {
    use Bio::Root::Test;
    test_begin(
	-tests => 40,
	-requires_modules => [qw( Bio::Root::Root XML::Simple)]
	);
    use_ok('Bio::DB::HIV::HIVQueryHelper');
}


# lanl-schema.xml characteristics as of $Date: 2008-12-11 08:05:24 -0500 (Thu, 11 Dec 2008) $
my ($naliases, $nfields, $ntables) = (165, 97, 14);
my ($Q, $r, $q);
# object tests
isa_ok(new HIVSchema(), "HIVSchema");
isa_ok($Q = new QRY(), "QRY");
isa_ok($r = new R(), "R");
isa_ok($q = new Q(), "Q");

#HIVSchema tests
my $tobj;
ok( $tobj = new HIVSchema(Bio::Root::IO->catfile(qw(Bio DB HIV lanl-schema.xml))), "schema load");

# methods
can_ok( $tobj, qw (
                   tables
                   columns
                   fields
                   options
                   aliases
                   ankh
                   tablepart
                   tbl
                   columnpart
                   col
                   primarykey
                   pk
                   foreignkey
                   fk
                   foreigntable
                   ftbl
                   loadSchema
                   _sfieldh
                   loadHIVSchema
                   ));

# lanl-schema.xml tests
is( scalar $tobj->fields, $nfields, "fields complete");
is( scalar $tobj->tables, $ntables, "tables complete");
is( scalar $tobj->aliases, $naliases, "aliases complete");
my ($tbl, $fld, $col, $als);
ok( ($fld) = grep /sequenceentry.se_sequence/, $tobj->fields, "test field present");
ok( $tbl = $tobj->tablepart($fld), "test field syntax ok");
ok( $col = $tobj->columnpart($fld), "test field syntax ok");
ok( $als = $tobj->aliases($fld), "test alias by field name");
is( $tobj->primarykey($tbl), 'sequenceentry.se_id', "correct primary key for SequenceEntry");
is( scalar $tobj->foreignkey('author'), 2, "correct number of foreign keys for AUthor");
is( $tobj->foreigntable(($tobj->foreignkey('author'))[1]), 'publication', "correct foreign table for au_pub_id");
is_deeply( {$tobj->ankh('seq_sample.ssam_second_receptor')}, 
{'seq_sample.ssam_second_receptor' => {'ankey'=>'coreceptor','antype'=>'Virus'}}, "correct annotation key hash");


#Query parser tests

#methods
can_ok($Q, qw(
              requests
              put_requests
              isnull
              A
              len
              clone
              ));
can_ok($r, qw(
              atoms
              put_atoms
              del_atoms
              fields
              isnull
              A
              len
              clone
              ));
can_ok($q, qw(
              fld
              dta
              A
              clone
              ));
# operations tests
ok($Q->isnull, "null QRY");
ok($r->isnull, "null R (request object)");
ok($q->isnull, "null Q (atomic query object)");

my ($R1, $R2);
ok($R1 = new R( new Q('X', 'a b c') ), "R obj create and init (1)");
ok($R2 = new R( new Q('X', 'a'), new Q('Y', 'u v w') ), "R obj create and init (2)");
ok(R::In($R2, $R1), "R::In");
ok(!R::In($R1, $R2), "!R::In");
ok(R::Eq($R1, $R1->clone), "R::Eq");
my ($Q1, $Q2, $Q3);
ok($Q1 = new QRY( $R1 ), "QRY obj create and init (1)");
ok($Q2 = new QRY( $R2 ), "QRY obj create and init (2)");
ok($Q3 = new QRY( new R( new Q('X', 'a'), new Q('Y', 'w v u'))), "QRY obj create and init (3)");

ok(QRY::Eq($Q1 | $Q1, $Q1), "QRY overload |");
ok(QRY::Eq($Q1 & $Q1,$Q1), "QRY overload &");
ok(QRY::Eq($Q1 & $Q2, $Q3), "QRY nontrivial &"); 

# parse and make query tests

my $pt;
ok($pt=QRY::_parse_q("('odds bodkins', a)[X] m[Y] u[Z] OR 'b'[X] {A B [C] [D]}"), "parse: ('odds bodkins', a)[X] m[Y] u[Z] OR 'b'[X] {A B [C] [D]} ");
is(scalar QRY::_make_q($pt), 2, "make: 2 queries returned");
is_deeply((QRY::_make_q($pt))[0]->{annot}, ['A','B','C','D'], "{annotation fields} parsed correctly");
ok($pt=QRY::_parse_q("('odds bodkins', a)[X] m[Y] u[Z] AND b[X] {A B [C] [D]}"), "parse: ('odds bodkins', a)[X] m[Y] u[Z] AND b[X] {A B [C] [D]} ");
is_deeply((QRY::_make_q($pt))[0],{}, "above query is null");


