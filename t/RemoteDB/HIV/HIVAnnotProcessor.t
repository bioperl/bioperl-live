# testing Bio::DB::HIVAnnotProcessor.pm
# $Id: HIVAnnotProcessor.t 231 2008-12-11 14:32:00Z maj $
use strict;
use warnings;

BEGIN {
    use Bio::Root::Test;
    test_begin(
	-tests => 11,
	);
    use_ok('Bio::Seq');
    use_ok('Bio::SeqIO');
    use_ok('Bio::DB::HIV::HIVAnnotProcessor');
}

my $tobj = new Bio::DB::HIV::HIVAnnotProcessor();

#object tests
isa_ok($tobj, 'Bio::DB::HIV::HIVAnnotProcessor');

#compliance tests
isa_ok($tobj, 'Bio::Root::Root');
can_ok($tobj, qw( source_stream next_seq write_seq ));

#methods
can_ok($tobj, qw( hiv_query ));

#exception tests
throws_ok {$tobj->hiv_query(bless({},"narb"))} qr/BadParameter/, "bad type set exception";

#stream tests
my $fas = test_output_file();
open my $FAS, '>', $fas or die "Could not write file '$fas': $!\n";
print $FAS ">goob\natcg\n";
close $FAS;
ok( $tobj->source_stream(new Bio::SeqIO(-file=>$fas, -format=>'fasta')), "attach stream");
throws_ok {$tobj->write_seq(new Bio::Seq(-sequence=>"atcg"))} qr/IOException/, "write exception";
ok( $tobj->next_seq, "access stream");
