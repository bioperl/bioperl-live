# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib 't/lib';
	use BioperlTest;
	
	test_begin(-tests => 15,
			   -requires_module => 'IO::String');
    
	use_ok('Cwd');
	use_ok('Bio::SearchIO');
	use_ok('Bio::Index::Blast');
}

END {  unlink qw( Wibbl Wibbl.pag Wibbl.dir ); }

my $index = Bio::Index::Blast->new(-filename => 'Wibbl',
				  -write_flag => 1);
ok($index);

$index->make_index(test_input_file('multi_blast.bls'));
($index->dbm_package eq 'SDBM_File') ? 
	(ok(-e "Wibbl.pag" && -e "Wibbl.dir")) :
	(ok(-e "Wibbl"));

foreach my $id ( qw(CATH_RAT PAPA_CARPA) ) {
	my $fh = $index->get_stream($id);
	ok($fh);
	ok( ! eof($fh) );
	my $report = Bio::SearchIO->new(-noclose => 1,
				   -format  => 'blast',
				   -fh      => $fh);
	my $result = $report->next_result;
	like($result->query_name, qr/$id/);
	ok( $result->next_hit);
	
	like( $index->fetch_report($id)->query_name, qr/$id/);
}
