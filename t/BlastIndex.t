# -*-Perl-*-
## $Id$

use strict;

BEGIN {
    use lib 't/lib';
	use BioperlTest;
	
	test_begin(-tests => 16,
			   -requires_modules => ['IO::String']);
    
	use_ok('Cwd');
	use_ok('Bio::SearchIO');
	use_ok('Bio::Index::Blast');
	use_ok('Bio::Root::IO');
}

END {  unlink qw( Wibbl Wibbl.pag Wibbl.dir ); }

my $index = Bio::Index::Blast->new(-filename => 'Wibbl',
				  -write_flag => 1);
ok($index);

$index->make_index(Bio::Root::IO->catfile(cwd,"t","data","multi_blast.bls"));
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

