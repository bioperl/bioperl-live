# -*-Perl-*- Test Harness script for Bioperl
# $Id: Blast.t 16293 2009-10-27 20:03:02Z cjfields $

use strict;

BEGIN {
	use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 27,
		-requires_module => 'IO::String');
    
	use_ok('Cwd');
	use_ok('Bio::SearchIO');
	use_ok('Bio::Index::BlastTable');
}

                       # -m 9                -m 8
my @test_cases = qw(multi.blast.m9      multi.blast.m8);

for my $file (@test_cases) {
    my $index = Bio::Index::BlastTable->new(-filename => 'Wibbl',
                                        -write_flag => 1);
    ok($index);
    $index->id_parser(\&my_id_parser);
    $index->make_index(test_input_file($file));
    ($index->dbm_package eq 'SDBM_File') ? 
        (ok(-e "Wibbl.pag" && -e "Wibbl.dir")) :
        (ok(-e "Wibbl"));
    
    foreach my $id ( qw(SP130_MOUSE IKZF1_MOUSE) ) {
        my $fh = $index->get_stream($id);
        ok($fh);
        ok( ! eof($fh) );
        my $report = Bio::SearchIO->new(-noclose => 1,
                       -format  => 'blasttable',
                       -fh      => $fh);
        my $result = $report->next_result;
        like($result->query_name, qr/$id/);
        ok( $result->next_hit);
        
        like( $index->fetch_report($id)->query_name, qr/$id/);
    }
    # ActivePerl will not allow deletion if the tie-hash is still active
    $index->DESTROY;
    unlink qw( Wibbl Wibbl.pag Wibbl.dir Wibbl.index);
}

# test id_parser
sub my_id_parser {
    if ($_[0] =~ /^\S+\|(\S+)/) {
        return $1;
    } else {
        return;
    }
}
