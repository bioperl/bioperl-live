# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;
use File::Spec;
use File::Temp;

BEGIN {
    use Bio::Root::Test;

    test_begin(-tests => 5,
               -requires_module => 'IO::String');

    use_ok('Cwd');
    use_ok('Bio::SearchIO');
    use_ok('Bio::Index::Blast');
}

sub test_results {
    my $index = shift;
    my @ids = @_;

    foreach my $id (@ids) {
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
}


subtest 'BLASTP' => sub {
    plan tests => 2 + 2*5;

    my $dir = File::Temp->newdir();
    my $basename = 'Wibbl';
    my $filepath = File::Spec->catfile($dir, $basename);

    my $index = Bio::Index::Blast->new(-filename => $filepath,
                                       -write_flag => 1);
    ok($index);

    $index->make_index(test_input_file('multi_blast.bls'));
    ($index->dbm_package eq 'SDBM_File') ?
        (ok(-e "$filepath.pag" && -e "$filepath.dir")) :
        (ok(-e "$filepath"));

    test_results($index, qw(CATH_RAT PAPA_CARPA));
};

subtest 'RPS-BLAST' => sub {
    plan tests => 1 + 2*5;

    my $dir = File::Temp->newdir();
    my $basename = 'Wibbl.index';
    my $filepath = File::Spec->catfile($dir, $basename);

    my $index = Bio::Index::Blast->new(-filename => $filepath,
                                       -write_flag => 1);
    ok($index);

    $index->make_index(test_input_file('rpsblast.bls'));

    test_results($index, qw(orf20 orf40));
};
