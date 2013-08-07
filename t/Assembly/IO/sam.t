use strict;
use warnings;
my %ASSEMBLY_TESTS;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin( -tests => 441+5,
                -requires_modules => [
                    'DB_File',
                    'Bio::DB::Sam',
                    'Bio::Tools::Run::Samtools',
                   ],
               );

    use_ok('Bio::Seq');
    use_ok('Bio::LocatableSeq');
    use_ok('Bio::Seq::Quality');
    use_ok('Bio::Assembly::IO');
    use_ok('Bio::Assembly::Singlet');
}

use Bio::Root::IO;

my ($aio, $assembly, @contig_seq_ids, @singlet_ids, @contig_ids, @all_seq_ids);
my $file = 'test.bam';
my $refdb = 'test.ref.fas';
ok $aio = Bio::Assembly::IO->new( -file => test_input_file($file),
                                  -refdb => test_input_file($refdb),
                                  -format => 'sam' ), "init sam IO object";
isa_ok($aio, 'Bio::Assembly::IO');
$aio->_current_refseq_id( ($aio->sam->seq_ids)[0] ); # kludge

while (my $contig = $aio->next_contig) {
    isa_ok($contig, 'Bio::Assembly::Contig');
}
ok $aio = Bio::Assembly::IO->new( -file => test_input_file($file),
                                  -refdb => test_input_file($refdb),
                                  -format => 'sam' ),"reopen";
ok $assembly = $aio->next_assembly, "get sam assy";
cmp_ok( $assembly->get_nof_contigs, '>=', 21, "got all contigs");


ok(@contig_seq_ids = $assembly->get_contig_seq_ids, "get_contig_seq_ids");
is(@contig_seq_ids, 334);
# trashing these for now; not much a test really anyway/maj
# for my $contig_seq_id (@contig_seq_ids) {
#     ok ($contig_seq_id =~ m/^SRR/i);
# }

ok(@contig_ids = $assembly->get_contig_ids, "get_contig_ids");
cmp_ok(@contig_ids, '>=', 21);
my $correct_ids = grep {$_ =~ m/sam_assy/i} @contig_ids;
is($correct_ids, scalar(@contig_ids));

ok(@singlet_ids = $assembly->get_singlet_ids, "get_singlet_ids");
is(@singlet_ids, 35);
# trashing these/maj
# for my $singlet_id (@singlet_ids) {
#     ok ($singlet_id =~ m/^SRR/i);
# }

ok(@all_seq_ids = $assembly->get_all_seq_ids, "get_all_seq_ids");
for my $seq_id (@all_seq_ids) {
    ok ($seq_id =~ m/^SRR/i);
}
is(@all_seq_ids, 369);

for my $f (qw(test.bam.bai test.ref.fas.fai)) {
    my $path = Bio::Root::IO->catfile('t', 'data', $f);
    unlink $path if -e $path
}
