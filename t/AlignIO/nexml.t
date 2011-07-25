#-*-perl-*-
# $Id$
use strict;

use Bio::Root::Test;
test_begin(-tests => 125,  -requires_modules => [qw(Bio::Phylo)]);
use_ok('Bio::AlignIO::nexml'); # checks that your module is there and loads ok
diag("WARNING: NeXML parsing for NeXML v0.9 is currently very experimental support");
SKIP: {
    skip("NeXML parsing for NeXML v0.9 is currently very experimental support", 124);

    #Read in Data
    ok( my $inAlnStream = Bio::AlignIO->new(-file => test_input_file("nexml","characters.nexml.xml"),
                                            -format => 'nexml'), 'make stream');
        
    #Read in aln objects
    ok( my $aln_obj = $inAlnStream->next_aln(), 'nexml matrix to aln' );
    isa_ok($aln_obj, 'Bio::SimpleAlign', 'obj ok');
    is ($aln_obj->id,   'DNA sequences', 'aln id');
    my $num =0;
    my @expected_seqs = ('ACGCTCGCATCGCATC', 'ACGCTCGCATCGCATC', 'ACGCTCGCATCGCATC');
    #checking sequence objects
    foreach my $seq_obj ($aln_obj->each_seq()) {
        $num++;
        
        is( $seq_obj->alphabet, 'dna', "alphabet" );
        
        TODO: {
            local $TODO = 'primary/display_id broken with NeXML 0.9';
            is( $seq_obj->display_id, "DNA sequences.row_$num", "display_id");
        }
        is( $seq_obj->seq, $expected_seqs[$num-1], "sequence correct");
    }
        
    #Write Data
    diag('Begin tests for write/read roundtrip');
    my $outdata = test_output_file();
    ok( my $outAlnStream = Bio::AlignIO->new(-file => ">$outdata", -format => 'nexml'), 'out stream ok');;
    ok( $outAlnStream->write_aln($aln_obj), 'write nexml');
    close($outdata);
    
    
    #Read in the written file to test roundtrip
    ok( my $inAlnStream2 = Bio::AlignIO->new(-file => $outdata, -format => 'nexml'), 'reopen');
        
    #Read in aln objects
    ok( my $aln_obj2 = $inAlnStream2->next_aln(),'get aln (rt)' );
    isa_ok($aln_obj2, 'Bio::SimpleAlign', 'aln obj (rt)');
    is ($aln_obj2->id, 'DNA sequences', "aln id (rt)");
    $num =0;
    @expected_seqs = ('ACGCTCGCATCGCATC', 'ACGCTCGCATCGCATC', 'ACGCTCGCATCGCATC');
    #checking sequence objects
    foreach my $seq_obj ($aln_obj2->each_seq()) {
        $num++;
        
        is( $seq_obj->alphabet, 'dna', "alphabet (rt)" );
        is( $seq_obj->display_id, "DNA sequences.row_$num", "display_id (rt)");
        is( $seq_obj->seq, $expected_seqs[$num-1], "sequence (rt)");
    }
    #check taxa object
    my %expected_taxa = ('DNA sequences.row_1' => 'Homo sapiens',
                         'DNA sequences.row_2' => 'Pan paniscus',
                         'DNA sequences.row_3' => 'Pan troglodytes');
    my @feats = $aln_obj2->get_all_SeqFeatures();
    foreach my $feat (@feats) {
        if ($feat->has_tag('taxa_id')){
            is ( ($feat->get_tag_values('taxa_id'))[0], 'taxa1', 'taxa id ok' );
            is ( ($feat->get_tag_values('taxa_label'))[0], 'Primary taxa block', 'taxa label ok');
            is ( $feat->get_tag_values('taxon'), 5, 'Number of taxa ok')
        }
        else{
            my $seq_num = ($feat->get_tag_values('id'))[0];
            is ( ($feat->get_tag_values('taxon'))[0], $expected_taxa{$seq_num}, "$seq_num taxon ok" )
        }
    }

}