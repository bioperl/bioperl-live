#-*-perl-*-
# $Id$

use strict;

use Bio::Root::Test;
test_begin(-tests => 126,
	    -requires_modules => [qw(Bio::Phylo)]);
use_ok( 'Bio::PrimarySeq' );
use_ok('Bio::SeqIO::nexml'); # checks that your module is there and loads ok
diag("WARNING: NeXML parsing for NeXML v0.9 is currently very experimental support");
SKIP: {
    skip("NeXML parsing for NeXML v0.9 is currently very experimental support", 124);

#Read Data
ok( my $SeqStream = Bio::SeqIO->new(
                    -file => test_input_file('nexml', 'characters.nexml.xml'),
                    -format => 'nexml'),'stream ok');

#checking first sequence object
ok( my $seq_obj = $SeqStream->next_seq(), 'seq obj' );
isa_ok($seq_obj, 'Bio::Seq');
is( $seq_obj->alphabet, 'dna', "alphabet" );
TODO: {
    local $TODO = 'primary/display_id broken with NeXML 0.9';
    is( $seq_obj->primary_id, 'DNA sequences.seq_1', "primary_id");
    is( $seq_obj->display_id, 'DNA sequences.seq_1', "display_id");
}
is( $seq_obj->seq, 'ACGCTCGCATCGCATC', "sequence");
#check taxa
my %expected_taxa = ('Homo sapiens' => 1,
                     'Pan paniscus' => 1,
                     'Pan troglodytes' => 1,
                     'Gorilla gorilla' => 1,
                     'Pongo pygmaeus' => 1);
my $feat = ($seq_obj->get_SeqFeatures())[0];
is( ($feat->get_tag_values('taxa_id'))[0], 'taxa1', 'taxa id');
is( ($feat->get_tag_values('taxa_label'))[0], 'Primary taxa block', 'taxa label');
is( ($feat->get_tag_values('my_taxon'))[0], 'Homo sapiens', "taxon ok" );
my @taxa = $feat->get_tag_values('taxon');
is( @taxa, 5, 'number of taxa');
foreach my $taxon (@taxa) {
    ok( $expected_taxa{$taxon}, 'taxon ok')	
}

#checking second sequence object
ok( $seq_obj = $SeqStream->next_seq() );
is( $seq_obj->alphabet, 'dna', "alphabet" );
TODO: {
    local $TODO = 'primary/display_id broken with NeXML 0.9';
    is( $seq_obj->primary_id, 'DNA sequences.seq_2', "primary_id");
    is( $seq_obj->display_id, 'DNA sequences.seq_2', "display_id");
}
is( $seq_obj->seq, 'ACGCTCGCATCGCATC', "sequence");
$SeqStream->next_seq();
$SeqStream->next_seq();

#checking fifth sequence object
ok( $seq_obj = $SeqStream->next_seq() );
is( $seq_obj->alphabet, 'rna', "alphabet" );
TODO: {
    local $TODO = 'primary/display_id broken with NeXML 0.9';
    is( $seq_obj->primary_id, 'RNA sequences.seq_2', "primary_id");
    is( $seq_obj->display_id, 'RNA sequences.seq_2', "display_id defaults to primary");
}
is( $seq_obj->seq, 'ACGCUCGCAUCGCAUC', "sequence");

#Write Data
diag('Begin tests for writing seq files');
my $outdata = test_output_file();
ok( my $outSeqStream = Bio::SeqIO->new(-file => ">$outdata",
                                       -format => 'nexml'), 'out stream ok');
ok $outSeqStream->write_seq($seq_obj), 'write nexml seq';
close($outdata);

diag("write_seq support for NeXML 0.9 NYI");

#Read in the out file to test roundtrip
my $inSeqStream = Bio::SeqIO->new(-file => $outdata, -format => 'nexml');

#checking fifth sequence object
ok( my $seq_obj2 = $inSeqStream->next_seq() );
is( $seq_obj2->alphabet, 'rna', "alphabet" );
is( $seq_obj2->primary_id, 'RNA sequences.seq_2', "primary_id");
is( $seq_obj2->display_id, 'RNA sequences.seq_2', "display_id defaults to primary");
is( $seq_obj2->seq, 'ACGCUCGCAUCGCAUC', "sequence");

#check taxa
my $feat1 = ($seq_obj2->get_SeqFeatures())[0];
is( ($feat1->get_tag_values('taxa_id'))[0], 'taxa1', 'taxa id');
is( ($feat1->get_tag_values('taxa_label'))[0], 'Primary taxa block', 'taxa label');
is( ($feat1->get_tag_values('my_taxon'))[0], 'Pan paniscus', "taxon ok" );
my @taxa2 = $feat1->get_tag_values('taxon');
is( @taxa2, 5, 'number of taxa');
foreach my $taxon (@taxa2) {
    ok( $expected_taxa{$taxon}, 'taxon ok')	
}

}
