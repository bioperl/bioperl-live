use strict;
use warnings;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    test_begin(-tests           => 115,
               -requires_module => 'DB_File');
    
    use_ok('Bio::Matrix::PSM::IO');
    use_ok('Bio::DB::TFBS');
    use_ok('Bio::DB::Taxonomy');
}

#*** need to test getting all ids of a certain kind, like $db->get_matrix_ids();
#    but hard to do without a complete tax dump

my $temp_dir = test_output_dir();
my $tax_db = Bio::DB::Taxonomy->new(-source => 'flatfile',
                                    -directory => $temp_dir,
                                    -nodesfile => test_input_file('taxdump', 'nodes.dmp'),
                                    -namesfile => test_input_file('taxdump', 'names.dmp'));

# test transfac pro (local flat files)
{
    ok my $db = Bio::DB::TFBS->new(-source => 'transfac_pro',
                                   -index_dir => $temp_dir,
                                   -dat_dir => test_input_file('transfac_pro'),
                                   -tax_db => $tax_db,
                                   -force => 1);
    
    # reference.dat
    {
        ok my ($ref_id) = $db->get_reference_ids(-pubmed => 16574738);
        is $ref_id, 'RE0047775';
        ok my $ref = $db->get_reference($ref_id);
        isa_ok $ref, 'Bio::Annotation::Reference';
        is $ref->primary_id, 16574738;
        is $ref->pubmed, $ref->primary_id;
        is $ref->database, 'PUBMED';
        is $ref->authors, '..Bet S . ,.u i rMeK ,,d. vWeWk KaS.ee.nyNk mJMMih. a, i P';
        is $ref->location, 'Mc (o0o.. 0n)lnir.do 2E:6l';
        is $ref->title, 'INDD VDGT C1AALEBEI.EIT IYIHLA6ITTE E ANV  ITSL MTRTANYE TM NISP TNBAUTPOIORSL I- NVTOD,MHIRRLINSDX TRPY NO CAELUAOA SNMMNT CED5CTH NII TERTOI2IMTVPEH3DSAI';
        
        my @sites = $db->get_site_ids(-reference => $ref_id);
        is join(' ', sort @sites), 'R19310 R19311 R19312 R19313 R19314 R19315 R19316';
        my @genes = $db->get_gene_ids(-reference => $ref_id);
        is "@genes", 'G036757';
        my @ref_ids = $db->get_reference_ids(-site => 'R19310');
        is "@ref_ids", $ref_id;
        @ref_ids = $db->get_reference_ids(-gene => 'G036757');
        is "@ref_ids", $ref_id;
        
        $ref_id = 'RE0047531';
        my @matrices = $db->get_matrix_ids(-reference => $ref_id);
        is join(' ', sort @matrices), 'M01123 M01124 M01125';
        my @factors = $db->get_factor_ids(-reference => $ref_id);
        like "@factors", qr/T08800/;
        @ref_ids = $db->get_reference_ids(-matrix => 'M01123');
        is join(' ', sort @ref_ids), "$ref_id RE0047626";
        @ref_ids = $db->get_reference_ids(-factor => 'T08800');
        is join(' ', sort @ref_ids), "$ref_id RE0047634 RE0047637 RE0047645";
		
		$ref_id = 'RE0023998';
		my %fragments = map { $_ => 1 } $db->get_fragment_ids(-reference => $ref_id);
		ok $fragments{'FR0002267'};
		@ref_ids = $db->get_reference_ids(-fragment => 'FR0002267');
		is "@ref_ids", $ref_id;
    }
    
    # gene.dat
    {
        ok my ($gene_id) = $db->get_gene_ids(-name => 'P5');
        is $gene_id, 'G000001';
		
		#*** get_genemap with ensembl lookup being fantastically slow
        #ok defined Bio::Map::Gene->set_from_db; # will try and do ensembl lookups for gene info
        #ok my $gene_map = $db->get_genemap($gene_id, 1000);
        #Bio::Tools::Run::Ensembl->_stats;
        #ok $gene_map->isa('Bio::Map::GeneMap');
        #ok $gene_map->unique_id, 'G000001';
        #ok $gene_map->universal_name, 'P5';
        #ok $gene_map->species->scientific_name, 'Adeno-associated virus';
        #my @factors = grep { $_->isa("Bio::Map::TranscriptionFactor") } $gene_map->get_elements;
        #ok @factors, 3;
        
        ($gene_id) = $db->get_gene_ids(-id => 'AAV$P5');
        is $gene_id, 'G000001';
        my @gene_ids = $db->get_gene_ids(-species => '9606');
        is @gene_ids, 5;
        is [sort @gene_ids]->[0], 'G000060'; # in real data this would be G000174, but since our taxdump doesn't have chicken in it, G000060 was changed to human
        ($gene_id) = $db->get_gene_ids(-site => 'R03174');
        is $gene_id, 'G000001';
        ($gene_id) = $db->get_gene_ids(-factor => 'T00267');
        is $gene_id, 'G000060';
		my %gene_ids = map { $_ => 1 } $db->get_gene_ids(-fragment => 'FR0002267');
		ok $gene_ids{'G020751'};
        # get_gene_ids(-reference => ...) already tested
        
        my @site_ids = $db->get_site_ids(-gene => 'G000001');
        is join(' ', sort @site_ids), 'R03174 R03175 R03176';
        my @factor_ids = $db->get_factor_ids(-gene => 'G000060');
        is join(' ', sort @factor_ids), 'T00267 T08293'; # only found for genes that encode factors
		my %fragment_ids = map { $_ => 1 } $db->get_fragment_ids(-gene => 'G020751');
		ok $fragment_ids{'FR0002267'};
        # get_reference_ids(-gene => ...) already tested
    }
    
    # site.dat
    {
        ok my ($site_id) = $db->get_site_ids(-id => 'HS$IFI616_01');
        is $site_id, 'R00001';
        ok my $seq = $db->get_seq($site_id);
        isa_ok $seq, 'Bio::Seq';
        is $seq->id, 'HS$IFI616_01';
        is $seq->accession_number, 'R00001';
        is $seq->seq, 'aGAGACATAAGTgA';
        my $annot = $seq->annotation;
        is [$annot->get_Annotations('relative_start')]->[0]->value, -172;
        is [$annot->get_Annotations('relative_end')]->[0]->value, -98;
        is [$annot->get_Annotations('relative_type')]->[0]->value, 'TSS';
        is [$annot->get_Annotations('relative_to')]->[0]->value, 'G000176';
        is $seq->species, 9606;
        
        my @site_ids = $db->get_site_ids(-species => '9606');
        is @site_ids, 14;
        is [sort @site_ids]->[0], 'R00001';
        # get_site_ids(-gene => ...) already tested
        ($site_id) = $db->get_site_ids(-matrix => 'M00972');
        is $site_id, 'R00001';
        my %site_ids = map { $_ => 1 } $db->get_site_ids(-factor => 'T00428');
        ok $site_ids{R00001};
        # get_site_ids(-reference => ...) already tested
        
        # get_gene_ids(-site => ...) already tested
        my @matrix_ids = $db->get_matrix_ids(-site => 'R00001');
        is "@matrix_ids", 'M00972';
        my @factor_ids = $db->get_factor_ids(-site => 'R00001');
        is "@factor_ids", 'T00428';
        # get_reference_ids(-site => ...) already tested
    }
    
    # matrix.dat
    {
        ok my ($matrix_id) = $db->get_matrix_ids(-id => 'V$E47_01');
        is $matrix_id, 'M00002';
        ok my $matrix = $db->get_matrix($matrix_id);
        isa_ok $matrix, 'Bio::Matrix::PSM::SiteMatrix';
        
        # detailed psm tests
        {
            # Lets try to compress and uncompress the frequencies, see if
            # there is no considerable loss of data.
            my $fA = $matrix->get_compressed_freq('A');
            my @check = Bio::Matrix::PSM::SiteMatrix::_uncompress_string($fA,1,1);
            my @A = $matrix->get_array('A');
            my ($var, $max) = (0, 0);
            for (my $i = 0; $i < @check; $i++) {
                my $diff = abs(abs($check[$i]) - abs($A[$i]));
                $var += $diff;
                $max = $diff if ($diff > $max);
            }
            my $avg = $var / @check;
            cmp_ok $avg, '<', 0.01; # Loss of data under 1 percent
            
            # SiteMatrixI methods
            is $matrix->id, 'V$E47_01';
            is $matrix->accession_number, $matrix_id;
            is $matrix->consensus, 'ATGCATGCATGC';
            is $matrix->IUPAC, 'NNNNNNNNNNNN';
            is $matrix->regexp, '\S\S\S\S\S\S\S\S\S\S\S\S';
            is $matrix->width, 12;
            is $matrix->sites, 5;
            ok ! $matrix->IC;
            ok ! $matrix->e_val;
        }
        
        ok my $aln = $db->get_aln($matrix_id);
        isa_ok $aln, 'Bio::SimpleAlign';
        is $aln->length, 12;
        is $aln->num_residues, 132;
        ok $aln->is_flush;
        is $aln->num_sequences, 11;
        my @ids = qw(R05108 R05109 R05110 R05111 R05112 R05113 R05114 R05115 R05116 R05117 R05118);
        foreach my $seq ($aln->each_alphabetically) {
            is $seq->id, shift(@ids);
        }
        is @ids, 0;
        ok ! $db->get_aln('M00001'); # no seqs in db
        ok $aln = $db->get_aln('M00001', 1); # force to find seqs, store in db
        ok $aln = $db->get_aln('M00001'); # seqs now in db
        is $aln->num_sequences, 5;
		
        ($matrix_id) = $db->get_matrix_ids(-name => 'MyoD');
        is $matrix_id, 'M00001';
        # get_matrix_ids(-site =>  ...) already tested
        my %matrix_ids = map { $_ => 1 } $db->get_matrix_ids(-factor => 'T00526');
        ok $matrix_ids{M00001};
        # get_matrix_ids(-reference => ...) already tested
        
        # get_site_ids(-matrix => ...) already tested
        my @factor_ids = $db->get_factor_ids(-matrix => 'M00001');
        is join(' ', sort @factor_ids), 'T00526 T09177';
        # get_reference_ids(-matrix => ...) already tested
    }
    
	# fragment.dat
	{
		ok my ($fragment_id) = $db->get_fragment_ids(-id => 'FR0002267');
        is $fragment_id, 'FR0002267'; # id and accession are the same for fragments
		ok my $seq = $db->get_fragment($fragment_id);
		isa_ok $seq, 'Bio::SeqI';
        is $seq->id, 'FR0002267';
        is $seq->seq, 'GTCTACAACACTCTTGCGGACGGAGAGCCGAAGAGCAAAGCGTCGCCGGGTAAGACGAACGCTCAAGGGGGTACGAGCAGCGTAACGACGGAAACGGTGACGCCCCGGGATTTGGGGCTCAGCTAGGGTCGCCGAGTAGGGGGCCGCGGGGACAACGGGGGCGACACGCCGCTTTCCCTGCGTCTGTGGAGCCTATGGTACGGCGTAACCGGTTGTGTGATGAACTGTCCAGACCGCACGTAGTCCCAGCGCAAGGTCTATGCCGCCTAGAGGCAAGACGGGCCGTCTCCTACTTAGTAGCCAGCTACGGGGCGTTGGTCCCCTCGGTAGTGCAACTATCCAGCCACGGCGTCCGCCGGGCTGAGCCTCAGCAGAGCTGGGGGGGTATCATTCCGACGCTGTTTAATTCGTCAGCAGGACCCACTACACGCTCTGTCATTCGCCTGAGCAGTTGTAAATTAGCGCGGCGATCTTGCAAGAGACAAGGAGGCGAACCTGGGGTCGGGACGTAAGGACGAACGGCAGTACAGACGCTGGGGGACGCCACGTGCCAGAACCTCTCACGACCGGAGGTTCAACGCTGATTGGGGCGCAACAGAGGGCGGAGCAGCGAGGTGGCGCTGGTGGGATGGGGCGAGACAAACCCAAGCTGACGCCGAAGGGCCCGCGTGGCCGGGCTGGGGCCCGTAGAACGAGGGAATTGTATGCGGCGCCTGAATGGGCGCACCACA';
		is $seq->species, 9606;
		
        # -id -species -gene -factor -reference
        my @fragment_ids = $db->get_fragment_ids(-species => '9606');
        is @fragment_ids, 2;
        is [sort @fragment_ids]->[0], 'FR0000001';
        my %fragment_ids = map { $_ => 1 } $db->get_fragment_ids(-factor => 'T03828');
        ok $fragment_ids{'FR0002267'};
        # get_fragment_ids(-gene => ...) already tested
        # get_fragment_ids(-reference => ...) already tested
        
        my ($factor_id) = $db->get_factor_ids(-fragment => 'FR0002267');
        is $factor_id, 'T03828';
        # get_gene_ids(-fragment => ...) already tested
        # get_reference_ids(-fragment => ...) already tested
	}
	
    # factor.dat
    {
        ok my ($factor_id) = $db->get_factor_ids(-id => 'T00001');
        is $factor_id, 'T00001'; # id and accession are the same for factors
        ok my $factor = $db->get_factor($factor_id);
        isa_ok $factor, 'Bio::Map::TranscriptionFactor';
        is $factor->id, 'T00001';
        is $factor->universal_name, 'AAF';
        is $factor->known_maps, 1;
        my @positions = $factor->get_positions;
        is @positions, 1;
        
        ($factor_id) = $db->get_factor_ids(-name => 'AAF');
        is $factor_id, 'T00001';
        my @factor_ids = $db->get_factor_ids(-species => '9606');
        is @factor_ids, 7;
        is [sort @factor_ids]->[0], 'T00001';
        @factor_ids = $db->get_factor_ids(-interactors => 'T03200');
        is [sort @factor_ids]->[0], 'T00002';
        # get_factor_ids(-gene => ...) already tested
        # get_factor_ids(-site => ...) already tested
        # get_factor_ids(-matrix => ...) already tested
        # get_factor_ids(-fragment => ...) already tested
        # get_factor_ids(-reference => ...) already tested
        
        # get_*_ids(-factor => ...) already tested
    }
}

# how to get something like ok $psmIO->release, '10.2--2006-06-30'; ?
# or all factors, all sites, all matrices, all genes etc.?
