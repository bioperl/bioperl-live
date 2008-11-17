# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_fasta.t 14995 2008-11-16 06:20:00Z cjfields $

use strict;

BEGIN {
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 261);
	
	use_ok('Bio::SearchIO');
	use_ok('Bio::SearchIO::Writer::HitTableWriter');
	use_ok('Bio::SearchIO::Writer::HTMLResultWriter');
}

my ($searchio, $result, $hit, $hsp);

# simple FASTA
$searchio = Bio::SearchIO->new(-format => 'fasta',
				 -file   => test_input_file('HUMBETGLOA.FASTA'));
$result = $searchio->next_result;
like($result->database_name, qr/dros_clones.2.5/);
is($result->database_letters, 112936249);
is($result->database_entries, 657);
is($result->algorithm, 'FASTN');
is($result->algorithm_version, '3.3t08');
is($result->query_name, "HUMBETGLOA");
is($result->query_description, "Human haplotype C4 beta-globin gene, complete cds.");
is($result->query_length, 3002);
is($result->get_parameter('gapopen'), -16);
is($result->get_parameter('gapext'), -4);
is($result->get_parameter('ktup'), 6);

is($result->get_statistic('lambda'), 0.0823);
is($result->get_statistic('dbletters'), 112936249);
is($result->get_statistic('dbentries'), 657);

my @valid = ( [ 'BACR21I23', 73982, 'BACR21I23', '0.017', 44.2],
	   [ 'BACR40P19', 73982, 'BACR40P19', '0.017', 44.2],
	   [ 'BACR30L17', 32481, 'BACR30L17', '0.018', 44.1]);
my $count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is($hit->significance, shift @$d );
    is($hit->raw_score, shift @$d );
    is($hit->rank, $count + 1);
    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 31);
            is($hsp->query->end, 289);
            is($hsp->query->strand, -1);
            is($hsp->hit->end, 65167);
            is($hsp->hit->start, 64902);
            is($hsp->hit->strand, 1);
            is($hsp->length('total'), 267);	    
            is($hsp->evalue , 0.017);
            is($hsp->score, 134.5);
            is($hsp->bits,44.2);
            is(sprintf("%.2f",$hsp->percent_identity), '57.30');
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.5907); 
            is(sprintf("%.4f",$hsp->frac_identical('hit')), 0.5752);
			# these are really UNGAPPED values not CONSERVED
			# otherwise ident and conserved would be identical for
			# nucleotide alignments
			is(sprintf("%.4f",$hsp->frac_conserved('total')), 0.5955); 
			is(sprintf("%.4f",$hsp->frac_conserved('query')), 0.6139); 
			is(sprintf("%.4f",$hsp->frac_conserved('hit')), 0.5977); 
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 0);
            TODO: {
                local $TODO = 'Is this correct for gap total?';
                is($hsp->gaps('total'), 9);
            }
            is($hsp->gaps('query'), 8);
            is($hsp->gaps('hit'),1);
            is($hsp->query_string, 'GATTAAAACCTTCTGGTAAGAAAAGAAAAAATATATATATATATATATGTGTATATGTACACACATACATATACATATATATGCATTCATTTGTTGTTGTTTTTCTTAATTTGCTCATGCATGCTA----ATAAATTATGTCTAAAAATAGAAT---AAATACAAATCAATGTGCTCTGTGCATTA-GTTACTTATTAGGTTTTGGGAAACAAGAGGTAAAAAACTAGAGACCTCTTAATGCAGTCAAAAATACAAATAAATAAAAAGTCACTTACAACCCAAAGTGTGACTATCAATGGGGTAATCAGTGGTGTCAAATAGGAGGT');
            is($hsp->hit_string, 'GATGTCCTTGGTGGATTATGGTGTTAGGGTATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATATAAAATATAATATAAAATATAATATAAAATAAAATATAAAATAAAATATAAAATAAAATATAAAATAAAATATAAAATAAAATAT-AATATAAAATATAAAATAAAATATAATATAAAATATAATATAAAATATAATATAAAATATAATATAAAATA');
            is($hsp->homology_string, '                              :::::::::::::::::: : ::::: :: : : ::: ::::: ::::::::  ::  :: : :   : : : : :  ::    : :: ::   ::    : ::: :::     :::::: :::   ::::: ::  :::  :    :    : ::   :::  : ::   : :   : : :: :   :: : : :: : :       ::  : : ::: ::: ::  ::::: ::: : :  :: ::   ::: : : : ::: ::   '.' 'x60);
            is(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '33 37 39 41 43 47-49 52 55 56 58 60 64 70 71 74 78 82 84 86 87 90-96 98 100 103 105 107 110-112 114 117 119 121-123 125 127-129 132 134 135 139 140 143 145-148 150-153 155 156 160 161 164 173 180-184 188 192 196-198 201 204 206-209 212 213 215 217 219 221 223-225 227 229 232 233 236 237 246 252 256 258 260 263 269 271');
            is(join(' ', $hsp->seq_inds('query', 'mismatch',1)), '33 37 41 43 47-49 52 55 56 58 60 64 70 71 74 78 82 84 86 87 90-96 98 100 103 105 107 110-112 114 117 119 121-123 125 127-129 132 134 135 139 140 143 145-148 150-153 155 156 160 161 164 173 180-184 188 192 196-198 201 204 206-209 212 213 215 217 219 221 223-225 227 229 232 233 236 237 246 252 256 258 260 263 269 271');
            is(join(' ', $hsp->seq_inds('query', 'conserved',1)), '31 32 34-36 38 40 42 44-46 50 51 53 54 57 59 61-63 65-69 72 73 75-77 79-81 83 85 88 89 97 99 101 102 104 106 108 109 113 115 116 118 120 124 126 130 131 133 136-138 141 142 144 149 154 157-159 162 163 165-172 174-179 185-187 189-191 193-195 199 200 202 203 205 210 211 214 216 218 220 222 226 228 230 231 234 235 238-245 247-251 253-255 257 259 261 262 264-268 270 272-289');
            is(join(' ', $hsp->seq_inds('hit', 'nomatch',1)), '64920 64922 64928 64931 64933 64935 64939 64945 64954 64955 64958 64959 64962 64964 64966-64968 64970 64972 64974 64976 64978 64979 64982-64985 64987 64990 64993-64995 64998-65001 65003 65007 65011-65015 65022 65026-65028 65034 65037 65038 65042 65043 65045-65048 65050-65053 65055 65058-65060 65064 65065 65067 65070-65072 65074 65076-65078 65080 65082 65085 65087-65089 65092 65094 65096 65099 65101 65103-65109 65112 65113 65115 65117 65121 65125 65128 65129 65135 65139 65141 65143 65144 65147 65150-65152 65156 65158 65161 65165');
            is(join(' ', $hsp->seq_inds('hit', 'mismatch',1)), '64920 64922 64928 64931 64933 64935 64939 64945 64954 64955 64958 64959 64962 64964 64966-64968 64970 64972 64974 64976 64978 64979 64982-64985 64987 64990 64993-64995 65003 65007 65011-65015 65022 65034 65037 65038 65042 65043 65045-65048 65050-65053 65055 65059 65060 65064 65065 65067 65070-65072 65074 65076-65078 65080 65082 65085 65087-65089 65092 65094 65096 65099 65101 65103-65109 65112 65113 65115 65117 65121 65125 65128 65129 65135 65139 65141 65143 65144 65147 65150-65152 65156 65158 65161 65165');
            is(join(' ', $hsp->seq_inds('hit', 'conserved',1)), '64902-64919 64921 64923-64927 64929 64930 64932 64934 64936-64938 64940-64944 64946-64953 64956 64957 64960 64961 64963 64965 64969 64971 64973 64975 64977 64980 64981 64986 64988 64989 64991 64992 64996 64997 65002 65004-65006 65008-65010 65016-65021 65023-65025 65029-65033 65035 65036 65039-65041 65044 65049 65054 65056 65057 65061-65063 65066 65068 65069 65073 65075 65079 65081 65083 65084 65086 65090 65091 65093 65095 65097 65098 65100 65102 65110 65111 65114 65116 65118-65120 65122-65124 65126 65127 65130-65134 65136-65138 65140 65142 65145 65146 65148 65149 65153-65155 65157 65159 65160 65162-65164 65166 65167');            
            is(join(' ', $hsp->seq_inds('query', 'gap',1)), '141 170 194');
            is($hsp->ambiguous_seq_inds, '');
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            is(sprintf("%.2f",$hsp->get_aln->percentage_identity()), '59.30');
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);

$searchio = Bio::SearchIO->new(-format => 'fasta',
				 -file   => test_input_file('cysprot1.FASTA'));
$result = $searchio->next_result;
like($result->database_name, qr/ecoli.aa/);
is($result->database_letters, 1358987);
is($result->database_entries, 4289);
is($result->algorithm, 'FASTP');
is($result->algorithm_version, '3.3t08');
is($result->query_name, 'CYS1_DICDI');
is($result->query_length, 343);
is($result->get_parameter('gapopen'), -12);
is($result->get_parameter('gapext'), -2);
is($result->get_parameter('ktup'), 2);

is($result->get_statistic('lambda'), 0.1456);
is($result->get_statistic('dbletters'), 1358987);
is($result->get_statistic('dbentries'), 4289);


@valid = ( [ 'gi|1787478|gb|AAC74309.1|', 512, 'AAC74309', 1787478, 1.2, 29.2],
	   [ 'gi|1790635|gb|AAC77148.1|', 251, 'AAC77148', 1790635, 2.1, 27.4],
	   [ 'gi|1786590|gb|AAC73494.1|', 94, 'AAC73494', 1786590, 2.1, 25.9]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is($hit->ncbi_gi, shift @$d);
    is($hit->significance, shift @$d );
    is($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 125);
            is($hsp->query->end, 305);
            is($hsp->query->strand, 0);
            is($hsp->hit->start, 2);
            is($hsp->hit->end, 181);
            is($hsp->hit->strand, 0);
            is($hsp->length('total'), 188);	    
            is($hsp->evalue , 1.2);
            is($hsp->score, 109.2);
            is($hsp->bits,29.2);
            is(sprintf("%.2f",$hsp->percent_identity), 23.94);
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.2486);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), '0.2500');
            is(sprintf("%.4f",$hsp->frac_conserved('query')), '0.2707');
            is(sprintf("%.4f",$hsp->frac_conserved('hit')), '0.2722');
	    # there is slight rounding different here so file says 26.012%
	    # but with the rounding this ends up as 0.2606
            is(sprintf("%.4f",$hsp->frac_conserved('total')), '0.2606');
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 0);
            is($hsp->gaps('query'), 7);
            is($hsp->gaps, 49);	    
            is($hsp->query_string, 'NKEAIFTDDLPVADYLDDEFINSIPTAFDWRTRGAVTPVKNQGQCGSCWSFSTT-GNV----EGQHFISQNKLVSLSEQNLVDCDHECME-YEGEEACDEGCNGGLQPNAYNYIIKNGGIQTESSYPYTAETGTQCNFNSANIGAKISNFTMIPKNETVMAGYIVSTGP-LAIAADAVEWQFYIGGVFDIPCNPNSLDHGILIVGYSAKNTIFRKNMPYWIVKNSWGADWGEQGYIYLRRGKNTCGVSNFVSTSII');
            is($hsp->hit_string, (' 'x29).'MKIRSQVGMVLNLDKCIGCHTCSVTCKNVWTSREGVEYAWFNNVETKPGQGF-PTDWENQEKYKGGWI--RKINGKLQPRMGNRAMLLGKIFANPHLPGIDDYYEPFDFDYQNLHTAPEG----SKSQPIARPRSLITGERMAKIEKGPNWEDDLGGEFDKLAKDKNFDN-IQKAMYSQFENTFMMYLPRLCEHCLNPACVATCPSGAIYKREEDGIVLIDQDKCRGWRMCITGCPYKKIYFNWKSGKSEKCIFCYPRIEAGQPTVCSETC');
            is($hsp->homology_string, '                              . :. :  : :  .: .: . :.:  ::    :: ..   :.. .   :..   : : .: :.:     .  :: :::   :  .  : : ..   :   .     .:.  :. .   .     :.. .     . ::  .:    . .:.  .:: ::   . ...:. :  . ::  .. :   .:                      '.' 'x60);
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            is(sprintf("%.2f",$hsp->get_aln->percentage_identity()), 26.01);
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);

is($result->hits, 8);
$searchio = Bio::SearchIO->new(-format => 'fasta',
				 -file   => test_input_file('cysprot_vs_gadfly.FASTA'));
$result = $searchio->next_result;
like($result->database_name, qr/gadflypep2/);
is($result->database_letters, 7177762);
is($result->database_entries, 14334);
is($result->algorithm, 'FASTP');
is($result->algorithm_version, '3.3t08');
is($result->query_name, 'cysprot.fa');
is($result->query_length, 2385);
is($result->get_parameter('gapopen'), -12);
is($result->get_parameter('gapext'), -2);
is($result->get_parameter('ktup'), 2);
is($result->get_parameter('matrix'), 'BL50');

is($result->get_statistic('lambda'), 0.1397);
is($result->get_statistic('dbletters'), 7177762 );
is($result->get_statistic('dbentries'), 14334);


@valid = ( [ 'Cp1|FBgn0013770|pp-CT20780|FBan0006692', 341, 
	     'FBan0006692', '3.1e-59', 227.8],
	   [ 'CG11459|FBgn0037396|pp-CT28891|FBan0011459', 336, 
	     'FBan0011459', '6.4e-41',  166.9],
	   [ 'CG4847|FBgn0034229|pp-CT15577|FBan0004847', 390, 
	     'FBan0004847',  '2.5e-40', 165.2]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;

    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    is($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 1373);
            is($hsp->query->end, 1706);
            is($hsp->query->strand, 0);
            is($hsp->hit->start, 5);
            is($hsp->hit->end, 341);
            is($hsp->hit->strand, 0);
            is($hsp->length('total'), 345);	    
            is(sprintf("%g",$hsp->evalue), sprintf("%g",'3.1e-59') );
            is($hsp->score, 1170.6);
            is($hsp->bits,227.8);
            is(sprintf("%.2f",$hsp->percent_identity), 53.04);
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.5479);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), '0.5430');
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 0);
            is($hsp->gaps('query'), 11);
            is($hsp->gaps, 194);
            is($hsp->hit_string, (' 'x26).'MRTAVLLPLLAL----LAVAQA-VSFADVVMEEWHTFKLEHRKNYQDETEERFRLKIFNENKHKIAKHNQRFAEGKVSFKLAVNKYADLLHHEFRQLMNGFNYTLHKQLRAADESFKGVTFISPAHVTLPKSVDWRTKGAVTAVKDQGHCGSCWAFSSTGALEGQHFRKSGVLVSLSEQNLVDCSTKYGNNGCNGGLMDNAFRYIKDNGGIDTEKSYPYEAIDDSCHFNKGTVGATDRGFTDIPQGDEKKMAEAVATVGPVSVAIDASHESFQFYSEGVYNEPQCDAQNLDHGVLVVGFGTDESGED---YWLVKNSWGTTWGDKGFIKMLRNKENQCGIASASSYPLV');
            is($hsp->query_string, 'SNWGNNGYFLIERGKNMCGLAACASYPIPQVMNPTLILAAFCLGIASATLTFDHSLEAQWTKWKAMHNRLY-GMNEEGWRRAVWEKNMKMIELHNQEYREGKHSFTMAMNAFGDMTSEEFRQVMNGFQ---NRKPR------KGKVFQEPLFYEAPRSVDWREKGYVTPVKNQGQCGSCWAFSATGALEGQMFRKTGRLISLSEQNLVDCSGPQGNEGCNGGLMDYAFQYVQDNGGLDSEESYPYEATEESCKYNPKYSVANDTGFVDIPK-QEKALMKAVATVGPISVAIDAGHESFLFYKEGIYFEPDCSSEDMDHGVLVVGYGFESTESDNNKYWLVKNSWGEEWGMGGYVKMAKDRRNHCGIASAASYPTVMTPLLLLAVLCLGTALATPKFDQTFNAQWHQWKSTHRRLYGTNEE');
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            is(sprintf("%.2f",$hsp->get_aln->percentage_identity()), 56.13);
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);
is($result->hits, 21);

# test on TFASTXY
$searchio = Bio::SearchIO->new(-format => 'fasta',
			      -file   => test_input_file('5X_1895.FASTXY'));
$result = $searchio->next_result;
like($result->database_name, qr/yeast_nrpep.fasta/);
is($result->database_letters, 4215311);
is($result->database_entries, 9190);
is($result->algorithm, 'FASTY');
is($result->algorithm_version, '3.4t07');
is($result->query_name, '5X_1895.fa');
is($result->query_length, 7972);
is($result->get_parameter('gapopen'), -14);
is($result->get_parameter('gapext'), -2);
is($result->get_parameter('ktup'), 2);
is($result->get_parameter('matrix'), 'BL50');

is($result->get_statistic('lambda'), 0.1711);
is($result->get_statistic('dbletters'), 4215311);
is($result->get_statistic('dbentries'), 9190);


@valid = ( [ 'NR_SC:SW-YNN2_YEAST', 1056, 'NR_SC:SW-YNN2_YEAST','1.6e-154', '547.0'],
	   [ 'NR_SC:SW-MPCP_YEAST', 311, 'NR_SC:SW-MPCP_YEAST', '1.3e-25', 117.1],
	   [ 'NR_SC:SW-YEO3_YEAST', 300, 'NR_SC:SW-YEO3_YEAST', '5.7e-05', 48.5]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;

    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    is($hit->raw_score, shift @$d );

    if ( $count == 1 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 6947);
            is($hsp->query->end, 7938);
            is($hsp->query->strand, -1);
            is($hsp->hit->start, 19);
            is($hsp->hit->end, 287);
            is($hsp->hit->strand, 0);
            is($hsp->length('total'), 330);
            
            is(sprintf("%g",$hsp->evalue), sprintf("%g",'1.3e-25'));
            is($hsp->score, 563.4);
            is($hsp->bits,'117.1');
            is(sprintf("%.2f",$hsp->percent_identity), 51.67);
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.5244);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), 0.5728);
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 0);
            is($hsp->gaps, 220);
            is($hsp->query_string, 'RFALAGALGCAVTHGALTPVDVVKTRIQLEPEVYNRVGRFFNSS*GF*EL*GVVLMSQT\KGMVASFRQIIAKEGAGALLTGFGPTAVGYAIQGAFKFGG*VMMSLQITA*SRANLPISYEFWKKKAIDLVGVDKARENRQAIYLGASAIAEFFADIALCPLEATRIRLVSQPSFANGLSGGFLRILREEGPAAFYAGFGPILFKQVPYTMAKFAV*VDRTA*\YQTFG*YYRSYEVAVEKILKATGKSKDSLTGGQLTGLNLTSGLIAGLAAAVISQPADTLLSKINKTKGAPGQSTTSRLVQMAGQLGVSGLFTGMTTRLVMIGTLTAGQ');
            is($hsp->hit_string, 'KFALAGAIGCGSTHSSMVPIDVVKTRIQLEPTVYN-------------------------KGMVGSFKQIIAGEGAGALLTGFGPTLLGYSIQGAFKFGG-------------------YEVFKKFFIDNLGYDTASRYKNSVYMGSAAMAEFLADIALCPLEATRIRLVSQPQFANGLVGGFSRILKEEGIGSFYSGFTPILFKQIPYNIAKFLVFERASEF-YYGFAG------------------PKEKLSSTSTTLLNLLSGLTAGLAAAIVSQPADTLLSKVNKTKKAPGQSTVGLLAQLAKQLGFFGSFAGLPTRLVMVGTLTSLQ');
            is($hsp->homology_string, '.::::::.::. ::....:.::::::::::: :::                         ::::.::.:::: ::::::::::::: .::.:::::::::                   :: .::  :: .: : : . ....:.:..:.:::.:::::::::::::::::::.::::: ::: :::.::: ..::.:: ::::::.::..::: :    .   :  :.                    :..:.. . : ::: ::: ::::::..::::::::::.:::: ::::::.. :.:.: :::  : :.:. :::::.::::. :');
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            is(sprintf("%.2f",$hsp->get_aln->overall_percentage_identity()),
               '54.27');
            is(sprintf("%.2f",$hsp->get_aln->average_percentage_identity()),
               '67.42');
            TODO: {
                local $TODO = 'calculating seq_inds for FASTA including frameshifts not supported yet';
                is(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '');
                is(join(' ', $hsp->seq_inds('query', 'mismatch',1)), '');
            }
            is(join(' ', $hsp->seq_inds('hit', 'nomatch',1)), '30 50 66 80 96 100 101 104 107 109 111 113 154 158 166 174 189 191-194 196 197 199 200 203 204 211 213 215 219 223 247 256 262 266 267 269 274 286');
            is(join(' ', $hsp->seq_inds('hit', 'mismatch',1)), '30 50 66 80 96 100 101 104 107 109 111 113 154 158 166 174 189 191-194 196 197 199 200 203 204 211 213 215 219 223 247 256 262 266 267 269 274 286');
            is(join(' ', $hsp->seq_inds('query', 'gap',1)), '');
            is(join(' ', $hsp->seq_inds('hit', 'gap',1)), '53 93 203');
            is($hsp->ambiguous_seq_inds, 'query');
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);
is($result->hits, 58);

# test FASTA v35.04, params encoding changed 
# test on TFASTXY
$searchio = Bio::SearchIO->new(-format => 'fasta',
			      -file   => test_input_file('BOSS_DROME.FASTP_v35_04'));
$result = $searchio->next_result;
like($result->database_name, qr/wormpep190/, 'TFASTXY');
is($result->database_letters, 10449259);
is($result->database_entries, 23771);
is($result->algorithm, 'FASTA');
is($result->algorithm_version, '35.04');
is($result->query_name, 'BOSS_DROME');
is($result->query_length, 896);
is($result->get_parameter('gapopen'), -10);
is($result->get_parameter('gapext'), -2);
is($result->get_parameter('ktup'), 2);
is($result->get_parameter('matrix'), 'BL50');
is($result->get_parameter('wordsize'), 16);
is($result->get_parameter('filter'), '15:-5');

is($result->get_statistic('lambda'), 0.122629);
is($result->get_statistic('dbletters'), 10449259);
is($result->get_statistic('dbentries'), 23771);
is($result->get_statistic('effectivespace'),23771);

# bug 2473 - fasta3.4 parsing with -U option

$searchio = Bio::SearchIO->new(-format => 'fasta',
							  -file   => test_input_file('bug2473.fasta'));

while(my $res = $searchio->next_result) {
    is($res->query_name, 'total:39860_L:12096_-3:12346_0:617_+3:14801');
    is($res->query_description, '');
    is($res->query_length, 22);
    is($res->algorithm, 'FASTN');
}

