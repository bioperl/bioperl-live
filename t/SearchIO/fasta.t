# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_fasta.t 14995 2008-11-16 06:20:00Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 299);
	
	use_ok('Bio::SearchIO');
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
    float_is($hit->significance, shift @$d);
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
            float_is($hsp->evalue, 0.017);
            is($hsp->score, 134.5);
            is($hsp->bits,44.2);
            is(sprintf("%.2f",$hsp->percent_identity), '57.30');
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.5907); 
            is(sprintf("%.4f",$hsp->frac_identical('hit')), 0.5752);
			# these are really UNGAPPED values not CONSERVED
			# otherwise ident and conserved would be identical for
			# nucleotide alignments
			is(sprintf("%.4f",$hsp->frac_conserved('total')), '0.5918'); 
			is(sprintf("%.4f",$hsp->frac_conserved('query')), '0.6100'); 
			is(sprintf("%.4f",$hsp->frac_conserved('hit')), '0.5940'); 
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 0);
            is($hsp->gaps('query'), 8);
            is($hsp->gaps('hit'),1);
            is($hsp->gaps('total'), 9);
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
            is(join(' ', $hsp->seq_inds('hit', 'frameshift')), '');
            is(join(' ', $hsp->seq_inds('query', 'frameshift')), '');
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
    float_is($hit->significance, shift @$d);
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
            float_is($hsp->evalue, 1.2);
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
            is($hsp->gaps('hit'), 8);
            is($hsp->gaps, 15);
            is($hsp->query_string, 'NKEAIFTDDLPVADYLDDEFINSIPTAFDWRTRGAVTPVKNQGQCGSCWSFSTT-GNV----EGQHFISQNKLVSLSEQNLVDCDHECME-YEGEEACDEGCNGGLQPNAYNYIIKNGGIQTESSYPYTAETGTQCNFNSANIGAKISNFTMIPKNETVMAGYIVSTGP-LAIAADAVEWQFYIGGVFDIPCNPNSLDHGILIVGYSAKNTIFRKNMPYWIVKNSWGADWGEQGYIYLRRGKNTCGVSNFVSTSII');
            is($hsp->hit_string, (' 'x29).'MKIRSQVGMVLNLDKCIGCHTCSVTCKNVWTSREGVEYAWFNNVETKPGQGF-PTDWENQEKYKGGWI--RKINGKLQPRMGNRAMLLGKIFANPHLPGIDDYYEPFDFDYQNLHTAPEG----SKSQPIARPRSLITGERMAKIEKGPNWEDDLGGEFDKLAKDKNFDN-IQKAMYSQFENTFMMYLPRLCEHCLNPACVATCPSGAIYKREEDGIVLIDQDKCRGWRMCITGCPYKKIYFNWKSGKSEKCIFCYPRIEAGQPTVCSETC');
            is($hsp->homology_string, '                              . :. :  : :  .: .: . :.:  ::    :: ..   :.. .   :..   : : .: :.:     .  :: :::   :  .  : : ..   :   .     .:.  :. .   .     :.. .     . ::  .:    . .:.  .:: ::   . ...:. :  . ::  .. :   .:                      '.' 'x60);
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            is(sprintf("%.2f",$hsp->get_aln->percentage_identity()), 26.01);
            is(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '126 129 131 132 134 136 137 140 143 145 149 154 157-159 163 165-167 171-173 175 177 183-187 189 190 193 197-199 201 202 204 205 207 209 212-214 216-218 220-224 228 229 232 234-236 238-242 246 248-252 254 257 260-263 265 269 270 274 277-279 281 287 289 290 292 295 296 299 301-303');
            is(join(' ', $hsp->seq_inds('query', 'mismatch',1)), '126 129 131 132 134 136 137 140 143 145 149 154 157-159 163 165-167 172 173 175 177 183-185 189 190 193 197-199 201 202 204 205 207 209 212-214 216-218 220-224 228 229 232 234-236 242 246 248-252 254 257 260-263 265 269 270 274 277-279 281 289 290 292 295 296 299 301-303');
            is(join(' ', $hsp->seq_inds('query', 'conserved',1)), '125 127 128 130 133 135 138 139 141 142 144 146-148 150-153 155 156 160-162 164 168-170 174 176 178-182 188 191 192 194-196 200 203 206 208 210 211 215 219 225-227 230 231 233 237 243-245 247 253 255 256 258 259 264 266-268 271-273 275 276 280 282-286 288 291 293 294 297 298 300 304 305');
            is(join(' ', $hsp->seq_inds('hit', 'nomatch',1)), '3 6 8 9 11 13 14 17 20 22 26 27 30-33 36 39-41 45 47-49 53 54 56 58 61 65-67 69 70 73 77-79 81 82 84 85 87 89 92-94 96-98 100-104 108 109 112 114-116 118 122 124-128 130 133 134 137-140 142 146 147 151 154-156 158 165 166 168 171 172 175 177-179');
            is(join(' ', $hsp->seq_inds('hit', 'mismatch',1)), '3 6 8 9 11 13 14 17 20 22 27 36 39-41 45 47-49 53 54 56 58 65-67 69 70 73 77-79 81 82 84 85 87 89 92-94 96-98 100-104 108 109 112 114-116 118 122 124-128 130 133 137-140 142 146 147 151 154-156 158 165 166 168 171 172 175 177-179');
            is(join(' ', $hsp->seq_inds('hit', 'conserved',1)), '2 4 5 7 10 12 15 16 18 19 21 23-25 28 29 34 35 37 38 42-44 46 50-52 55 57 59 60 62-64 68 71 72 74-76 80 83 86 88 90 91 95 99 105-107 110 111 113 117 119-121 123 129 131 132 135 136 141 143-145 148-150 152 153 157 159-164 167 169 170 173 174 176 180 181');
            is(join(' ', $hsp->seq_inds('query', 'gap',1)), '148 151 179 257');
            is(join(' ', $hsp->seq_inds('hit', 'frameshift')), '');
            is(join(' ', $hsp->seq_inds('query', 'frameshift')), '');
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
    float_is($hit->significance, shift @$d);
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
            float_is($hsp->evalue, 3.1e-59);
            is($hsp->score, 1170.6);
            is($hsp->bits,227.8);
            is(sprintf("%.2f",$hsp->percent_identity), 53.04);
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.5479);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), '0.5430');
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 0);
            is($hsp->gaps('query'), 11);
            is($hsp->gaps('hit'), 8);
            is($hsp->gaps, 19);
            is($hsp->hit_string, (' 'x26).'MRTAVLLPLLAL----LAVAQA-VSFADVVMEEWHTFKLEHRKNYQDETEERFRLKIFNENKHKIAKHNQRFAEGKVSFKLAVNKYADLLHHEFRQLMNGFNYTLHKQLRAADESFKGVTFISPAHVTLPKSVDWRTKGAVTAVKDQGHCGSCWAFSSTGALEGQHFRKSGVLVSLSEQNLVDCSTKYGNNGCNGGLMDNAFRYIKDNGGIDTEKSYPYEAIDDSCHFNKGTVGATDRGFTDIPQGDEKKMAEAVATVGPVSVAIDASHESFQFYSEGVYNEPQCDAQNLDHGVLVVGFGTDESGED---YWLVKNSWGTTWGDKGFIKMLRNKENQCGIASASSYPLV');
            is($hsp->query_string, 'SNWGNNGYFLIERGKNMCGLAACASYPIPQVMNPTLILAAFCLGIASATLTFDHSLEAQWTKWKAMHNRLY-GMNEEGWRRAVWEKNMKMIELHNQEYREGKHSFTMAMNAFGDMTSEEFRQVMNGFQ---NRKPR------KGKVFQEPLFYEAPRSVDWREKGYVTPVKNQGQCGSCWAFSATGALEGQMFRKTGRLISLSEQNLVDCSGPQGNEGCNGGLMDYAFQYVQDNGGLDSEESYPYEATEESCKYNPKYSVANDTGFVDIPK-QEKALMKAVATVGPISVAIDAGHESFLFYKEGIYFEPDCSSEDMDHGVLVVGYGFESTESDNNKYWLVKNSWGEEWGMGGYVKMAKDRRNHCGIASAASYPTVMTPLLLLAVLCLGTALATPKFDQTFNAQWHQWKSTHRRLYGTNEE');
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            is(sprintf("%.2f",$hsp->get_aln->percentage_identity()), 56.13);
            is(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '1375 1377 1379 1381-1384 1391 1395-1397 1399 1400 1403 1404 1407 1408 1410 1412 1414 1415 1419 1422 1423 1429 1431 1433 1434 1440 1444 1447 1452 1457 1458 1473 1477 1480 1481 1483-1487 1495 1498 1501 1524 1530 1544-1546 1558 1580 1588-1592 1596 1607 1609 1630 1638 1658 1660 1662 1663 1665-1667 1677 1678 1681 1682 1688 1705');
            is(join(' ', $hsp->seq_inds('query', 'mismatch',1)), '1375 1377 1379 1395-1397 1399 1400 1403 1404 1407 1408 1410 1412 1414 1415 1419 1422 1423 1429 1431 1433 1434 1440 1444 1447 1452 1457 1458 1473 1477 1480 1481 1483-1487 1495 1498 1501 1524 1530 1544-1546 1558 1580 1588-1592 1596 1607 1609 1630 1638 1658 1660 1662 1663 1677 1678 1681 1682 1688 1705');
            is(join(' ', $hsp->seq_inds('query', 'conserved',1)), '1373 1374 1376 1378 1380 1385-1390 1392-1394 1398 1401 1402 1405 1406 1409 1411 1413 1416-1418 1420 1421 1424-1428 1430 1432 1435-1439 1441-1443 1445 1446 1448-1451 1453-1456 1459-1472 1474-1476 1478 1479 1482 1488-1494 1496 1497 1499 1500 1502-1523 1525-1529 1531-1543 1547-1557 1559-1579 1581-1587 1593-1595 1597-1606 1608 1610-1629 1631-1637 1639-1657 1659 1661 1664 1668-1676 1679 1680 1683-1687 1689-1704 1706');
            is(join(' ', $hsp->seq_inds('hit', 'nomatch',1)), '7 9 11 22-24 26 27 30 31 34 35 37 39 41-43 47 50 51 57 59 61 62 68 72 75 80 85 86 98-100 104 106-111 114 117 118 120-124 132 135 138 161 167 181-183 195 217 225-229 233 241 245 247 268 276 296 298 300 301 312 313 316 317 323 340');
            is(join(' ', $hsp->seq_inds('hit', 'mismatch',1)), '7 9 11 22-24 26 27 30 31 34 35 37 39 42 43 47 50 51 57 59 61 62 68 72 75 80 85 86 104 114 117 118 120-124 132 135 138 161 167 181-183 195 217 225-229 233 245 247 268 276 296 298 300 301 312 313 316 317 323 340');
            is(join(' ', $hsp->seq_inds('hit', 'conserved',1)), '5 6 8 10 12-21 25 28 29 32 33 36 38 40 44-46 48 49 52-56 58 60 63-67 69-71 73 74 76-79 81-84 87-97 101-103 105 112 113 115 116 119 125-131 133 134 136 137 139-160 162-166 168-180 184-194 196-216 218-224 230-232 234-240 242-244 246 248-267 269-275 277-295 297 299 302-311 314 315 318-322 324-339 341');
            is(join(' ', $hsp->seq_inds('query', 'gap',1)), '1413 1469 1474 1603');
            is(join(' ', $hsp->seq_inds('hit', 'frameshift')), '');
            is(join(' ', $hsp->seq_inds('query', 'frameshift')), '');
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
    float_is($hit->significance, shift @$d);
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
            
            float_is($hsp->evalue, 1.3e-25);
            is($hsp->score, 563.4);
            is($hsp->bits,'117.1');
            is(sprintf("%.2f",$hsp->percent_identity), 54.24);
            is(sprintf("%.4f",$hsp->frac_identical('query')), '0.5413');
            is(sprintf("%.4f",$hsp->frac_identical('hit')), '0.6654');
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 0);
            is($hsp->gaps('query'), 0);
            is($hsp->gaps('hit'), 61);
            is($hsp->gaps, 61);
            is($hsp->query_string, 'RFALAGALGCAVTHGALTPVDVVKTRIQLEPEVYNRVGRFFNSS*GF*EL*GVVLMSQT\KGMVASFRQIIAKEGAGALLTGFGPTAVGYAIQGAFKFGG*VMMSLQITA*SRANLPISYEFWKKKAIDLVGVDKARENRQAIYLGASAIAEFFADIALCPLEATRIRLVSQPSFANGLSGGFLRILREEGPAAFYAGFGPILFKQVPYTMAKFAV*VDRTA*\YQTFG*YYRSYEVAVEKILKATGKSKDSLTGGQLTGLNLTSGLIAGLAAAVISQPADTLLSKINKTKGAPGQSTTSRLVQMAGQLGVSGLFTGMTTRLVMIGTLTAGQ');
            is($hsp->hit_string, 'KFALAGAIGCGSTHSSMVPIDVVKTRIQLEPTVYN-------------------------KGMVGSFKQIIAGEGAGALLTGFGPTLLGYSIQGAFKFGG-------------------YEVFKKFFIDNLGYDTASRYKNSVYMGSAAMAEFLADIALCPLEATRIRLVSQPQFANGLVGGFSRILKEEGIGSFYSGFTPILFKQIPYNIAKFLVFERASEF-YYGFAG------------------PKEKLSSTSTTLLNLLSGLTAGLAAAIVSQPADTLLSKVNKTKKAPGQSTVGLLAQLAKQLGFFGSFAGLPTRLVMVGTLTSLQ');
            is($hsp->homology_string, '.::::::.::. ::....:.::::::::::: :::                         ::::.::.:::: ::::::::::::: .::.:::::::::                   :: .::  :: .: : : . ....:.:..:.:::.:::::::::::::::::::.::::: ::: :::.::: ..::.:: ::::::.::..::: :    .   :  :.                    :..:.. . : ::: ::: ::::::..::::::::::.:::: ::::::.. :.:.: :::  : :.:. :::::.::::. :');
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            is(sprintf("%.2f",$hsp->get_aln->overall_percentage_identity()),
               '53.92');
            is(sprintf("%.2f",$hsp->get_aln->average_percentage_identity()),
               '67.29');
            is(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '6950-6952 6986-6988 7001-7003 7007-7012 7022-7024 7040-7042 7067-7069 7139-7141 7151-7153 7163-7165 7169-7171 7175-7177 7196-7255 7262-7267 7272-7277 7281-7292 7296-7298 7341-7343 7365-7367 7389-7391 7401-7403 7524-7526 7530-7532 7536-7538 7542-7544 7551-7553 7560-7565 7575-7577 7584-7640 7680-7682 7722-7724 7762-7833 7843-7845 7903-7905');
            is(join(' ', $hsp->seq_inds('query', 'mismatch',1)), '6950-6952 6986-6988 7001-7003 7007-7012 7022-7024 7040-7042 7067-7069 7139-7141 7151-7153 7163-7165 7169-7171 7175-7177 7196-7198 7253-7255 7262-7267 7272-7277 7281-7292 7296-7298 7341-7343 7365-7367 7389-7391 7401-7403 7524-7526 7530-7532 7536-7538 7542-7544 7551-7553 7560-7565 7575-7577 7680-7682 7722-7724 7843-7845 7903-7905');
            is(join(' ', $hsp->seq_inds('query', 'conserved',1)), '6947-6949 6953-6985 6989-7000 7004-7006 7013-7021 7025-7039 7043-7066 7070-7138 7142-7150 7154-7162 7166-7168 7172-7174 7178-7195 7256-7261 7269-7272 7278-7280 7293-7295 7299-7340 7344-7364 7368-7388 7392-7400 7404-7523 7527-7529 7533-7535 7539-7541 7545-7550 7554-7559 7566-7574 7578-7583 7641-7679 7683-7721 7725-7757 7759-7762 7834-7842 7846-7902 7906-7938');
            is(join(' ', $hsp->seq_inds('query', 'gap',1)), '');
            is(join(' ', $hsp->seq_inds('hit', 'nomatch',1)), '30 50 66 80 96 100 101 104 107 109 111 113 154 158 166 174 189 191-194 196 197 199 200 203 204 211 213 215 219 223 247 256 262 266 267 269 274 286');
            is(join(' ', $hsp->seq_inds('hit', 'mismatch',1)), '30 50 66 80 96 100 101 104 107 109 111 113 154 158 166 174 189 191-194 196 197 199 200 203 204 211 213 215 219 223 247 256 262 266 267 269 274 286');
            is(join(' ', $hsp->seq_inds('hit', 'conserved',1)), '19-29 31-49 51-65 67-79 81-95 97-99 102 103 105 106 108 110 112 114-153 155-157 159-165 167-173 175-188 190 195 198 201 202 205-210 212 214 216-218 220-222 224-246 248-255 257-261 263-265 268 270-273 275-285 287');
            is(join(' ', $hsp->seq_inds('hit', 'gap',1)), '53 93 203');
            is(join(' ', $hsp->seq_inds('query', 'frameshift', 1)), '7271-7273 7761-7763');
            is(join(' ', $hsp->seq_inds('hit', 'frameshift')), '');
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
is($result->algorithm, 'FASTP');
is($result->algorithm_version, '35.04');
is($result->query_name, 'BOSS_DROME');
is($result->query_length, 896);
is($result->get_parameter('gapopen'), -10);
is($result->get_parameter('gapext'), -2);
is($result->get_parameter('ktup'), 2);
is($result->get_parameter('matrix'), 'BL50');

# wordsize is the same as ktup, not opt width, as we used to parse
# is($result->get_parameter('wordsize'), 16);

# this is the range of the scoring matrix, not a filter (which is meant
# to capture whether xS seg filtering used)
# is($result->get_parameter('filter'), '15:-5');

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

# bug 2537

$searchio = Bio::SearchIO->new(
   -format => 'fasta',
   -file   => test_input_file('in.fasta'),
);

while ( my $result = $searchio->next_result() ) {
   while( my $hit = $result->next_hit ) {
      # process the Bio::Search::Hit::HitI object
      while( my $hsp = $hit->next_hsp ) {
         # process the Bio::Search::HSP::HSPI object
         is($hsp->num_identical, 30, "num_identical()");
         is($hsp->num_conserved, 30, "num_conserved()");
      }
   }
}


# bug 2937 and version 35.04
$searchio = Bio::SearchIO->new(
   -format => 'fasta',
   -file   => test_input_file('bug2937.fasta'),
);
ok($searchio, 'bug 2937 and FASTA version 3.5');

while ( my $result = $searchio->next_result() ) {
   is ($result->algorithm_version, '35.04', 'algorithm version');
   is($result->query_name, 'ILTV-miR1', 'query name');
   is($result->query_description, '', 'query description');
   is($result->query_length, 70, 'query length');
   is($result->algorithm, 'FASTN', 'algorithm');
   
   while( my $hit = $result->next_hit ) {
      # process the Bio::Search::Hit::HitI object
      while( my $hsp = $hit->next_hsp ) {
         # process the Bio::Search::HSP::HSPI object
         is($hsp->num_identical, 26, "num_identical()");
         is($hsp->num_conserved, 26, "num_conserved()");
         is($hsp->strand('hit'), '1', 'hsp->strand(hit)');
		 is($hsp->hit->strand, '1', 'hsp->hit->strand');
		 is($hsp->strand('query'), '-1', 'hsp->strand(query)');
		 is($hsp->query->strand, '-1', 'hsp->query->strand');
      }
   }
}
