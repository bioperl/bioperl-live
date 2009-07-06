use strict;
use warnings;

use Bio::Root::Test;
use Bio::AlignIO;

test_begin( -tests => 6 ,
            -requires_modules => [qw(IO::String)],
          );

# CLUSTAL
{
    my $io = Bio::AlignIO->new( -file => test_input_file("testaln.aln") );
    my $aln = $io->next_aln();
    isa_ok( $aln, 'Bio::SimpleAlign' );
    my $consensus = <<EOU;
MNEGEHQIKLDELFEKLLRARKIFKNKDVLRHSWEPKDLPHRHEQIEALAQILV
PVLRGETMKIIFCGHHACELGEDRGTKGFVIDELKDVDEDRNGKVDVIEINCEH
MDTHYRVLPNIAKLFDDCTGIGVPMHGGPTDEVTAKLKQVIDMKERFVIIVLDE
IDKLVKKSGDEVLYSLTRINTELKRAKVSVIGISNDLKFKEYLDPRVLSSLSEE
EVVFPPYDANQLRDILTQRAEEAFYPGVLDEGVIPLCAALAAREHGDARKALDL
LRVAGEIAEREGASKVTEKHVWKAQEKIEQDMMEEVIKTLPLQSKVLLYAIVLL
DENGDLPANTGDVYAVYRELCEYIDLEPLTQRRISDLINELDMLGIINAKVVSK
GRYGRTKEIRLMVTSYKIRNVLRYDYSIQPLLTISLKSEQRRLI
EOU
    $consensus =~ s/\n//g;

    is( $aln->consensus_string, $consensus, 'consensus string looks ok' );

    is( aln2str( $aln => 'pfam' ), <<EOA, 'looks like correct unmasked alignment (from clustalw)' );
P84139/1-420              MNEGEHQIKLDELFEKLLRARKIFKNKDVLRHSYTPKDLPLRHEQIETLAQILVPVLRGETPSNIFVYG-KTGTGKTVTVK-FVTEELKRISEKYNIPVDVIYINCEIVDTHYRVLANIVNYFKDETGIGVPMVGWPTDEVYAKLKQVIDMKERFVIIVLDEIDKLVKKSGDEVLYSLTRINTELKRAKVSVIGISNDLKFKEYLDPRVLSSLSEEEVVFPPYDANQLRDILTQRAEEAFYPGVLDEGVIPLCAALAAREHGDARKALDLLRVAGEIAEREGASKVTEKHVWKAQEKIEQDMMEEVIKTRPLQSKVLLYAIVLLDENGDLPANTGDVYAVYRELCEYIDLEPLTQRRISDLINELDMLGIINAKVVSKGRYGRTKEIRLNVTSYKIRNVLRYDYSIQPLLTISLKSEQRRLI
P814153/1-420             MNEGMHQIKLDVLFEKLLRARKIFKNKDVLRHSYTPKDLPHRHEQIETLAQILVPVLRGETPSNIFVYG-KTGTGKTVTVK-FVTEELKRISEKYNIPVDVIYINCEIVDTHYRVLANIVNYFKDETGIEVPMVGWPTDEVYAKLKQVIDMKERFVIIVLDEIDKLVKKSGDEVLYSLTRINTELKRAKVSVIGISNDLKFKEYLDPRVLSSLSEEEVVFPPYDANQLRDILTQRAEEAFYPGVLDEGVIPLCAALAAREHGDARKALDLLRVAGEIAEREGASKVTEKHVWKAQEKIEQDMMEEVIKTLPLQSKVLLYAIVLLDENGDLPANTGDVYAVYRELCEYIDLEPLTQRRISDLINELDMLGIINAKVVSKGRYGRTKEIRLMVTSYKIRNVLRYDYSIQPLLTISLKSEQRRLI
P851414/1-60              -------------------------------------------------------------MKIVWCGH-ACFLVEDRGTK-ILIDPYPDVDEDRIGKVDYILQTHEHMD-HYGKTPLIAKLSD----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
P841414/1-60              -------------------------------------------------------------MKIVWCGH-ACFLVEDRGTK-ILIDPYPDVDEDRIGKVDYILVTHEHMD-HYGKTPLIAKLSD----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
BAB68554/1-141            --------------------MLTEDDKQLIQHVWEKVLEHQEDFGAEALERMFIVYPSTKTYFPHFDLHHDSEQIRHHGKK-VVGALGDAVKHIDNLSATLSELSNLHCY-NLRVDPVNFKLLSHCFQVVLGAHLG--REYTPQVQVAYDKFLAAVSAVLAEKYR-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
gb|443893|124775/1-331    -MRFRFGVVVPPAVAGARPELLVVGSRPELG-RWEPRGAVRLRPAGTAAGDGALALQEPGLWLGEVELA-AEEAAQDGAEPGRVDTFWYKFLKREPGGELSWEGNGPHHDRCCTYNENNLVDGVYCLPIG---HWGEATGHTNEMKHTTDFYFNIAGHQAMHYSRILPNIWLGSCPRQVEHVTIKLKHELGITAVMN-FQTEWDIVQNSSGCNRYPEPMTPDTMIKLYREEGLAYIWMP-TPDMSTEGRVQMLPQAVCLLHALLEKGHIVY-----VHCNAGVGRSTAAVCGWLQYVMGWNLRKVQYFLMAKRPAVYIDEEALARAQEDFFQKFGKVRSSVCSL------------------------------------------------------------------------------
EOA

    my $newaln = $aln->mask_columns(12,20,'?');
    is( aln2str( $newaln, 'pfam' ), <<EOA, 'looks like correct masked alignment (from clustalw)' );
P84139/1-420              MNEGEHQIKLD?????????RKIFKNKDVLRHSYTPKDLPLRHEQIETLAQILVPVLRGETPSNIFVYG-KTGTGKTVTVK-FVTEELKRISEKYNIPVDVIYINCEIVDTHYRVLANIVNYFKDETGIGVPMVGWPTDEVYAKLKQVIDMKERFVIIVLDEIDKLVKKSGDEVLYSLTRINTELKRAKVSVIGISNDLKFKEYLDPRVLSSLSEEEVVFPPYDANQLRDILTQRAEEAFYPGVLDEGVIPLCAALAAREHGDARKALDLLRVAGEIAEREGASKVTEKHVWKAQEKIEQDMMEEVIKTRPLQSKVLLYAIVLLDENGDLPANTGDVYAVYRELCEYIDLEPLTQRRISDLINELDMLGIINAKVVSKGRYGRTKEIRLNVTSYKIRNVLRYDYSIQPLLTISLKSEQRRLI
P814153/1-420             MNEGMHQIKLD?????????RKIFKNKDVLRHSYTPKDLPHRHEQIETLAQILVPVLRGETPSNIFVYG-KTGTGKTVTVK-FVTEELKRISEKYNIPVDVIYINCEIVDTHYRVLANIVNYFKDETGIEVPMVGWPTDEVYAKLKQVIDMKERFVIIVLDEIDKLVKKSGDEVLYSLTRINTELKRAKVSVIGISNDLKFKEYLDPRVLSSLSEEEVVFPPYDANQLRDILTQRAEEAFYPGVLDEGVIPLCAALAAREHGDARKALDLLRVAGEIAEREGASKVTEKHVWKAQEKIEQDMMEEVIKTLPLQSKVLLYAIVLLDENGDLPANTGDVYAVYRELCEYIDLEPLTQRRISDLINELDMLGIINAKVVSKGRYGRTKEIRLMVTSYKIRNVLRYDYSIQPLLTISLKSEQRRLI
P851414/1-69              -----------?????????-----------------------------------------MKIVWCGH-ACFLVEDRGTK-ILIDPYPDVDEDRIGKVDYILQTHEHMD-HYGKTPLIAKLSD----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
P841414/1-69              -----------?????????-----------------------------------------MKIVWCGH-ACFLVEDRGTK-ILIDPYPDVDEDRIGKVDYILVTHEHMD-HYGKTPLIAKLSD----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
BAB68554/1-150            -----------?????????MLTEDDKQLIQHVWEKVLEHQEDFGAEALERMFIVYPSTKTYFPHFDLHHDSEQIRHHGKK-VVGALGDAVKHIDNLSATLSELSNLHCY-NLRVDPVNFKLLSHCFQVVLGAHLG--REYTPQVQVAYDKFLAAVSAVLAEKYR-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
gb|443893|124775/1-331    -MRFRFGVVVP?????????LLVVGSRPELG-RWEPRGAVRLRPAGTAAGDGALALQEPGLWLGEVELA-AEEAAQDGAEPGRVDTFWYKFLKREPGGELSWEGNGPHHDRCCTYNENNLVDGVYCLPIG---HWGEATGHTNEMKHTTDFYFNIAGHQAMHYSRILPNIWLGSCPRQVEHVTIKLKHELGITAVMN-FQTEWDIVQNSSGCNRYPEPMTPDTMIKLYREEGLAYIWMP-TPDMSTEGRVQMLPQAVCLLHALLEKGHIVY-----VHCNAGVGRSTAAVCGWLQYVMGWNLRKVQYFLMAKRPAVYIDEEALARAQEDFFQKFGKVRSSVCSL------------------------------------------------------------------------------
EOA
}


# my $outfile = test_output_file();
# my $strout = Bio::AlignIO->new
#     ( -file   => ">$outfile",
#       -format => 'clustalw'
#     );
# $status = $strout->write_aln($aln);
# is $status, 1, "clustalw (.aln) output test";
# undef $strout;
# $str = Bio::AlignIO->new(
#    '-file'=> $outfile, 
# 			   '-format' => 'clustalw');
# $aln = $str->next_aln($aln);
# isa_ok($aln,'Bio::Align::AlignI');
# is $aln->get_seq_by_pos(1)->get_nse, 'P84139/1-420', "clustalw (.aln) input test";


###### test with phylip

{
    my $phy_fh = IO::String->new( <<EOF );
 3   37
seq1        AAAATGGGGG TGGT------ GGTACCT--- -------
seq2        -----GGCGG TGGTGNNNNG GGTTCCCTNN NNNNNNN
new         AAAATGGNGG TGGTN----N GGTNCCNTNN NNNNNNN
EOF

    my $in = Bio::AlignIO->new( -fh => $phy_fh, -format => 'phylip' );

    my $aln = $in->next_aln();
    is( aln2str( $aln, 'phylip' ), <<EOU );
 3 37
seq1         AAAATGGGGG TGGT------ GGTACCT--- ------- 
seq2         -----GGCGG TGGTGNNNNG GGTTCCCTNN NNNNNNN 
new          AAAATGGNGG TGGTN----N GGTNCCNTNN NNNNNNN 

EOU

    my $newaln = $aln->mask_columns(15,20,'?');
    #$out->write_aln($newaln);
    is( aln2str( $newaln,'phylip' ), <<EOU, 'align after looks ok' );
 3 37
seq1         AAAATGGGGG TGGT?????? GGTACCT--- ------- 
seq2         -----GGCGG TGGT?????? GGTTCCCTNN NNNNNNN 
new          AAAATGGNGG TGGT?????? GGTNCCNTNN NNNNNNN 

EOU
}

######## SUBROUTINES

sub aln2str {
    my ( $aln, $fmt ) = @_;
    my $out;
    my $out_fh = IO::String->new( $out );
    my $alignio_out = Bio::AlignIO->new(-fh => $out_fh, -format => $fmt);
    $alignio_out->write_aln( $aln );
    return $out;
}
