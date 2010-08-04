#!/usr/bin/perl

# PROGRAM  : psw.pl
# PURPOSE  : Simple driver for Bio::Tools::pSW
# AUTHOR   : Ewan Birney birney@sanger.ac.uk 
# CREATED  : Tue Oct 27 1998
#
# INSTALLATION
#  
# you almost certainly have to have installed bioperl 
# from the makefile system for this to work. This is
# because this module use XS extensions (C source code
# 'compiled into' perl)
#
# The lib system below is just so that I (ewan) can test it
# on site... 
#

use lib "/nfs/disk100/pubseq/wise/PerlMod/";

# 
# This is a simple example script. We are going
# to make 3 sequences directly from memory and
# then align them once using blosum matrix and once
# using a gonnet matrix. These matrices should
# in the examples directory.
#

use Bio::Tools::pSW;

# redundant, as Bio::Tools::pSW uses them, but useful to say
# precisely what we are using ;)

use Bio::Seq; 
use Bio::SimpleAlign;
use Bio::AlignIO;

# for legibility - write with newlines and then strip them!

$tseq = 'SKSESPKEPEQLRKLFIGGLSFETTDESLRSHFEQWGTLTDCVVMRDPNTKRSRGFGFVT
YATVEEVDAAMNARPHKVDGRVVEPKRAVSREDSQRPGAHLTVKKIFVGGIKEDTEEHHL
RDYFEQYGKIEVIEIMTDRGSGKKRGFAFVTFDDHDSVDKIVIQKYHTVNGHNCEVRKAL
SKQEMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSG
DGYNGFGNDGGYGGGGPGYSGGSRGYGSGGQGYGNQGSGYGGSGSYDSYNNGGGRGFGGG
SGSNFGGGGSYNDFGNYNNQSSNFGPMKGGNFGGRSSGPYGGGGQYFAKPRNQGGYGGSS
SSSSYGSGRRF';

$tseq =~ s/[^A-Z]//g;

$seq1 = Bio::Seq->new(-id=>'roa1_human',-seq=>$tseq);

$tseq = 'MVNSNQNQNGNSNGHDDDFPQDSITEPEHMRKLFIGGLDYRTTDENLKAHFEKWGNIVDV
VVMKDPRTKRSRGFGFITYSHSSMIDEAQKSRPHKIDGRVVEPKRAVPRQDIDSPNAGAT
VKKLFVGALKDDHDEQSIRDYFQHFGNIVDINIVIDKETGKKRGFAFVEFDDYDPVDKVV
LQKQHQLNGKMVDVKKALPKQNDQQGGGGGRGGPGGRAGGNRGNMGGGNYGNQNGGGNWN
NGGNNWGNNRGGNDNWGNNSFGGGGGGGGGYGGGNNSWGNNNPWDNGNGGGNFGGGGNNW
NNGGNDFGGYQQNYGGGPQRGGGNFNNNRMQPYQGGGGFKAGGGNQGNYGGNNQGFNNGG
NNRRY';

$tseq =~ s/[^A-Z]//g;

$seq2 = Bio::Seq->new(-id=>'roa1_drome',-seq=>$tseq);

$tseq = 'MHKSEAPNEPEQLRKLFIGGLSFETTDESLREHFEQWGTLTDCVVMRDPNSKRSRGFGFV
TYLSTDEVDAAMTARPHKVDGRVVEPKRAVSREDSSRPGAHLTVKKIFVGGIKEDTEEDH
LREYFEQYGKIEVIEIMTDRGSGKKRGFAFVTFEDHDSVDKIVIQKYHTVNNHNSQVRKA
LSKQEMASVSGSQRERGGSGNYGSRGGFGNDNFGGRGGNFGGNRGGGGGFGNRGYGGDGY
NGDGQLWWQPSLLGWNRGYGAGQGGGYGAGQGGGYGGGGQGGGYGGNGGYDGYNGGGSGF
SGSGGNFGSSGGYNDFGNYNSQSSSNFGPMKGGNYGGGRNSGPYGGGYGGGSASSSSGYG
GGRRF';

$tseq =~ s/[^A-Z]//g;
$seq3 = Bio::Seq->new(-id=>'roa1_xenla',-seq=>$tseq);


#
# Now make an Alignment Factory with blosum62 as a matrix
# gap -12 and ext -2
#

$fac = Bio::Tools::pSW->new(-matrix => 'blosum62.bla',-gap => 12, -ext => 2);


#
# run seq1 vs seq2 and seq1 vs seq3 and write the output direct 
# to stdout using the 'pretty' method
#

$fac->align_and_show($seq1,$seq2,STDOUT);
print "Next alignment\n";
$fac->align_and_show($seq1,$seq3,STDOUT);


#
# a different factory, using gonnet, and now make a simple align and 
# provide MSF format
#

$fac = Bio::Tools::pSW->new(-matrix => 'gon250.bla',-gap => 12, -ext => 2);


# switch on reporting this time and change the amount of memory it is allowed

print STDOUT "Doing the next calculation in limited memory, with a progress report\n";

$fac->report(1);
$fac->kbyte(100);

$al = $fac->pairwise_alignment($seq1,$seq2);


# write out a MSF file
my $out = Bio::AlignIO->newFh('-fh'=> \*STDOUT,  '-format' => 'msf');
my $status = print $out $al;
#$al->write_MSF(\*STDOUT);



