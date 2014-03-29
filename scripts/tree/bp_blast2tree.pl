#!/usr/bin/perl
# Author:  Jason Stajich <jason@bioperl.org>
# Purpose: Blast Report -> MSA -> Tree 


# This needs lots more error checking, cmdline input of parameters
# and support for other treebuilding -- only Phylip Neighbor for 
# protein data is supported

# Also proper pulling in of the correct sequence data from the
# alignment - multiple hits on different parts of a protein aren't
# going to work properly right now.  So this is mostly and example
# starting point which needs a lot more work to be made robust.

use strict;
use warnings;
use Bio::AlignIO;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Tools::Run::Phylo::Phylip::ProtDist;
use Bio::Tools::Run::Phylo::Phylip::Neighbor;
use Bio::Tools::Run::Phylo::Molphy::ProtML;
use Bio::Tools::Run::Phylo::Phylip::ProtPars;
use Bio::SearchIO;
use Bio::SimpleAlign;
use Bio::PrimarySeq;
use Bio::TreeIO;
use Getopt::Long;

my $IDLENGTH = 12;


# we could in fact pull the tree out of the guide tree calculated
# by Clustalw in the alignment, but I believe that is UPGMA 
# which would *NOT* be a good thing to be giving people.

my $aln_factory = Bio::Tools::Run::Alignment::Clustalw->new
    ('ktuple' => 2, "quiet" => 1,
     'matrix' => 'BLOSUM');
my ($report,$format,$tree_method,$cutoff,$keepall);

$format = 'blast';
$tree_method = 'neighbor';
$cutoff = '0.01';
GetOptions(
	   'h|help'       => sub { exec('perldoc', $0);
				   exit(0); },
	   'i|input:s'    => \$report,
	   'f|format:s'   => \$format,
	   'm|method:s'   => \$tree_method,
	   'e|evalue:s'   => \$cutoff,
	   'k|keepall:s'  => \$keepall, # keep all HSPs not just best
	   );


unless( $format =~ /blast|fasta|hmmer/i ) {
    die("Must request a valid search report format (fasta,blast,hmmer)");
}

unless ( $tree_method =~ /nj|neighbor/i || $tree_method =~ /protml|ml/i ) {
    die("Must request a valid tree building method (neighbor,protml)");
} 

my (@alns,@seqs);

my $in = new Bio::SearchIO(-file   => $report,
			   -format => $format);
while( my $r = $in->next_result ) {
    # Let's build the simplest system first
    die("Cannot process report that does not contain protein sequence") 
	unless ($r->algorithm =~ /HMMER|BLASTP|FASTP/i );
    my @seqs;
    while( my $hit = $r->next_hit ) {
	next if $hit->significance > $cutoff;
	while( my $hsp = $hit->next_hsp ) {
	    next if $hsp->evalue > $cutoff;
	    my $seq = $hsp->get_aln->get_seq_by_pos(2)->seq();
	    push @seqs, new Bio::PrimarySeq(-seq => $seq,
					    -id  => $hsp->hit->seq_id,
					    -desc => $r->algorithm . " best align to ". $hsp->query->seq_id );
	    last unless $keepall;
	}
    }
    push @alns, $aln_factory->align(\@seqs);    
}

my $out = new Bio::AlignIO(-format => 'phylip',
			   -interleaved => 1,
			   -idlength => $IDLENGTH,
			   -file => ">alignfile.phy");
$out->write_aln(@alns);
$out = undef;

# these need to be parameterized for cmdline arguments
my @params = ('idlength'=>$IDLENGTH,
	      'model'=>'cat',
	      'gencode'=>'U',
	      'category'=>'H',
	      'probchange'=>'0.4',
	      'trans'=>'5',
	      'freq'=>'0.25,0.5,0.125,0.125');
my $dist_factory = Bio::Tools::Run::Phylo::Phylip::ProtDist->new(@params);
$dist_factory->quiet(1);
@params = ('type'=>'NJ',
	   'outgroup'=>1,
	   'upptri'=>1,
	   'jumble'=>17);

my $tree_factory = Bio::Tools::Run::Phylo::Phylip::Neighbor->new(@params);
$tree_factory->quiet(1);
my $count = 1;
my $outtrees = new Bio::TreeIO(-file => ">trees.tre",
			       -format => 'newick');
foreach my $aln ( @alns ) {
# NOTE NOTE NOTE 
# This interface is probably going to change from create_tree to
# next_tree per some discussions I'm having with Shawn - we may need
# to tweak any scripts before you publish

# also may move the create_distance_matrix method around some
# and need to write in the switched support for Molphy's ProtML
    
    my $matrix = $dist_factory->create_distance_matrix($aln);
    my @seqnames = keys %$matrix; 
    open my $MATRIX, '>', "Group$count.dist" or die "Could not write file 'Group$count.dist': $!\n";
    printf $MATRIX "%4d\n",scalar @seqnames;
    for(my $i =0; $i< (scalar @seqnames  - 1); $i++ ) {
	printf $MATRIX "%-12s     ", $seqnames[$i];
	for( my $j = $i+1; $j < scalar @seqnames; $j++ ) {
	    print $MATRIX $matrix->{$seqnames[$i]}->{$seqnames[$j]},"  ";
	}
	print $MATRIX "\n";
    }
    close $MATRIX;

    my $tree = $tree_factory->create_tree("Group$count.dist");
    $outtrees->write_tree($tree);   
    $count++;
}

=head1 NAME

tree_from_seqsearch - builds a phylogenetic tree based on a sequence
search (FastA,BLAST,HMMER)

=head1 DESCRIPTION

This script requires that the bioperl-run pkg be also installed.

=cut
