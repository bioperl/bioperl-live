#!perl -w
# $Id: pairwise_kaks.PLS,v 1.2 2003-01-10 20:04:12 jason Exp $
use strict;
# Author Jason Stajich <jason@bioperl.org>

=head1 NAME

pairwise_kaks - script to calculate pairwise Ka,Ks for a set of sequences

=head1 SYNOPSIS

pairwise_kaks.PLS -i t/data/worm_fam_2785.cdna [-f fasta/genbank/embl...] [-msa tcoffee/clustal] [-kaks yn00/codeml]

=head1 DESCRIPTION 

  This script will take as input a dataset of cDNA sequences verify
 that they contain no stop codons, align them in protein space,
 project the alignment back into cDNA and estimate the Ka
 (non-synonymous) and Ks (synonymous) substitutions based on the ML
 method of Yang with the PAML package.

 Requires:
 * bioperl-run package
 * PAML program codeml or yn00
 * Multiple sequence alignment programs Clustalw OR T-Coffee

 Often there are specific specific parameters you want to run when you
 a computing Ka/Ks ratios so consider this script a starting point and
 do not rely it on for every situation.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

 http://bugzilla.bioperl.org/

=head1 AUTHOR

 Jason Stajich <jason@bioperl.org>

=cut

# Ka/Ks estimators
use Bio::Tools::Run::Phylo::PAML::Codeml;
use Bio::Tools::Run::Phylo::PAML::Yn00;

# Multiple Sequence Alignment programs
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Tools::Run::Alignment::TCoffee;

# for projecting alignments from protein to R/DNA space
use Bio::Align::Utilities qw(aa_to_dna_aln);

# for input of the sequence data
use Bio::SeqIO;
use Bio::AlignIO;

# for the command line argument parsing
use Getopt::Long;

my ($aln_prog, $kaks_prog,$format,
    $cdna,$verbose,$help) = qw(clustalw codeml fasta);

GetOptions(
	   'i|input:s'      => \$cdna,
	   'f|format:s'     => \$format,
	   'msa:s'          => \$aln_prog,
	   'kaks:s'         => \$kaks_prog,
	   'v|verbose'      => \$verbose,
	   'h|help'         => \$help,
	   );

if( $help ) {
    exec('perldoc',$0);
    exit(0);
}
$verbose = -1 unless $verbose;
my ($aln_factory,$kaks_factory);
if( $aln_prog =~ /clus/i ) {
    $aln_factory = Bio::Tools::Run::Alignment::Clustalw->new(-verbose => $verbose);
} elsif( $aln_prog =~ /t\_?cof/i ) {
    $aln_factory = Bio::Tools::Run::Alignment::TCoffee->new(-verbose => $verbose);
} else { 
    warn("Did not provide either 'clustalw' or 'tcoffee' as alignment program names");
    exit(0);
}
unless( $aln_factory->executable ) {
    warn("Could not find the executable for $aln_prog, make sure you have installed it and have either set ".\U$aln_prog."DIR or it is in your PATH");
    exit(0);
}


if( $kaks_prog =~ /yn00/i ) {
    $kaks_factory = Bio::Tools::Run::Phylo::PAML::Yn00->new(-verbose => $verbose);
} elsif( $kaks_prog =~ /codeml/i ) {
    # change the parameters here if you want to tweak your Codeml running!
    $kaks_factory = Bio::Tools::Run::Phylo::PAML::Codeml->new
	(-verbose => $verbose,
	 -params => { 'runmode' => -2,
		      'seqtype' => 1,
		  }
	 );
}
unless ( $kaks_factory->executable ) {
    warn("Could not find the executable for $kaks_prog, make sure you have installed it and you have defined PAMLDIR or it is in your PATH");
    exit(0);
}

unless ( $cdna && -f $cdna && -r $cdna &&  ! -z $cdna ) {
    warn("Did not specify a valid cDNA sequence file as input"); 
}

my $seqin = new Bio::SeqIO(-file   => $cdna, 
			   -format => $format);

my %seqs;
my @prots;
while( my $seq = $seqin->next_seq ) {
    $seqs{$seq->display_id} = $seq;
    my $protein = $seq->translate();
    my $pseq = $protein->seq();
    if( $pseq =~ /\*/ && 
	$pseq !~ /\*$/ ) {
	warn("provided a cDNA sequence with a stop codon, PAML will choke!");
	exit(0);
    }
    # Tcoffee can't handle '*'
    $pseq =~ s/\*//g;
    $protein->seq($pseq);
    push @prots, $protein;
}

my $aa_aln = $aln_factory->align(\@prots);
my $dna_aln = &aa_to_dna_aln($aa_aln, \%seqs);

$kaks_factory->alignment($dna_aln);
my ($rc,$parser) = $kaks_factory->run();
my $result = $parser->next_result;
my $MLmatrix = $result->get_MLmatrix();

my @otus = $result->get_seqs();
print join("\t", qw(SEQ1 SEQ2 Ka Ks Ka/Ks)), "\n";
for( my $i = 0; $i < (scalar @otus -1) ; $i++) {
    for( my $j = $i+1; $j < scalar @otus; $j++ ) {
	print join("\t",  
		   $otus[$i]->display_id,
		   $otus[$j]->display_id,$MLmatrix->[$i]->[$j]->{'dN'},
		   $MLmatrix->[$i]->[$j]->{'dS'},
		   $MLmatrix->[$i]->[$j]->{'omega'}), "\n";
    }
}
