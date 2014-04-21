#!perl
use strict;
use warnings;
# Author Jason Stajich <jason-at-bioperl-dot-org>

=head1 NAME

bp_pairwise_kaks - script to calculate pairwise Ka,Ks for a set of sequences

=head1 SYNOPSIS

bp_pairwise_kaks.PLS -i t/data/worm_fam_2785.cdna [-f fasta/genbank/embl...] [-msa tcoffee/clustal] [-kaks yn00/codeml]

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

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

 Jason Stajich jason-at-bioperl-dot-org

=cut

eval {
    # Ka/Ks estimators
    require Bio::Tools::Run::Phylo::PAML::Codeml;
    require Bio::Tools::Run::Phylo::PAML::Yn00;
    
    # Multiple Sequence Alignment programs
    require Bio::Tools::Run::Alignment::Clustalw;
    require Bio::Tools::Run::Alignment::TCoffee;
};
if( $@ ) {
    die("Must have bioperl-run pkg installed to run this script");
}
# for projecting alignments from protein to R/DNA space
use Bio::Align::Utilities qw(aa_to_dna_aln);

# for input of the sequence data
use Bio::SeqIO;
use Bio::AlignIO;

# for the command line argument parsing
use Getopt::Long;

my ($aln_prog, $kaks_prog,$format, $output,
    $cdna,$verbose,$help) = qw(clustalw codeml fasta);

GetOptions(
	   'i|input:s'      => \$cdna,
	   'o|output:s'     => \$output,
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
    warn("Could not find the executable for $aln_prog, make sure you have installed it and have either set ".uc($aln_prog)."DIR or it is in your PATH");
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
    exit(0);
}

my $seqin = new Bio::SeqIO(-file   => $cdna, 
			   -format => $format);

my %seqs;
my @prots;
while( my $seq = $seqin->next_seq ) {
    $seqs{$seq->display_id} = $seq;
    my $protein = $seq->translate();
    my $pseq = $protein->seq();
    
    $pseq =~ s/\*$//;
    if( $pseq =~ /\*/ ) {
	warn("provided a cDNA (".$seq->display_id.") sequence with a stop codon, PAML will choke!");
	exit(0);
    }
    # Tcoffee can't handle '*'
    $pseq =~ s/\*//g;
    $protein->seq($pseq);
    push @prots, $protein;
}
if( @prots < 2 ) {
    warn("Need at least 2 cDNA sequences to proceed");
    exit(0);
}

local * OUT;
if( $output ) {
    open(OUT, ">$output") || die("cannot open output $output for writing");
} else { 
    *OUT = *STDOUT;
}

my $aa_aln = $aln_factory->align(\@prots);
my $dna_aln = &aa_to_dna_aln($aa_aln, \%seqs);

my @each = $dna_aln->each_seq();


$kaks_factory->alignment($dna_aln);

my ($rc,$parser) = $kaks_factory->run();
if( $rc <= 0 ) { 
    warn($kaks_factory->error_string,"\n");
    exit;
}
my $result = $parser->next_result;

if ($result->version =~ m/3\.15/) {
	warn("This script does not work with v3.15 of PAML!  Please use 3.14 instead.");
	exit(0);
}

my $MLmatrix = $result->get_MLmatrix();

my @otus = $result->get_seqs();

my @pos = map { 
    my $c= 1;
    foreach my $s ( @each ) {
	last if( $s->display_id eq $_->display_id );
	$c++;
    }
    $c; 
} @otus; 

print OUT join("\t", qw(SEQ1 SEQ2 Ka Ks Ka/Ks PROT_PERCENTID CDNA_PERCENTID)), "\n";
for( my $i = 0; $i < (scalar @otus -1) ; $i++) {
    for( my $j = $i+1; $j < (scalar @otus); $j++ ) {
	my $sub_aa_aln = $aa_aln->select_noncont($pos[$i],$pos[$j]);
	my $sub_dna_aln = $dna_aln->select_noncont($pos[$i],$pos[$j]);
	print OUT join("\t",  
		       $otus[$i]->display_id,
		       $otus[$j]->display_id,$MLmatrix->[$i]->[$j]->{'dN'},
		       $MLmatrix->[$i]->[$j]->{'dS'},
		       $MLmatrix->[$i]->[$j]->{'omega'},
		       sprintf("%.2f",$sub_aa_aln->percentage_identity),
		       sprintf("%.2f",$sub_dna_aln->percentage_identity),
		       ), "\n";
    }
}
