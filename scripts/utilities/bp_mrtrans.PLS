#!/usr/bin/perl -w
# Author:      Jason Stajich <jason@bioperl.org>
# Description: Perl implementation of Bill Pearson's mrtrans 
#              to project protein alignment back into cDNA coordinates
#
use strict;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::AlignIO;
use Bio::SeqIO;
use Getopt::Long;

# TODO - finish documentation, 
#      - add support for extra options in output alignment formats
#        such as idnewline in phylip out to support Molphy input files

my ($iformat,$seqformat,$oformat,$seqdb,$input,$output) = ('clustalw','fasta',
							   'phylip');
my ($help,$usage);

$usage = "usage: bp_mrtrans.pl -i prot_alignment -if align_format -o out_dna_align -of output_format -s cDNA_seqdb -sf cDNA_seqdb\n".
"defaults: -if clustalw
           -of phylip
           -sf fasta\n";

GetOptions(
	   'if|iformat:s'  => \$iformat,
	   'i|input:s'     => \$input,
	   'o|output:s'    => \$output,
	   'of|outformat:s'=> \$oformat,
	   's|seqdb|db:s'  => \$seqdb,
	   'sf|seqformat:s'=> \$seqformat,
	   'h|help'        => sub{ die($usage) },
	   );


$input ||= shift;
$seqdb ||= shift;
$output ||= shift;
if( ! defined $seqdb ) { 
    die("cannot proceed without a valid seqdb\n$usage");
}
if( ! defined $input ) { 
    die("cannot proceed without an input file\n$usage");
}
my $indb = new Bio::SeqIO(-file => $seqdb,
			  -format=>$seqformat);
my %seqs;
while( my $seq = $indb->next_seq ) {
    $seqs{$seq->id} = $seq;
}

my $in = new Bio::AlignIO(-format => $iformat,
			  -file   => $input);
my $out = new Bio::AlignIO(-format => $oformat,
			   -idlength => 22,
			   -interleaved => 0,
			   defined $output ? (-file   => ">$output") : () );
			   
while( my $aln = $in->next_aln ) {
    my $dnaaln = aa_to_dna_aln($aln,\%seqs);
    $dnaaln->set_displayname_flat(1);
    $out->write_aln($dnaaln);
}
