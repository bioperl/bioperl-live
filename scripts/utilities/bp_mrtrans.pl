#!perl
use strict;

# Author:      Jason Stajich <jason-at-bioperl-dot-org>
# Description: Perl implementation of Bill Pearson's mrtrans
#              to project protein alignment back into cDNA coordinates
#

=head1 NAME

bp_mrtrans - implement a transformer of alignments from protein to mrna coordinates

=head1 SYNOPSIS

Usage:
  bp_mrtrans -i inputfile -o outputfile [-if input format] [-of output format] [-s cDNA sequence database]  [-sf cDNA sequence format] [-h]

=head1 DESCRIPTION

This script will convert a protein alignment back into a cDNA.  Loosely
based on Bill Pearson's mrtrans.

The options are:

   -o filename          - the output filename [default STDOUT]
   -of format           - output sequence format
                          (multiple sequence alignment)
                          [default phylip]
   -i filename          - the input filename [required]
   -if format           - input sequence format
                          (multiple sequence alignment)
                          [ default clustalw]
   -s --seqdb filename  - the cDNA sequence database file
   -sf --seqformat      - the cDNA seq db format (flatfile sequence format)
   -h                   - this help menu

=head1 AUTHOR

Jason Stajich, jason-at-bioperl-dot-org

=cut

use strict;
use warnings;
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

$usage = "usage: bp_mrtrans.pl -i prot_alignment -if align_format -o out_dna_align -of output_format -s cDNA_seqdb -sf fasta\n".
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
	   'h|help'        => sub{ exec('perldoc',$0);
				   exit(0)
				   },
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

__END__
