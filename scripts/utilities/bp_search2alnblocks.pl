#!perl
# Author:      Jason Stajich <jason-at-bioperl-dot-org>
# Description: Turn SearchIO parseable report(s) into a set of Aligned blocks
#

=head1 NAME

bp_search2alnblocks - Turn SearchIO parseable reports(s) into a set of aligned blocks

=head1 SYNOPSIS

  bp_search2alnblocks --minid PERCENTID --minlen LEN --minevalue EVALUE file1.
  blast file2.blast ...> out.fas

=head1 DESCRIPTION

This script will parse and filter BLAST (or other formats
L<Bio::SearchIO> can parse) output and format the alignment as blocks
of alignments based on the HSPs.  Note this can only work if the input
file parsed contains the necessary.

Typically this can be used to turn BLAST output into a FASTA alignment format for input into the QRNA comparative gene finder for RNA genes (E.Rivas).

=head1 OPTIONS

 --maxevalue   Maximum E-value for an HSP
 --minevalue   Minimum E-value for an HSP 
 --minlen      Minimum length of an HSP [default 0] 
 --maxid       Maximum Percent Id [default 100]
               (to help remove sequences which are really close)
 --minid       Minimum Percent Identity for an HSP [default 0]
 -i/--input    An optional input filename (expects input on STDIN by default)
 -o/--output   An optional output filename (exports to STDOUT by default)
 -f/--format   Specify a different Search Alignment format- 
               {fasta, axt, waba, blast, blastxml} are all permitted
               although the format must have actual alignment 
               sequence for this script to work
               See L<Bio::SearchIO> for more information.
 -of/--outformat Output format for the alignment blocks, anything
               L<Bio::AlignIO> supports.
 -v/--verbose  Turn on debugging

=head1 AUTHOR - Jason Stajich

Jason Stajich, jason-at-bioperl-dot-org.

=cut


use strict;
use warnings;

use Bio::SearchIO;
use Bio::AlignIO;
use Math::BigFloat;
use Getopt::Long;

my $Usage = 'search2alnblocks --minid PERCENTID --minlen LEN --maxevalue EVALUE file1.blast file2.blast ... > blocks.fas';

my ($min_id,$min_len,$max_id,$max_len,$max_evalue,$min_evalue,$format,
    $outformat,$verbose,$input,$output);
$min_id  = 0;
$max_evalue = 0;
$min_evalue = undef;
$min_len = 0;
$format = 'blast';
$outformat= 'fasta';
GetOptions(
	   'minid:s'      => \$min_id,
	   'maxid:s'      => \$max_id,
	   'minlen:s'     => \$min_len,
	   'maxlen:s'     => \$max_len,
	   'minevalue:s'  => \$min_evalue,
	   'maxevalue:s'  => \$max_evalue,
	   'f|format:s'   => \$format,
	   'i|input:s'    => \$input,
	   'o|output:s'   => \$output,
	   'of|oformat:s' => \$outformat,
	   'v|verbose'    => \$verbose,
	   'h|help'       => sub { system('perldoc', $0);
				   exit(0) },
	   );
$max_evalue =~ s/^e/1e/;

# if $input is undef then will read from STDIN
my $parser =  new Bio::SearchIO(-format => $format,
				-file   => $input,
				-verbose=> $verbose);
my $out;

if( $output ) {
    $out = new Bio::AlignIO(-format => $outformat,
			    -file   => ">$output");
} else { 
    $out = new Bio::AlignIO(-format => $outformat);
}

my $id = 1;
while( my $r = $parser->next_result ) {
    while( my $hit = $r->next_hit ) {
	while( my $hsp = $hit->next_hsp ) {
	    my $hsplen = $hsp->length('total');
	    next if( $min_len && $hsplen < $min_len);
	    my $pid = $hsp->percent_identity;
	    next if( ($min_id && $pid < $min_id) || 
		     ($max_id && $pid > $max_id ) );
	    next if( defined $min_evalue && 
		     $hsp->evalue > $min_evalue ); 
	    next if( $max_evalue && 
		     $hsp->evalue < $max_evalue);	    
	    $verbose && $hsp->verbose($verbose);	    
	    my $aln = $hsp->get_aln();
	    my @seqs;
	    foreach my $seq ( $aln->each_seq ) {
		if( $seq->display_id =~ /(\S+)[\/\_](\d+)\-(\d+)/ ) {
		    $seq->display_id($1);
		    $seq->start($seq->start + $2 - 1);
		    $seq->end($seq->end + $2 - 1);
		}
		$seq->description(sprintf("PID=%.2f LEN=%d HSP=%d",
					  $pid,$hsplen,$id));
		push @seqs, $seq;
	    }
	    $aln = new Bio::SimpleAlign();
	    $aln->add_seq(shift @seqs);
	    $aln->add_seq(shift @seqs);
	    
	    $id++;
	    $out->write_aln($aln);
	}
    }
}
