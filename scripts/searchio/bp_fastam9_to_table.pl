#!/usr/bin/perl

=head1 NAME 

fastm9_to_table  - turn FASTA -m 9 output into NCBI -m 9 tabular output

=head1 SYNOPSIS

 fastm9_to_table [-e evaluefilter] [-b bitscorefilter] [--header] [-o outfile] inputfile1 inputfile2 ... 

=head1 DESCRIPTION

Command line options:
  -e/--evalue   evalue    -- filter by evalue
  -b/--bitscore bitscore  -- filter by bitscore
  --header                -- boolean flag to print column header 
  -o/--out                -- optional outputfile to write data, 
                             otherwise will write to STDOUT
  -h/--help               -- show this documentation

Not technically a SearchIO script as this doesn't use any Bioperl
components but is a useful and fast.  The output is tabular output
with the standard NCBI -m9 columns.

 queryname
 hit name
 percent identity
 alignment length
 number mismatches 
 number gaps
 query start  (if on rev-strand start > end)
 query end 
 hit start (if on rev-strand start > end)
 hit end 
 evalue
 bit score

Additionally 3 more columns are provided
 fasta score
 sw-score
 percent similar
 query length
 hit length
 query gaps
 hit gaps

=head1 AUTHOR - Jason Stajich

Jason Stajich jason_at_bioperl-dot-org

=cut

use strict;
use warnings;
use Getopt::Long;
my $hitsection = 0;
my %data;

my ($evalue,$bitscore,$header,$outfile) = ( 10, 0);
GetOptions(
	   'b|bitscore|bits:f'   => \$bitscore,
	   'e|evalue:f'          => \$evalue,
	   'header'              => \$header,
	   'o|out|outfile:s'     => \$outfile,
	   'h|help'              => sub { exec('perldoc',$0); exit; }
	   );

my $outfh;
if( $outfile ) { 
    open $outfh, '>', $outfile or die "Could not write file '$outfile': $!\n";
} else { 
    $outfh = \*STDOUT; 
}
	   
# query start -- an0
# query en    -- ax0
# hit start   -- an1
# hit end     -- ax1

my @fields = qw(qname hname percid alen mmcount gapcount 
		qstart qend hstart hend evalue score bits fs sw-score 
		percsim qlen hlen qgap hgap);
print $outfh "#",uc(join("", map{ sprintf("%-10s",$_) } @fields)), "\n" if $header;

while(<>) {
    my $linestr = $_;
    if( /^\s*\d+>>>(\S+).+/ ) {	
	$data{'qname'} = $1;
	if( /\-?\s+(\d+)\s+(aa|nt)\s+$/ ){
	    $data{'qlen'} = $1;
	}
    } elsif( $hitsection && (/^>>>\Q$data{'qname'}/ || /^>>>/) ) {	
	$hitsection = 0;
    } elsif( /^The best scores are:/ ) {	       
	$hitsection = 1;
    } elsif( /^\s+$/ ) {
    } elsif( $hitsection ) {
	if( s/^(\S+)\s+(.+)\(\s*(\d+)\)\s+// ) {
	    my ($hit, $desc,$hitlen) = ($1,$2,$3);
	    my ($dir) = ( s/^\[(r|f)\]\s+// );
	    my @line = split(/\s+/,$_);
	    $data{'hname'}    = $hit;
	    $data{'hlen'}     = $hitlen;
	    $data{'score'}    = shift @line;
	    $data{'bits'}     = shift @line;
	    $data{'evalue'}   = shift @line;
	    $data{'percid'}   = shift @line;
	    
	    $data{'percsim'}  = shift @line;
	    $data{'sw-score'} = shift @line;
	    $data{'alen'}     = shift @line;
	    $data{'qstart'}   = shift @line;
	    $data{'qend'}     = shift @line;
	    $data{'pn0'}      = shift @line; # pn0
	    $data{'px0'}      = shift @line; # px0
	    $data{'hstart'}   = shift @line; # an1 
	    $data{'hend'}     = shift @line; # ax1
	    $data{'pn1'}      = shift @line; # pn1
	    $data{'px1'}      = shift @line; # px1
	    # query + hit gaps
	    $data{'qgap'}     = shift @line;
	    $data{'hgap'}     = shift @line;
	    $data{'gapcount'} = $data{'qgap'} + $data{'hgap'};
	    $data{'fs'}       = shift @line;
	    
	    $data{'mmcount'} = $data{'alen'} - ( int($data{'percid'} * $data{'alen'}) + $data{'gapcount'});
	    #next if( $data{'evalue'} > $evalue || 
	#	     $data{'bits'} < $bitscore );
	    
	    for ( $data{'percid'}, $data{'percsim'} ) {
		$_ = sprintf("%.2f",$_*100);
	    }
	    print $outfh join( "\t",map { $data{$_} } @fields),"\n";
	} else { 
	    # print STDERR "unrecognized line \n$linestr";
	}
    } else { 
	# warn("skipping a line like this: $_");
    }
}
