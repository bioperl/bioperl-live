#!perl
use strict;
use warnings;

=head1 NAME 

bp_hmmer_to_table  - turn HMMER output into tabular format

=head1 SYNOPSIS

 bp_hmmer_to_table [-e evaluefilter] [-b bitscorefilter] [--header] [-o outfile] inputfile1 inputfile2 ... 

=head1 DESCRIPTION

Command line options:
  -e/--evalue   evalue    -- filter by evalue
  -b/--bitscore bitscore  -- filter by bitscore
  --header                -- boolean flag to print column header 
  -o/--out                -- optional outputfile to write data, 
                             otherwise will write to STDOUT
  -h/--help               -- show this documentation

Not technically a SearchIO script as this doesn't use any Bioperl
components but is a useful and fast.  The output is tabular output.

  query sequence/domain (these are flip-flopped for hmmsearch / hmmpfam)  
  query start
  query end
  domain/sequence name or PFAM accession 
  hit start
  hit end
  score
  e-value
  domain/sequence name (these are flip-flopped for hmmsearch / hmmpfam)  

=head1 AUTHOR - Jason Stajich

Jason Stajich jason_at_bioperl-dot-org

=cut

use Getopt::Long;

my ($evalue,$bitscore,$header,$outfile);
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

my @fields = qw(QNAME QSTART QEND HACCESSION HSTART HEND SCORE EVALUE HNAME);
if( $header ) {
    print $outfh join("\t", @fields), "\n";
}
my %dat;
while(<>) {
    if( s/^Query(\s+(sequence|HMM))?:\s+// ) {	
	s/\s+$//;
	$dat{'Query'} = $_;
    } elsif( /^Parsed for domains:/ ) {
	my $ready = 0;
	while(<>) {
	    if(/^Model|Sequence\s+Domain/ ) { $ready = 1; }
	    elsif( $ready && /^\-\-/) { $ready = 2; }
	    elsif( /^Alignments of/ ) { undef %dat; last; }
	    elsif( $ready == 2 ) {
		if( my ($n,$domainnum,$domainct, @vals) = 
		    (m!^(\S+)\s+      # domain name
		     (\d+)\/(\d+)\s+  # num/num (ie 1 of 2) 
		     (\d+)\s+(\d+).+? # sequence start and end
		     (\d+)\s+(\d+)\s+ # hmm start and end
		     \S+\s+           # []
		     (\S+)\s+         # score
		     (\S+)            # evalue
		     \s*$!ox)  ) {
		    next if( defined $bitscore && $vals[4] < $bitscore );
		    next if (defined $evalue && $vals[5] > $evalue);
		    print $outfh join("\t",
				      $dat{'Query'},
				      $vals[0], $vals[1],
				      $n,
				      $vals[2],$vals[3],
				      $vals[4],$vals[5],
				      $n),"\n";
		}
	    }
	}
    }
}
