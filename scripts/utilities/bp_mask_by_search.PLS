#!perl
# Author:  Jason Stajich <jason-at-bioperl-dot-org>


=head1 NAME

bp_mask_by_search - mask sequence(s) based on its alignment results

=head1 SYNOPSIS 

  bp_mask_by_search.pl -f blast genomefile blastfile.bls > maskedgenome.fa

=head1 DESCRIPTION

Mask sequence based on significant alignments of another sequence.
You need to provide the report file and the entire sequence data which
you want to mask.  By default this will assume you have done a TBLASTN
(or TFASTY) and try and mask the hit sequence assuming you've provided
the sequence file for the hit database.  If you would like to do the
reverse and mask the query sequence specify the -t/--type query flag.

This is going to read in the whole sequence file into memory so for
large genomes this may fall over.  I'm using DB_File to prevent
keeping everything in memory, one solution is to split the genome into
pieces (BEFORE you run the DB search though, you want to use the exact
file you BLASTed with as input to this program).

Below the double dash (--) options are of the form
--format=fasta or --format fasta
or you can just say
-f fasta

By -f/--format I mean either are acceptable options.  The =s or =n
or =c specify these arguments expect a 'string'

Options:
    -f/--format=s    Search report format (fasta,blast,axt,hmmer,etc)
    -sf/--sformat=s  Sequence format (fasta,genbank,embl,swissprot)
    --hardmask       (booelean) Hard mask the sequence
                     with the maskchar [default is lowercase mask]
    --maskchar=c     Character to mask with [default is N], change 
                     to 'X' for protein sequences
    -e/--evalue=n    Evalue cutoff for HSPs and Hits, only 
                     mask sequence if alignment has specified evalue 
                     or better
    -o/--out/
    --outfile=file   Output file to save the masked sequence to.
    -t/--type=s      Alignment seq type you want to mask, the 
                     'hit' or the 'query' sequence. [default is 'hit']
    --minlen=n       Minimum length of an HSP for it to be used 
                     in masking [default 0]
    -h/--help        See this help information

=head1 AUTHOR - Jason Stajich

Jason Stajich, jason-at-bioperl-dot-org.

=cut 


use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;
use Getopt::Long;
use Bio::Seq;
use DB_File;
# assuming tblastn or tfasty type alignment

my $format = 'blast';
my $sformat= undef;
my $evalue = undef;
my $type   = 'hit';
my $minlen = 50;
my $hardmask = 0; # mask with $maskchar instead of lowercase
my $maskchar = 'N'; # if we hard mask, mask with this cahr
my $outfile;
GetOptions(
	   'f|format:s'  => \$format,
	   'sf|sformat:s'=> \$sformat,
	   'hardmask'    => \$hardmask,
	   'maskchar:s'  => \$maskchar,
	   'e|evalue:s'  => \$evalue,
	   'o|out|outfile:s' => \$outfile,
	   't|type:s'    => \$type,
	   'minlen:s'    => \$minlen,
	   'h|help'      => sub { system('perldoc', $0);
				  exit; },
	   );
if( $type !~ /^(hit|query)/i ) {
    die("type must be query or hit[default] not $type") ;
}
$type = lc($type);

if(length($maskchar) > 1 ) {
    die("expected a mask character, not a string (you gave $maskchar)");
}
my $genomefile = shift || die('need a file containing the genome');
my $reportfile = shift;

# this could be problem for large genomes, figure out a 
# better way to do this later on
# or force people to split it up
my $genomeparser = new Bio::SeqIO(-file  => $genomefile,
				  -format=> $sformat);
my %seqs; 
unlink('/tmp/genome.idx');
tie(%seqs,'DB_File','/tmp/genome.idx');
while( my $seq = $genomeparser->next_seq ) {
    # should we pre-force to upper case?
    $seqs{$seq->display_id} = $seq->seq();
}

my $parser = new Bio::SearchIO(-file   => $reportfile,
			       -format => $format);

while( my $r = $parser->next_result ) {
    while( my $h = $r->next_hit ) {
	last if( defined $evalue && $h->significance > $evalue );
	my $hname = $h->name;
	if( ! $seqs{$hname} ) { 
	    die("Cannot find sequence $hname in genome seq");
	}
	while( my $hsp = $h->next_hsp ) {
	    last if( defined $evalue && $hsp->evalue > $evalue );
	    next if( $hsp->length('total') < $minlen);
	    my ($s,$len) = ( $hsp->$type()->start,
			     $hsp->$type()->length);
	    
	    if( $hardmask ) { 
		substr($seqs{$hname}, $s,$len, $maskchar x $len);
	    } else { 
		substr($seqs{$hname}, $s,$len, 
		       lc(substr($seqs{$hname}, $s,$len)));
	    }
	}
    }
}

my $out;
if( $outfile ) { 
    $out = new Bio::SeqIO(-file   => ">$outfile",
			  -format => $sformat);
} else { 
    $out = new Bio::SeqIO(-format => $sformat);
}

while( my ($seqname,$seq) = each %seqs ) {
    $out->write_seq(Bio::Seq->new(-seq        => $seq,
				  -display_id => $seqname,
				  -description=> 'MASKED'));
}
END { 
    unlink('/tmp/genome.idx');
}
