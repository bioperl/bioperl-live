#!/usr/bin/perl -w

use strict;
use vars qw($USAGE %VALIDALIGN $CODONSIZE);
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Getopt::Long;
use Data::Dumper;

BEGIN {
    $CODONSIZE = 3; # parametrize everything like a good little programmer
    if( ! defined $ENV{'CLUSTALDIR'} ) { 
	$ENV{'CLUSTALDIR'} = '/usr/local/bin';
    } 
    if( ! defined $ENV{'TCOFFEEDIR'} ) { 
	$ENV{'TCOFFEEDIR'} = '/usr/local/bin';
    }
    $USAGE = 
qq{align_on_codons.pl < file.fa
-f/--frame            Translation Frame (0,1,2) are valid (defaults to '0')
-ct/--codontable      Codon table to use (defaults to '0')
                        see perldoc Bio::PrimarySeq for more information
-i/--input            Input Filename (defaults to STDIN)
-o/--output           Output Filename (defaults to STDOUT)
-sf/--seqformat       Input format (defaults to FASTA/Pearson)
-af/--alignformat     Alignment output format (clustal,fasta,nexus,phylip,
		      msf,pfam,mase,meme,prodom,selex,stockholm)
-ap/--alignprog       ClustalW, TCoffee (currently only support 
					 local execution) 
};

    %VALIDALIGN = ('clustalw' => 'Bio::Tools::Run::Alignment::Clustalw',
		   'tcoffee' => 'Bio::Tools::Run::Alignment::TCoffee');
}

my ($help, $input, $output);

my ($alignprog, $sformat, $aformat,$frame,
    $codontable) = ('clustalw', 'fasta', 'clustalw', 
		    0,0);

GetOptions( 'h|help'            => \$help,
	    'i|input:s'         => \$input,
	    'o|output:s'        => \$output,
	    'sf|seqformat:s'    => \$sformat,
	    'af|alignformat:s'  => \$aformat,
	    'ap|alignprog:s'    => \$alignprog,
	    
	    # for translate
	    'f|frame:s'         => \$frame,
	    'ct|codontable:s'   => \$codontable,
	    
	    );

if( ! $alignprog || !defined $VALIDALIGN{$alignprog} ) {
    die("Cannot use $alignprog as 'alignprog' parameter");  
} else {
    my $modname = $VALIDALIGN{$alignprog} .".pm";
    $modname =~ s/::/\//g;
    require $modname;
}

my $alignout;
if( $output ) {
    $alignout = new Bio::AlignIO('-format' => $aformat,
				 '-file'   => ">$output");
} else { 
    $alignout = new Bio::AlignIO('-format' => $aformat);
}

my (@nucseqs,@protseqs);
my $seqio;

if( $input ) {
    $seqio = new Bio::SeqIO('-format' => $sformat,
			    '-file'   => $input);
} else { 
    $seqio = new Bio::SeqIO('-format' => $sformat,
			    '-fh'     => \*STDIN);
}

while( my $seq = $seqio->next_seq ) {
    push @nucseqs, $seq;
    push @protseqs, $seq->translate(undef,undef,$frame,$codontable,0,0);
}

if( @nucseqs <= 1 ) {
    die("Must specify > 1 sequence for alignment on codons");
}

# allow these to be tweaked by cmdline parameters at some point
my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM'); 

my $alignengine = $VALIDALIGN{$alignprog}->new(@params);

my $aln = $alignengine->align(\@protseqs);

my $dnaalign = new Bio::SimpleAlign;
my $seqorder = 0;
my $alnlen = $aln->length;
foreach my $seq ( $aln->each_seq ) {    
    my $newseq;
    
    foreach my $pos ( 1..$alnlen ) {
	my $loc = $seq->location_from_column($pos);
	my $dna; 
	if( !defined $loc || $loc->location_type ne 'EXACT' ) {
	    $dna = '---';
	} else {
	    # to readjust to codon boundaries
	    # end needs to be +1 so we can just multiply by CODONSIZE 
	    # to get this
	    $dna = $nucseqs[$seqorder]->subseq
		((($loc->start - 1)*$CODONSIZE) +1,
		 ($loc->end)*$CODONSIZE);
	}
	$newseq .= $dna;
    }
    $seqorder++;
    # funky looking math is to readjust to codon boundaries and deal
    # with fact that sequence start with 1
    my $newdna = new Bio::LocatableSeq(-display_id  => $seq->id(),
				       -start => (($seq->start - 1) * 
						   $CODONSIZE) + 1, 
					-end   => ($seq->end * $CODONSIZE),
					-strand => $seq->strand,
					-seq   => $newseq);    
    $dnaalign->add_seq($newdna);
}

$alignout->write_aln($dnaalign);

