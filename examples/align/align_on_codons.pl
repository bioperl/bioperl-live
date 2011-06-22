#!/usr/bin/perl

use strict;
use vars qw($USAGE %VALIDALIGN $CODONSIZE);
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Getopt::Long;
use Bio::Tools::CodonTable;
use Carp;

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
-h/--help                 See this information
-f/--frame            Translation Frame (0,1,2) are valid (defaults to '0')
-ct/--codontable      Codon table to use (defaults to '1')
                        see perldoc Bio::PrimarySeq for more information
-i/--input            Input Filename (defaults to STDIN)
-o/--output           Output Filename (defaults to STDOUT)
-sf/--seqformat       Input format (defaults to FASTA/Pearson)
-af/--alignformat     Alignment output format (clustal,fasta,nexus,phylip,
		      msf,pfam,mase,meme,prodom,selex,stockholm)
-ap/--alignprog       ClustalW, TCoffee (currently only support 
					 local execution) 
-v/--verbose          Run in verbose mode
};

    %VALIDALIGN = ('clustalw' => 'Bio::Tools::Run::Alignment::Clustalw',
		   'tcoffee' => 'Bio::Tools::Run::Alignment::TCoffee');
}

my ($help, $input, $output);

my ($alignprog, $sformat, $aformat, $frame, $codontable, $verbose) 
  = ('clustalw', 'fasta', 'clustalw', 0, 1, 0);

GetOptions( 'h|help'            => \$help,
				'i|input:s'         => \$input,
				'o|output:s'        => \$output,
				'sf|seqformat:s'    => \$sformat,
				'af|alignformat:s'  => \$aformat,
				'ap|alignprog:s'    => \$alignprog,
				# for translate
				'f|frame:s'         => \$frame,
				'ct|codontable:s'   => \$codontable,
				'v|verbose'         => \$verbose,
			 );

if( $help ) { 
    die($USAGE);
}
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

my $table = new Bio::Tools::CodonTable();
while( my $seq = $seqio->next_seq ) {
    
	#    if( $frame == 0 && $alignprog eq 'tcoffee' ) {
	#	print "last codon is ",$seq->subseq($seq->length() -2,
	#					    $seq->length()), "\n";
	#	if( $table->is_ter_codon($seq->subseq($seq->length() -2,
	#					      $seq->length())) ) {
	#	    $seq->seq($seq->subseq(1,$seq->length() - 3));
	#	}
	#    }

	push @nucseqs, $seq;    
	push @protseqs, $seq->translate(-frame => $frame,
											   -codontable_id => $codontable );
}

if( @nucseqs <= 1 ) {
	die("Must specify > 1 sequence for alignment on codons");
}

# allow these to be tweaked by cmdline parameters at some point
my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM'); 

my $alignengine = $VALIDALIGN{$alignprog}->new('-verbose' => $verbose,
					       @params);

my $aln = $alignengine->align(\@protseqs);

my $dnaalign = new Bio::SimpleAlign;
my $seqorder = 0;
my $alnlen = $aln->length;
foreach my $seq ( $aln->each_seq ) {    
	my $newseq;
    
	foreach my $pos ( 1..$alnlen ) {
		my $loc = $seq->location_from_column($pos);
		my $dna = ''; 
		if( !defined $loc || $loc->location_type ne 'EXACT' ) {
			$dna = '---';
		} else {
			# to readjust to codon boundaries
			# end needs to be +1 so we can just multiply by CODONSIZE 
			# to get this
			my ($start,$end) = ((($loc->start - 1)*$CODONSIZE) +1,
									  ($loc->end)*$CODONSIZE);
			if( $start <=0 || $end > $nucseqs[$seqorder]->length() ) {
				print "start is ", $loc->start, " end is ", $loc->end, "\n";
				warn("codons don't seem to be matching up for $start,$end");
				$dna = '---';
			} else {
				$dna = $nucseqs[$seqorder]->subseq($start,$end);
			}
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

