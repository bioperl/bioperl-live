#!/usr/bin/perl -w 

#  $Id$
# Jason Stajich <jason@bioperl.org>
#
# Convert an input mRNA/cDNA sequence into protein
# 
# -f/--frame           Specifies frame [0,1,2]
#
# One can also specify 

# -t/--terminator      Stop Codon character (defaults to '*')
# -u/--unknown         Unknown Protein character (defaults to 'X')
# -ct/--codontable     Codon table to use (defaults to '0')
#                      see perldoc Bio::PrimarySeq for more information
# -cds/--fullcds       Expected Full CDS (with start and Stop codon)
# -throwOnError        Throw an error if no Full CDS
# -i/--input           Input Filename (defaults to STDIN)
# -l/--format          Input format (defaults to FASTA/Pearson)
# -o/--output          Output Filename (defaults to STDOUT)

use strict;
use Bio::SeqIO;
use Getopt::Long;

use vars qw($USAGE);

BEGIN { 
    $USAGE = 
qq{make_mrna_protein.pl < file.fa > file.prots
-f/--frame            Translation Frame (0,1,2) are valid (defaults to '0')
-t/--terminator	      Stop Codon Character ('*' by default)
-u/--unknown          Unknown Protein character (defaults to 'X')
-ct/--codontable      Codon table to use (defaults to '0')
                        see perldoc Bio::PrimarySeq for more information
-cds/--fullcds        Expected Full CDS (with start and Stop codon)
-throwOnError         Throw an error if no Full CDS
-i/--input            Input Filename (defaults to STDIN)
-l/--format           Input format (defaults to FASTA/Pearson)
-o/--output           Output Filename (defaults to STDOUT)
};	

}
my ($format,$frame,$termchar,$unknownProt,$codontable,$fullCDS,
    $throw_on_Incomp_CDS,$help) = ('fasta',0,undef, undef,0,0,0);
my ($input,$output);

GetOptions('f|frame:s'       => \$frame,
	   't|terminator:s'  => \$termchar,
	   'u|unknown:s'     => \$unknownProt,
	   'ct|codontable:s' => \$codontable,
	   'cds|fullcds'     => \$fullCDS,
	   'throwOnError'    => \$throw_on_Incomp_CDS,
	   'h|help'          => \$help,
	   'i|input:s'       => \$input,
	   'l|format:s'        => \$format,
	   'o|output:s'      => \$output,
	   );

die $USAGE if( $help );

my ($in,$out);
if( $input ) { 
    $in = new Bio::SeqIO('-format' => $format, '-file' => $input);
} else { 
    $in = new Bio::SeqIO('-format' => $format, '-fh' => \*STDIN);
}

if( $output ) { 
    $out = new Bio::SeqIO('-format' => $format, '-file' => ">$output" ); 
} else { 
    $out = new Bio::SeqIO('-format' => $format ); 
}

while( my $seq = $in->next_seq ) { 
    my $protseq = $seq->translate($termchar,$unknownProt, $frame,$codontable,
				  $fullCDS, $throw_on_Incomp_CDS);
    $out->write_seq($protseq);
}
