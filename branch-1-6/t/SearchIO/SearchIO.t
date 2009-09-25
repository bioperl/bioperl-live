# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 19);
	
	use_ok('Bio::SearchIO');
}

# only methods defined in Bio::SearchIO should be tested here; all
# parsers should have their own separate SearchIO_*.t test file

# Let's test if _guess_format is doing its job correctly
my %pair = ( 'filename.blast'  => 'blast',
	     'filename.bls'    => 'blast',
	     'f.blx'           => 'blast',
	     'f.tblx'          => 'blast',
	     'fast.bls'        => 'blast',
	     'f.fasta'         => 'fasta',
	     'f.fa'            => 'fasta',
	     'f.fx'            => 'fasta',
	     'f.fy'            => 'fasta',
	     'f.ssearch'       => 'fasta',
	     'f.SSEARCH.m9'    => 'fasta',
	     'f.m9'            => 'fasta',
	     'f.psearch'       => 'fasta',
	     'f.osearch'       => 'fasta',
	     'f.exon'          => 'exonerate',
	     'f.exonerate'     => 'exonerate',
	     'f.blastxml'      => 'blastxml',
	     'f.xml'           => 'blastxml');
while( my ($file,$expformat) = each %pair ) {
    is(Bio::SearchIO->_guess_format($file),$expformat, "$expformat for $file");
}
