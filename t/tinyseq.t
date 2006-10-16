# -*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error;
use strict;
use vars qw($DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'};


use strict;

BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    use vars qw($TESTCOUNT);
    $TESTCOUNT = 15;
    plan tests => $TESTCOUNT;
    
    $error  = 0;
    eval { require XML::Parser::PerlSAX; };
    if( $@ ) {
	print STDERR "XML::Parser::PerlSAX not loaded. This means game test cannot be executed. Skipping\n";
	foreach ( $Test::ntest..$TESTCOUNT ) {
	    skip('XML::Parser::PerlSAX installed',1);
	}
	$error = 1;
    } 
    # make sure we can load it, assuming that the prerequisites are really met

    if( $error == 0 ) {
	eval { require Bio::SeqIO::tinyseq; };
	if( $@ ) {
	    print STDERR "tinyseq.pm not loaded. This means tinyseq test cannot be executed. Skipping\n";
	    foreach ( $Test::ntest..$TESTCOUNT ) {
		skip('tinyseq.pm not loaded because XML::Writer not loaded',1);
	    }
	    $error = 1;
	} 
    }
}

if( $error == 1 ) {
    exit(0);
}

use Bio::SeqIO;
use Bio::Seq;

my $file = File::Spec->catfile(qw(t data NM_002253.tseq));
my $outfile = 'tinyseqout.xml';

my $instream = Bio::SeqIO->new( -file 		=> $file,
				-format		=> 'tinyseq' );

my $outstream = Bio::SeqIO->new( -file		=> ">$outfile",
				 -format	=> 'tinyseq' );

my $seq = $instream->next_seq;
ok(defined $seq);
ok(defined $seq->seq);
ok($seq->length, 5830);
ok($seq->accession_number,'NM_002253');
ok($seq->species);
ok($seq->species->binomial, 'Homo sapiens');
ok($seq->species->ncbi_taxid, 9606);   
$outstream->write_seq($seq);
undef $outstream;
#$outstream->close_writer;
 
ok(-e $outfile);

my $reread = Bio::SeqIO->new( -file 		=> $outfile,
			      -format		=> 'tinyseq' );

my $seq2 = $reread->next_seq;

ok($seq2);
ok($seq2->seq);
ok($seq2->length, 5830);
ok($seq2->accession_number, 'NM_002253');
ok($seq2->species);
ok($seq2->species->binomial, 'Homo sapiens');
ok($seq2->species->ncbi_taxid, 9606);   

unlink $outfile;
