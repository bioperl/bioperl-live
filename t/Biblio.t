# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl Biblio.t'

use strict;
use vars qw($NUMTESTS);

my $error;

BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    $error = 0;
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 24;
}

my $testnum;
my $verbose = 0;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my $serror = 0; my $serror2 = 0;
my $ferror = 0; my $ferror2 = 0;
my $xerror = 0;

my $format = ($ENV{'TEST_DETAILS'} ? '%-25s' : '');

unless (eval "require SOAP::Lite; 1;") {
    print STDERR "SOAP::Lite not installed. Skipping some tests.\n";
    $serror = 1;
}

unless (eval "require IO::String; 1;") {
    print STDERR "IO::String not installed. Skipping some tests.\n";
    $serror2 = 1;
}

unless (eval "require XML::Parser; 1;") {
    print STDERR "XML::Parser not installed. Skipping some tests.\n";
    $xerror = 1;
}

use Bio::Root::IO;
my $testfile = Bio::Root::IO->catfile ('t','data','stress_test_medline.xml');
unless (-e $testfile) {
    print STDERR "Cannot find testing data '$testfile'. Skipping some tests.\n";
    $ferror = 1;
}
my $testfile2 = Bio::Root::IO->catfile ('t','data','stress_test_pubmed.xml');
unless (-e $testfile2) {
    print STDERR "Cannot find testing data '$testfile2'. Skipping some tests.\n";
    $ferror2 = 1;
}


# check 'use ...'
eval { require Bio::Biblio };
print sprintf ($format, 'use Bio::Biblio'); ok (%Bio::Biblio::);
print $@ if $@;

# check 'new...'
my $biblio;
eval { $biblio = new Bio::Biblio (-location => 'http://localhost:4567'); };
print sprintf ($format, "new Bio::Biblio "); skip ($serror, defined $biblio);
print $@ if $@;


# check 'use ...IO...'
eval { require Bio::Biblio::IO };
print sprintf ($format, "use Bio::Biblio::IO "); ok (%Bio::Biblio::IO::);

my $io;

# check MEDLINE XML parser
print sprintf ($format, "new Bio::Biblio::IO (1)");
skip ($ferror || $xerror,
      defined (eval { $io = new Bio::Biblio::IO ('-format' => 'medlinexml',
						 '-file'   => $testfile,
						 '-result' => 'raw') }));
print $@ if $@;

print "Reading and parsing MEDLINE XML file...\n";
print sprintf ($format, "    citation 1 "); skip ($ferror || $xerror, eval { $io->next_bibref->{'medlineID'} }, 'Text1');
print sprintf ($format, "    citation 2 "); skip ($ferror || $xerror, eval { $io->next_bibref->{'medlineID'} }, 'Text248');
print sprintf ($format, "    citation 3 "); skip ($ferror || $xerror, eval { $io->next_bibref->{'medlineID'} }, 'Text495');

print "Getting citations using callback...\n";
my (@ids) = ('Text1', 'Text248', 'Text495');
my $callback_used = 'no';
if ($ferror || $xerror) {
    foreach my $i (1..3) {
	print sprintf ($format, "    citation $i "); skip (1,1);
    }
} else {
    $io = new Bio::Biblio::IO ('-format'   => 'medlinexml',
			       '-file'     => $testfile,
#			       '-result'   => 'medline2ref',  # this is default
			       '-callback' => \&callback);
}
print sprintf ($format, "    calling callback "); skip ($ferror || $xerror, $callback_used, 'yes');

sub callback {
    my $citation = shift;
    $callback_used = 'yes';
    print sprintf ($format, '    citation ' . (@ids+0) . ' '); skip ($ferror, $citation->{'_identifier'}, shift @ids);
}

print "Reading and parsing XML string...\n";
if ($xerror) {
    print sprintf ($format, "    citation 1 "); skip (1, 1);
    print sprintf ($format, "    citation 2 "); skip (1, 1);
} else {
    $io = new Bio::Biblio::IO ('-format'   => 'medlinexml',
			       '-data'     => <<XMLDATA,
<MedlineCitationSet>
<MedlineCitation>
<MedlineID>12345678</MedlineID>
<Article><Journal></Journal></Article>
</MedlineCitation>
<MedlineCitation>
<MedlineID>abcdefgh</MedlineID>
<Article><Journal></Journal></Article>
</MedlineCitation>
</MedlineCitationSet>
XMLDATA
			       '-result'   => 'medline2ref',
			       );
    print sprintf ($format, "    citation 1 "); ok ($io->next_bibref->{'_identifier'}, '12345678');
    print sprintf ($format, "    citation 2 "); ok ($io->next_bibref->{'_identifier'}, 'abcdefgh');
}

print "Reading and parsing XML string handle...\n";
#use IO::String;
if ($xerror || $serror2) {
    print sprintf ($format, "    citation 1 "); skip (1,1);
    print sprintf ($format, "    citation 2 "); skip (1,1);
} else {
    my $data = <<XMLDATA;
<MedlineCitationSet>
<MedlineCitation>
<MedlineID>87654321</MedlineID>
<Article><Journal></Journal></Article>
</MedlineCitation>
<MedlineCitation>
<MedlineID>hgfedcba</MedlineID>
<Article><Journal></Journal></Article>
</MedlineCitation>
</MedlineCitationSet>
XMLDATA

    $io = new Bio::Biblio::IO ('-format' => 'medlinexml',
			       '-fh'     => IO::String->new ($data),
			       );
    print sprintf ($format, "    citation 1 "); ok (eval { $io->next_bibref->identifier }, '87654321');
    print sprintf ($format, "    citation 2 "); ok (eval { $io->next_bibref->identifier }, 'hgfedcba');
}

# check PUBMED XML parser
print sprintf ($format, "new Bio::Biblio::IO (2)");
skip ($ferror2 || $xerror,
      defined (eval { $io = new Bio::Biblio::IO ('-format' => 'pubmedxml',
						 '-file'   => $testfile2,
						 '-result' => 'pubmed2ref') }));
print "Reading and parsing PUBMED XML file...\n";
if ($xerror) {
    foreach my $i (1..4) {
	print sprintf ($format, "    citation $i "); skip (1,"Can't read citation from PUBMED XML");
    }
} else {
    print sprintf ($format, "    citation 1 "); skip ($ferror2, eval { $io->next_bibref->identifier }, '11223344');
    print sprintf ($format, "    citation 2 "); skip ($ferror2, eval { $io->next_bibref->identifier }, '21583752');
    print sprintf ($format, "    citation 3 "); skip ($ferror2, eval { $io->next_bibref->identifier }, '21465135');
    print sprintf ($format, "    citation 4 "); skip ($ferror2, eval { $io->next_bibref->identifier }, '21138228');
}

# test for FH
my $fh;
my @expvals = qw(11223344 21583752 21465135 21138228);
print "Testing FH\n";
eval { 
    $fh = Bio::Biblio::IO->newFh('-format' => 'pubmedxml',
				  '-file'   => $testfile2,
				  '-result' => 'pubmed2ref');
    while(<$fh>) {
	ok($_->identifier,shift @expvals);
    }
};
if( $@) {
    foreach ( 1..4 ) { skip(1,"unable to use pubmedxml"); }
}
__END__
