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
    $NUMTESTS = 15;
    plan tests => $NUMTESTS;
}

my $testnum;
my $verbose = 0;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my $serror = 0;
my $ferror = 0;
my $format = '%-25s';

use Bio::Root::IO;
my $testfile = Bio::Root::IO->catfile ('t','data','stress_test_medline.xml');
unless (-e $testfile) {
    print STDERR "Cannot find testing data '$testfile'. Skipping some tests.\n";
    $ferror = 1;
}

# check 'use ...'
eval { require Bio::Biblio };
print sprintf ($format, 'use Bio::Biblio'); skip ($error, %Bio::Biblio::);
print $@ if $@;

# check 'new...'
my $biblio;
print sprintf ($format, "new Bio::Biblio "); skip ($serror,
						   defined ($biblio = new Bio::Biblio (-location => 'http://localhost:4567')));

# check MEDLINE XML parser
eval { require Bio::Biblio::IO };
print sprintf ($format, "use Bio::Biblio::IO "); skip ($error, %Bio::Biblio::IO::);

my $io;
print sprintf ($format, "new Bio::Biblio::IO ");
skip ($ferror,
      defined (eval { $io = new Bio::Biblio::IO ('-format' => 'medlinexml',
						 '-file'   => $testfile,
						 '-result' => 'raw') }));
print "Reading and parsing XML file...\n";
print sprintf ($format, "    citation 1 "); skip ($ferror, eval { $io->next_bibref->{'medlineID'} }, 'Text1');
print sprintf ($format, "    citation 2 "); skip ($ferror, eval { $io->next_bibref->{'medlineID'} }, 'Text248');
print sprintf ($format, "    citation 3 "); skip ($ferror, eval { $io->next_bibref->{'medlineID'} }, 'Text495');

print "Getting citations using callback...\n";
my (@ids) = ('Text1', 'Text248', 'Text495');
my $callback_used = 'no';
unless ($ferror) {
    $io = new Bio::Biblio::IO ('-format'   => 'medlinexml',
			       '-file'     => $testfile,
#			       '-result'   => 'medline2ref',  # this is default
			       '-callback' => \&callback);
}
print sprintf ($format, "    calling callback "); skip ($ferror, $callback_used, 'yes');

sub callback {
    my $citation = shift;
    $callback_used = 'yes';
    print sprintf ($format, '    citation ' . (@ids+0) . ' '); skip ($ferror, $citation->{'_identifier'}, shift @ids);
}

print "Reading and parsing XML string...\n";
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

print "Reading and parsing XML string handle...\n";
use IO::String;
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
print sprintf ($format, "    citation 1 "); ok ($io->next_bibref->{'_identifier'}, '87654321');
print sprintf ($format, "    citation 2 "); ok ($io->next_bibref->{'_identifier'}, 'hgfedcba');

__END__
