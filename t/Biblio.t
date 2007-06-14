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
    eval { require Test::More; };
    $error = 0;
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;
    plan tests => 27;
    use_ok('Bio::Root::IO');
}

my $testnum;
my $verbose = 0;

## End of black magic.

map {$_ = 0} my ($serror, $serror2, $ferror, $ferror2, $xerror);

my $format = ($ENV{'TEST_DETAILS'} ? '%-25s' : '');

unless (eval "require SOAP::Lite; 1;") {
    $serror = "SOAP::Lite not installed. Skipping some tests.";
}

unless (eval "require IO::String; 1;") {
    $serror2 = "IO::String not installed. Skipping some tests.";
}

unless (eval "require XML::Parser; 1;") {
    $xerror = "XML::Parser not installed. Skipping some tests.";
}

# these are always present in CVS; commenting out irrelevant bits
my $testfile = Bio::Root::IO->catfile ('t','data','stress_test_medline.xml');
#unless (-e $testfile) {
#    print STDERR "Cannot find testing data '$testfile'. Skipping some tests.\n";
#    $ferror = 1;
#}
my $testfile2 = Bio::Root::IO->catfile ('t','data','stress_test_pubmed.xml');
#unless (-e $testfile2) {
#    print STDERR "Cannot find testing data '$testfile2'. Skipping some tests.\n";
#    $ferror2 = 1;
#}

# check 'use ...'
require_ok('Bio::Biblio');
print sprintf ($format, 'use Bio::Biblio');
# I'm puzzled about the reasoning for the following test... cjf 3/7/2007
ok (%Bio::Biblio::);

# check 'new...'
my $biblio = Bio::Biblio->new(-location => 'http://localhost:4567');
print sprintf ($format, "Bio::Biblio->new() ");
SKIP: {
    skip($serror,1) if $serror;
    ok (defined $biblio);
}

# check 'use ...IO...'
require_ok('Bio::Biblio::IO');
print sprintf ($format, "use Bio::Biblio::IO ");
# I'm puzzled about the reasoning for the following test... cjf
ok (%Bio::Biblio::IO::);

my $io;

# check MEDLINE XML parser
print sprintf ($format, "Bio::Biblio::IO->new(1)");
SKIP:{
    skip($ferror || $xerror,4) if $ferror || $xerror;
    ok defined ($io = Bio::Biblio::IO->new('-format' => 'medlinexml',
						 '-file'   => $testfile,
						 '-result' => 'raw'));
    print $@ if $@;    
    print "Reading and parsing MEDLINE XML file...\n";
    print sprintf ($format, "    citation 1 ");
    is ($io->next_bibref->{'medlineID'}, 'Text1');
    print sprintf ($format, "    citation 2 ");
    is ($io->next_bibref->{'medlineID'}, 'Text248');
    print sprintf ($format, "    citation 3 ");
    is ($io->next_bibref->{'medlineID'}, 'Text495');
    print "Getting citations using callback...\n";
    my (@ids) = ('Text1', 'Text248', 'Text495');
    my $callback_used = 'no';
    $io = Bio::Biblio::IO->new('-format'   => 'medlinexml',
                   '-file'     => $testfile,
#			       '-result'   => 'medline2ref',  # this is default
                   '-callback' => \&callback);
    
    print sprintf ($format, "    calling callback ");
    is ( $callback_used, 'yes');

    sub callback {
        my $citation = shift;
        $callback_used = 'yes';
        print sprintf ($format, '    citation ' . (@ids+0) . ' ');
        is ($citation->{'_identifier'}, shift @ids);
    }

    print "Reading and parsing XML string...\n";
    if ($xerror) {
        skip($xerror,2);
    } else {
        $io = Bio::Biblio::IO->new('-format'   => 'medlinexml',
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
        print sprintf ($format, "    citation 1 ");
        is ($io->next_bibref->{'_identifier'}, '12345678');
        print sprintf ($format, "    citation 2 ");
        is ($io->next_bibref->{'_identifier'}, 'abcdefgh');
    }
    
    print "Reading and parsing XML string handle...\n";
    #use IO::String;
    if ($xerror || $serror2) {
        skip ($xerror || $serror2, 2);
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
    
        $io = Bio::Biblio::IO->new('-format' => 'medlinexml',
                       '-fh'     => IO::String->new ($data),
                       );
        print sprintf ($format, "    citation 1 ");
        is (eval { $io->next_bibref->identifier }, '87654321');
        print sprintf ($format, "    citation 2 ");
        is (eval { $io->next_bibref->identifier }, 'hgfedcba');
    }
    
    # check PUBMED XML parser
    print sprintf ($format, "Bio::Biblio::IO->new(2)");
    skip ($ferror2 || $xerror, 1) if $ferror2 || $xerror;
    ok defined (eval { $io = Bio::Biblio::IO->new('-format' => 'pubmedxml',
                             '-file'   => $testfile2,
                             '-result' => 'pubmed2ref') });
    print "Reading and parsing PUBMED XML file...\n";
    if ($xerror) {
        skip ("Can't read citation from PUBMED XML, $xerror", 4);
    } else {
        print sprintf ($format, "    citation 1 ");
        is ($io->next_bibref->identifier, '11223344');
        print sprintf ($format, "    citation 2 ");
        is ($io->next_bibref->identifier, '21583752');
        print sprintf ($format, "    citation 3 ");
        is ($io->next_bibref->identifier, '21465135');
        print sprintf ($format, "    citation 4 ");
        is ($io->next_bibref->identifier, '21138228');
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
            is($_->identifier,shift @expvals);
        }
    };
    if( $@) {
        skip("unable to use pubmedxml",4);
    }

}

__END__
