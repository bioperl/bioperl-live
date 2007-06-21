# This is -*-Perl-*- code
# $Id$

use strict;

BEGIN { 
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 25);
	
    use_ok('Bio::Root::IO');
	use_ok('Bio::Biblio');
	use_ok('Bio::Biblio::IO');
}

my $testnum;
my $verbose = test_debug();

my $testfile = Bio::Root::IO->catfile ('t','data','stress_test_medline.xml');
my $testfile2 = Bio::Root::IO->catfile ('t','data','stress_test_pubmed.xml');

# check 'new...'
SKIP: {
    test_skip(-tests => 1, -requires_modules => ['SOAP::Lite']);
	ok my $biblio = Bio::Biblio->new(-location => 'http://localhost:4567');
}

# check MEDLINE XML parser
my $io;
SKIP: {
	test_skip(-tests => 4, -requires_modules => ['XML::Parser']);
    
    ok defined ($io = Bio::Biblio::IO->new('-format' => 'medlinexml',
						 '-file'   => $testfile,
						 '-result' => 'raw'));
	
    print "Reading and parsing MEDLINE XML file...\n" if $verbose;
    is ($io->next_bibref->{'medlineID'}, 'Text1', 'citation 1');
    is ($io->next_bibref->{'medlineID'}, 'Text248', 'citation 2');
    is ($io->next_bibref->{'medlineID'}, 'Text495', 'citation 3');
}


print "Getting citations using callback...\n" if $verbose;
my (@ids) = ('Text1', 'Text248', 'Text495');
my $callback_used = 'no';
$io = Bio::Biblio::IO->new('-format'   => 'medlinexml',
			   '-file'     => $testfile,
			  #'-result'   => 'medline2ref',  # this is default
			   '-callback' => \&callback);

is ( $callback_used, 'yes', 'calling callback');

sub callback {
	my $citation = shift;
	$callback_used = 'yes';
	is ($citation->{'_identifier'}, shift @ids);
}

SKIP: {
	test_skip(-tests => 2, -requires_modules => ['XML::Parser']);
	
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
							   '-result'   => 'medline2ref');
	
	is ($io->next_bibref->{'_identifier'}, '12345678', 'citation 1');
	is ($io->next_bibref->{'_identifier'}, 'abcdefgh', 'citation 2');
}

SKIP: {
	test_skip(-tests => 2, -requires_modules => ['XML::Parser', 'IO::String']);
	
    print "Reading and parsing XML string handle...\n" if $verbose;
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
							   '-fh'     => IO::String->new ($data));
	is (eval { $io->next_bibref->identifier }, '87654321', 'citation 1');
	is (eval { $io->next_bibref->identifier }, 'hgfedcba', 'citation 2');
}

SKIP: {
	test_skip(-tests => 5, -requires_modules => ['XML::Parser']);
	
    # check PUBMED XML parser
    ok defined (eval { $io = Bio::Biblio::IO->new('-format' => 'pubmedxml',
                             '-file'   => $testfile2,
                             '-result' => 'pubmed2ref') });
	
    print "Reading and parsing PUBMED XML file...\n" if $verbose;
    
	is ($io->next_bibref->identifier, '11223344', 'citation 1');
	is ($io->next_bibref->identifier, '21583752', 'citation 2');
	is ($io->next_bibref->identifier, '21465135', 'citation 3');
	is ($io->next_bibref->identifier, '21138228', 'citation 4');
}

SKIP: {
    # test for FH
    my $fh;
    my @expvals = qw(11223344 21583752 21465135 21138228);
    print "Testing FH\n" if $verbose;
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

