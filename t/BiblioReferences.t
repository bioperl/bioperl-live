# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

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

    plan tests => 23 }

use Bio::Biblio::TechReport;
ok(1);

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my($citation, $citation2, $me, $you, @persons, $link1, $link2);

ok $citation = new Bio::Biblio::TechReport (-id => '123');
ok $citation->id, '123';
ok $citation->id('123'), '123';

my @string_methods = qw(type title rights language format date spatial_location 
		     temporal_period last_modified repository_subset abstract 
		     abstract_type toc toc_type 
		     ); 

foreach my $method (@string_methods) {
    ok $citation->$method('string'), 'string';
    ok $citation->$method, 'string';
}

use Bio::Biblio::Person;
ok 1;
ok $me = new Bio::Biblio::Person (-name => 'me');
ok $you = new Bio::Biblio::Person (-name => 'you');

ok $citation->publisher($me);
ok $citation->publisher->name(), 'me';

ok $citation->add_author($me);
$citation->add_author($you);

@persons = $citation->each_author;
ok $persons[1]->name, 'you';


ok $citation->add_contributor($me);
$citation->add_contributor($you);

@persons = $citation->each_contributor;
ok $persons[1]->name, 'you';

use Bio::Annotation::DBLink;
$link1 = new Bio::Annotation::DBLink(-database => 'here',
				     -primary_id => '001'
				     );
$link2 = new Bio::Annotation::DBLink(-database => 'there',
				     -primary_id => '002'
				     );

ok $citation->add_DBLink($link1);
$citation->add_DBLink($link2);

@persons = $citation->each_DBLink;
ok $persons[1]->primary_id, '002';

# all Ref methods tested

use Bio::Biblio::Proceeding;
ok 1;
ok $citation = new Bio::Biblio::Proceeding (-title => 'this');

use Bio::Biblio::Thesis;
ok 1;
ok $citation = new Bio::Biblio::Thesis (-language => 'English');

use Bio::Biblio::WebResource;
ok 1;
ok $citation = new Bio::Biblio::WebResource (-format => 'XML');


ok $citation->url('http://my.own.domain/'), 'http://my.own.domain/';
ok $citation->url, 'http://my.own.domain/';

ok $citation->size(12345), 12345;
ok $citation->size, 12345;


ok $citation->cost('1 euro'), '1 euro';
ok $citation->cost, '1 euro';



use Bio::Biblio::Patent;
ok 1;
ok $citation = new Bio::Biblio::Patent (-rights => 'no');


ok $citation->doc_number('123456789'), '123456789';
ok $citation->doc_number, '123456789';

ok $citation->doc_office('EPO'), 'EPO';
ok $citation->doc_office, 'EPO';

ok $citation->doc_type('plain'), 'plain';
ok $citation->doc_type, 'plain';

ok $citation->add_applicant($me);
$citation->add_applicant($you);

@persons = $citation->each_applicant;
ok $persons[1]->name, 'you';


use Bio::Biblio::Book;
ok 1;
ok $citation = new Bio::Biblio::Book (-series => 'Masterpieces');
ok $citation->series, 'Masterpieces';
ok $citation->series('Rubbish'), 'Rubbish';

ok $citation->isbn('123-123-123'), '123-123-123';;
ok $citation->isbn(), '123-123-123';

ok $citation->edition('first'), 'first';
$citation->edition(), 'first';

ok $citation->volume('XX');
ok $citation->volume(), 'XX';

ok $citation->editor($me);
ok $citation->editor->name(), 'me';


use Bio::Biblio::Article;
ok 1;
ok $citation2 = new Bio::Biblio::Article (-start => 222);
ok $citation2->start, 222;
ok $citation2->start(333), 333;

ok $citation2->end(444), 444;;
ok $citation2->end(), 444;

use Bio::Biblio::BookArticle;
ok 1;
ok $citation2 = new Bio::Biblio::BookArticle (-start => 222);
ok $citation2->start, 222;
ok $citation2->start(333), 333;

ok $citation2->book($citation);
ok $citation2->book->edition, 'first';


use Bio::Biblio::Journal;
ok 1;
my $journal;
ok $journal = new Bio::Biblio::Journal (-name => 'Nature');
ok $journal->name, 'Nature';
ok $journal->name('Science'), 'Science';


use Bio::Biblio::JournalArticle;
ok 1;
ok $citation = new Bio::Biblio::JournalArticle (-id => '9999999');
ok $citation->id, '9999999';
ok $citation->id('000000'), '000000';

ok $citation->volume(123), 123;;
ok $citation->volume(), '123';

ok $citation->issue('12b'), '12b';
ok $citation->issue(), '12b';

ok $citation->issue_supplement('c'), 'c';
ok $citation->issue_supplement(), 'c';


ok $citation->volume('XX'); ok $citation->volume(), 'XX';


ok $citation->journal($journal);
ok $citation->journal->name(), 'Science';
