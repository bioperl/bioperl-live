# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

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
    $NUMTESTS = 537;
    plan tests => $NUMTESTS;
}

## End of black magic.

my ($biblio, $count, $str, @args);
my ($citation, $provider);

my $format = ($ENV{'TEST_DETAILS'} ? "\t%-45s" : '');

print "Testing 'use Bio::Biblio:: ...'\n";

eval { require Bio::Biblio::Article };
print sprintf ($format, "use Bio::Biblio::Article"); ok (%Bio::Biblio::Article::);
print $@ if $@;

eval { require Bio::Biblio::Book };
print sprintf ($format, "use Bio::Biblio::Book"); ok (%Bio::Biblio::Book::);
print $@ if $@;

eval { require Bio::Biblio::BookArticle };
print sprintf ($format, "use Bio::Biblio::BookArticle"); ok (%Bio::Biblio::BookArticle::);
print $@ if $@;

eval { require Bio::Biblio::Journal };
print sprintf ($format, "use Bio::Biblio::Journal"); ok (%Bio::Biblio::Journal::);
print $@ if $@;

eval { require Bio::Biblio::JournalArticle };
print sprintf ($format, "use Bio::Biblio::JournalArticle"); ok (%Bio::Biblio::JournalArticle::);
print $@ if $@;

eval { require Bio::Biblio::MedlineArticle };
print sprintf ($format, "use Bio::Biblio::MedlineArticle"); ok (%Bio::Biblio::MedlineArticle::);
print $@ if $@;

eval { require Bio::Biblio::MedlineBook };
print sprintf ($format, "use Bio::Biblio::MedlineBook"); ok (%Bio::Biblio::MedlineBook::);
print $@ if $@;

eval { require Bio::Biblio::MedlineBookArticle };
print sprintf ($format, "use Bio::Biblio::MedlineBookArticle"); ok (%Bio::Biblio::MedlineBookArticle::);
print $@ if $@;

eval { require Bio::Biblio::MedlineJournal };
print sprintf ($format, "use Bio::Biblio::MedlineJournal"); ok (%Bio::Biblio::MedlineJournal::);
print $@ if $@;

eval { require Bio::Biblio::MedlineJournalArticle };
print sprintf ($format, "use Bio::Biblio::MedlineJournalArticle"); ok (%Bio::Biblio::MedlineJournalArticle::);
print $@ if $@;

eval { require Bio::Biblio::Organisation };
print sprintf ($format, "use Bio::Biblio::Organisation"); ok (%Bio::Biblio::Organisation::);
print $@ if $@;

eval { require Bio::Biblio::Patent };
print sprintf ($format, "use Bio::Biblio::Patent"); ok (%Bio::Biblio::Patent::);
print $@ if $@;

eval { require Bio::Biblio::Person };
print sprintf ($format, "use Bio::Biblio::Person"); ok (%Bio::Biblio::Person::);
print $@ if $@;

eval { require Bio::Biblio::Proceeding };
print sprintf ($format, "use Bio::Biblio::Proceeding"); ok (%Bio::Biblio::Proceeding::);
print $@ if $@;

eval { require Bio::Biblio::Provider };
print sprintf ($format, "use Bio::Biblio::Provider"); ok (%Bio::Biblio::Provider::);
print $@ if $@;

eval { require Bio::Biblio::Ref };
print sprintf ($format, "use Bio::Biblio::Ref"); ok (%Bio::Biblio::Ref::);
print $@ if $@;

eval { require Bio::Biblio::Service };
print sprintf ($format, "use Bio::Biblio::Service"); ok (%Bio::Biblio::Service::);
print $@ if $@;

eval { require Bio::Biblio::TechReport };
print sprintf ($format, "use Bio::Biblio::TechReport"); ok (%Bio::Biblio::TechReport::);
print $@ if $@;

eval { require Bio::Biblio::Thesis };
print sprintf ($format, "use Bio::Biblio::Thesis"); ok (%Bio::Biblio::Thesis::);
print $@ if $@;

eval { require Bio::Biblio::WebResource };
print sprintf ($format, "use Bio::Biblio::WebResource"); ok (%Bio::Biblio::WebResource::);
print $@ if $@;

eval { require Bio::Biblio::PubmedArticle };
print sprintf ($format, "use Bio::Biblio::PubmedArticle"); ok (%Bio::Biblio::PubmedArticle::);
print $@ if $@;

eval { require Bio::Biblio::PubmedBookArticle };
print sprintf ($format, "use Bio::Biblio::PubmedBookArticle"); ok (%Bio::Biblio::PubmedBookArticle::);
print $@ if $@;

eval { require Bio::Biblio::PubmedJournalArticle };
print sprintf ($format, "use Bio::Biblio::PubmedJournalArticle"); ok (%Bio::Biblio::PubmedJournalArticle::);
print $@ if $@;

print "Testing 'Bio::Biblio::->new() ...'\n";
foreach my $object (
		    qw(
		     Bio::Biblio::Article
		     Bio::Biblio::Book
		     Bio::Biblio::BookArticle
		     Bio::Biblio::Journal
		     Bio::Biblio::JournalArticle
		     Bio::Biblio::MedlineArticle
		     Bio::Biblio::MedlineBook
		     Bio::Biblio::MedlineBookArticle
		     Bio::Biblio::MedlineJournal
		     Bio::Biblio::MedlineJournalArticle
		     Bio::Biblio::Organisation
		     Bio::Biblio::Patent
		     Bio::Biblio::Person
		     Bio::Biblio::Proceeding
		     Bio::Biblio::Provider
		     Bio::Biblio::Ref
		     Bio::Biblio::Service
		     Bio::Biblio::TechReport
		     Bio::Biblio::Thesis
		     Bio::Biblio::WebResource
		     Bio::Biblio::PubmedArticle
		     Bio::Biblio::PubmedBookArticle
		     Bio::Biblio::PubmedJournalArticle
		       )) {
    print sprintf ($format, "new $object");
	ok defined ($biblio = new $object);
}

my @scalar_methods_for_ref =
    qw(
     abstract
     abstract_language
     abstract_type
     author_list_complete
     cross_references_list_complete
     date
     date_completed
     date_created
     date_revised
     format
     identifier
     language
     last_modified_date
     repository_subset
     rights
     spatial_location
     subject_headings_source
     temporal_period
     title
     toc
     toc_type
     type
     );
my @other_methods_for_ref =
    qw(
     authors
     cross_references
     codes
     contributors
     keywords
     publisher
     subject_headings
    );

my @scalar_methods_for_book =
    qw(
     edition
     isbn
     series
     volume
     );
my @other_methods_for_book =
    qw(
     editor
     );

my @scalar_methods_for_bookarticle =
    qw(
     );
my @other_methods_for_bookarticle =
    qw(
     book
     );

my @scalar_methods_for_article =
    qw(
     first_page
     last_page
     );
my @other_methods_for_article =
    qw(
     );

my @scalar_methods_for_journalarticle =
    qw(
     issue
     issue_supplement
     volume
     );
my @other_methods_for_journalarticle =
    qw(
     journal
     );

my @scalar_methods_for_medlinearticle =
    qw(
     affiliation
     citation_owner
     date_of_electronic_publication
     gene_symbols
     grant_list_complete
     medline_date
     medline_id
     medline_page
     number_of_references
     other_languages
     pmid
     season
     status
     vernacular_title
     );
my @other_methods_for_medlinearticle =
    qw(
     chemicals
     comment_ins
     comment_ons
     erratum_fors
     erratum_ins
     general_notes
     grants
     mesh_headings
     original_report_ins
     other_abstracts
     other_ids
     republished_froms
     republished_ins
     retraction_ins
     retraction_ofs
     summary_for_patients_ins
     update_ins
     update_ofs
     );

my @scalar_methods_for_medlinejournalarticle =
    qw(
     );
my @other_methods_for_medlinejournalarticle =
    qw(
     journal
     );

my @scalar_methods_for_medlinebookarticle =
    qw(
     );
my @other_methods_for_medlinebookarticle =
    qw(
     book
     );

my @scalar_methods_for_medlinebook =
    qw(
     );
my @other_methods_for_medlinebook =
    qw(
     );



my @scalar_methods_for_pubmedarticle =
    qw(
     pubmed_status
     pubmed_provider_id
     );
my @other_methods_for_pubmedarticle =
    qw(
     pubmed_history_list
     pubmed_article_id_list
     pubmed_url_list
     );


my @scalar_methods_for_journal =
    qw(
     abbreviation
     issn
     name
     );
my @other_methods_for_journal =
    qw(
     );

my @scalar_methods_for_medlinejournal =
    qw(
     coden
     country
     medline_code
     medline_ta
     nlm_unique_id
     );
my @other_methods_for_medlinejournal =
    qw(
     );

my @scalar_methods_for_patent =
    qw(
     doc_number
     doc_office
     doc_type
     );
my @other_methods_for_patent =
    qw(
     applicants
     );

my @scalar_methods_for_webresource =
    qw(
     url
     estimated_size
     cost
     );
my @other_methods_for_webresource =
    qw(
     );

my @scalar_methods_for_provider =
    qw(
     type
     );

my @scalar_methods_for_person =
    qw(
     affiliation
     email
     firstname
     forename
     initials
     lastname
     middlename
     postal_address
     suffix
     );

my @scalar_methods_for_organisation =
    qw(
     name
     );

my @scalar_methods_for_service =
    qw(
     name
     );

#
# Bio::Biblio::MedlineJournalArticle
#
print "Testing Bio::Biblio::MedlineJournalArticle ...\n";
$citation = Bio::Biblio::MedlineJournalArticle->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_ref,
                    @scalar_methods_for_article,
                    @scalar_methods_for_journalarticle,
                    @scalar_methods_for_medlinearticle,
                    @scalar_methods_for_medlinejournalarticle) {
    $str = 'string' . ($count++);
    print sprintf ($format, "set '$method' ");
	is $citation->$method ($str), $str;
    print sprintf ($format, "get '$method' ");
	is $citation->$method(), $str;
    push (@args, ("-$method" => $str));
}
print sprintf ($format, "set all attributes in a constructor");
ok defined ($biblio = Bio::Biblio::MedlineJournalArticle->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
    print sprintf ($format, "   $method");
	is $citation->$method(), $args[$i+1];
}
foreach my $method (@other_methods_for_ref,
                    @other_methods_for_article,
                    @other_methods_for_journalarticle,
                    @other_methods_for_medlinearticle,
                    @other_methods_for_medlinejournalarticle) {
    print sprintf ($format, "get '$method' ");
	is $citation->$method(), undef;
}
my ($me) = Bio::Biblio::Person->new(-lastname => 'me');
my ($you) = Bio::Biblio::Person->new(-lastname => 'you');
print sprintf ($format, "add_author 1");
ok $citation->add_author ($me);
print sprintf ($format, "add_author 2");
ok $citation->add_author ($you);
print sprintf ($format, "get authors");
is ${ $citation->authors }[1]->lastname, 'you';

print sprintf ($format, "add_contributor 1");
ok $citation->add_contributor ($me);
print sprintf ($format, "add_contributor 2");
ok $citation->add_contributor ($you);
print sprintf ($format, "get contributors");
is ${ $citation->contributors }[1]->lastname, 'you';

use Bio::Annotation::DBLink;
my $link1 = Bio::Annotation::DBLink->new(-database => 'here',
				        -primary_id => '001'
				        );
my $link2 = Bio::Annotation::DBLink->new(-database => 'there',
				        -primary_id => '002'
				        );
print sprintf ($format, "add_cross_reference 1");
ok $citation->add_cross_reference ($link1);
print sprintf ($format, "add_cross_reference 2");
ok $citation->add_cross_reference ($link2);
print sprintf ($format, "get cross_references");
is ${ $citation->cross_references }[0]->database, 'here';
print sprintf ($format, "get cross_references");
is ${ $citation->cross_references }[1]->primary_id, '002';


#
# Bio::Biblio::MedlineBookArticle
#
print "Testing Bio::Biblio::MedlineBookArticle ...\n";
$citation = Bio::Biblio::MedlineBookArticle->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_ref,
                    @scalar_methods_for_article,
                    @scalar_methods_for_bookarticle,
                    @scalar_methods_for_medlinearticle,
                    @scalar_methods_for_medlinebookarticle) {
    $str = 'string' . ($count++);
    print sprintf ($format, "set '$method' ");
	is $citation->$method ($str), $str;
    print sprintf ($format, "get '$method' ");
	is $citation->$method(), $str;
    push (@args, ("-$method" => $str));
}
print sprintf ($format, "set all attributes in a constructor");
ok defined ($biblio = Bio::Biblio::MedlineBookArticle->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
    print sprintf ($format, "   $method");
	is $citation->$method(), $args[$i+1];
}
foreach my $method (@other_methods_for_ref,
                    @other_methods_for_article,
                    @other_methods_for_bookarticle,
                    @other_methods_for_medlinearticle,
                    @other_methods_for_medlinebookarticle) {
    print sprintf ($format, "get '$method' ");
	is $citation->$method(), undef;
}


#
# Bio::Biblio::MedlineBook
#
print "Testing Bio::Biblio::MedlineBook ...\n";
$citation = Bio::Biblio::MedlineBook->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_ref,
                    @scalar_methods_for_book,
                    @scalar_methods_for_medlinebook) {
    $str = 'string' . ($count++);
    print sprintf ($format, "set '$method' ");
	is $citation->$method ($str), $str;
    print sprintf ($format, "get '$method' ");
	is $citation->$method(), $str;
    push (@args, ("-$method" => $str));
}
print sprintf ($format, "set all attributes in a constructor");
ok defined ($biblio = Bio::Biblio::MedlineBook->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
    print sprintf ($format, "   $method");
	is $citation->$method(), $args[$i+1];
}
foreach my $method (@other_methods_for_ref,
                    @other_methods_for_book,
                    @other_methods_for_medlinebook) {
    print sprintf ($format, "get '$method' ");
	is $citation->$method(), undef;
}

#
# Bio::Biblio::MedlineJournal
#
print "Testing Bio::Biblio::MedlineJournal ...\n";
$citation = Bio::Biblio::MedlineJournal->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_journal,
                    @scalar_methods_for_medlinejournal) {
    $str = 'string' . ($count++);
    print sprintf ($format, "set '$method' ");
	is $citation->$method ($str), $str;
    print sprintf ($format, "get '$method' ");
	is $citation->$method(), $str;
    push (@args, ("-$method" => $str));
}
print sprintf ($format, "set all attributes in a constructor");
ok defined ($biblio = Bio::Biblio::MedlineJournal->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
    print sprintf ($format, "   $method");
	is $citation->$method(), $args[$i+1];
}
foreach my $method (@other_methods_for_journal,
                    @other_methods_for_medlinejournal) {
    print sprintf ($format, "get '$method' ");
	is $citation->$method(), undef;
}

#
# Bio::Biblio::Patent
#
print "Testing Bio::Biblio::Patent ...\n";
$citation = Bio::Biblio::Patent->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_patent) {
    $str = 'string' . ($count++);
    print sprintf ($format, "set '$method' ");
	is $citation->$method ($str), $str;
    print sprintf ($format, "get '$method' ");
	is $citation->$method(), $str;
    push (@args, ("-$method" => $str));
}
print sprintf ($format, "set all attributes in a constructor");
ok defined ($biblio = Bio::Biblio::Patent->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
    print sprintf ($format, "   $method");
	is $citation->$method(), $args[$i+1];
}
foreach my $method (@other_methods_for_patent) {
    print sprintf ($format, "get '$method' ");
	is $citation->$method(), undef;
}

#
# Bio::Biblio::WebResource
#
print "Testing Bio::Biblio::WebResource ...\n";
$citation = Bio::Biblio::WebResource->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_webresource) {
    $str = 'string' . ($count++);
    print sprintf ($format, "set '$method' ");
	is $citation->$method ($str), $str;
    print sprintf ($format, "get '$method' ");
	is $citation->$method(), $str;
    push (@args, ("-$method" => $str));
}
print sprintf ($format, "set all attributes in a constructor");
ok defined ($biblio = Bio::Biblio::WebResource->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
    print sprintf ($format, "   $method");
	is $citation->$method(), $args[$i+1];
}
foreach my $method (@other_methods_for_webresource) {
    print sprintf ($format, "get '$method' ");
	is $citation->$method(), undef;
}


#
# Bio::Biblio::Person
#
print "Testing Bio::Biblio::Person ...\n";
$provider = Bio::Biblio::Person->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_provider,
                    @scalar_methods_for_person) {
    $str = 'string' . ($count++);
    print sprintf ($format, "set '$method' ");
	is $provider->$method ($str), $str;
    print sprintf ($format, "get '$method' ");
	is $provider->$method(), $str;
    push (@args, ("-$method" => $str));
}
print sprintf ($format, "set all attributes in a constructor");
ok defined ($biblio = Bio::Biblio::Person->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
    print sprintf ($format, "   $method");
	is $provider->$method(), $args[$i+1];
}

#
# Bio::Biblio::Organisation
#
print "Testing Bio::Biblio::Organisation ...\n";
$provider = Bio::Biblio::Organisation->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_provider,
                    @scalar_methods_for_organisation) {
    $str = 'string' . ($count++);
    print sprintf ($format, "set '$method' ");
	is $provider->$method ($str), $str;
    print sprintf ($format, "get '$method' ");
	is $provider->$method(), $str;
    push (@args, ("-$method" => $str));
}
print sprintf ($format, "set all attributes in a constructor");
ok defined ($biblio = Bio::Biblio::Organisation->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
    print sprintf ($format, "   $method");
	is $provider->$method(), $args[$i+1];
}

#
# Bio::Biblio::Service
#
print "Testing Bio::Biblio::Service ...\n";
$provider = Bio::Biblio::Service->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_provider,
                    @scalar_methods_for_organisation) {
    $str = 'string' . ($count++);
    print sprintf ($format, "set '$method' ");
	is $provider->$method ($str), $str;
    print sprintf ($format, "get '$method' ");
	is $provider->$method(), $str;
    push (@args, ("-$method" => $str));
}
print sprintf ($format, "set all attributes in a constructor");
ok defined ($biblio = Bio::Biblio::Service->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
    print sprintf ($format, "   $method");
	is $provider->$method(), $args[$i+1];
}

#
# Bio::Biblio::PubmedJournalArticle
#
print "Testing Bio::Biblio::PubmedJournalArticle ...\n";
$citation = Bio::Biblio::PubmedJournalArticle->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_pubmedarticle) {
    $str = 'string' . ($count++);
    print sprintf ($format, "set '$method' ");
	is $citation->$method ($str), $str;
    print sprintf ($format, "get '$method' ");
	is $citation->$method(), $str;
    push (@args, ("-$method" => $str));
}
print sprintf ($format, "set all attributes in a constructor");
ok defined ($biblio = Bio::Biblio::PubmedJournalArticle->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
    print sprintf ($format, "   $method");
	is $citation->$method(), $args[$i+1];
}
foreach my $method (@other_methods_for_pubmedarticle) {
    print sprintf ($format, "get '$method' ");
	is $citation->$method(), undef;
}
__END__
