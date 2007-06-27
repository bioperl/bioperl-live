# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

our @modules;

BEGIN {
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 537);
	
	@modules = qw(
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
		Bio::Biblio::PubmedJournalArticle);
	
	foreach my $module (@modules) {
		use_ok($module);
	}
}

my $verbose = test_debug();

my ($biblio, $count, $str, @args);
my ($citation, $provider);

print "Testing 'Bio::Biblio::->new() ...'\n" if $verbose;
foreach my $object (@modules) {
	ok defined ($biblio = $object->new());
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
print "Testing Bio::Biblio::MedlineJournalArticle ...\n" if $verbose;
$citation = Bio::Biblio::MedlineJournalArticle->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_ref,
                    @scalar_methods_for_article,
                    @scalar_methods_for_journalarticle,
                    @scalar_methods_for_medlinearticle,
                    @scalar_methods_for_medlinejournalarticle) {
    $str = 'string' . ($count++);
	is $citation->$method ($str), $str, "set '$method'";
	is $citation->$method(), $str, "get '$method'";
    push (@args, ("-$method" => $str));
}

ok defined ($biblio = Bio::Biblio::MedlineJournalArticle->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
	is $citation->$method(), $args[$i+1], $method;
}
foreach my $method (@other_methods_for_ref,
                    @other_methods_for_article,
                    @other_methods_for_journalarticle,
                    @other_methods_for_medlinearticle,
                    @other_methods_for_medlinejournalarticle) {
	is $citation->$method(), undef, "get '$method'";
}
my ($me) = Bio::Biblio::Person->new(-lastname => 'me');
my ($you) = Bio::Biblio::Person->new(-lastname => 'you');
ok $citation->add_author ($me), "add_author 1";
ok $citation->add_author ($you), "add_author 2";
is ${ $citation->authors }[1]->lastname, 'you', "get authors";

ok $citation->add_contributor ($me), "add_contributor 1";
ok $citation->add_contributor ($you), "add_contributor 2";
is ${ $citation->contributors }[1]->lastname, 'you', "get contributors";

use Bio::Annotation::DBLink;
my $link1 = Bio::Annotation::DBLink->new(-database => 'here',
				        -primary_id => '001'
				        );
my $link2 = Bio::Annotation::DBLink->new(-database => 'there',
				        -primary_id => '002'
				        );

ok $citation->add_cross_reference ($link1), "add_cross_reference 1";
ok $citation->add_cross_reference ($link2), "add_cross_reference 2";
is ${ $citation->cross_references }[0]->database, 'here', "get cross_references";
is ${ $citation->cross_references }[1]->primary_id, '002', "get cross_references";


#
# Bio::Biblio::MedlineBookArticle
#
print "Testing Bio::Biblio::MedlineBookArticle ...\n" if $verbose;
$citation = Bio::Biblio::MedlineBookArticle->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_ref,
                    @scalar_methods_for_article,
                    @scalar_methods_for_bookarticle,
                    @scalar_methods_for_medlinearticle,
                    @scalar_methods_for_medlinebookarticle) {
    $str = 'string' . ($count++);
	is $citation->$method ($str), $str, "set '$method'";
	is $citation->$method(), $str, "get '$method'";
    push (@args, ("-$method" => $str));
}

ok defined ($biblio = Bio::Biblio::MedlineBookArticle->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
	is $citation->$method(), $args[$i+1], $method;
}
foreach my $method (@other_methods_for_ref,
                    @other_methods_for_article,
                    @other_methods_for_bookarticle,
                    @other_methods_for_medlinearticle,
                    @other_methods_for_medlinebookarticle) {
	is $citation->$method(), undef, "get '$method'";
}


#
# Bio::Biblio::MedlineBook
#
print "Testing Bio::Biblio::MedlineBook ...\n" if $verbose;
$citation = Bio::Biblio::MedlineBook->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_ref,
                    @scalar_methods_for_book,
                    @scalar_methods_for_medlinebook) {
    $str = 'string' . ($count++);
	is $citation->$method ($str), $str, "set '$method'";
	is $citation->$method(), $str, "get '$method'";
    push (@args, ("-$method" => $str));
}

ok defined ($biblio = Bio::Biblio::MedlineBook->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
	is $citation->$method(), $args[$i+1], $method;
}
foreach my $method (@other_methods_for_ref,
                    @other_methods_for_book,
                    @other_methods_for_medlinebook) {
	is $citation->$method(), undef, "get '$method'";
}

#
# Bio::Biblio::MedlineJournal
#
print "Testing Bio::Biblio::MedlineJournal ...\n" if $verbose;
$citation = Bio::Biblio::MedlineJournal->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_journal,
                    @scalar_methods_for_medlinejournal) {
    $str = 'string' . ($count++);
	is $citation->$method ($str), $str, "set '$method'";
	is $citation->$method(), $str, "get '$method'";
    push (@args, ("-$method" => $str));
}

ok defined ($biblio = Bio::Biblio::MedlineJournal->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
	is $citation->$method(), $args[$i+1], $method;
}
foreach my $method (@other_methods_for_journal,
                    @other_methods_for_medlinejournal) {
	is $citation->$method(), undef, "get '$method'";
}

#
# Bio::Biblio::Patent
#
print "Testing Bio::Biblio::Patent ...\n" if $verbose;
$citation = Bio::Biblio::Patent->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_patent) {
    $str = 'string' . ($count++);
	is $citation->$method ($str), $str, "set '$method'";
	is $citation->$method(), $str, "get '$method'";
    push (@args, ("-$method" => $str));
}

ok defined ($biblio = Bio::Biblio::Patent->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
	is $citation->$method(), $args[$i+1], $method;
}
foreach my $method (@other_methods_for_patent) {
	is $citation->$method(), undef, "get '$method'";
}

#
# Bio::Biblio::WebResource
#
print "Testing Bio::Biblio::WebResource ...\n" if $verbose;
$citation = Bio::Biblio::WebResource->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_webresource) {
    $str = 'string' . ($count++);
	is $citation->$method ($str), $str, "set '$method'";
	is $citation->$method(), $str, "get '$method'";
    push (@args, ("-$method" => $str));
}

ok defined ($biblio = Bio::Biblio::WebResource->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
	is $citation->$method(), $args[$i+1], $method;
}
foreach my $method (@other_methods_for_webresource) {
	is $citation->$method(), undef, "get '$method'";
}


#
# Bio::Biblio::Person
#
print "Testing Bio::Biblio::Person ...\n" if $verbose;
$provider = Bio::Biblio::Person->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_provider,
                    @scalar_methods_for_person) {
    $str = 'string' . ($count++);
	is $provider->$method ($str), $str, "set '$method'";
	is $provider->$method(), $str, "get '$method'";
    push (@args, ("-$method" => $str));
}

ok defined ($biblio = Bio::Biblio::Person->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
	is $provider->$method(), $args[$i+1], $method;
}

#
# Bio::Biblio::Organisation
#
print "Testing Bio::Biblio::Organisation ...\n" if $verbose;
$provider = Bio::Biblio::Organisation->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_provider,
                    @scalar_methods_for_organisation) {
    $str = 'string' . ($count++);
	is $provider->$method ($str), $str, "set '$method'";
	is $provider->$method(), $str, "get '$method'";
    push (@args, ("-$method" => $str));
}

ok defined ($biblio = Bio::Biblio::Organisation->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
	is $provider->$method(), $args[$i+1], $method;
}

#
# Bio::Biblio::Service
#
print "Testing Bio::Biblio::Service ...\n" if $verbose;
$provider = Bio::Biblio::Service->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_provider,
                    @scalar_methods_for_organisation) {
    $str = 'string' . ($count++);
	is $provider->$method ($str), $str, "set '$method'";
	is $provider->$method(), $str, "get '$method'";
    push (@args, ("-$method" => $str));
}

ok defined ($biblio = Bio::Biblio::Service->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
	is $provider->$method(), $args[$i+1], $method;
}

#
# Bio::Biblio::PubmedJournalArticle
#
print "Testing Bio::Biblio::PubmedJournalArticle ...\n" if $verbose;
$citation = Bio::Biblio::PubmedJournalArticle->new();
@args = ();
$count = 1;
foreach my $method (@scalar_methods_for_pubmedarticle) {
    $str = 'string' . ($count++);
	is $citation->$method ($str), $str, "set '$method'";
	is $citation->$method(), $str, "get '$method'";
    push (@args, ("-$method" => $str));
}

ok defined ($biblio = Bio::Biblio::PubmedJournalArticle->new(@args));
for (my $i = 0; $i < @args; $i += 2) {
    my $method = substr ($args[$i], 1);
	is $citation->$method(), $args[$i+1], $method;
}
foreach my $method (@other_methods_for_pubmedarticle) {
	is $citation->$method(), undef, "get '$method'";
}
