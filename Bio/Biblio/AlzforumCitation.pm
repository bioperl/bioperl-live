package Bio::Biblio::AlzforumCitation;
use vars qw(@ISA @AUTO_ATTRIBUTES @OTHER_ATTRIBUTES %SYNONYMS);
use strict;
use Bio::Root::AutoClass;
@ISA = qw(Bio::Root::AutoClass);
BEGIN {
  @AUTO_ATTRIBUTES=qw(Country
		      DerivedPubDate
		      langID
		      ReaderChoice
		      CopyrightInfo
		      MileStone
		      PubDay
		      ARFRec
		      Affiliation
		      WFStatus
		      PubMonth
		      Initials
		      lastName
		      Pagination
		      Volume
		      Issue
		      citID
		      PubYear
		      IpID
		      IpDateCreated
		      PubStatus
		      SourceTypeID
		      MedlineTA
		      DateChanged
		      ISSN
		      DateCreated
		      ArticleTitle
		     );
  @OTHER_ATTRIBUTES=qw();
  %SYNONYMS=(pmid=>'IpID',
	     alzid=>'citID',
	     medline_page=>'Pagination',
	     medline_ta=>'MedlineTA',
	     pubmed_status=>'PubStatus',
	     title=>'ArticleTitle',
	     date=>'DerivedPubDate');
  Bio::Root::AutoClass::declare(__PACKAGE__,\@AUTO_ATTRIBUTES,\%SYNONYMS,'lower');
}



1;
