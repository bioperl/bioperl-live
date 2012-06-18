#
# BioPerl module Bio::Biblio::IO::medline2ref.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::IO::medline2ref - A converter of a raw hash to MEDLINE citations

=head1 SYNOPSIS

 # to be written

=head1 DESCRIPTION

 # to be written

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR

Martin Senger (senger@ebi.ac.uk)

=head1 COPYRIGHT

Copyright (c) 2002 European Bioinformatics Institute. All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

Here is the rest of the object methods.  Internal methods are preceded
with an underscore _.

=cut


# Let the code begin...


package Bio::Biblio::IO::medline2ref;

use strict;

use Bio::Biblio::MedlineJournal;
use Bio::Biblio::MedlineBook;
use Bio::Biblio::Provider;
use Bio::Biblio::Person;
use Bio::Biblio::Organisation;

use base qw(Bio::Root::Root);

# -----------------------------------------------------------------------------
sub new {
    my ($caller, @args) = @_;
    my $class = ref ($caller) || $caller;

    # object creation and blessing    
    my ($self) = $class->SUPER::new (@args);	
    
    # make a hashtable from @args
    my %param = @args;
    @param { map { lc $_ } keys %param } = values %param; # lowercase keys

    # copy all @args into this object (overwriting what may already be
    # there) - changing '-key' into '_key', and making keys lowercase
    my $new_key;
    foreach my $key (keys %param) {
	($new_key = $key) =~ s/^-/_/;
	$self->{ lc $new_key } = $param { $key };
    }

    # done
    return $self;
}

# ---------------------------------------------------------------------
#
#   Here is the core...
#
# ---------------------------------------------------------------------

sub _load_instance {
    my ($self, $source) = @_;

    #
    # MEDLINE has only JournalArticles and BookArticles
    # but we may create a general Ref if there is no attribute 'article'
    #
    my $result;
    my $article = $$source{'article'};
    if (defined $article) {
	if (defined $$article{'journal'}) {
	    $result = $self->_new_instance ('Bio::Biblio::MedlineJournalArticle');
	    $result->type ('JournalArticle');
	} elsif (defined $$article{'book'}) {
	    $result = $self->_new_instance ('Bio::Biblio::MedlineBookArticle');
	    $result->type ('BookArticle');
	} else {
	    $result->type ('MedlineArticle');
	}
    }
    $result = $self->_new_instance ('Bio::Biblio::Ref') unless defined $result;
    return $result;
}

sub convert {
   my ($self, $source) = @_;
   my $result = $self->_load_instance ($source);

   if (defined $result->type) {
       if ($result->type eq 'JournalArticle') {
	   &_convert_journal_article ($result, $source);
       } elsif ($result->type eq 'BookArticle') {
	   &_convert_book_article ($result, $source);
       } elsif ($result->type eq 'Article') {
	   &_convert_article ($result, $source);
       }
   }

   #
   # now do the attributes which are the same for all resource types
   #

   # ...identification is now by MedlineID but the trend is to replace
   # it by PMID (I have heard) theefore we keep both also separately
   # from the 'identifier'
   if (defined $$source{'medlineID'}) {
       $result->identifier ($$source{'medlineID'});
   } else {
       $result->identifier ($$source{'PMID'});
   }
   $result->pmid ($$source{'PMID'}) if defined $$source{'PMID'};
   $result->medline_id ($$source{'medlineID'}) if defined $$source{'medlineID'};

   # ...few others
   $result->citation_owner ($$source{'owner'}) if defined $$source{'owner'};
   $result->status ($$source{'status'}) if defined $$source{'status'};
   $result->number_of_references ($$source{'numberOfReferences'}) if defined $$source{'numberOfReferences'};

   # ...entry status of the citation in the repository
   my $date;
   if (defined $$source{'dateRevised'}) {
       $result->last_modified_date (&_convert_date ($$source{'dateRevised'}));
       $date = &_convert_date ($$source{'dateCreated'});
       $result->date_created ($date) if defined $date;
       $date = &_convert_date ($$source{'dateCompleted'});
       $result->date_completed ($date) if defined $date;
   } elsif (defined $$source{'dateCompleted'}) {
       $result->last_modified_date (&_convert_date ($$source{'dateCompleted'}));
       $date = &_convert_date ($$source{'dateCreated'});
       $result->date_created ($date) if defined $date;
   } elsif (defined $$source{'dateCreated'}) {
       $result->last_modified_date (&_convert_date ($$source{'dateCreated'}));
   }

   # ...put citation subsets in a comma-separated string
   if (defined $$source{'citationSubsets'}) {
       $result->repository_subset (join (',', @{ $$source{'citationSubsets'} }));
   }

   # ...MEDLINE's Comments & Corrections will be arrays of hashes
   if (defined $$source{'commentsCorrections'}) {
       my $corr = $$source{'commentsCorrections'};
       $result->comment_ons ($$corr{'commentOns'}) if defined $$corr{'commentOns'};
       $result->comment_ins ($$corr{'commentIns'}) if defined $$corr{'commentIns'};
       $result->erratum_ins ($$corr{'erratumIns'}) if defined $$corr{'erratumIns'};
       $result->erratum_fors ($$corr{'erratumFors'}) if defined $$corr{'erratumFors'};
       $result->original_report_ins ($$corr{'originalReportIns'}) if defined $$corr{'originalReportIns'};
       $result->republished_froms ($$corr{'republishedFroms'}) if defined $$corr{'republishedFroms'};
       $result->republished_ins ($$corr{'republishedIns'}) if defined $$corr{'republishedIns'};
       $result->retraction_ofs ($$corr{'retractionOfs'}) if defined $$corr{'retractionOfs'};
       $result->retraction_ins ($$corr{'retractionIns'}) if defined $$corr{'retractionIns'};
       $result->summary_for_patients_ins ($$corr{'summaryForPatientsIns'}) if defined $$corr{'summaryForPatientsIns'};
       $result->update_ins ($$corr{'updateIns'}) if defined $$corr{'updateIns'};
       $result->update_ofs ($$corr{'updateOfs'}) if defined $$corr{'updateOfs'};
   }

   # ...MEDLINE's GeneSymbols are put in a comma-separated string
   if (defined $$source{'geneSymbols'}) {
       $result->gene_symbols (join (',', @{ $$source{'geneSymbols'} }));
   }

   # ...MEDLINE's GeneralNotes into an array of hashtables, each one
   # having keys for the 'owner' and the 'note'
   $result->general_notes ($$source{'generalNotes'}) if defined $$source{'generalNotes'};

   # ...MEDLINE's PersonalNameSubjects into contributors (TBD: is that correct?)
   if (defined $$source{'personalNameSubjects'}) {
       my @contributors;
       foreach my $person ( @{ $$source{'personalNameSubjects'} } ) {
	   push (@contributors, &_convert_personal_name ($person));
       }
       $result->contributors (\@contributors);
   }

   # ...MEDLINE's OtherAbstract into an array of hashtables, each one
   # having keys for the 'type', 'AbstractText' and the 'copyright'
   $result->other_abstracts ($$source{'otherAbstracts'}) if defined $$source{'otherAbstracts'};
#   if (defined $$source{'otherAbstracts'}) {
#	my @other_abstracts = ();
#	foreach my $oa ( @{ $$source{'otherAbstracts'} } ) {
#	    if (defined $$oa{'abstractText'}) {
#		my $abstract = $$oa{'abstractText'};
#		delete $$oa{'abstractText'};
#		$$oa{'abstract'} = $$abstract{'abstractText'};
#		$$oa{'rights'} = $$abstract{'copyrightInformation'} if defined $$abstract{'copyrightInformation'};
#		push (@other_abstracts, $oa);
#	    }
#	}
#	$result->other_abstracts (\@other_abstracts);
#    }

   # ...MEDLINE's OtherIDs into an array of hashtables, each one
   # having keys for the 'id', and 'source'
   $result->other_ids ($$source{'otherIDs'}) if defined $$source{'otherIDs'};

   # ...MEDLINE's Chemicals - store them as an array of hashtables
   # (each one for each Chemical)
   $result->chemicals ($$source{'chemicals'}) if defined $$source{'chemicals'};

   # MeshHeadings are put on two places:
   # - a complete information in a property called "MeshHeadings", and
   # - only descriptors in the hashtable "subject_headings", together
   #   with the word "MeSH" in "subject_headings_source"
   if (defined $$source{'meshHeadings'}) {
       $result->mesh_headings ($$source{'meshHeadings'});
       my %subject_headings;
       foreach my $mesh ( @{ $$source{'meshHeadings'} } ) {
	   $subject_headings{ $$mesh{'descriptorName'} } = 1 if defined $$mesh{'descriptorName'};
       }
       if (%subject_headings) {
	   $result->subject_headings (\%subject_headings);
	   $result->subject_headings_source ('Mesh');
       }
   }

   # ...MEDLINE's keyword lists are merger all together (this may not
   # be good idea - but again the keywords are better accessible
   # -TBD?)
   if (defined $$source{'keywordLists'}) {
       my %keywords;
       foreach my $keywords ( @{ $$source{'keywordLists'} } ) {
	   if ($$keywords{'keywords'}) {
	       foreach my $keyword ( @{ $$keywords{'keywords'} } ) {
		   $keywords{$keyword} = 1;
	       }
	   }
       }
       $result->keywords (\%keywords) if %keywords;
   }

   # Done!
   return $result;
}

# load a module (given as a real module name, e.g. 'Bio::Biblio::MedlineJournalArticle'),
# call new() method on it, and return the instance returned by the new() method
sub _new_instance {
    my ($self, $module) = @_;
    my ($filename);
    ($filename = $module . '.pm') =~ s|\:\:|/|g;
    eval { require $filename; };
    $self->throw ("Loading error when trying '$filename'. $@\n") if $@;
    return $module->new;
}

#
# see OpenBQS specification (http://www.ebi.ac.uk/~senger/openbqs/) how
# a date should be coded;
# TBD: this can be improved - checking is missing, timezones,
#      converting to UTC...
# Also note that this routine does not convert 'medline_date' - it
# is stored in a separate attribute without ant conversion.
#
sub _convert_date {
    my ($date) = @_;
    return unless
	exists $$date{'year'} or
	    exists $$date{'month'} or
		exists $$date{'day'} or
		    exists $$date{'hour'} or
			exists $$date{'minute'} or
			    exists $$date{'second'};


    my $converted = (exists $$date{'year'} ? $$date{'year'} : '0000');

    if (exists $$date{'month'}) {
	$converted .= '-' . $$date{'month'};
    } elsif (exists $$date{'day'}) {
	$converted .= '-00';
    }

    if (exists $$date{'day'}) {
	$converted .= '-' . $$date{'day'};
    } elsif (exists $$date{'hour'}) {
	$converted .= '-00';
    }

    if (exists $$date{'hour'}) {
	$converted .= 'T' . $$date{'hour'} .
	    ':' . (exists $$date{'minute'} ? $$date{'minute'} : '00') .
		':' . (exists $$date{'second'} ? $$date{'second'} : '00') . 'Z';
    }
    return $converted;
}

# $person is a hash with persons attributes - we need to create and
# return a Bio::Biblio::Person object
sub _convert_personal_name {
    my ($person) = @_;
    foreach my $key (keys %$person) {
	$$person{"_$key"} = $$person{$key};
	delete $$person{$key};
    }
    Bio::Biblio::Person->new(%$person);
}

#
# takes journal article related attributes from $article and convert
# them into $result and at the end call _convert_article (which is
# shared with book article)
#
sub _convert_journal_article {
    my ($result, $source) = @_;
    my $article = $$source{'article'};

    # create and populate both a Journal and a resulting Article objects
    my $from_journal = $$article{'journal'};
    my $journal = Bio::Biblio::MedlineJournal->new();
    $journal->name ($$from_journal{'title'}) if defined $$from_journal{'title'};
    $journal->issn ($$from_journal{'iSSN'}) if defined $$from_journal{'iSSN'};
    $journal->abbreviation ($$from_journal{'iSOAbbreviation'}) if defined $$from_journal{'iSOAbbreviation'};
    $journal->coden ($$from_journal{'coden'}) if defined $$from_journal{'coden'};
    if (defined $$from_journal{'journalIssue'}) {
	my $issue = $$from_journal{'journalIssue'};
	$result->volume ($$issue{'volume'}) if defined $$issue{'volume'};
	$result->issue ($$issue{'issue'}) if defined $$issue{'issue'};

	if (defined $$issue{'pubDate'}) {
	    my $pub_date = $$issue{'pubDate'};
	    my $converted = &_convert_date ($pub_date);
	    $result->date ($converted) if defined $converted;

	    # Some parts of a MEDLINE date are stored just as properties
	    # because they have almost non-parseable format :-).
	    $result->medline_date ($$pub_date{'medlineDate'}) if defined $$pub_date{'medlineDate'};
	    $result->season ($$pub_date{'season'}) if defined $$pub_date{'season'};
	}
    }

    # ...some attributes are in journalInfo (which is outside of the article)
    my $journal_info = $$source{'journalInfo'};
    if (defined $journal_info) {
	$journal->country ($$journal_info{'country'}) if defined $$journal_info{'country'};
	$journal->medline_ta ($$journal_info{'medlineTA'}) if defined $$journal_info{'medlineTA'};
	$journal->medline_code ($$journal_info{'medlineCode'}) if defined $$journal_info{'medlineCode'};
	$journal->nlm_unique_id ($$journal_info{'nlmUniqueID'}) if defined $$journal_info{'nlmUniqueID'};
    }

    $result->journal ($journal);
    &_convert_article ($result, $source);
}

#
# takes book article related attributes from $article and convert
# them into $result and at the end call _convert_article (which is
# shared with journal article)
#
sub _convert_book_article {
    my ($result, $source) = @_;
    my $article = $$source{'article'};

    # create and populate both book and resulting article objects
    my $from_book = $$article{'book'};
    my $book = Bio::Biblio::MedlineBook->new();
    $book->title ($$from_book{'title'}) if defined $$from_book{'title'};
    $book->volume ($$from_book{'volume'}) if defined $$from_book{'volume'};
    $book->series ($$from_book{'collectionTitle'}) if defined $$from_book{'collectionTitle'};

    if (defined $$from_book{'pubDate'}) {
	my $pub_date = $$from_book{'pubDate'};
	my $converted = &_convert_date ($pub_date);
	$result->pub_date ($converted) if defined $converted;

	# Some parts of a MEDLINE date are stored just as properties
	# because they have almost non-parseable format :-).
	$result->medline_date ($$pub_date{'medlineDate'}) if defined $$pub_date{'medlineDate'};
	$result->season ($$pub_date{'season'}) if defined $$pub_date{'season'};
    }

    if (defined $$from_book{'publisher'}) {
	my $publisher = Bio::Biblio::Organisation->new();
	$publisher->name ($$from_book{'publisher'});
        $book->publisher ($publisher);
    }

    my @authors = &_convert_providers ($$from_book{'authors'});
    $book->authors (\@authors) if @authors;

    $result->book ($book);
    &_convert_article ($result, $source);
}

#
# takes from $source article related attributes and convert them into
# $article (these attributes are the same both for journal and book
# articles
#
sub _convert_article {
    my ($article, $source) = @_;
    my $from_article = $$source{'article'};

    $article->title ($$from_article{'articleTitle'}) if defined $$from_article{'articleTitle'};
    $article->affiliation ($$from_article{'affiliation'}) if defined $$from_article{'affiliation'};
    $article->vernacular_title ($$from_article{'vernacularTitle'}) if defined $$from_article{'vernacularTitle'};
    $article->date_of_electronic_publication
	($$from_article{'dateOfElectronicPublication'}) if defined $$from_article{'dateOfElectronicPublication'};

    if (defined $$from_article{'pagination'}) {
	my $pagination = $$from_article{'pagination'};
	$article->first_page ($$pagination{'startPage'}) if defined $$pagination{'startPage'};
	$article->last_page ($$pagination{'endPage'}) if defined $$pagination{'endPage'};
	$article->medline_page ($$pagination{'medlinePgn'}) if defined $$pagination{'medlinePgn'};
    }

    if (defined $$from_article{'abstract'}) {
	my $abstract = $$from_article{'abstract'};
	$article->abstract ($$abstract{'abstractText'}) if defined $$abstract{'abstractText'};
	$article->abstract_type ('text/plain');
	$article->rights ($$abstract{'copyrightInformation'}) if defined $$abstract{'copyrightInformation'};
    }

    if (defined $$from_article{'languages'}) {
	my $languages = $$from_article{'languages'};  # ref-array
	if ( @{ $languages } > 0) {
	    $article->language ( $$languages[0] );
	}
	if ( @{ $languages } > 1) {
	    $article->other_languages (join (',', @{ $languages }));
	}
    }

    my @authors = &_convert_providers ($$from_article{'authors'});
    if (@authors) {
	$article->authors (\@authors);
	$article->author_list_complete
	    ($$from_article{'authorListComplete'}) if defined $$from_article{'authorListComplete'};
    }

    # references to database entries are prefixed with database name
    # (separated by a slash)
    use Bio::Annotation::DBLink;
    if (defined $$from_article{'dataBanks'}) {
	my $databanks = $$from_article{'dataBanks'};  # a ref-array
	my @references;
	foreach my $bank ( @{ $databanks } ) {
	    my $db_name = $$bank{'dataBankName'};
	    if (defined $$bank{'accessionNumbers'}) {
		foreach my $accn ( @{ $$bank{'accessionNumbers'} } ) {
		    my $dblink = Bio::Annotation::DBLink->new(-primary_id => $accn);
		    $dblink->database ($db_name);   # it does not matter if it is undef
		    push (@references, $dblink);
		}
	    }
	}
	if (@references) {
	    $article->cross_references (\@references);
	    $article->cross_references_list_complete
		($$from_article{'dataBankListComplete'}) if defined $$from_article{'dataBankListComplete'};
	}
    }

    # grants are stored in an array of hashtables (each of the
    # hashtables has keys agency, grantID and acronym)
    $article->grants ($$from_article{'grants'}) if defined $$from_article{'grants'};
    $article->grant_list_complete
	    ($$from_article{'grantListComplete'}) if defined $$from_article{'grandListComplete'};

}

#
# takes a ref-array of providers - they can be persons or
# organisations, and returns an array of converted providers
#
sub _convert_providers {
    my ($providers) = @_;
    return () unless defined $providers;

    my @results;
    foreach my $provider ( @{ $providers } ) {
	if (defined $$provider{'personalName'}) {
	    my $converted = &_convert_personal_name ($$provider{'personalName'});
	    push (@results, $converted) if defined $converted;
	} elsif (defined $$provider{'collectiveName'}) {
	    push (@results, Bio::Biblio::Organisation->new(-name => $$provider{'collectiveName'}));
	} else {
            Bio::Biblio::Provider->new();
	}
    }
    return () unless @results;
    return @results;
}

1;
__END__
