#
# BioPerl module Bio::Biblio::IO::medlinexml.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::IO::medlinexml - A converter of XML files with MEDLINE citations

=head1 SYNOPSIS

Do not use this object directly, it is recommended to access it and use
it through the I<Bio::Biblio::IO> module:

  use Bio::Biblio::IO;
  my $io = Bio::Biblio::IO->new(-format => 'medlinexml');

=head1 DESCRIPTION

This object reads bibliographic citations in XML/MEDLINE format and
converts them into I<Bio::Biblio::RefI> objects. It is an
implementation of methods defined in I<Bio::Biblio::IO>.

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

The main documentation details are to be found in
L<Bio::Biblio::IO>.

Here is the rest of the object methods.  Internal methods are preceded
with an underscore _.

=cut


# Let the code begin...


package Bio::Biblio::IO::medlinexml;
use vars qw(@Citations $Callback $Convert @ObjectStack @PCDataStack);
use vars qw(%PCDATA_NAMES %SIMPLE_TREATMENT %POP_DATA_AND_PEEK_OBJ %POP_OBJ_AND_PEEK_OBJ);
use vars qw(%POP_AND_ADD_ELEMENT %POP_AND_ADD_DATA_ELEMENT);

use strict;

use XML::Parser;

use base qw(Bio::Biblio::IO);

# -----------------------------------------------------------------------------

sub _initialize {
    my ($self, @args) = @_;
    
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

    # find the format for output - and put it into a global $Convert
    # because it will be used by the event handler who knows nothing
    # about this object
    my $result = $self->{'_result'} || 'medline2ref';
    $result = "\L$result";	# normalize capitalization to lower case

    # a special case is 'raw' when no converting module is loaded
    # and citations will be returned as a hashtable (the one which
    # is created during parsing XML file/stream)
    unless ($result eq 'raw') {

	# load module with output converter - as defined in $result
	if (defined &Bio::Biblio::IO::_load_format_module ($result)) {
	    $Convert = "Bio::Biblio::IO::$result"->new (@args);
	}
    }

    # create an instance of the XML parser
    # (unless it is already there...)
    $self->{'_xml_parser'} = new XML::Parser (Handlers => {Init  => \&handle_doc_start,
							   Start => \&handle_start,
							   End   => \&handle_end,
							   Char  => \&handle_char,
							   Final => \&handle_doc_end})
	unless $self->{'_xml_parser'};

    # if there is an argument '-callback' then start parsing at once -
    # the registered event handlers will use 'callback' to report
    # back after each citation
    #
    # we need to remember this situation also in a global variable
    # because the event handler subroutines know nothing about this
    # object (unfortunately)
    if ($Callback = $self->{'_callback'}) {
	$self->_parse;
    }
}

# -----------------------------------------------------------------------------

sub _parse {
    my ($self) = shift;


    if (defined $self->{'_file'}) {
	$self->{'_xml_parser'}->parsefile ($self->{'_file'});
    } elsif (defined $self->{'_fh'}) {
	my $fh = $self->{'_fh'};
	if (ref ($fh) and UNIVERSAL::isa ($fh, 'IO::Handler')) {
	    $self->{'_xml_parser'}->parse ($fh);
	} else {
	    my $data;
	    $data .= $_ while <$fh>;
	    $self->{'_xml_parser'}->parse ($data);
	}
    } elsif ($self->{'_data'}) {
	$self->{'_xml_parser'}->parse ($self->{'_data'});
    } else {
	$self->throw ("XML source to be parsed is unknown. Should be given in the new().");
    }

    # when parsing is done all citations have already been delivered
    # to the caller using her callbacks - and nothing to be stored
    # here, or parser put all citations into global @Cittaions where
    # we want to copy there into this instance - so any caller can
    # start parsing other XML input without overwriting already read
    # citations from the first parser
    if (@Citations) {
	$self->{'_citations'} = [];
	foreach my $cit (@Citations) {
	    push (@{ $self->{'_citations'} }, $cit);
	    undef $cit;
	}
	undef @Citations;
    }
}

# ---------------------------------------------------------------------
#
#   Here is an implementation of Bio::Biblio::IO methods
#
# ---------------------------------------------------------------------

# global variables used by the XML event handlers
# TBD: make them accessible at least ONLY from this module...
@Citations = ();
$Callback = undef;
$Convert = undef;
@ObjectStack = ();   # it has Hash-Ref elements
@PCDataStack = ();   # it has String elements

sub next_bibref {
   my ($self) = @_;
   $self->throw ("Method 'next_bibref' should not be called when a '-callback' argument given.")
       if $self->{'_callback'};

   # parse the whole input into memory (global @Citations)
   # and then copy it into this object
   $self->_parse unless $self->{'_citations'};

   # return the next citation (and forget it here)
   shift (@{ $self->{'_citations'} });
}

# ---------------------------------------------------------------------
#
#   Here are the event handlers (they do the real job!)
#
# Note that these methods do not know anything about the object they
# are part of - they are called as subroutines. not as methods.
# It also means that they need to use global variables to store and
# exchnage intermediate results.
#
# ---------------------------------------------------------------------

#
# This is a list of #PCDATA elements.
#
%PCDATA_NAMES = (
		 'AbstractText' => 1,
		 'AccessionNumber' => 1,
		 'Acronym' => 1,
		 'Affiliation' => 1,
		 'Agency' => 1,
		 'ArticleTitle' => 1,
		 'CASRegistryNumber' => 1,
		 'CitationSubset' => 1,
		 'Coden' => 1,
		 'CollectionTitle' => 1,
		 'CollectiveName' => 1,
		 'CopyrightInformation' => 1,
		 'Country' => 1,
		 'DataBankName' => 1,
		 'DateOfElectronicPublication' => 1,
		 'Day' => 1,
		 'Descriptor' => 1,
		 'DescriptorName' => 1,
		 'EndPage' => 1,
		 'FirstName' => 1,
		 'ForeName' => 1,
		 'GeneralNote' => 1,
		 'GeneSymbol' => 1,
		 'GrantID' => 1,
		 'Hour' => 1,
		 'ISOAbbreviation' => 1,
		 'ISSN' => 1,
		 'Initials' => 1,
		 'Issue' => 1,
		 'Keyword' => 1,
		 'Language' => 1,
		 'LastName' => 1,
		 'MedlineCode' => 1,
		 'MedlineDate' => 1,
		 'MedlineID' => 1,
		 'MedlinePgn' => 1,
		 'MedlineTA' => 1,
		 'MiddleName' => 1,
		 'Minute' => 1,
		 'Month' => 1,
		 'NameOfSubstance' => 1,
		 'NlmUniqueID' => 1,
		 'Note' => 1,
		 'NumberOfReferences' => 1,
		 'OtherID' => 1,
		 'PMID' => 1,
		 'PublicationType' => 1,
		 'Publisher' => 1,
		 'QualifierName' => 1,
		 'RefSource' => 1,
		 'RegistryNumber' => 1,
		 'Season' => 1,
		 'Second' => 1,
		 'SpaceFlightMission' => 1,
		 'StartPage' => 1,
		 'SubHeading' => 1,
		 'Suffix' => 1,
		 'Title' => 1,
		 'VernacularTitle' => 1,
		 'Volume' => 1,
		 'Year' => 1,
		 );

%SIMPLE_TREATMENT = (
		     'MeshHeading' => 1,
		     'Author' => 1,
		     'Article' => 1,
		     'Book' => 1,
		     'Investigator' => 1,
		     'Chemical' => 1,
		     'Pagination' => 1,
		     'MedlineJournalInfo' => 1,
		     'JournalIssue' => 1,
		     'Journal' => 1,
		     'DateCreated' => 1,
		     'DateCompleted' => 1,
		     'DateRevised' => 1,
		     'PubDate' => 1,
		     'Abstract' => 1,
		     'Grant' => 1,
		     'CommentsCorrections' => 1,
		     'CommentOn' => 1,
		     'CommentIn' => 1,
		     'ErratumFor' => 1,
		     'ErratumIn' => 1,
		     'OriginalReportIn' => 1,
		     'RepublishedFrom' => 1,
		     'RepublishedIn' => 1,
		     'RetractionOf' => 1,
		     'RetractionIn' => 1,
		     'SummaryForPatientsIn' => 1,
		     'UpdateIn' => 1,
		     'UpdateOf' => 1,
		     'DataBank' => 1,
		     'KeywordList' => 1,
		     'DeleteCitation' => 1,
		     );

%POP_DATA_AND_PEEK_OBJ = (
			  'Descriptor' => 1,
			  'DescriptorName' => 1,
			  'Year' => 1,
			  'Month' => 1,
			  'Day' => 1,
			  'LastName' => 1,
			  'Initials' => 1,
			  'FirstName' => 1,
			  'ForeName' => 1,
			  'NameOfSubstance' => 1,
			  'RegistryNumber' => 1,
			  'CASRegistryNumber' => 1,
			  'MiddleName' => 1,
			  'NlmUniqueID' => 1,
			  'MedlineTA' => 1,
			  'MedlinePgn' => 1,
			  'MedlineCode' => 1,
			  'Country' => 1,
			  'ISSN' => 1,
			  'ArticleTitle' => 1,
			  'Issue' => 1,
			  'AbstractText' => 1,
			  'VernacularTitle' => 1,
			  'GrantID' => 1,
			  'Agency' => 1,
			  'Acronym' => 1,
			  'MedlineDate' => 1,
			  'NumberOfReferences' => 1,
			  'RefSource' => 1,
			  'DataBankName' => 1,
			  'CopyrightInformation' => 1,
			  'Suffix' => 1,
			  'Note' => 1,
			  'CollectiveName' => 1,
			  'Hour' => 1,
			  'Minute' => 1,
			  'Second' => 1,
			  'Season' => 1,
			  'Coden' => 1,
			  'ISOAbbreviation' => 1,
			  'Publisher' => 1,
			  'CollectionTitle' => 1,
			  'DateOfElectronicPublication' => 1,
			  'StartPage' => 1,
			  'EndPage' => 1,
			  'Volume' => 1,
			  'Title' => 1,
			  );

%POP_OBJ_AND_PEEK_OBJ = (
			 'Pagination' => 1,
			 'JournalIssue' => 1,
			 'Journal' => 1,
			 'DateCreated' => 1,
			 'Article' => 1,
			 'DateCompleted' => 1,
			 'DateRevised' => 1,
			 'CommentsCorrections' => 1,
			 'Book' => 1,
			 'PubDate' => 1,
			 'Abstract' => 1,
			 );

%POP_AND_ADD_DATA_ELEMENT = (
			     'Keyword' => 'keywords',
			     'PublicationType' => 'publicationTypes',
			     'CitationSubset' => 'citationSubsets',
			     'Language' => 'languages',
			     'AccessionNumber' => 'accessionNumbers',
			     'GeneSymbol' => 'geneSymbols',
			     'SpaceFlightMission' => 'spaceFlightMissions',
			     );


%POP_AND_ADD_ELEMENT = (
			'OtherAbstract' => 'otherAbstracts',
			'Chemical' => 'chemicals',
			'KeywordList' => 'keywordLists',
			'Grant' => 'grants',
			'UpdateIn' => 'updateIns',
			'CommentOn' => 'commentOns',
			'CommentIn' => 'commentIns',
			'DataBank' => 'dataBanks',
			'PersonalNameSubject' => 'personalNameSubjects',
			'ErratumFor' => 'erratumFors',
			'ErratumIn' => 'erratumIns',
			'RepublishedFrom' => 'republishedFroms',
			'RepublishedIn' => 'republishedIns',
			'RetractionOf' => 'retractionOfs',
			'RetractionIn' => 'retractionIns',
			'UpdateOf' => 'updateOfs',
			'OriginalReportIn' => 'originalReportIns',
			'SummaryForPatientsIn' => 'summaryForPatientsIns',
			'MeshHeading' => 'meshHeadings',
			);

sub handle_doc_start {
    @Citations = ();
    @ObjectStack = ();
    @PCDataStack = ();
}

sub handle_doc_end {
    undef @ObjectStack;
    undef @PCDataStack;
}

sub handle_char {
    my ($expat, $str) = @_;

    # this may happen with whitespaces between tags;
    # but because I have not created an entry for data on the stack
    # I can also ignore such data, can't I
    return if $#PCDataStack < 0;

    $PCDataStack [$#PCDataStack] .= $str;
}




=head2 VERSION and Revision

 Usage   : print $Bio::Biblio::IO::medlinexml::VERSION;
           print $Bio::Biblio::IO::medlinexml::Revision;

=cut


sub handle_start {
    my ($expat, $e, %attrs) = @_; 
#    &_debug_object_stack ("START", $e);

    #
    # The #PCDATA elements which have an attribute list must
    # be first here - because for them I create entries both on
    # the @PCDataStack _and_ on @ObjectStack.
    #
    if ($e eq 'QualifierName' or
	$e eq 'SubHeading') {
	my %p = ();
	$p{'majorTopic'} = $attrs{'MajorTopicYN'} if $attrs{'MajorTopicYN'};
	push (@ObjectStack, \%p);
    }

    if ($e eq 'GeneralNote') {
	my %p = ();
	$p{'owner'} = $attrs{'Owner'} if $attrs{'Owner'};
	push (@ObjectStack, \%p);
    }

    if ($e eq 'OtherID') {
	my %p = ();
	$p{'source'} = $attrs{'Source'};
	push (@ObjectStack, \%p);
    }

    #
    # A special treatment is for attributes for personal name.
    # Because there is no XML element 'PersonalName' I need to
    # to put yet another object on @ObjectStack unless there is
    # already one.
    #
    if ($e eq 'LastName' or
	$e eq 'FirstName' or
	$e eq 'MidleName' or
	$e eq 'Initials' or
	$e eq 'ForeName' or
	$e eq 'Suffix') {
	my $peek = $ObjectStack[$#ObjectStack];
	push (@ObjectStack, {'type' => 'PersonalName'})
	    unless (ref $peek and &_eq_hash_elem ($peek, 'type', 'PersonalName'));
    }

    #
    # Then we have #PCDATA elements without an attribute list.
    # For them I create an entry on @PCDataStack.
    #
    if (exists $PCDATA_NAMES{$e}) {
	push (@PCDataStack, '');

    #
    # And finally, all non-PCDATA elements go to the objectStack
    #
    } elsif (exists $SIMPLE_TREATMENT{$e}) {
	push (@ObjectStack, {});

    } elsif ($e eq 'PersonalNameSubject') {
	push (@ObjectStack, {'type' => 'PersonalName'});

    } elsif ($e eq 'DescriptorName' or
	     $e eq 'Descriptor') {
	if (&_eq_hash_elem (\%attrs, 'MajorTopicYN', "Y")) {
	    my $peek = $ObjectStack[$#ObjectStack];
	    $$peek{'descriptorMajorTopic'} = "Y";
	}
	    
    } elsif ($e eq 'MedlineCitation' ||
	     $e eq 'NCBIArticle') {
	my %p = ( 'type' => 'MedlineCitation' );
	$p{'owner'} = $attrs{'Owner'} if $attrs{'Owner'};
	$p{'status'} = $attrs{'Status'} if $attrs{'Status'};
	push (@ObjectStack, \%p);

    } elsif ($e eq 'GrantList') {
	if (&_eq_hash_elem (\%attrs, 'CompleteYN', "N")) {
	    my $peek = $ObjectStack[$#ObjectStack];
	    $$peek{'grantListComplete'} = "N";
	}

    } elsif ($e eq 'DataBankList') {
	if (&_eq_hash_elem (\%attrs, 'CompleteYN', "N")) {
	    my $peek = $ObjectStack[$#ObjectStack];
	    $$peek{'dataBankListComplete'} = "N";
	}

    } elsif ($e eq 'AuthorList') {
	if (&_eq_hash_elem (\%attrs, 'CompleteYN', "N")) {
	    my $peek = $ObjectStack[$#ObjectStack];
	    $$peek{'authorListComplete'} = "N";
	}

    } elsif ($e eq 'OtherAbstract') {
	my %p = ();
	$p{'type'} = $attrs{'Type'} if $attrs{'Type'};
	push (@ObjectStack, \%p);
#	push (@ObjectStack, { 'type' => 'Abstract' });
	      
    }
}

sub handle_end {
    my ($expat, $e) = @_;
    #
    # First I have to deal with those elements which are both PCDATA
    # (and therefore they are on the pcdataStack) and which have an
    # attribute list (therefore they are also known as a separate
    # p-object on the objectStack.
    #
    if ($e eq 'QualifierName' or
	$e eq 'SubHeading') {
	my $p = pop @ObjectStack;   # pSubHeading
        $$p{'subHeading'} = pop @PCDataStack;
	&_add_element ('subHeadings', $p);  # adding to pMeshHeadings
#	&_debug_object_stack ("END", $e);
	return;

    } elsif ($e eq 'GeneralNote') {
	my $p = pop @ObjectStack;  # pGeneralNote
        $$p{'generalNote'} = pop @PCDataStack;
	&_add_element ('generalNotes', $p);  # adding to pMedlineCitation
#	&_debug_object_stack ("END", $e);
	return;

    } elsif ($e eq 'OtherID') {
	my $p = pop @ObjectStack;  # pOtherID
        $$p{'otherID'} = pop @PCDataStack;
	&_add_element ('otherIDs', $p);  # adding to pMedlineCitation
#	&_debug_object_stack ("END", $e);
	return;
    }

    #
    # both object and pcdata stacks elements mixed here together
    # (the element names appear in the order of frequency in the
    # medline data set)
    #

    if (exists $POP_DATA_AND_PEEK_OBJ{$e}) {
	&_data2obj ("\l$e");

    } elsif (exists $POP_OBJ_AND_PEEK_OBJ{$e}) {
	&_obj2obj ("\l$e");

    } elsif (exists $POP_AND_ADD_ELEMENT{$e}) {
	&_add_element ($POP_AND_ADD_ELEMENT{$e}, pop @ObjectStack);

    } elsif (exists $POP_AND_ADD_DATA_ELEMENT{$e}) {
	&_add_element ($POP_AND_ADD_DATA_ELEMENT{$e});

    } elsif ($e eq 'Author' or
	     $e eq 'Investigator') {
	my $pAuthor;
	my $p = pop @ObjectStack;  # pPersonalName or pAuthor
	if (&_eq_hash_elem ($p, 'type', 'PersonalName')) {
	    $pAuthor = pop @ObjectStack;
	    $$pAuthor{'personalName'} = $p;
	} else {
	    $pAuthor = $p;
	}
	my $peek = $ObjectStack[$#ObjectStack];   # pMedlineCitation, pArticle or pBook
	if (&_eq_hash_elem ($peek, 'type', 'MedlineCitation')) {
	    &_add_element ('investigators', $pAuthor);
	} else {
	    &_add_element ('authors', $pAuthor);
	}

    } elsif ($e eq 'MedlineJournalInfo') {
	&_obj2obj ('journalInfo');

    } elsif ($e eq 'PMID') {
	my $peek = $ObjectStack[$#ObjectStack];   # pMedlineCitation, pReference or pDeleteCitation
	if (&_eq_hash_elem ($peek, 'type', 'DeleteCitation')) {
	    &_add_element ('PMIDs');
	} else {
	    $$peek{'PMID'} = pop @PCDataStack;
	}

    } elsif ($e eq 'MedlineID') {
	my $peek = $ObjectStack[$#ObjectStack];   # pMedlineCitation, pReference or pDeleteCitation
	if (&_eq_hash_elem ($peek, 'type', 'DeleteCitation')) {
	    &_add_element ('MedlineIDs');
	} else {
	    $$peek{'medlineID'} = pop @PCDataStack;
	}

#    } elsif ($e eq 'OtherAbstract') {
#	my $pAbstract = pop @ObjectStack;
#	my $pOtherAbstract = pop @ObjectStack;
#	$$pOtherAbstract{'abstract'} = $pAbstract
#	    &_add_element ('otherAbstracts', $pOtherAbstract);

    } elsif ($e eq 'Affiliation') {
	my $peek = $ObjectStack[$#ObjectStack];
	if (&_eq_hash_elem ($peek, 'type', 'PersonalName')) {
	    my $peek2 = $ObjectStack[$#ObjectStack - 1];
	    $$peek2{'affiliation'} = pop @PCDataStack;
	} else {
	    $$peek{'affiliation'} = pop @PCDataStack;
	}

    } elsif ($e eq 'DeleteCitation') {
	pop @ObjectStack;
###	warn ("'DeleteCitation' tag found. Not known what to do with it.");   # silently ignored

    } elsif ($e eq 'MedlineCitation') {

	#
	# Here we finally have the whole citation ready.
	#
	&_process_citation (pop @ObjectStack);

    #
    # ERROR: if we are here, there was an unexpected element
    #
    } elsif (exists $PCDATA_NAMES{$e}) {
	pop @PCDataStack;
	warn ("An unexpected element found: $e");
    }
#    &_debug_object_stack ("END", $e);

}

# what to do when we have the whole $citation ready
sub _process_citation {
    my ($citation) = @_;
    $citation = $Convert->convert ($citation) if defined $Convert;

    if ($Callback) {
	&$Callback ($citation);
    } else {
	push (@Citations, $citation);
    }
}

# add $element into an array named $key to the top object at @ObjectStack;
# if $element is empty, take it from @PCDataStack
sub _add_element {
    my ($key, $element) = @_;
    my $peek = $ObjectStack[$#ObjectStack];
    $$peek{$key} = [] unless $$peek{$key};
    push (@{ $$peek{$key} }, (defined $element ? $element : pop @PCDataStack));
}

# remove top of @PCDataStack and put it into top object at @ObjectStack under name $key
sub _data2obj {
    my ($key) = @_;
    my $peek = $ObjectStack[$#ObjectStack];
    $$peek{$key} = pop @PCDataStack;
}

# remove top of @ObjectStack and put it into now-top at @ObjectStack under name $key
sub _obj2obj {
    my ($key) = @_;
    my $p = pop @ObjectStack;
    my $peek = $ObjectStack[$#ObjectStack];
    $$peek{$key} = $p;
}

# check if a $key exists in a ref-hash $rh and if it is equal to $value
sub _eq_hash_elem {
    my ($rh, $key, $value) = @_;
    return (defined $$rh{$key} and $$rh{$key} eq $value);
}

#
# --- only for debugging
#
use vars qw(%DEBUGSTACK);
%DEBUGSTACK = ();
sub _debug_object_stack {
    my ($action, $element) = @_;
    if ($action =~ /^START/o) {
	$DEBUGSTACK{$element} = (@ObjectStack+0);
    } else {
	return if $element eq 'LastName';
	print "Element $element starts on " .
	    $DEBUGSTACK{$element} . 'and ends on ' . (@ObjectStack+0) . "\n"
		if $DEBUGSTACK{$element} != (@ObjectStack+0);
    }
}

1;
__END__
