#
# BioPerl module Bio::Biblio::IO::pubmed2ref.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::IO::pubmed2ref - A converter of a raw hash to PUBMED citations

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


package Bio::Biblio::IO::pubmed2ref;

use strict;

use base qw(Bio::Biblio::IO::medline2ref);


# ---------------------------------------------------------------------
#
#   Here is the core...
#
# ---------------------------------------------------------------------

sub _load_instance {
    my ($self, $source) = @_;

    my $result;
    my $article = $$source{'article'};
    if (defined $article) {
	if (defined $$article{'journal'}) {
	    $result = $self->_new_instance ('Bio::Biblio::PubmedJournalArticle');
	    $result->type ('JournalArticle');
	} elsif (defined $$article{'book'}) {
	    $result = $self->_new_instance ('Bio::Biblio::PubmedBookArticle');
	    $result->type ('BookArticle');
	} else {
	    $result->type ('PubmedArticle');
	}
    }
    $result = $self->_new_instance ('Bio::Biblio::Ref') unless defined $result;
    return $result;
}

sub convert {
    my ($self, $source) = @_;
    my $result = $self->SUPER::convert ($source->{'Citation'});	

    # here we do PUBMED's specific stuff
    my $pubmed_data = $$source{'PubmedData'};
    if (defined $pubmed_data) {

	# ... just take it (perhaps rename it)
	$result->pubmed_status ($$pubmed_data{'publicationStatus'}) if defined $$pubmed_data{'publicationStatus'};
	$result->pubmed_provider_id ($$pubmed_data{'providerId'}) if defined $$pubmed_data{'providerId'};
	$result->pubmed_article_id_list ($$pubmed_data{'pubmedArticleIds'}) if defined $$pubmed_data{'pubmedArticleIds'};
	$result->pubmed_url_list ($$pubmed_data{'pubmedURLs'}) if defined $$pubmed_data{'pubmedURLs'};

	# ... put all dates from all 'histories' into one array
	if (defined $$pubmed_data{'histories'}) {
	    my @history_list;
	    foreach my $history ( @{ $$pubmed_data{'histories'} } ) {
		my $ra_pub_dates = $$history{'pubDates'};
		foreach my $pub_date ( @{ $ra_pub_dates } ) {
		    my %history = ();
		    my $converted_date = &Bio::Biblio::IO::medline2ref::_convert_date ($pub_date);
		    $history{'date'} = $converted_date if defined $converted_date;
		    $history{'pub_status'} = $$pub_date{'pubStatus'} if defined $$pub_date{'pubStatus'};
		    push (@history_list, \%history);
		}
	    }
	    $result->pubmed_history_list (\@history_list);
	}
    }

    # Done!
    return $result;
}

1;
__END__
