# $Id$
#
# BioPerl module for Bio::Biblio::PubmedJournalArticle
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::PubmedJournalArticle - Representation of a PUBMED journal article

=head1 SYNOPSIS

#

=head1 DESCRIPTION

#

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHORS

Heikki Lehvaslaiho (heikki@ebi.ac.uk)
Martin Senger (senger@ebi.ac.uk)

=head1 COPYRIGHT

Copyright (c) 2002 European Bioinformatics Institute. All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut


# Let the code begin...


package Bio::Biblio::PubmedJournalArticle;
use strict;
use vars qw(@ISA);

use Bio::Biblio::MedlineJournalArticle;

@ISA = qw(Bio::Biblio::MedlineJournalArticle);

#
# a closure with a list of allowed attribute names (these names
# correspond with the allowed 'get' and 'set' methods); each name also
# keep what type the attribute should be (use 'undef' if it is a
# simple scalar)
#
{
    my %_allowed =
	(
	 _pubmed_status => undef,
	 _pubmed_provider_id => undef,
	 _pubmed_history_list => 'ARRAY',
	 _pubmed_article_id_list => 'ARRAY',
	 _pubmed_url_list => 'ARRAY',
	 );

    # return 1 if $attr is allowed to be set/get in this class
    sub _accessible {
	my ($self, $attr) = @_;
	return 1 if exists $_allowed{$attr};
        foreach my $parent (@ISA) {
	    return 1 if $parent->_accessible ($attr);
	}
    }

    # return an expected type of given $attr
    sub _attr_type {
	my ($self, $attr) = @_;
	if (exists $_allowed{$attr}) {
	    return $_allowed{$attr};
	} else {
	    foreach my $parent (@ISA) {
		if ($parent->_accessible ($attr)) {
		    return $parent->_attr_type ($attr);
		}
	    }
	}
	return 'unknown';
    }
}


1;
__END__
