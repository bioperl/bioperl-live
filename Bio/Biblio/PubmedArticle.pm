# $Id$
#
# BioPerl module for Bio::Biblio::PubmedArticle
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::PubmedArticle - Representation of a PUBMED article

=head1 SYNOPSIS

    $obj = new Bio::Biblio::PubmedArticle
                  (-pubmed_history_list =>
                       [ { 'pub_status' => 'pubmed',
			   'date' => '2001-12-1T10:0:00Z' },
			 { 'pub_status' => 'medline',
			   'date' => '2002-1-5T10:1:00Z' } ],
                   -pubmed_status => 'ppublish');
 --- OR ---

    $obj = new Bio::Biblio::PubmedArticle;
    $obj->pubmed_status ('ppublish');

=head1 DESCRIPTION

A storage object for a general PUBMED article.
See its place in the class hierarchy in
http://industry.ebi.ac.uk/openBQS/images/bibobjects_perl.gif

=head2 Attributes

The following attributes are specific to this class
(however, you can also set and get all attributes defined in the parent classes):

  pubmed_status
  pubmed_provider_id
  pubmed_history_list       type: array ref of hashes
  pubmed_article_id_list    type: array ref of hashes
  pubmed_url_list           type: array ref of hashes

=head1 SEE ALSO

=over

=item *

OpenBQS home page: http://industry.ebi.ac.uk/openBQS

=item *

Comments to the Perl client: http://industry.ebi.ac.uk/openBQS/Client_perl.html

=back

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
  http://bugzilla.bioperl.org/

=head1 AUTHOR

Martin Senger (senger@ebi.ac.uk)

=head1 COPYRIGHT

Copyright (c) 2002 European Bioinformatics Institute. All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut


# Let the code begin...


package Bio::Biblio::PubmedArticle;
use strict;
use vars qw(@ISA);

use Bio::Biblio::MedlineArticle;
@ISA = qw(Bio::Biblio::MedlineArticle);

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
	 _linkout => 'Bio::Biblio::PubmedLinkout', # NG 03-06-18
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
