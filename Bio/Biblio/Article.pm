# $Id$
#
# BioPerl module for Bio::Biblio::Article
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::Article - Representation of a general article

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


package Bio::Biblio::Article;
use strict;
use vars qw(@ISA);

use Bio::Biblio::Ref;

@ISA = qw( Bio::Biblio::Ref);

#
# a closure with a list of allowed attribute names (these names
# correspond with the allowed 'get' and 'set' methods); each name also
# keep what type the attribute should be (use 'undef' if it is a
# simple scalar)
#
{
    my %_allowed =
	(
	 _abstract => undef,
	 _abstract_type => undef,
	 _affiliation => undef,
	 _author_list_complete => undef,
	 _authors => 'ARRAY',
	 _cross_references => 'ARRAY',
	 _cross_references_list_complete => undef,
	 _date_of_electronic_publication => undef,
	 _first_page => undef,
	 _grant_list_complete => undef,
	 _grants => 'ARRAY',
	 _language => undef,
	 _last_page => undef,
	 _medline_page => undef,
	 _other_abstracts => 'ARRAY',
	 _other_languages => undef,
	 _rights => undef,
	 _title => undef,
	 _vernacular_title => undef,
	 );

    # return 1 if $attr is allowed to be set/get in this class
    sub _accessible {
	my ($self, $attr) = @_;
	exists $_allowed{$attr} or $self->SUPER::_accessible ($attr);
    }

    # return an expected type of given $attr
    sub _attr_type {
	my ($self, $attr) = @_;
	if (exists $_allowed{$attr}) {
	    return $_allowed{$attr};
	} else {
	    return $self->SUPER::_attr_type ($attr);
	}
    }
}


1;
__END__
