# $Id$
#
# BioPerl module for Bio::Biblio::MedlineArticle
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::edlineArticle - Representation of a general MEDLINE article

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


package Bio::Biblio::MedlineArticle;
use strict;
use vars qw(@ISA);

use Bio::Biblio::Article;

@ISA = qw(Bio::Biblio::Article);

#
# a closure with a list of allowed attribute names (these names
# correspond with the allowed 'get' and 'set' methods); each name also
# keep what type the attribute should be (use 'undef' if it is a
# simple scalar)
#
{
    my %_allowed =
	(
	 _affiliation => undef,
	 _chemicals => 'ARRAY',
	 _citation_owner => undef,
	 _comment_ins => 'ARRAY',
	 _comment_ons => 'ARRAY',
	 _date_of_electronic_publication => undef,
	 _erratum_fors => 'ARRAY',
	 _erratum_ins => 'ARRAY',
	 _general_notes => 'ARRAY',
	 _grant_list_complete => undef,
	 _grants => 'ARRAY',
	 _medline_date => undef,
	 _medline_id => undef,
	 _medline_page => undef,
	 _mesh_headings => 'ARRAY',
	 _number_of_references => undef,
	 _original_report_ins => 'ARRAY',
	 _other_abstracts => 'ARRAY',
	 _other_ids => 'ARRAY',
	 _other_languages => undef,
	 _pmid => undef,
	 _republished_froms => 'ARRAY',
	 _republished_ins => 'ARRAY',
	 _retraction_ins => 'ARRAY',
	 _retraction_ofs => 'ARRAY',
	 _season => undef,
	 _status => undef,
	 _summary_for_patients_ins => 'ARRAY',
	 _update_ins => 'ARRAY',
	 _update_ofs => 'ARRAY',
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
