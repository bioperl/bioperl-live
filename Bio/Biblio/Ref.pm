# $Id$
#
# BioPerl module for Bio::Biblio::Ref
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::Ref - Representation of a bibliographic reference

=head1 SYNOPSIS

  # to be written

=head1 DESCRIPTION

Super class and interface class for bibliographic references. The
central class of the Bio::Biblio name space.

The class names and attributes comes from the Martin Senger's java implementation
of Biblio objects - the I<OpenBQS> project pages are at
http://industry.ebi.ac.uk/openBQS/.

See Martin's the UML diagram at 
http://industry.ebi.ac.uk/openBQS/images/bibobjects_java.jpg
for an overview.


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

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are preceded with a _

=cut


# Let the code begin...


package Bio::Biblio::Ref;
use strict;
use vars qw(@ISA $AUTOLOAD);

use Bio::Biblio::BiblioBase;

@ISA = qw(Bio::Biblio::BiblioBase);

#
# a closure with a list of allowed attribute names (these names
# correspond with the allowed 'get' and 'set' methods); each name also
# keep what type the attribute should be (use 'undef' if it is a
# simple scalar)
#
{
    my %_allowed =
	(
	 _author_list_complete => undef,
	 _authors => 'ARRAY',
	 _cross_references => 'ARRAY',
	 _cross_references_list_complete => undef,
	 _abstract => undef,
	 _abstract_language => undef,
	 _abstract_type => undef,
	 _codes => 'HASH',
	 _contributors => 'ARRAY',
	 _date => undef,
	 _date_completed => undef,
	 _date_created => undef,
	 _date_revised => undef,
	 _format => undef,
	 _identifier => undef,
	 _keywords => 'HASH',
	 _language => undef,
	 _last_modified_date => undef,
	 _publisher => 'Bio::Biblio::Provider',
	 _repository_subset => undef,
	 _rights => undef,
	 _spatial_location => undef,
	 _subject_headings => 'HASH',
	 _subject_headings_source => undef,
	 _temporal_period => undef,
	 _title => undef,
	 _toc => undef,
	 _toc_type => undef,
	 _type => undef,
	 );

    # return 1 if $attr is allowed to be set/get in this class
    sub _accessible {
	my ($self, $attr) = @_;
	exists $_allowed{$attr};
    }

    # return an expected type of given $attr
    sub _attr_type {
	my ($self, $attr) = @_;
	$_allowed{$attr};
    }
}
1;
__END__
