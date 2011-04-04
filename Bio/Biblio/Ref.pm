#
# BioPerl module for Bio::Biblio::Ref
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::Ref - Representation of a bibliographic reference

=head1 SYNOPSIS

    $obj = Bio::Biblio::Ref->new(-type  => 'Letter',
                                 -title => 'Onegin to Tatiana');
  #--- OR ---

    $obj = Bio::Biblio::Ref->new();
    $obj->type ('Letter');

=head1 DESCRIPTION

A storage object for a general bibliographic reference (a citation).
See its place in the class hierarchy in
http://www.ebi.ac.uk/~senger/openbqs/images/bibobjects_perl.gif

=head2 Attributes

The following attributes are specific to this class,
and they are inherited by all citation types.

  author_list_complete            values: 'Y'  (default) or 'N'
  authors                         type:   array ref of Bio::Biblio::Provider's
  cross_references                type:   array ref of Bio::Annotation::DBLink's
  cross_references_list_complete  values: 'Y' (default) or 'N'
  abstract
  abstract_language
  abstract_type
  codes                           type:   hash ref
  contributors                    type:   array ref of Bio::Biblio::Provider's
  date
  date_completed
  date_created
  date_revised
  format
  identifier
  keywords
  language
  last_modified_date
  publisher                       type:   Bio::Biblio::Provider
  repository_subset
  rights
  spatial_location
  subject_headings                type:   hash ref
  subject_headings_source
  temporal_period
  title
  toc
  toc_type
  type

=head1 SEE ALSO

=over 4

=item *

OpenBQS home page: http://www.ebi.ac.uk/~senger/openbqs/

=item *

Comments to the Perl client: http://www.ebi.ac.uk/~senger/openbqs/Client_perl.html

=back

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

=head1 AUTHORS

Heikki Lehvaslaiho (heikki-at-bioperl-dot-org),
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
use vars qw($AUTOLOAD);

use Bio::Annotation::DBLink;

use base qw(Bio::Biblio::BiblioBase);

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
	 _authors => 'ARRAY',  # of Bio::Biblio::Provider
	 _cross_references => 'ARRAY',   # of Bio::Annotation::DBLink
	 _cross_references_list_complete => undef,
	 _abstract => undef,
	 _abstract_language => undef,
	 _abstract_type => undef,
	 _codes => 'HASH',
	 _contributors => 'ARRAY',  # of Bio::Biblio::Provider
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


=head2 add_cross_reference

 Usage   : $self->add_cross_reference
               (Bio::Annotation::DBLink->new(-database   => 'EMBL',
					     -primary_id => 'V00808');
 Function: adding a link to a database entry
 Returns : new value of 'cross_references'
 Args    : an object of type Bio::Annotation::DBLink

=cut

sub add_cross_reference {
    my ($self, $value) = @_;
    $self->throw ($self->_wrong_type_msg (ref $value, 'Bio::Annotation::DBLink'))
	unless (UNIVERSAL::isa ($value, 'Bio::Annotation::DBLink'));
    (defined $self->cross_references) ?
	push (@{ $self->cross_references }, $value) :
	    return $self->cross_references ( [$value] );
    return $self->cross_references;
}


=head2 add_author

 Usage   : $self->add_author (Bio::Biblio::Person->new(-lastname => 'Novak');
 Function: adding an author to a list of authors
 Returns : new value of 'authors' (a full list)
 Args    : an object of type Bio::Biblio::Provider

=cut


sub add_author {
    my ($self, $value) = @_;
    $self->throw ($self->_wrong_type_msg (ref $value, 'Bio::Biblio::Provider'))
	unless (UNIVERSAL::isa ($value, 'Bio::Biblio::Provider'));
    (defined $self->authors) ?
	push (@{ $self->authors }, $value) :
	    return $self->authors ( [$value] );
    return $self->authors;
}

=head2 add_contributor

 Usage   : $self->add_contributor (Bio::Biblio::Person->new(-lastname => 'Novak');
 Function: adding a contributor to a list of contributors
 Returns : new value of 'contributors' (a full list)
 Args    : an object of type Bio::Biblio::Provider

=cut

sub add_contributor {
    my ($self, $value) = @_;
    $self->throw ($self->_wrong_type_msg (ref $value, 'Bio::Biblio::Provider'))
	unless (UNIVERSAL::isa ($value, 'Bio::Biblio::Provider'));
    (defined $self->contributors) ?
	push (@{ $self->contributors }, $value) :
	    return $self->contributors ( [$value] );
    return $self->contributors;
}


1;
__END__
