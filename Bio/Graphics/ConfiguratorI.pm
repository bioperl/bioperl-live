# $Id$
#
# BioPerl module for Bio::Graphics::ConfiguratorI
#
# Cared for by Robert Hubley <rhubley@systemsbiology.org>
#
# Copyright Robert Hubley
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Graphics::ConfiguratorI - A sectioned map of configuration
options (a map of maps), with a default section.  Intended to augment
existing tag-E<gt>value semantics (ie. of Bio::AnnotationCollectionI) for
object-representation information (eg. foreground color), and for
general interface preferences (eg. image width in gbrowse).

=head1 SYNOPSIS

    # get a ConfiguratorI somehow
    my $fg_color = $configurator->get('fgcolor');

=head1 DESCRIPTION

This object contains various configuration parameters.  It is divided
up into sections and tags.  This is essentially a multi-level map
(section-E<gt>tag-E<gt>value).  There is also the concept of a default
section which is referenced when no section is passed to the
ConfiguratorI methods.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Robert Hubley

Email rhubley@systemsbiology.org

=head1 CONTRIBUTORS

Paul Edlefsen, pedlefsen@systemsbiology.org
Lincoln Stein, lstein@cshl.org
Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Graphics::ConfiguratorI;
use strict;
use Carp;

use base qw(Bio::Root::RootI);

=head2 get_sections

 Title   : get_sections
 Usage   : my @values = $configurator->get_sections();
 Function: Returns a list of the valid sections except
           the default or undef.
 Returns : A list of the sections which can be queried.
 Args    : (optional section as string, tag as string)

=cut

sub get_sections {
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 get_tags

 Title   : get_tags
 Usage   : my @values = $configurator->get_tags();
           or
           my @values = $configurator->get_tags('dna');
 Function: Returns a list of tags for a given section
           or only the default tags section if no section
           is given.
 Returns : A scalar list of tags
 Args    :

=cut

sub get_tags {
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 get

 Title   : get
 Usage   : my $value = $configurator->get('height');
           or
           my $value = $configurator->get('dna','height');
 Function: Returns a tag value from a configurator from the
           either the default "_general" section or from
           a specified section or undef.
 Returns : A scalar value for the tag
 Args    : (optional section as string, tag as string)

=cut

sub get {
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 set

 Title   : set
 Usage   : $configurator->set('fgcolor','chartreuse');
           or
           $configurator->set('EST','fgcolor','chartreuse');
 Function: Set a value for a tag
 Returns : The old value of the tag
 Args    : (optional section as string, tag as string, value as scalar)

=cut

sub set {
   my ($self) = @_;
   $self->throw_not_implemented();
}


=head2 get_and_eval

 Title   : get_and_eval
 Usage   : my $value = $configurator->get_and_eval('height');
           or
           my $value = $configurator->get_and_eval('dna','height');
 Function: This works like get() except that it is
           also able to evaluate code references.  These are
           options whose values begin with the characters
           "sub {".  In this case the value will be passed to
           an eval() and the resulting codereference returned.
 Returns : A value of the tag or undef.
 Args    : (optional section as string, tag as string)

=cut

sub get_and_eval {
   my ($self) = @_;
   $self->throw_not_implemented();
}

1;
