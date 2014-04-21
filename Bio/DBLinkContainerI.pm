#
# BioPerl module for Bio::DBLinkContainerI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DBLinkContainerI - Abstract interface for any object wanting to use  
                        database cross references

=head1 SYNOPSIS

    # get an objects containing database cross reference

        foreach $obj ( @objs ) {
                if( $obj->isa('Bio::DBLinkContainerI') ) {
                        foreach $dblink ( $obj->each_DBLink() ) {
                                # do stuff
                        }
                }
        }

=head1 DESCRIPTION

This interface defines the functions one can expect for any object
wanting to use database cross-references. This class does not actually
provide any implementation, it just provides the definitions of what
methods one can call.

The database cross-references are implemented as L<Bio::Annotation::DBLink>
objects.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DBLinkContainerI;
use strict;

use Carp;

use base qw(Bio::Root::RootI);

=head2 each_DBLink

 Title   : each_DBLink
 Usage   : foreach $ref ( $self->each_DBlink() )
 Function: gets an array of DBlink of objects
 Example :
 Returns : an array of Bio::Annotation::DBLink objects
 Args    : none


=cut

sub each_DBLink{
   my ($self) = @_;
   my $class = ref($self) || $self;
   $self->throw("Class $class did not define method 'each_DBLink' for interface DBLinkContainerI");
}

1;

