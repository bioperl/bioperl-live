# $Id$
#
# BioPerl module for Bio::Tree::NodeI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::NodeI - DESCRIPTION of Interface

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the interface here

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

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tree::NodeI;
use vars qw(@ISA);
use strict;
use Carp;

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  confess "Abstract method '$caller' defined in interfaceBio::Tree::NodeI not implemented by pacakge $package. Not your fault - author of $package should be blamed!";
}

=head2 get_parent

 Title   : get_parent
 Usage   : my $node = $node->parent;
 Function: Gets a Node\'s parent node
 Returns : Null if this is top level node
 Args    : none

=cut

sub get_parent{
   my ($self,@args) = @_;
   $self->_abstractDeath();
}

=head2 is_leaf

 Title   : is_leaf
 Usage   : if( $node->is_leaf ) 
 Function: Get/Set Leaf status
 Returns : boolean
 Args    : (optional) boolean

=cut

sub is_leaf{
   my ($self) = @_;
   $self->_abstractDeath();
}


1;
