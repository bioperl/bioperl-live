#------------------------------------------------------------------
# $Id$
#
# BioPerl module for Bio::Factory::GenericObjectFactory
#
# Cared for by Steve Chervitz <sac@bioperl.org> and Jason Stajich <jason@bioperl.org>
#
# Copyright Steve Chervitz, Jason Stajich
#
# You may distribute this module under the same terms as perl itself
#------------------------------------------------------------------

# POD documentation - main docs before the code

=head1 NAME

Bio::Factory::GenericObjectFactory - A generic implementation of ObjectFactoryI for creating any type of object

=head1 SYNOPSIS

    use Bio::Factory::GenericObjectFactory;
    my $factory = Bio::Factory::GenericObjectFactory->new(-type=>'Bio::Root::Root');
    my $obj = $factory->create(@args);

=head1 DESCRIPTION

This is a generic implementation of Bio::Factory::ObjectFactoryI.
It can be used for creating objects of any type.

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

=head1 AUTHOR - Steve Chervitz, Jason Stajich

Based on Bioperl factory modules created by Jason Stajich.

Email sac@bioperl.org
Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Factory::GenericObjectFactory;
use vars qw(@ISA $DEFAULT_TYPE);
use strict;

use Bio::Root::Root;
use Bio::Factory::ObjectFactoryI;

@ISA = qw(Bio::Root::Root Bio::Factory::ObjectFactoryI );

BEGIN { 
    $DEFAULT_TYPE = 'Bio::Root::Root';
}

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Factory::GenericObjectFactory(%args);
 Function: Builds a new Bio::Factory::GenericObjectFactory object 
 Returns : Bio::Factory::GenericObjectFactory
 Args    : -type => $module_name  (optional)
           If -type is not supplied, a default type is specified by
           $Bio::Factory::GenericObjectFactory::DEFAULT_TYPE

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($type) = $self->_rearrange([qw(TYPE)],@args);
  $self->type($type) if defined $type;
  return $self;
}

=head2 create

 Title   : create
 Usage   : $factory->create(%args)
 Function: Create a new object 
 Returns : Object of type defined by type()
 Args    : hash of initialization parameters

See Also: L<type()>

=cut

sub create{
   my ($self,@args) = @_;
   my $type = $self->type;
   eval { $self->_load_module($type) };
   if( $@ ) { $self->throw("Unable to load module $type: $@"); }
   return $type->new(@args);
}


=head2 type

 Title   : type
 Usage   : $factory->type( $module_name );
 Function: Get/Set the Result creation type
 Returns : string
 Args    : [optional] string to set the type

=cut

sub type{
    my ($self,$type) = @_;
   if( defined $type ) { 
       # redundancy with the create method which also calls _load_module
       # I know - but this is not a highly called object so I am going 
       # to leave it in
       eval {$self->_load_module($type) };
       if( $@ ){ $self->warn("Cannot find module $type, unable to set type"); }
       else { $self->{'_type'} = $type; }
   } 
    return $self->{'_type'} || $DEFAULT_TYPE;
}

1;
