#
# BioPerl module for Bio::Search::Hit::HitFactory
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Hit::HitFactory - A factory to create Bio::Search::Hit::HitI objects 

=head1 SYNOPSIS

    use Bio::Search::Hit::HitFactory;
    my $factory = Bio::Search::Hit::HitFactory->new();
    my $resultobj = $factory->create(@args);

=head1 DESCRIPTION

This is a general way of hiding the object creation process so that we
can dynamically change the objects that are created by the SearchIO
parser depending on what format report we are parsing.

This object is for creating new Hits.

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

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::Hit::HitFactory;
use vars qw($DEFAULT_TYPE);
use strict;


use base qw(Bio::Root::Root Bio::Factory::ObjectFactoryI);

BEGIN { 
    $DEFAULT_TYPE = 'Bio::Search::Hit::GenericHit'; 
}

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::Hit::HitFactory->new();
 Function: Builds a new Bio::Search::Hit::HitFactory object 
 Returns : Bio::Search::Hit::HitFactory
 Args    :


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
 Function: Create a new L<Bio::Search::Hit::HitI> object  
 Returns : L<Bio::Search::Hit::HitI>
 Args    : hash of initialization parameters


=cut

sub create{
   my ($self,@args) = @_;
   my $type = $self->type;
   eval { $self->_load_module($type) };   
   if( $@ ) { $self->throw("Unable to load module $type"); }
   return $type->new(@args);
}


=head2 type

 Title   : type
 Usage   : $factory->type('Bio::Search::Hit::GenericHit');
 Function: Get/Set the Hit creation type
 Returns : string
 Args    : [optional] string to set 


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
