#
# BioPerl module for Bio::Search::Result::ResultFactory
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

Bio::Search::Result::ResultFactory - A factory to create Bio::Search::Result::ResultI objects 

=head1 SYNOPSIS

    use Bio::Search::Result::ResultFactory;
    my $factory = Bio::Search::Result::ResultFactory->new();
    my $resultobj = $factory->create(@args);

=head1 DESCRIPTION

This is a general way of hiding the object creation process so that we
can dynamically change the objects that are created by the SearchIO
parser depending on what format report we are parsing.

This object is for creating new Results.

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


package Bio::Search::Result::ResultFactory;
use vars qw($DEFAULT_TYPE);
use strict;


use base qw(Bio::Root::Root Bio::Factory::ObjectFactoryI);

BEGIN { 
    $DEFAULT_TYPE = 'Bio::Search::Result::GenericResult'; 
}

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::Result::ResultFactory->new();
 Function: Builds a new Bio::Search::Result::ResultFactory object 
 Returns : Bio::Search::Result::ResultFactory
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
 Function: Create a new L<Bio::Search::Result::ResultI> object  
 Returns : L<Bio::Search::Result::ResultI>
 Args    : hash of initialization parameters


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
 Usage   : $factory->type('Bio::Search::Result::GenericResult');
 Function: Get/Set the Result creation type
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
