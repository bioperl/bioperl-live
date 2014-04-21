#
# BioPerl module for Bio::Factory::ObjectFactory
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

#
# (c) Hilmar Lapp, hlapp at gmx.net, 2003.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2003.
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
# 
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#

# POD documentation - main docs before the code

=head1 NAME

Bio::Factory::ObjectFactory - Instantiates a new Bio::Root::RootI (or derived class) through a factory

=head1 SYNOPSIS

    use Bio::Factory::ObjectFactory;

    my $factory = Bio::Factory::ObjectFactory->new(-type => 'Bio::Ontology::GOterm');
    my $term = $factory->create_object(-name => 'peroxisome',
                                       -ontology => 'Gene Factory',
                                       -identifier => 'GO:0005777');


=head1 DESCRIPTION

This object will build L<Bio::Root::RootI> objects generically.

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

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net


=head1 CONTRIBUTORS

This is mostly copy-and-paste with subsequent adaptation from
Bio::Seq::SeqFactory by Jason Stajich. Most credits should in fact go
to him.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Factory::ObjectFactory;
use strict;


use base qw(Bio::Root::Root Bio::Factory::ObjectFactoryI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Factory::ObjectFactory->new();
 Function: Builds a new Bio::Factory::ObjectFactory object 
 Returns : Bio::Factory::ObjectFactory
 Args    : -type      => string, name of a L<Bio::Root::RootI> derived class.
                         There is no default.
           -interface => string, name of the interface or class any type
                         specified needs to at least implement.
                         The default is Bio::Root::RootI.

=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
  
    my ($type,$interface) = $self->_rearrange([qw(TYPE INTERFACE)], @args);

    $self->{'_loaded_types'} = {};
    $self->interface($interface || "Bio::Root::RootI");
    $self->type($type) if $type;

    return $self;
}


=head2 create_object

 Title   : create_object
 Usage   : my $seq = $factory->create_object(<named parameters>);
 Function: Instantiates a new object of the previously set type.

           This object allows us to genericize the instantiation of
           objects.

           You must have provided -type at instantiation, or have
           called type($mytype) before you can call this method.

 Returns : an object of the type returned by type()

           The return type is configurable using new(-type =>"..."),
           or by calling $self->type("My::Fancy::Class").
 Args    : Initialization parameters specific to the type of
           object we want. Check the POD of the class you set as type.

=cut

sub create_object {
   my ($self,@args) = @_;

   my $type = $self->type(); # type has already been loaded upon set
   return $type->new(-verbose => $self->verbose, @args);
}

=head2 type

 Title   : type
 Usage   : $obj->type($newval)
 Function: Get/set the type of object to be created.

           This may be changed at any time during the lifetime of this
           factory.

 Returns : value of type (a string)
 Args    : newvalue (optional, a string)


=cut

sub type{
    my $self = shift;

    if(@_) {
	my $type = shift;
	if($type && (! $self->{'_loaded_types'}->{$type})) {
	    eval {
		$self->_load_module($type);
	    };
	    if( $@ ) {
		$self->throw("module for '$type' failed to load: ".
			     $@);
	    }
	    my $o = bless {},$type;
	    if(!$self->_validate_type($o)) { # this may throw an exception
		$self->throw("'$type' is not valid for factory ".ref($self));
	    }
	    $self->{'_loaded_types'}->{$type} = 1;
	}
	return $self->{'type'} = $type;
    }
    return $self->{'type'};
}

=head2 interface

 Title   : interface
 Usage   : $obj->interface($newval)
 Function: Get/set the interface or base class that supplied types
           must at least implement (inherit from).
 Example : 
 Returns : value of interface (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub interface{
    my $self = shift;
    my $interface = shift;

    if($interface) {
	return $self->{'interface'} = $interface;
    }
    return $self->{'interface'};
}

=head2 _validate_type

 Title   : _validate_type
 Usage   : $factory->_validate_type($object)
 Function: Called to let derived factories validate the type set
           via type().

           The default implementation here checks whether the supplied
           object skeleton implements the interface set via -interface
           upon factory instantiation.

 Example :
 Returns : TRUE if the type is to be considered valid, and FALSE otherwise.
           Instead of returning FALSE this method may also just throw
           an informative exception.

           The default implementation here will throw an exception
           if the supplied object does not inherit from the interface
           provided by the interface() method.

 Args    : A hash reference blessed into the specified type, allowing
           queries like isa().


=cut

sub _validate_type{
    my ($self,$obj) = @_;

    if(! $obj->isa($self->interface())) {
	$self->throw("invalid type: '".ref($obj).
		     "' does not implement '".$self->interface()."'");
    }
    return 1;
}

#####################################################################
# aliases for naming consistency or other reasons                   #
#####################################################################

*create = \&create_object;

1;
