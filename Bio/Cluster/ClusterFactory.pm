# $Id$
#
# BioPerl module for Bio::Cluster::ClusterFactory
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

#
# (c) Hilmar Lapp, hlapp at gmx.net, 2002.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2002.
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

Bio::Cluster::ClusterFactory - Instantiates a new Bio::ClusterI (or derived class) through a factory

=head1 SYNOPSIS

    use Bio::Cluster::ClusterFactory;
    # if you don't provide a default type, the factory will try
    # some guesswork based on display_id and namespace
    my $factory = new Bio::Cluster::ClusterFactory(-type => 'Bio::Cluster::UniGene');
    my $clu = $factory->create_object(-description => 'NAT',
                                      -display_id  => 'Hs.2');


=head1 DESCRIPTION

This object will build L<Bio::ClusterI> objects generically.

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


package Bio::Cluster::ClusterFactory;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::Factory::ObjectFactoryI;

@ISA = qw(Bio::Root::Root Bio::Factory::ObjectFactoryI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Cluster::ClusterFactory();
 Function: Builds a new Bio::Cluster::ClusterFactory object 
 Returns : Bio::Cluster::ClusterFactory
 Args    : -type => string, name of a ClusterI derived class.
                    If not provided, the factory will have to guess
                    from ID and namespace, which may or may not be
                    successful.

=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
  
    my ($type) = $self->_rearrange([qw(TYPE)], @args);

    $self->{'_loaded_types'} = {};
    $self->type($type) if $type;

    return $self;
}


=head2 create_object

 Title   : create_object
 Usage   : my $seq = $factory->create_object(<named parameters>);
 Function: Instantiates new Bio::ClusterI (or one of its child classes)

           This object allows us to genericize the instantiation of
           cluster objects.

 Returns : L<Bio::ClusterI> compliant object
           The return type is configurable using new(-type =>"...").
 Args    : initialization parameters specific to the type of cluster
           object we want.  Typically 
           -display_id  => $name
           -description => description of the cluster
           -members     => arrayref, members of the cluster

=cut

sub create_object {
   my ($self,@args) = @_;

   my $type = $self->type();
   if(! $type) {
       # we need to guess this
       $type = $self->_guess_type(@args);
       if(! $type) {
	   $self->throw("No cluster type set and unable to guess.");
       }
       # load dynamically if it hasn't been loaded yet
       if(! $self->{'_loaded_types'}->{$type}) {
	   eval {
	       $self->_load_module($type);
	       $self->{'_loaded_types'}->{$type} = 1;
	   };
	   if($@) {
	       $self->throw("Bio::ClusterI implementation $type ".
			    "failed to load: ".$@);
	   }
       }
   }
   return $type->new(-verbose => $self->verbose, @args);
}

=head2 type

 Title   : type
 Usage   : $obj->type($newval)
 Function: 
 Returns : value of type
 Args    : newvalue (optional)


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
		$self->throw("Cluster implementation '$type' failed to load: ".
			     $@);
	    }
	    my $a = bless {},$type;
	    if( ! $a->isa('Bio::ClusterI') ) {
		$self->throw("'$type' does not implement Bio::ClusterI. ".
			     "Too bad.");
	    }
	    $self->{'_loaded_types'}->{$type} = 1;
	}
	return $self->{'type'} = $type;
    }
    return $self->{'type'};
}

=head2 _guess_type

 Title   : _guess_type
 Usage   :
 Function: Guesses the right type of L<Bio::ClusterI> implementation
           based on initialization parameters for the prospective
           object.
 Example :
 Returns : the type (a string, the module name)
 Args    : initialization parameters to be passed to the prospective
           cluster object


=cut

sub _guess_type{
    my ($self,@args) = @_;
    my $type;

    # we can only guess from a certain number of arguments
    my ($dispid, $ns, $members) =
	$self->_rearrange([qw(DISPLAY_ID
			      NAMESPACE
			      MEMBERS
			      )], @args);
    # Unigene namespace or ID?
    if($ns && (lc($ns) eq "unigene")) {
	$type = 'Bio::Cluster::UniGene';
    } elsif($dispid && ($dispid =~ /^Hs\.[0-9]/)) {
	$type = 'Bio::Cluster::UniGene';
    }
    # what else could we look for?
    return $type;
}

#####################################################################
# aliases for naming consistency or other reasons                   #
#####################################################################

*create = \&create_object;

1;
