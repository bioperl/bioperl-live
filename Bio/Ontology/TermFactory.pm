# $Id$
#
# BioPerl module for Bio::Ontology::TermFactory
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

Bio::Ontology::TermFactory - Instantiates a new Bio::Ontology::TermI (or derived class) through a factory

=head1 SYNOPSIS

    use Bio::Ontology::TermFactory;
    # the default type is Bio::Ontology::Term
    my $factory = Bio::Ontology::TermFactory->new(-type => 'Bio::Ontology::GOterm');
    my $clu = $factory->create_object(-name => 'peroxisome',
                                      -category => 'cellular component',
                                      -identifier => 'GO:0005777');


=head1 DESCRIPTION

This object will build L<Bio::Ontology::TermI> objects generically.

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


package Bio::Ontology::TermFactory;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::Factory::ObjectFactoryI;

@ISA = qw(Bio::Root::Root Bio::Factory::ObjectFactoryI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Ontology::TermFactory();
 Function: Builds a new Bio::Ontology::TermFactory object 
 Returns : Bio::Ontology::TermFactory
 Args    : -type => string, name of a L<Bio::Ontology::TermI> derived class.
                    The default is L<Bio::Ontology::Term>.

=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
  
    my ($type) = $self->_rearrange([qw(TYPE)], @args);

    $self->{'_loaded_types'} = {};
    $type = "Bio::Ontology::Term" unless $type;
    $self->type($type) if $type;

    return $self;
}


=head2 create_object

 Title   : create_object
 Usage   : my $seq = $factory->create_object(<named parameters>);
 Function: Instantiates new Bio::Ontology::TermI (or one of its child classes)

           This object allows us to genericize the instantiation of
           cluster objects.

 Returns : L<Bio::Ontology::TermI> compliant object
           The return type is configurable using new(-type =>"...").
 Args    : initialization parameters specific to the type of term
           object we want.  Typically 
           -name  => $name
           -identifier => identifier for the term
           -category     => category for the term

=cut

sub create_object {
   my ($self,@args) = @_;

   my $type = $self->type(); # type has already been loaded upon set
   return $type->new(-verbose => $self->verbose, @args);
}

=head2 type

 Title   : type
 Usage   : $obj->type($newval)
 Function: Get/set the type of L<Bio::AnnotationI> object to be created.

           This may be changed at any time during the lifetime of this
           factory.

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
		$self->throw("Term implementation '$type' failed to load: ".
			     $@);
	    }
	    my $a = bless {},$type;
	    if( ! $a->isa('Bio::Ontology::TermI') ) {
		$self->throw("'$type' does not implement Bio::Ontology::TermI. ".
			     "Too bad.");
	    }
	    $self->{'_loaded_types'}->{$type} = 1;
	}
	return $self->{'type'} = $type;
    }
    return $self->{'type'};
}

#####################################################################
# aliases for naming consistency or other reasons                   #
#####################################################################

*create = \&create_object;

1;
