#
# BioPerl module for Bio::Annotation::TypeManager
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::TypeManager - Manages types for annotation collections

=head1 SYNOPSIS

    # default type manager

    $tm = Bio::Annotation::TypeManager->new();

    # $key is a string or a Bio::Ontology::TermI compliant object
    print "The type for $key is ",$tm->type_for_key($key),"\n";

    if( !$tm->is_valid($key,$object) ) {
        $self->throw("Invalid object for key $key");
    }

=head1 DESCRIPTION

Manages types for annotation collections.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.
Bug reports can be submitted via the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Annotation::TypeManager;
use strict;

# Object preamble - inherits from Bio::Root::Root



use base qw(Bio::Root::Root);
# new() can be inherited from Bio::Root::Root

=head2 new

 Title   : new
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub new{
   my ($class,@args) = @_;

   my $self = $class->SUPER::new(@args);

   $self->{'_type'} = {};

   $self->_add_type_map('reference',"Bio::Annotation::Reference");
   $self->_add_type_map('comment',"Bio::Annotation::Comment");
   $self->_add_type_map('dblink',"Bio::Annotation::DBLink");

   return $self;
}


=head2 type_for_key

 Title   : type_for_key
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub type_for_key{
   my ($self,$key) = @_;

   $key = $key->name() if ref($key) && $key->isa("Bio::Ontology::TermI");
   return $self->{'_type'}->{$key};
}


=head2 is_valid

 Title   : is_valid
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub is_valid{
   my ($self,$key,$object) = @_;

   if( !defined $object || !ref $object ) {
       $self->throw("Cannot type an object [$object]!");
   }

   if( !$object->isa($self->type_for_key($key)) ) {
       return 0;
   } else {
       return 1;
   }
}


=head2 _add_type_map

 Title   : _add_type_map
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub _add_type_map{
   my ($self,$key,$type) = @_;

   $key = $key->name() if ref($key) && $key->isa("Bio::Ontology::TermI");
   $self->{'_type'}->{$key} = $type;
}

1;
