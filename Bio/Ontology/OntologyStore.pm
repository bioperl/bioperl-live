#
# BioPerl module for Bio::Ontology::OntologyStore
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Ontology::OntologyStore - A repository of ontologies

=head1 SYNOPSIS

  #----------
  #SCENARIO 1
  #----------

  #make an ontology object manually. via OntologyIO
  my $io = Bio::OntologyIO->new(
                                #params to fetch Cell Ontology here
                               );
  my $cell_ontology = $io->next_ontology;

  #this is a singleton that caches the fact that you've created
  #a 'Cell Ontology' intance...
  my $store = Bio::Ontology::OntologyStore->get_instance();

  #...and it can hand you back a copy of it at any time.
  my $cell_ontology_copy = $store->get_ontology('Cell Ontology');


  #----------
  #SCENARIO 2
  #----------

  my $store = Bio::Ontology::OntologyStore->get_instance();
  #this use case allows the construction of an ontology on
  #demand just by supplying the name.
  my $ontology = $store->get_ontology('Sequence Ontology');


=head1 DESCRIPTION

The primary purpose of this module is that of a singleton repository
of L<Bio::Ontology::OntologyI> instances from which an Ontology
instance can be retrieved by name or identifier. This enables TermI
implementations to return their corresponding OntologyI through using
this singleton store instead of storing a direct reference to the
Ontology object. The latter would almost inevitably lead to memory
cycles, and would therefore potentially blow up an application.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp

  Hilmar Lapp E<lt>hlapp@gmx.netE<gt>
  Allen Day E<lt>allenday@ucla.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Ontology::OntologyStore;
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Ontology::DocumentRegistry;
use Bio::OntologyIO;
use FileHandle;
use File::Spec::Functions;


use base qw(Bio::Root::Root);

# these are the static ontology stores by name and by identifier - there is
# only one of each in any application
my %ont_store_by_name = ();
my %ont_store_by_id = ();
my %ont_aliases = (
                   'Gene Ontology' => 'Gene_Ontology'
                    );
# also, this is really meant as a singleton object, so we try to enforce it
my $instance = undef;

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Ontology::OntologyStore->new();
 Function: Returns the Bio::Ontology::OntologyStore object.

           Unlike usual implementations of new, this implementation
           will try to return a previously instantiated store, if
           there is any. It is just a synonym for get_instance. In
           order to avoid ambiguities in your code, you may rather
           want to call rather get_instance explicitly, which also
           usually is better associated with this kind of behaviour.

 Returns : an instance of Bio::Ontology::OntologyStore
 Args    :

=cut

sub new {
    return shift->get_instance(@_);
}

=head2 get_instance

 Title   : get_instance
 Usage   :
 Function: Get an instance of this class for perusal.

           Since by design this class is meant to be used as a
           singleton, the implementation will return a previously
           instantianted store if there is one, and instantiate a new
           one otherwise. In order to use this class by means of an
           instance, call this method for added code clarity, not
           new().

 Example :
 Returns : an instance of this class
 Args    : named parameters, if any (currently, there are no 
           class-specific parameters other than those accepted by
           Bio::Root::Root.

See L<Bio::Root::Root>.

=cut

sub get_instance{
   my ($self,@args) = @_;

   if(! $instance) {
       $instance = $self->SUPER::new(@args);
   }
   return $instance;
}

=head2 get_ontology

 Title   : get_ontology
 Usage   :
 Function: Get a previously instantiated and registered instance of
           this class by name or by identifier. 

           One of the main purposes of this class is to enable TermI
           implementations to return their respective ontology without
           keeping a strong reference to the respective ontology
           object. Only objects previously registered objects can be
           retrieved.

           This is a class method, hence you can call it on the class
           name, without dereferencing an object.

 Example :
 Returns : a Bio::Ontology::OntologyI implementing object, or undef
           if the query could not be satisfied
 Args    : Named parameters specifying the query. The following parameters
           are recognized:
              -name   query the store for an ontology with the given name
              -id     query for an ontology with the given identifier
           If both are specified, an implicit AND logical operator is
           assumed.

See L<Bio::Ontology::OntologyI>.

=cut

sub get_ontology{
  my ($self,@args) = @_;
  my $ont;

  my ($name,$id) = $self->_rearrange([qw(NAME ID)], @args);
  if($id) {
    $ont = $ont_store_by_id{$id};
    return unless $ont; # no AND can be satisfied in this case
  }

  if($name) {
    my $o = $ont_store_by_name{$name};

    if(!$o){
      my $doc_registry = Bio::Ontology::DocumentRegistry->get_instance();
      my($url,$def,$fmt) = $doc_registry->documents($name);

      if(ref($url) eq 'ARRAY'){
        my $io = Bio::OntologyIO->new(-url      => $url,
                                      -defs_url => $def,
                                      -format   => $fmt,
                                     );

        $o = $io->next_ontology();
        $ont_store_by_name{$name} = $o;
      } elsif($url){
        my $io = Bio::OntologyIO->new(-url      => $url,
                                      -defs_url => $def,
                                      -format   => $fmt,
                                     );
        $o = $io->next_ontology;
        $ont_store_by_name{$name} = $o;
      }
    }

    if((! $ont) || ($ont->identifier() eq $o->identifier())) {
      $ont = $o;
    } else {
      $ont = undef;
    }
  }
  
  return $ont;
}

=head2 register_ontology

 Title   : register_ontology
 Usage   :
 Function: Registers the given Ontology object for later retrieval
           by name and identifier.

 Example :
 Returns : TRUE on success and FALSE otherwise
 Args    : the Bio::Ontology::OntologyI object(s) to register

See L<Bio::Ontology::OntologyI>.

=cut

sub register_ontology {
  my ($self,@args) = @_;
  my $ret = 1;
  foreach my $ont (@args) {
    if(ref($ont) && $ont->isa('Bio::Ontology::OntologyI')){
      $ont_store_by_name{$ont->name()} = $ont if $ont->name;
      next;
    }

	if(! (ref($ont) && $ont->isa("Bio::Ontology::OntologyI"))) {
      $self->throw((ref($ont) ? ref($ont) : $ont)." does not implement ".
                   "Bio::Ontology::OntologyI or is not an object");
	}
	if($self->get_ontology(-name => $ont->name())) {
      $self->warn("ontology with name \"".$ont->name().
                  "\" already exists in the store, ignoring new one");
      $ret = 0;
      next;
	}
	if($self->get_ontology(-id => $ont->identifier())) {
      $self->warn("ontology with id \"".$ont->identifier().
                  "\" already exists in the store, ignoring new one");
      $ret = 0;
      next;
	}
	$ont_store_by_name{$ont->name()} = $ont;
	$ont_store_by_id{$ont->identifier()} = $ont;
  }
  return $ret;
}

=head2 remove_ontology

 Title   : remove_ontology
 Usage   :
 Function: Remove the specified ontology from the store.
 Example :
 Returns : TRUE on success and FALSE otherwise
 Args    : the Bio::Ontology::OntologyI implementing object(s)
           to be removed from the store

See L<Bio::Ontology::OntologyI>.

=cut

sub remove_ontology{
    my $self = shift;
    my $ret = 1;

    foreach my $ont (@_) {
	$self->throw(ref($ont)." does not implement Bio::Ontology::OntologyI")
	    unless $ont && ref($ont) && $ont->isa("Bio::Ontology::OntologyI");
	# remove it from both the id hash and the name hash
	delete $ont_store_by_id{$ont->identifier()};
	delete $ont_store_by_name{$ont->name()} if $ont->name();
    }
    return 1;
}

=head2 guess_ontology()

 Usage   : my $ontology = 
           Bio::Ontology::OntologyStore->guess_ontology('GO:0000001');
 Function: tries to guess which ontology a term identifier comes from, 
           loads it as necessary,
           and returns it as a Bio::Ontology::Ontology object.
 Example :
 Returns : a Bio::Ontology::Ontology object, or warns and returns undef
 Args    : an ontology term identifier in XXXX:DDDDDDD format.  
           Guessing is based on the XXXX string before the colon.

=cut

sub guess_ontology {
  my ($self,$id) = @_;

  my($prefix) = $id =~ /^(.+?):.+$/;

  my %prefix = (
                SO => 'Sequence Ontology',
                SOFA => 'Sequence Ontology Feature Annotation',
                GO => 'Gene Ontology',
               );

  return $prefix{$prefix} || undef;
}

1;
