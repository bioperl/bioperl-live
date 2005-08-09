# $Id$
#
# BioPerl module for Bio::Taxonomy::Node
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Taxonomy::Node - A node in a represented taxonomy

=head1 SYNOPSIS

  use Bio::Taxonomy::Node;
  # typically you will get a Node from a Bio::DB::Taxonomy object
  # but here is how you initialize one
  my $node = new Bio::Taxonomy::Node(-name      => $name,
                                     -object_id => $oid,
                                     -parent_id => $pid,
                                     -rank   => $rank,
                                     -division  => $div,
                                     -dbh       => $dbh);

  my $dbh = new Bio::DB::Taxonomy(-source   => 'flatfile',
                                  -directory=> '/tmp',
                                  -nodesfile=> '/path/to/nodes.dmp',
                                  -namesfile=> '/path/to/names.dmp');
  my $hum_node = $dbh->get_Taxonomy_Node(-name => 'Homo sapiens');
  my $hum_node2= $dbh->get_Taxonomy_Node(-taxonid => '9606');

  print "rank is ", $hum_node->rank, "\n";
  print "classification is ", join(" ", $hum_node->classification),"\n"; 
  print "division is ", $node->division, "\n";

=head1 DESCRIPTION

This is the next generation (for Bioperl) of representing Taxonomy
information.  Previously all information was managed by a single
object called Bio::Species.  This new implementation allows
representation of the intermediete nodes not just the species nodes
and can relate their connections.

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
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 CONTRIBUTORS

Juguang Xiao, juguang@tll.org.sg

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::Taxonomy::Node;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::IdentifiableI;
use Bio::DB::Taxonomy;

@ISA = qw(Bio::Root::Root Bio::IdentifiableI  );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Taxonomy::Node();
 Function: Builds a new Bio::Taxonomy::Node object 
 Returns : an instance of Bio::Taxonomy::Node
 Args    : -dbh       => a reference to a Bio::DB::Taxonomy object
           -name      => a string representing the node name
           -object_id => unique identifier - typically NCBI Taxid
           -parent_id => parent id (unique identifier for parent)
           -rank      => node rank (one of 'species', 'genus', etc)
           -division  => 'primates', 'rodents', or 'inv', etc
           -factory   => Bio::Taxonomy::FactoryI object for creating new nodes
           -classification => Arrayref of full lineage classification
           -ncbi_taxid  => alias for -object_id, only use one of these
           -common_name => genbank common name for an organism
           -organelle   => organelle type where appropriate
           -sub_species => if this is a subspecies, provide that as well,
                           most relavent for virii and others
           -variant     => provide variant info
           -genetic_code=> genetic code table number
           -mito_genetic_code => mitochondrial genetic code table number
           -create_date       => date created (where available)
           -update_date       => date last updated
           -pub_date          => date published 

=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($name,$uniqueid,$parentid,
	$rank,$div,$dbh, 
	$factory,
	$classification, $ncbitaxid, $commonname,
	$organelle, $sub_species,$variant,$gcode,$mitocode,
	$createdate,$updatedate,$pubdate) = 
	    $self->_rearrange([ qw ( NAME OBJECT_ID PARENT_ID RANK DIVISION
				     DBH 
				     FACTORY
				     CLASSIFICATION
				     NCBI_TAXID
				     COMMON_NAME
				     ORGANELLE
				     SUB_SPECIES
				     VARIANT
				     GENETIC_CODE
				     MITO_GENETIC_CODE
				     CREATE_DATE
				     UPDATE_DATE
				     PUB_DATE
				     )],
			      @args);
    if( defined $ncbitaxid && defined $uniqueid && $ncbitaxid ne $uniqueid  ) {
	$self->warn("Only provide one of -object_id or -ncbi_taxid, using $uniqueid\n");
    } elsif( ! defined $uniqueid && defined $ncbitaxid ) { 
	$uniqueid = $ncbitaxid;
    }
    $uniqueid && $self->object_id($uniqueid);

    if( defined $name && defined $commonname && $commonname ne $name ) { 
	$self->warn("Only provide one of -name or -common_name to $class, using -name value\n");
    } elsif( ! defined $name && defined $commonname ) { 
	$name = $commonname;
    }

    defined $name && $self->node_name($name);
    defined $rank && $self->rank($rank);
    defined $div  && $self->division($div);

    defined $gcode    && $self->genetic_code($gcode);
    defined $mitocode && $self->mitochondrial_genetic_code($mitocode);
    defined $createdate    && $self->create_date($createdate);
    defined $updatedate    && $self->update_date($updatedate);
    defined $pubdate       && $self->pub_date($pubdate);

#unless(defined $factory ){
    $self->db_handle($dbh || Bio::DB::Taxonomy->new(-source => 'entrez'));
#  } else { 
#      $self->factory($factory);
#  }

    defined $parentid && $self->parent_id($parentid);
    if( defined $classification ) {
	if( ref($classification) !~ /ARRAY/ ) {
	    $self->warn("Classification can only be initialize with an arrayref\n");
	}
	$self->classification(@$classification);
    }
    defined $organelle && $self->name('organelle', $organelle);
    defined $sub_species && $self->name('sub_species', $sub_species);

    defined $variant && $self->name('variant',$variant);


    return $self;
}

=head1 Bio::IdentifiableI interface 

Also see L<Bio::IdentifiableI>

=head2 version

 Title   : version
 Usage   : $obj->version($newval)
 Function: 
 Example : 
 Returns : value of version (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub version{
    my $self = shift;

    return $self->{'version'} = shift if @_;
    return $self->{'version'};
}

=head2 authority

 Title   : authority
 Usage   : $obj->authority($newval)
 Function: 
 Example : 
 Returns : value of authority (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub authority{
    my $self = shift;

    return $self->{'authority'} = shift if @_;
    return $self->{'authority'};
}


=head2 namespace

 Title   : namespace
 Usage   : $obj->namespace($newval)
 Function: 
 Example : 
 Returns : value of namespace (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub namespace{
    my $self = shift;

    return $self->{'namespace'} = shift if @_;
    return $self->{'namespace'};
}

=head1 Bio::Taxonomy::Node implementation

=head2 db_handle

 Title   : db_handle
 Usage   : $obj->db_handle($newval)
 Function: Get/Set Bio::DB::Taxonomy Handle
 Returns : value of db_handle (a scalar) (Bio::DB::Taxonomy object)
 Args    : on set, new value (a scalar or undef, optional) Bio::DB::Taxonomy object

Also see L<Bio::DB::Taxonomy>

=cut

sub db_handle {
    my $self = shift;
    if( @_ ) {
	my $v = shift;
	# until we establish some other higher level TaxonomyDB interface
	if( ! ref($v) || ! $v->isa('Bio::DB::Taxonomy') ) {
	    $self->throw("Must have provided a valid Bio::DB::Taxonomy object");
	}
	$self->{'db_handle'} = $v;
    }
    return $self->{'db_handle'};
}

=head2 factory

  Title:    factory
  Usage:    $factory->factory($newval);
  Function: Get/Set Bio::Taxonomy::FactoryI implementation
  Returns:  Bio:;Taxonomy::FactoryI
  Args:     Bio::Taxonomy::FactoryI

Also see L<Bio::Taxonomy::FactoryI>

=cut

sub factory {
    my $self = shift;
    if(@_){
        my $v = shift;
        unless(ref($v) || $v->isa('Bio::Taxonomy::FactoryI')){
            $self->throw('A Bio::Taxonomy::FactoryI object required');
        }
        $self->{_factory} = $v;
    }
    return $self->{_factory};
}

=head2 rank

 Title   : rank
 Usage   : $obj->rank($newval)
 Function: Get/set rank of this Node, 'species', 'genus', 'order', etc...
 Returns : value of rank (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub rank{
    my $self = shift;

    return $self->{'rank'} = shift if @_;
    return $self->{'rank'};
}

=head2 object_id

 Title   : object_id
 Usage   : $obj->object_id($newval)
 Function: Get/Set object id (NCBI Taxonomy ID in most cases)
 Returns : value of object_id (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub object_id {
    my $self = shift;

    return $self->{'_object_id'} = shift if @_;
    return $self->{'_object_id'};
}

*ncbi_taxid = \&object_id;

=head2 parent_id

 Title   : parent_id
 Usage   : $obj->parent_id($newval)
 Function: Get/Set parent ID
 Returns : value of parent_id (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub parent_id{
    my $self = shift;

    return $self->{'parent_id'} = shift if @_;
    return $self->{'parent_id'};
}

=head2 genetic_code

 Title   : genetic_code
 Usage   : $obj->genetic_code($newval)
 Function: Get/set genetic code table
 Returns : value of genetic_code (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub genetic_code{
    my $self = shift;

    return $self->{'genetic_code'} = shift if @_;
    return $self->{'genetic_code'};
}

=head2 mitochondrial_genetic_code

 Title   : mitochondrial_genetic_code
 Usage   : $obj->mitochondrial_genetic_code($newval)
 Function: Get/set mitochondrial genetic code table
 Returns : value of mitochondrial_genetic_code (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub mitochondrial_genetic_code{
    my $self = shift;

    return $self->{'mitochondrial_genetic_code'} = shift if @_;
    return $self->{'mitochondrial_genetic_code'};
}


=head2 create_date

 Title   : create_date
 Usage   : $obj->create_date($newval)
 Function: Get/Set Date this node was created (in the database)
 Returns : value of create_date (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub create_date{
    my $self = shift;

    return $self->{'create_date'} = shift if @_;
    return $self->{'create_date'};
}


=head2 update_date

 Title   : update_date
 Usage   : $obj->update_date($newval)
 Function: Get/Set Date this node was updated (in the database)
 Returns : value of update_date (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub update_date{
    my $self = shift;

    return $self->{'update_date'} = shift if @_;
    return $self->{'update_date'};
}


=head2 pub_date

 Title   : pub_date
 Usage   : $obj->pub_date($newval)
 Function: Get/Set Date this node was pubd (in the database)
 Returns : value of pub_date (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub pub_date{
    my $self = shift;

    return $self->{'pub_date'} = shift if @_;
    return $self->{'pub_date'};
}

=head2 get_Parent_Node

 Title   : get_Parent_Node
 Usage   : my $parentnode = $node->get_Parent_Node()
 Function: Retrieve the full Parent node from the database
 Returns : Bio::Taxonomy::Node
 Args    : none


=cut

sub get_Parent_Node {
   my ($self) = @_;

   if( ! $self->db_handle ||
       ! defined $self->parent_id ) {
       $self->warn("Cannot get the parent node for ".$self->node_name.
		   " because parent_id or db handle is not defined\n");
       return undef;
   }


   my $node = $self->db_handle->get_Taxonomy_Node(-taxonid => $self->parent_id);
   unless ( defined $node ) {
       $self->warn("Could not find node for parent id ". $self->parent_id);
       return undef;
   }
   return $node;
}


=head2 get_Children_Nodes

 Title   : get_Children_Nodes
 Usage   : my @nodes = $node->get_Children_Nodes();
 Function: Get the children of a node as L<Bio::Taxonomy::Node> objects
 Returns : Array of L<Bio::Taxonomy::Node> objects
 Args    : none


=cut

sub get_Children_Nodes{
   my ($self) = @_;
   if( ! $self->db_handle ||
       ! defined $self->parent_id ) {
       $self->warn("Cannot get the children nodes for ".$self->name('common').
		   " because the db handle is not defined\n");
       return undef;
   }
   my @nodes = $self->db_handle->get_Children_Taxids($self->object_id);
   my @children;
   for my $n ( @nodes ) {
       my $node = $self->db_handle->get_Taxonomy_Node(-taxonid => $n);
       unless( defined $node ) {
	   $self->warn("Could not find node for id $n child  of ".
		       $self->object_id);
       } else {
	   push @children, $node;
       }
   }
   return @children;
}

=head2 node_name

 Title   : node_name
 Usage   : $obj->node_name($newval)
 Function: 
 Example : 
 Returns : value of node_name (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub node_name{
    my $self = shift;
    my @v = @{$self->name('common',@_) || []};
    return ( wantarray ) ? @v : shift @v;
}

*common_name = \&node_name;

=head2 classification

 Title   : classification
 Usage   : $self->classification(@class_array);
           @classification = $self->classification();
 Function: Fills Returns the classification list in
           the object.  The array provided must be in
           the order SPECIES, GENUS ---> KINGDOM.
           Checks are made that species is in lower case,
           and all other elements are in title case.
 Returns : Classification array
 Args    : [optional] array of vals to set classification to

=cut

sub classification {
   my ($self,@vals) = @_;
   my $p;
   if( @vals ) {
       $self->{'_classification'} = [@vals];
       return @vals;
   } elsif ( defined $self->{'_classification'} ) {
       return @{$self->{'_classification'}};
   } else { 
       if( defined($p = $self->get_Parent_Node()) &&
	   $p->object_id != 1  ) {
	   # okay this won't really work - need to do proper recursion
	   push @vals, $p->classification;
       }
       if( $self->show_all || $self->rank ne 'no rank') {
	   push @vals,$self->node_name();
       }
       $self->{'_classification'} = \@vals;
   }
   return @vals;
}


=head2 division

 Title   : division
 Usage   : $obj->division($newval)
 Function: 
 Example : 
 Returns : value of division (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub division {
    my $self = shift;
    my @v = @{$self->name('division',@_) || []};
    return wantarray ? @v : shift @v;
}

=head2 show_all

 Title   : show_all
 Usage   : $obj->show_all($newval)
 Function: Boolean flag whether or not we should show all intermediete
           nodes that do not have actual ranks.
 Returns : value of show_all (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub show_all{
    my $self = shift;
    return $self->{'show_all'} = shift if @_;
    return $self->{'show_all'};
}

=head2 name

  Title:    name
  Usage:    $obj->name('scientific', 'sapiens');
            $obj->name('common', 'human', 'man');
            my @names = @{$obj->name('common')};
  Function: Get and set the names
  Returns:  names (a array reference)
  Args:     Arg1 => the name_class. You can assign any text, but the words
                'scientific' and 'common' have the special meaning, as
                scientific name and common name, respectively.
            Arg2 .. => the names

=cut

sub name {
    my ($self, $name_class, @names) = @_;
    $self->throw('No name class specified') unless defined $name_class;
    # scientific name should be special, because of its uniqueness.
    return [$self->scientific_name(@names)] if $name_class =~ /scientific/i;
    $self->{'_names_hash'} = {} unless exists $self->{'_names_hash'};
    if(@names){
        $self->{'_names_hash'}->{$name_class} = [] 
            unless exists $self->{'_names_hash'}->{$name_class};
        push @{$self->{'_names_hash'}->{$name_class}}, @names;
    }
    return $self->{'_names_hash'}->{$name_class};
}

=head2 scientific_name

  Title:    scientific_name
  Usage:    my $new_val = $obj->scientific_name($newval);
  Function: Get/Set the scientific name
  Returns:  a scalar text value
  Args:     a scalar text value

=cut

sub scientific_name {
    my $self = shift;
    if( @_ ) { 
	$self->warn("Scientfic name must be set through classification");
    }
    $self->binomial;
}

=head2 parent_taxon_id

  Title   : parent_taxon_id
  Usage   : $self->parent_taxon_id($newval);
            $val = $self->parent_taxon_id;
  Function: Get/Set for parent_taxon_id
  Return  : 
  Args    :    

=cut

sub parent_taxon_id {
    my $self = shift;
    return $self->{'_parent_taxon_id'} = shift if @_;
    return $self->{'_parent_taxon_id'};
}


=head2 species

 Title   : species
 Usage   : $self->species( $species );
           $species = $self->species();
 Function: Get or set the scientific species name.  The species
           name must be in lower case.
 Example : $self->species( 'sapiens' );
 Returns : Scientific species name as string
 Args    : Scientific species name as string

=cut


sub species {
    my($self, $species) = @_;

    if (defined $species) {
        $self->validate_species_name( $species );
        $self->{'_classification'}[0] = $species;
    }
    return $self->{'_classification'}[0];
}

=head2 genus

 Title   : genus
 Usage   : $self->genus( $genus );
           $genus = $self->genus();
 Function: Get or set the scientific genus name.  The genus
           must be in title case.
 Example : $self->genus( 'Homo' );
 Returns : Scientific genus name as string
 Args    : Scientific genus name as string

=cut


sub genus {
    my($self, $genus) = @_;
    if (defined $genus) {
        $self->{'_classification'}[1] = $genus;
    }
    return $self->{'_classification'}[1];
}

=head2 sub_species

 Title   : sub_species
 Usage   : $obj->sub_species($newval)
 Function:
 Returns : value of sub_species
 Args    : newvalue (optional)


=cut

sub sub_species {
    my $self = shift;
    my @v = @{$self->name('sub_species', @_) || []};
    return wantarray ? @v : shift @v;
}

=head2 binomial

 Title   : binomial
 Usage   : $binomial = $self->binomial();
           $binomial = $self->binomial('FULL');
 Function: Returns a string "Genus species", or "Genus species subspecies",
           if the first argument is 'FULL' (and the species has a subspecies).
 Args    : Optionally the string 'FULL' to get the full name including
           the subspecies.

=cut


sub binomial {
    my( $self, $full ) = @_;
    my $rank = $self->rank;
    if( $rank ne 'species' ) {
	$self->warn("Asking for binomial for a $rank node which is not a species, you may not get what you are expecting (genus, species is not available)\n");
    }
    my( $species, $genus ) = $self->classification();
    unless( defined $species) {
	$species = 'sp.';
	$self->warn("requested binomial but classification was not set");
    }
    $genus = ''   unless( defined $genus);
    my $bi = "$genus $species";
    if (defined($full) && $full =~ /full/i) { 
	my $ssp = $self->sub_species;
        $bi .= " $ssp" if $ssp;
    }
    return $bi;
}


=head2 organelle

 Title   : organelle
 Usage   : $self->organelle( $organelle );
           $organelle = $self->organelle();
 Function: Get or set the organelle name
 Example : $self->organelle('Chloroplast')
 Returns : The organelle name in a string
 Args    : String, which is the organelle name

=cut

sub organelle {
    my $self = shift;
    my @v = @{$self->name('organelle', @_) || []};
    return wantarray ? @v : shift @v;
}

=head2 variant

 Title   : variant
 Usage   : $obj->variant($newval)
 Function: Get/set variant information for this species object (strain,
           isolate, etc).
 Example : 
 Returns : value of variant (a scalar)
 Args    : new value (a scalar or undef, optional)


=cut

sub variant{
    my $self = shift;
    my @v = @{$self->name('variant',@_) || []};
    return wantarray ? @v : shift @v;
}

sub validate_species_name {
    my( $self, $string ) = @_;

    return 1 if $string eq "sp.";
    return 1 if $string =~ /^[a-z][\w\s]+$/i;
    $self->throw("Invalid species name '$string'");
}

sub validate_name {
    return 1; # checking is disabled as there is really not much we can
              # enforce HL 2002/10/03
#     my( $self, $string ) = @_;

#     return 1 if $string =~ /^[\w\s\-\,\.]+$/ or
#         $self->throw("Invalid name '$string'");
}


1;
