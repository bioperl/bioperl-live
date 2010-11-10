#
# BioPerl module for Bio::TreeIO::TreeEventDBBuilder
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <lapp@bioperl.org>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::TreeEventDBBuilder - Create tree and node records in a persistent storage database from parser events.

=head1 SYNOPSIS

# internal use only

=head1 DESCRIPTION

This object will take events fired by a Bio::TreeIO parser and
populate a a persistent storage database with it.

At present this event handler makes the assumption that tree nodes are
encountered in nested containment (depth first) order, i.e., a node
that is complete is the direct child of the most recent node that is
not yet complete. Formats such as newick work in this way. Formats
that store the hierarchy as an unordered list of edges or parent-child
relationships will at present not work with this handler.

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

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Hilmar Lapp

Email lapp-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::TreeEventDBBuilder;
use strict;

use base qw(Bio::Root::Root Bio::Event::EventHandlerI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::TreeIO::TreeEventDBBuilder->new();
 Function: Builds a new Bio::TreeIO::TreeEventDBBuilder object 
 Returns : Bio::TreeIO::TreeEventDBBuilder
 Args    : Named argument-value pairs:

                  -store   the Bio::DB::Tree::Store object to store the
                           trees in


=cut

sub new {
    my($class,@args) = @_;
    
    my $self = $class->SUPER::new(@args);
    my ($store,$treetype, $nodetype) = 
          $self->_rearrange([qw(STORE
                                )], @args);
    
    if (! (ref($store) && $store->is_a("Bio::DB::Tree::Store"))) {
        $self->throw("Need a Bio::DB::Tree::Store object for -store, not "
                     .ref($store));
    }

    return $self;
}

=head2 store

 Title   : store
 Usage   : $handler->store($newval)
 Function: Gets (or sets) the persistent storage database interface
           object to which the data are to be serialized.

 Example : 
 Returns : value of store (a Bio::DB::Tree::Store object)
 Args    : on set, new value (a Bio::DB::Tree::Store object or undef, optional)

=cut

sub store{
    my $self = shift;

    return $self->{'store'} = shift if @_;
    return $self->{'store'};
}

=head2 SAX methods

=cut

=head2 start_document

 Title   : start_document
 Usage   : $handler->start_document
 Function: Begins a Tree event cycle
 Returns : none 
 Args    : none

=cut

sub start_document {
    my $self = shift;
    $self->{'_keymap'} = {};
    $self->{'_nodestack'} = [];
    $self->{'_elementstack'} = [];
    $self->{'_attrstack'} = [];
    $self->{'_rootnode'} = undef;
    return;
}

=head2 end_document

 Title   : end_document
 Usage   : my @trees = $parser->end_document
 Function: Finishes a document, usually an input file. A document may
           contain multiple trees.
 Returns : An array  Bio::Tree::TreeI
 Args    : none

=cut

sub end_document {
    my $self = shift;
    my $label = shift;

    my $tree = $self->store->get_tree_by_root($self->{'_rootnode'});
    return $tree;       
}

=head2 start_element

 Title   : start_element
 Usage   :
 Function: See Bio::Event::EventHandlerI documentation.
 Example :
 Returns : 

 Args :    hashref with key 'Name' and value being the name of the
           element being started

=cut

sub start_element{
    my $self = shift;
    my $elem = shift->{'Name'};
    
    if( $elem eq 'node' ) {
        push @{$self->{'_nodestack'}}, {};
        push @{$self->{'_elemstack'}}, $elem;
    } elsif ( $elem eq 'tree' ) {
        push @{$self->{'_elemstack'}}, $elem;
    } else {
        push @{$self->{'_attrstack'}}, $elem;
    }
}

=head2 end_element

 Title   : end_element
 Usage   : 
 Function: See Bio::Event::EventHandlerI documentation.
 Returns : none
 Args    : hashref with key 'Name' and value being the name of the
           element that ended

=cut

sub end_element{
    my $self = shift;
    my $elem = shift->{'Name'};

    $self->debug("end of element: $elem\n");

    if( $elem eq 'node' ) {
        my $nodeh = pop @{$self->{'_nodestack'}};
        delete $nodeh->{'-NHXtagname'}; # cleanup - was an internal helper
        # if we have a lineage of parents, and one or more of them is
        # not in the database yet, we need to create them now so we
        # can establish parent-child relationships to create it now
        my $parenth = $self->{'_nodestack'}->[-1] if @{$self->{'_nodestack'}};
        my $parentkey; # will remain undef for the root
        if ($parenth) {
            if (!exists($parenth->{'dbkey'})) {
                $self->_create_parents($self->{'_nodestack'});
                $parentkey = $self->{'_nodestack'}->[-1]->{'dbkey'};
            } else {
                $parentkey = $parenth->{'dbkey'};
            }
        }
        # if we created the node before, we need to update it now;
        # otherwise create it
        if (exists($nodeh->{'dbkey'})) {
            $self->store->update_node($nodeh->{'dbkey'}, $nodeh, $parentkey);
        } else {
            my $dbkey = $self->store->insert_node($nodeh, $parentkey);
            $nodeh->{'dbkey'} = $dbkey;
        }
        # for nested containment, the last element to go off the stack
        # should be the root node
        $self->{'_rootnode'} = $nodeh if (! @{$self->{'_nodestack'}});
        pop @{$self->{'_elemstack'}};
    } elsif ( $elem eq 'tree' ) { 
        pop @{$self->{'_elemstack'}};
    } else {
        pop @{$self->{'_attrstack'}}
    }
}

=head2 in_element

 Title   : in_element
 Usage   :
 Function: See Bio::Event::EventHandlerI documentation.
 Example :
 Returns : boolean
 Args    :


=cut

sub in_element{
    my $self = shift;
    my $e = shift;

    return ($e eq $self->{'_attrstack'}->[-1]);
}

=head2 within_element

 Title   : within_element
 Usage   : 
 Function: See Bio::Event::EventHandlerI documentation.
 Example :
 Returns : boolean
 Args    : Either "node" or "tree".


=cut

sub within_element{
    my $self = shift;
    my $e = shift;
    return ($e eq $self->{'_elemstack'}->[-1]);
}

=head2 characters

 Title   : characters
 Usage   : $handler->characters($text);
 Function: See Bio::Event::EventHandlerI documentation.
 Returns : none
 Args    : text string


=cut

sub characters{
    my $self = shift;
    my $ch = shift;
    if( $self->within_element('node') ) {
       my $nodeh = $self->{'_nodestack'}->[-1];
       if( $self->in_element('bootstrap') ) {
	   # leading/trailing Whitespace-B-Gone
	   $ch =~ s/^\s+//; $ch =~ s/\s+$//;  
	   $nodeh->{'-bootstrap'} = $ch;
       } elsif( $self->in_element('branch_length') ) {
	   # leading/trailing Whitespace-B-Gone
	   $ch =~ s/^\s+//; $ch =~ s/\s+$//;
	   $nodeh->{'-branch_length'} = $ch;
       } elsif( $self->in_element('id')  ) {
	   $nodeh->{'-id'} = $ch;
       } elsif( $self->in_element('description') ) {
	   $nodeh->{'-desc'} = $ch;
       } elsif ( $self->in_element('tag_name') ) {
	   $nodeh->{'-NHXtagname'} = $ch;
       } elsif ( $self->in_element('tag_value') ) {
	   $nodeh->{'-nhx'}->{$nodeh->{'-NHXtagname'}} = $ch;
       } elsif( $self->in_element('leaf') ) {
	   $nodeh->{'-leaf'} = $ch;
       }
   }
   $self->debug("chars: $ch\n");
}

=head1 Private methods

Private methods are used by this module internally, but should not be
called from outside, except perhaps inheriting modules.

=head2 _create_parents

 Title   : _create_parents
 Usage   :
 Function: Takes an array of node hashes, assumed to be stacked in
           nested containment order (i.e., the node at index i is the
           parent of the node at index i+1, and thus recursively
           contains all nodes at indices greater than i), and makes
           sure that they are all created in the database, with proper
           parent-child relationships.

 Example :
 Returns : 
 Args    : A reference to the array of node hashes 

=cut

sub _create_parents{
    my $self = shift;
    my $parents = shift;

    # Assume nested containment - i.e., stacking depth representing
    # containment- and start at the first (highest up) node that isn't
    # in the database yet. Retain the last primary key as the next
    # parent key.
    my $parentkey;
    foreach my $parent (@$parents) {
        my $dbkey = $parent->{'dbkey'};
        # skip to next if key exists
        if ($dbkey) {
            $parentkey = $dbkey;
            next;
        } else {
            # insert as presumably it doesn't exist
            $dbkey = $self->store->insert_node($parent,$parentkey);
            $parent->{'dbkey'} = $dbkey;
            $parentkey = $dbkey;
        }
    }
}

1;
