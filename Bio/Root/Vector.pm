#-----------------------------------------------------------------------------
# PACKAGE : Bio::Root::Vector.pm
# AUTHOR  : Steve Chervitz (sac@bioperl.org)
# CREATED : 15 April 1997
# REVISION: $Id$
# STATUS  : Alpha
#
# WARNING: This is considered an experimental module.
#
# For documentation, run this module through pod2html
# (preferably from Perl v5.004 or better).
#
# MODIFIED:
#    sac --- Fri Nov  6 14:24:48 1998
#       * Added destroy() method (experimental).
#    0.023, 20 Jul 1998, sac:
#      * Improved memory management (_destroy_master()).
#
#   Copyright (c) 1997 Steve Chervitz. All Rights Reserved.
#             This module is free software; you can redistribute it and/or
#             modify it under the same terms as Perl itself.
#-----------------------------------------------------------------------------

package Bio::Root::Vector;

use Bio::Root::Global qw(:devel);
use Bio::Root::Object ();

# @ISA = qw(Bio::Root::Object);  # Eventually perhaps...

use vars qw($ID $VERSION);
$ID = 'Bio::Root::Vector';
$VERSION = 0.04;

use strict;
my @SORT_BY = ('rank','name');

## POD Documentation:

=head1 NAME

Bio::Root::Vector - Interface for managing linked lists of Perl5 objects.

=head1 SYNOPSIS

=head2 Object Creation

B<At present, Vector objects cannot be instantiated.> This
package is currently designed to be inherited along with another class
that provides a constructor (e.g., B<Bio::Root::Object.pm>).
The Vector provides a set of methods that can then be used for managing
sets of objects.

See L<the USAGE section | USAGE>.

=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.


=head1 DESCRIPTION

Bio::Root::Vector.pm  provides an interface for creating and manipulating
dynamic sets (linked lists) of Perl5 objects. This is an abstract class (ie.,
there is no constructor) and as such is expected to be inherited along with
some other class (see note above).

Vectors are handy when, for example, an object may contain one or more
other objects of a certain class. The container object knows only
that is has at least one such object; the multiplex nature of the contained
object is managed by the contained object via its Vector interface.
The methods for adding, removing, counting, listing, and sorting all objects
are bundled together in Vector.pm.

Thus, the current Bio::Root::Vector class is somewhat of a cross between an
interator and a composite design pattern. At present, a number of classes
utilize Bio::Root::Vector's composite-like behavior to implement a composite
pattern (Bio::SeqManager.pm, for example).
This is not necessarily ideal and is expected to change.

=head1 USAGE

For a usage demo of Bio::Root::Vector.pm see the scripts in the 
examples/root_object/vector directory.


=head1 DEPENDENCIES

Bio::Root::Vector.pm does not directly inherit from B<Bio::Root::Object.pm> but
creates an manager object which does.

=head1 BUGS/FEATURES

By default, all Vectors are doubly-linked lists. This relieves one from
the burden of worrying about whether a given Vector is single- or doubly-linked.
However, when generating lots of Vectors or extremely large vectors, memory
becomes an issue. In particular, signaling the GC to free the memory for
an object when you want to remove it. B<Until this memory issue is resolved,
the use of Vector.pm is not recommended for large projects.>

Although it is not required, all objects within a vector are expected
to derive from the same class (package). Support for heterogeneous
Vectors has not been explicitly added (but they should work fine, as long
as you know what you're doing).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other Bioperl modules.
Send your comments and suggestions preferably to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs and
their resolution. Bug reports can be submitted via email or the web:

    bioperl-bugs@bio.perl.org
    http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR 

Steve Chervitz  E<lt>sac@bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports and comments.

=head1 VERSION

Bio::Root::Vector.pm, 0.04

=head1 TODO

=over 4

=item * (Maybe) create an container class version of this module

to permit Vectors to be instantiated. Thus, instead of inherited
from both Object.pm and Vector.pm, you could create a Vector.pm object.

=item * Improve documentation.

=back

=head1 SEE ALSO

  Bio::Root::Object.pm    - Core object
  Bio::Root::Err.pm       - Error/Exception object
  Bio::Root::Global.pm    - Manages global variables/constants

  http://bio.perl.org/Projects/modules.html  - Online module documentation
  http://bio.perl.org/                       - Bioperl Project Homepage

=head1 ACKNOWLEDGEMENTS

This module was developed under the auspices of the Saccharomyces Genome
Database:
    http://genome-www.stanford.edu/Saccharomyces

=head1 COPYRIGHT

Copyright (c) 1996-98 Steve Chervitz. All Rights Reserved.
This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=cut


#'
##
###
#### END of main POD documentation.
###
##
#

=head1 APPENDIX

Methods beginning with a leading underscore are considered private
and are intended for internal use by this module. They are
B<not> considered part of the public interface and are described here
for documentation purposes only.

=cut

########################################################
#                CONSTRUCTOR                           #
########################################################

## No constructor. See _set_master() for construction of {Master} data member.

## Destructor: Use remove_all() or remove().

# New Idea for destructor
#-------------
sub destroy {
#-------------
    my $self = shift;
    local($^W) = 0;
    undef $self->{'_prev'};
    undef $self->{'_next'};
    undef $self->{'_rank'};
    undef $self->{'_master'};
}

#####################################################################################
##                                  ACCESSORS                                      ##
#####################################################################################


=head2 set_rank

 Purpose  : To set an object's rank to an arbitrary numeric
          : value to be used for sorting the objects of the Vector.
 Usage    : $self->set_rank(-RANK    =>numeric_ranking_data,
	  :     	    -RANK_BY =>ranking_criterion_string);
          : or without the named parameters as (rank, rank_by).
 Throws   : warning (if rank is set without also setting RANK_BY)
 Comments : It is redundant for every object to store RANK_BY data.
          : For this reason, the RANK_BY data is stored with the master
          : object associated with the vector.

See Also  : L<rank>(), L<rank_by>()

=cut

#-------------'
sub set_rank {
#-------------
    my( $self, %param) = @_;

    $self->_set_master($self) unless $self->{'_master'}->{'_set'};

    my($rank, $rank_by) = $self->master->_rearrange([qw(RANK RANK_BY)], %param);

    $DEBUG==1 and do{ print STDERR "$ID:set_rank() = $rank; Criteria: $rank_by."; <STDIN>; };

    $self->{'_rank'} = ($rank || undef);
    $self->{'_master'}->{'_rankBy'} = ($rank_by || undef);
    if( defined $self->{'_rank'} and not defined $self->{'_master'}->{'_rankBy'} ) {
	return $self->master->warn('Rank defined without ranking criteria.');
    }	
    1;
}

sub _set_rank_by {
    my( $self, $arg) = @_;
    $self->{'_master'}->{'_rankBy'} = $arg || 'unknown';
}

sub _set_master {
    ## A vector does not need a master object unless it needs to grow.
    my($self,$obj) = @_;

#    print "$ID: _set_master() new Master object for ${\$obj->name}."; <STDIN>;

    require Bio::Root::Object;
    my $master = {};
    bless $master, 'Bio::Root::Object';

    $master->{'_set'}  = 1;  ## Special member indicating that this method has been called.
                            ## Necessary since perl will generate an anonymous {Master}
                            ## hash ref on the fly. This ref will not be blessed however.
    $master->{'_first'} = $obj;
    $master->{'_last'}  = $obj;
    $master->{'_size'}  = 1;
    $master->{'_index'}->{$obj->name()} = $obj;
    $self->{'_master'} = $master;

    $self->{'_rank'} = 1;
    $self->{'_prev'} = undef;
    $self->{'_next'} = undef;
#    $self->{'_master'}->{'_rankBy'} = undef;  # Default rank is the order of addition to Vector.
}

sub _destroy_master {
# This is called when the last object in the vector is being remove()d
    my $self = shift;

    return if !$self->master or !$self->master->{'_set'};

    my $master = $self->master;

    ## Get rid of the Vector master object.
    ref $master->{'_first'} and (%{$master->{'_first'}} = (), undef $master->{'_first'});
    ref $master->{'_last'}  and (%{$master->{'_last'}} = (), undef $master->{'_last'});
    ref $master->{'_index'} and (%{$master->{'_index'}} = (), undef $master->{'_index'});
    %{$master} = ();
    undef $master;
}


=head2 clone_vector

 Purpose  : Call this method to clone the whole vector.
          : NOT calling this method will extract the vector element.
 Usage    : $self->clone_vector();
 Throws   : Exception if argument is not an object reference.
 Comments : This method is usually called from within a module's
          : _set_clone() method for modules that inherit from
          : Bio::Root::Vector.pm.

=cut

#-----------------'
sub clone_vector {
#-----------------
    my($self, $obj) = @_;

    ref($obj) || throw($self, "Can't clone $ID object: Not an object ref ($obj)");

    $self->{'_prev'} = $obj->{'_prev'};
    $self->{'_next'} = $obj->{'_next'};
    $self->{'_rank'} = $obj->{'_rank'};
    $self->{'_master'} = $obj->{'_master'};
}


=head2 prev

 Purpose  : Returns the previous object in the Vector or undef
          : if on first object.
 Usage    : $self->prev

=cut

#--------
sub prev { my $self = shift; $self->{'_prev'}; }
#--------



=head2 next

 Purpose  : Returns the next object in the Vector or undef
          : if on last object.
 Usage    : $self->next

=cut

#--------
sub next { my $self = shift; $self->{'_next'}; }
#--------



=head2 first

 Purpose  : Returns the first object in the Vector or $self
          : if Vector size = 1.
 Usage    : $self->first

=cut

#----------
sub first  {
#----------
    my $self = shift;
    defined $self->{'_master'} ? $self->{'_master'}->{'_first'} : $self;
}


=head2 last

 Purpose  : Returns the last object in the Vector or
          : $self if Vector size = 1.
 Usage    : $self->last

=cut

#-------
sub last   {
#-------
    my $self = shift;
    defined $self->{'_master'} ? $self->{'_master'}->{'_last'} : $self;
}



=head2 rank

 Purpose  : Returns the rank of the current object or 1
          : if rank is not defined.
 Usage    : $self->rank

See Also  : L<set_rank>()

=cut

#---------
sub rank { my $self = shift; $self->{'_rank'} || 1; }
#---------



=head2 rank_by

 Purpose  : Returns the ranking criterion or the string 'order of addition'
          : if rankBy has not been explicitly set.
 Usage    : $self->rank_by

See Also  : L<set_rank>()

=cut

#-----------
sub rank_by {
#-----------
    my $self = shift;
    defined $self->{'_master'} ? ($self->{'_master'}->{'_rankBy'}||'order of addition')
	: 'unranked';
}



=head2 size

 Purpose  : Returns the number of objects currently in the Vector.
 Usage    : $self->size

=cut

#---------
sub size {
#---------
    my $self = shift;
    defined $self->{'_master'} ? $self->{'_master'}->{'_size'} : 1;
}


=head2 master

 Purpose  : Returns the master object associated with the Vector.
          : (should probably be a private method).
 Usage    : $self->master

=cut

#-----------
sub master { my $self = shift; $self->{'_master'}; }
#-----------


## Not sure what these potentially dangerous methods are used for.
## Should be unnecessary and probably can be removed.
sub set_prev { my($self,$obj) = @_; $self->{'_prev'} = $obj;  }
sub set_next { my($self,$obj) = @_; $self->{'_next'} = $obj; }

#############################################################################
#                           INSTANCE METHODS                               ##
#############################################################################


=head2 is_first

 Purpose  : Test whether the current object is the first in the Vector.
 Usage    : $self->is_first

=cut

#------------
sub is_first { my($self) = shift; return not defined $self->{'_prev'}; }
#------------


=head2 is_last

 Purpose  : Test whether the current object is the last in the Vector.
 Usage    : $self->is_last

=cut

#------------
sub is_last { my($self) = shift; return not defined $self->{'_next'}; }
#------------




=head2 get

 Purpose  : Retrives an object from the Vector given its name.
          : Returns undef if the object cannot be found.
 Usage    : $self->get(object_name)
 Examples : $self->get($obj->name)

See Also  : L<list>()

=cut

#--------
sub get {
#--------
    my($self,$name) = @_;

    my ($obj);
#    print "$ID get(): getting $name\n";

    if($self->{'_master'}->{'_set'}) {
#	my @names = keys %{$self->{'_master'}->{'_index'}};
#	print "$ID: names in hash:\n@names";<STDIN>;
#	print "  returning $self->{'_master'}->{'_index'}->{$name}\n";
	local($^W) = 0;
	$obj = $self->{'_master'}->{'_index'}->{$name};
    }

    elsif($self->name =~ /$name/i) {
#	print "  returning self\n";
	$obj = $self;
    }

    if(not ref $obj) {
	$self->throw("Can't get object named \"$name\": object not set or name undefined.");
    }
    $obj;
}

## Former strategy: hunt through the list for the object.
## No longer needed since master indexes all objects.
#	 do{
#	     if($obj->name eq $name) { return $obj; }
#	
#	 } while($obj = $current->prev());




=head2 add

 Purpose  : Add an object to the end of a Vector.
 Usage    : $self->add(object_reference)

See also  : L<insert>(), L<remove>()

=cut

#--------
sub add {
#--------
    my($self,$new,$index) = @_;

    $self->_set_master($self) unless $self->{'_master'}->{'_set'};

#    print "\n\nADDING TO VECTOR ${\ref $self} ${\$self->name}\nFOR PARENT: ${\ref $self->parent} ${\$self->parent->name}\n\n";

    $self->{'_next'} = $new;
    $new->{'_prev'} = $self;
    $self->{'_master'}->{'_last'} = $new;
    $self->{'_master'}->{'_size'}++;
    $new->{'_master'} = $self->{'_master'};
    $new->_incrementRank();
    $new->Bio::Root::Vector::_index();

#    printf "NEW CONTENTS: (n=%s)\n", $self->size;
#    my $obj = $self->first;
#    my $count=0;
#    do { print "\n","Object #",++$count,"\n";
#	       $obj->display;
#     } while($obj=$obj->next);
#    <STDIN>;
}


sub _index {
    my($self) = @_;
    my $name = $self->name;

    # Generate unique name, if necessary, for indexing purposes.
    if( not $name or $name =~ /anonymous/) {
	$name ||= '';
	$name .= $self->size();
    }
#    print "$ID: _index() called for $name\n";

    $self->{'_master'}->{'_index'}->{$name} = $self;
}

sub _incrementRank {
    my $self = shift;
    return if not defined $self->{'_prev'};
    $self->{'_rank'} = $self->{'_prev'}->rank() + 1;
}


=head2 remove

 Purpose  : Remove the current object from the Vector.
 Usage    : $self->remove([-RET=>'first'|'last'|'next'|'prev'], [-UPDATE => 0|1])
 Returns  : The first, last, next, or previous object in the Vector
          : depending on the value of the -RET parameter.
          : Default = next.
 Argument : Named parameters: (TAGS CAN BE ALL UPPER OR ALL LOWER CASE)
          :    -RET    => string: 'first', 'last', 'prev' 'next'
          :               THis indicates the object to be returned.
          :    -UPDATE => boolean, if non-zero all objects in the vector
          :               will be re-ranked.
 Comments : The -UPDATE parameter should be set to true to force a re-updating
          : of the rank data for each object. (default = 0, no update).

See Also  : L<add>(), L<insert>(), L<shift>(), L<chop>()

=cut

#-----------
sub remove {
#-----------
    my($self,%param) = @_;
    my $updateRank = $param{-UPDATE} || $param{'-update'}  || 0;
    my $ret = $param{-RET} || $param{'-ret'} || 'next';

    $DEBUG==2 && do{ print STDERR "$ID: removing ${\$self->name}; ret = $ret";<STDIN>; };

    ## This set of conditionals involves primarily pointer shuffling.
    ## The special case of destroying a vector of size 1 is handled.

    if($self->is_first()) {
	$DEBUG==2 && print STDERR "---> removing first object: ${\$self->name()}.\n";
	if($self->is_last) {
#	    print "Removing only object in vector: ${\$self->name}.\n";
	    $self->_destroy_master();
	    return $self->destroy;
	} else {
	    undef ($self->{'_next'}->{'_prev'});
	    $self->_update_first($self->{'_next'});
	}

    } elsif($self->is_last()) {
	$DEBUG==2 && print STDERR "---> removing last object: ${\$self->name()}.\n";
	undef ($self->{'_prev'}->{'_next'});
	$self->_update_last($self->{'_prev'});

    } else {
	$DEBUG==2 && print STDERR "---> removing internal object.\n";
	$self->{'_prev'}->{'_next'} = $self->{'_next'};
	$self->{'_next'}->{'_prev'} = $self->{'_prev'};
    }

    $updateRank && $self->_update_rank();
    $self->{'_master'}->{'_size'}--;

#    print "new vector size = ",$self->size,"\n"; <STDIN>;

    my($retObj);

    if( $self->size) {
	if($ret eq 'first') { $retObj = $self->first(); }
	elsif($ret eq 'last') { $retObj = $self->last(); }
	elsif($ret eq 'next') { $retObj = $self->next(); }
	elsif($ret eq 'prev') { $retObj = $self->prev(); }
    }

    ## Destroy the object.
#    $self->destroy;

    $DEBUG && do{ print STDERR "$ID: returning ${\$retObj->name}";<STDIN>; };

    $retObj;
}

sub _update_first {
    my($self,$first) = @_;
    $DEBUG && print STDERR "Updating first.\n";
    undef ($first->{'_prev'});
    $self->{'_master'}->{'_first'} = $first;
}

sub _update_last {
    my($self,$last) = @_;
    $DEBUG && print STDERR "Updating last.\n";
    undef ($last->{'_next'});
    $self->{'_master'}->{'_last'} = $last;
}


=head2 remove_all

 Purpose  : Remove all objects currently in the Vector.
 Usage    : $self->remove_all

See Also  : L<remove>(), L<shift>(), L<chop>()

=cut

#---------------
sub remove_all {
#---------------
    my($self,%param) = @_;

    $DEBUG==2 && print STDERR "DESTROYING VECTOR $self ${\$self->name}";

#    print "$ID Removing all.";

    $self = $self->first();

    while(ref $self) {
#	print "$ID: removing ${\$self->name}\n";
	$self = $self->remove(-RET=>'next');
    }
}


=head2 shift

 Purpose  : Remove the first object from the Vector.
          : This is a wrapper for remove().
 Usage    : $self->shift([-RET=>'first'|'last'|'next'|'prev'])
 Returns  : The object returned by remove().

See Also  : L<remove>(), L<chop>()

=cut

#---------
sub shift {
#---------
    my($self,%param) = @_;
    $self = $self->first();
    $self = $self->remove(%param);
}


=head2 chop

 Purpose  : Remove the last object from the Vector.
          : This is a wrapper for remove().
 Usage    : $self->chop([-RET=>'first'|'last'|'next'|'prev'])
 Returns  : The object returned by remove().

See Also  : L<remove>(), L<shift>()

=cut

#----------
sub chop {
#----------
    my($self,%param) = @_;
    $self = $self->last();
    $self = $self->remove(%param);
}



=head2 insert

 Purpose  : Insert a new object into the vector relative to the current object.
 Usage    : $self->insert(object_ref, ['before'|'after'])
 Examples : $self->insert($obj)  # Default insert after $self
          : $self->insert($obj,'before')
 Returns  : The new number of objects in the vector (int).
 Throws   : exception if the first argument is not a reference.

See Also  : L<add>(), L<remove>()

=cut

#-----------
sub insert {
#-----------
    my($self,$object,$where) = @_;
    my($first);
    $where ||= 'after';

    $self->_set_master($self) unless $self->{'_master'}->{'_set'};

    ref($object) || return $self->master->throw("Can't insert. Not an object: $object");

    if($where eq 'before') {
	$object->{'_next'} = $self;
	$object->{'_prev'} = $self->{'_prev'};
	$object->{'_master'} = $self->{'_master'};
	$self->{'_prev'}->{'_next'} = $object;
	$self->{'_prev'} = $object;
    } else {
	$object->{'_prev'} = $self;
	$object->{'_next'} = $self->{'_next'};
	$object->{'_master'} = $self->{'_master'};
	$self->{'_next'}->{'_prev'} = $object;
	$self->{'_next'} = $object;
    }
    $self->{'_master'}->{'_size'}++;
    $object->Bio::Root::Vector::_index();  ##Fully qualified to disambiguate a potentially common method name.
    $self->_update_rank();
}

sub _update_rank {
    my($self) = @_;
    my $current = $self->first();
    my $count = 0;
    $DEBUG && print STDERR "$ID: Updating rank.\n";
    do{
	$count++;
	$current->{'_rank'} = $count;

    } while($current = $current->next());
}


=head2 list

 Purpose  : Returns objects in the Vector as an array or array slice.
 Usage    : $self->list([starting_object|'first'] [,ending_object|'last'])
 Examples : $self->list
          : $self->list('first',$self->prev)

See Also  : L<get>()

=cut

#----------
sub list {
#----------
    my($self,$start,$stop) = @_;
    my(@list);

    $start ||= 1;
    $stop  ||= 'last';

    if( $start =~ /first|beg|start/i or $start <= 1 ) {
	$start = $self->first();
    }

    if( $stop =~ /last|end|stop/i ) {
	$stop = $self->last();
    }

    ref($start) || ($start = $self->first());
    ref($stop)  || ($stop = $self->last());

    my $obj = $start;
    my $fini = 0;
    do{
	push @list, $obj;
	if($obj eq $stop) { $fini = 1; }
    } while( $obj = $obj->next() and !$fini);

    @list;
}


=head2 sort

 Purpose  : Sort the objects in the Vector.
 Usage    : $self->sort(['rank'|'name'], [reverse])
 Returns  : The last object of the sorted Vector.
 Argument : First argument can be 'name' or 'rank' to sort on
          : the object's name or rank data member, respectively.
          : If reverse is non-zero, sort will be in reverse order.
 Example  : $self->sort()   #  Default sort by rank, not reverse.
          : $self->sort('name','reverse')

=cut

#---------'
sub sort {
#---------
    my ($self,$sortBy,$reverse) = @_;
    my (@unsortedList,@sortedList);

    $sortBy ||= 'rank';
    my $rankBy = $self->rank_by;

    ### Build the initial unsorted list.
    my $obj = $self->first();
    do{
	push @unsortedList, $obj;
    } while( $obj = $obj->next());

#    print "UNSORTED LIST:\n";
#    foreach(@unsortedList) {print $_->name().' '};<STDIN>;

    ### Sort it.
    if( $sortBy =~ /rank/i) {
#	print "sorting by rank";
	if($reverse) {
#	    print " (reverse).\n";
	    @sortedList = reverse sort _sort_by_rank @unsortedList;
	} else {
	    @sortedList = sort _sort_by_rank @unsortedList;
	}
    } elsif( $sortBy =~ /name/i) {
#	print "sorting by name";
	if($reverse) {
#	    print "(reverse).\n";
	    @sortedList = reverse sort _sort_by_name @unsortedList;
	} else {
	    @sortedList = sort _sort_by_name @unsortedList;
	}
    } else {
#	print "unknown sort criteria: $sortBy\n";
	$self->warn("Invalid sorting criteria: $sortBy.",
		    "Sorting by rank.");
	@sortedList = sort _sort_by_rank @unsortedList;
    }


#    if($reverse) { @sortedList = reverse sort @sortedList;  }

#    print "SORTED LIST:\n";
#    foreach(@sortedList) {print $_->name().' '};<STDIN>;

    ### Re-load the Vector with the sorted list.
    my $count=0;

    $self = $sortedList[0];
    $self->_set_master($self);
    $self->_set_rank_by($rankBy);

    my($i);
    my $current = $self;
    for($i=1; $i<@sortedList; $current=$sortedList[$i], $i++) {
	$current->add($sortedList[$i]);
	if($i==$#sortedList) { $sortedList[$i]->{'_next'} = undef;}
    }

    $self->last();
}

sub _sort_by_rank { my $aRank = $a->rank(); my $bRank = $b->rank(); $aRank <=> $bRank; }

sub _sort_by_name { my $aName = $a->name(); my $bName = $b->name(); $aName cmp $bName; }



=head2 valid_any

 Purpose  : Determine if at least one object in the Vector is valid.
 Usage    : $self->valid_any
 Status   : Deprecated.
 Comments : A non-valid object should throw an exception that must be
          : be caught an dealt with on the spot.

See Also  : B<Bio::Root::Object::valid()>

=cut

#-------------
sub valid_any {
#-------------
    my $self = &shift(@_);

    my $obj = $self->first();
    do{
	return 1 if $obj->valid();
    } while( $obj = $obj->next());

   return undef;
}


=head2 valid_all

 Purpose  : Determine if all objects in the Vector are valid.
 Usage    : $self->valid_all
 Comments : A non-valid object should throw an exception that must be
          : be caught an dealt with on the spot.

See Also  : B<Bio::Root::Object::valid()>

=cut

#--------------
sub valid_all {
#--------------
    my $self = &shift(@_);

    my $obj = $self->first();
    do{
	return  unless $obj->valid();
    } while( $obj = $obj->next());

   return 1;
}

sub _display_stats {
# This could be fleshed out a bit...

    my( $self, $OUT ) = @_;

    printf ( $OUT "%-11s %s\n","RANK:", $self->rank());
    printf ( $OUT "%-11s %s\n","RANK BY:", $self->rank_by());
}

1;
__END__

#####################################################################################
#                                 END OF CLASS                                      #
#####################################################################################

=head1 FOR DEVELOPERS ONLY

=head2 Data Members

Information about the various data members of this module is provided for those
wishing to modify or understand the code. Two things to bear in mind:

=over 4

=item 1 Do NOT rely on these in any code outside of this module.

All data members are prefixed with an underscore to signify that they are private.
Always use accessor methods. If the accessor doesn't exist or is inadequate,
create or modify an accessor (and let me know, too!).

=item 2 This documentation may be incomplete and out of date.

It is easy for this documentation to become obsolete as this module is still evolving.
Always double check this info and search for members not described here.

=back

Bio::Root::Vector.pm objects currently cannot be instantiated. Vector.pm must be inherited
along with Bio::Root::Object.pm (or an object that provides a constructor).
Vector.pm defines the following fields:

 FIELD          VALUE
 ------------------------------------------------------------------------
  _prev         Reference to the previous object in the Vector.

  _next         Reference to the next object in the Vector.

  _rank         Rank relative to other objects in the Vector.
  	        Default rank = chronological order of addition to the Vector.

  _master       A reference to an Bio::Root::Object that acts as a manager for
 	        the given Vector. There is only one master per Vector.
 	        A master object is only needed when the Vector size is >1.
 	        The master object manages the following Vector data:

 	        _first  - Reference to the first object in the Vector.
 	        _last   - Reference to the last object in the Vector.
 	        _size   - Total number of objects in the Vector.
 	        _rankBy - Criteria used to rank the object.
 	                 Default: chronological order of addition.
 	        _index  - Hash reference for quick access to any object
 	                 based on its name.
 	        Bio::Root::Object{'_err'} - Holds any errors affecting the
 	                             Vector as a whole.

=cut

1;
