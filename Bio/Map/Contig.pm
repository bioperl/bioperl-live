#
# BioPerl module for Bio::Map::Contig
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Gaurav Gupta
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::Contig - A MapI implementation handling the contigs of a
Physical Map (such as FPC)

=head1 SYNOPSIS

    # get the contig object of $contig from the Bio::Map::Physical
    my $ctgobj = $physical->get_contigobj($contig);

    # acquire all the markers that lie in this contig
    foreach my $marker ($ctgobj->each_markerid()) {
	print "   +++$marker\n";
    }

    # find the group of this contig
    print "Group: ",$ctgobj->group(),"\n";

    # find the range of this contig
    print "RANGE: start:",$ctgobj->range()->start(),"\tend: ",
           $ctgobj->range()->end(),"\n";

    # find the position of this contig in $group (chromosome)
    print "Position in Group $group"," = ",$ctgobj->position($group),"\n";


=head1 DESCRIPTION

This is an implementation of Bio::Map::MapI.  It handles the
essential storage of name, species, type, and units as well as in
memory representation of the elements of a map.

Bio::Map::Contig has been tailored to work for FPC physical maps, but
could probably be used for others as well (with the appropriate MapIO
module).

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

=head1 AUTHOR - Gaurav Gupta

Email gaurav@genome.arizona.edu

=head1 CONTRIBUTORS

Sendu Bala  bix@sendu.me.uk

=head1 PROJECT LEADERS

Jamie Hatfield      jamie@genome.arizona.edu
Dr. Cari Soderlund  cari@genome.arizona.edu

=head1 PROJECT DESCRIPTION

The project was done in Arizona Genomics Computational Laboratory (AGCoL)
at University of Arizona.

This work was funded by USDA-IFAFS grant #11180 titled "Web Resources for 
the Computation and Display of Physical Mapping Data".

For more information on this project, please refer: 
  http://www.genome.arizona.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::Contig;
use vars qw($MAPCOUNT);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Range;

use base qw(Bio::Map::SimpleMap);
BEGIN { $MAPCOUNT = 1; }

=head2 new

 Title   : new
 Usage   : my $clone = Bio::Map::Contig->new
                      (
		       -name    => $name,
		       -chr_remark   => $cremark,
		       -user_remark  => $uremark,
		       -trace_remark => $tremark,
		       -group   => $group,
		       -subgroup=> $subgroup,
		       -anchor  => $anchor,
		       -markers => \%markers,
		       -clones  => \%clones,
		       -position => $pos
		       -range    => Bio::Range->new(-start =>$s,-end=>$e),
		       );

 Function: Initialize a new Bio::Map::Contig object
           Most people will not use this directly but get Markers
           through L<Bio::MapIO::fpc>
 Returns : L<Bio::Map::Contig> object
 Args    : ( -name    => name string,
	     -chr_remark   => chr remark string,
	     -user_remark  => userremark string,
	     -trace_remark => tremark string,
	     -group   => group string,
	     -subgroup=> subgroup string,
	     -anchor  => boolean if this is anchored or not,
	     -markers => hashref of contained markers,
	     -clones  => hashref of contained clones,
	     -position => position
	     -range    => L<Bio::Range>

=cut

sub new {
   my ($class,@args) = @_;
   my $self = $class->SUPER::new(@args);

   my ($name,$cremark,$uremark,$tremark,
       $group,$subgroup, $anchor,$markers, $clones,
       $position,$range) = $self->_rearrange([qw(NAME CHR_REMARK USER_REMARK
						 TRACE_REMARK GROUP SUBGROUP
						 ANCHOR MARKERS CLONES
						 POSITION RANGE)],@args);

   $self->name($name)                  if defined $name;
   $self->chr_remark($cremark)         if defined $cremark;
   $self->user_remark($uremark)        if defined $uremark;
   $self->trace_remark($tremark)       if defined $tremark;
   $self->group($group)                if defined $group;
   $self->subgroup($group)             if defined $subgroup;
   $self->anchor($anchor)              if defined $anchor;

   $self->set_markers($markers)        if defined $markers;
   $self->set_clones($clones)          if defined $clones;
   $self->range($range)                if defined $range;
   $self->position($position)          if defined $position;

   return $self;
}

=head2 Modifier methods

All methods present in L<Bio::Map::SimpleMap> are implemented by this class.
Most of the methods are inherited from SimpleMap.  The following methods
have been modified to reflect the needs of physical maps.

=head2 chr_remark

 Title   : chr_remark
 Usage   : my $chrremark = $contigobj->chr_remark();
 Function: Get/set the group remark for this contig
 Returns : scalar representing the current group_remark of this contig
 Args    : none to get, OR string to set

=cut

sub chr_remark {
    my ($self) = shift;
    $self->{'_cremark'} = shift if @_;
    return defined $self->{'_cremark'} ? $self->{'_cremark'} : '';
}

=head2 user_remark

 Title   : user_remark
 Usage   : my $userremark = $contigobj->user_remark();
 Function: Get/set the user remark for this contig
 Returns : scalar representing the current user_remark of this contig
 Args    : none to get, OR string to set

=cut

sub user_remark {
    my ($self) = shift;
    $self->{'_uremark'} = shift if @_;
    return defined $self->{'_uremark'} ? $self->{'_uremark'} : '';
}

=head2 trace_remark

 Title   : trace_remark
 Usage   : my $traceremark = $contigobj->trace_remark();
 Function: Get/set the trace remark for this contig
 Returns : scalar representing the current trace_remark of this contig
 Args    : none to get, OR string to set

=cut

sub trace_remark {
    my ($self) = shift;
    $self->{'_tremark'} = shift if @_;
    return defined $self->{'_tremark'} ? $self->{'_tremark'} : '';
}

=head2 range

 Title   : range
 Usage   : my $range = $contigobj->range();
 Function: Get/set the range for this Contig
 Returns : Bio::Range representing the current range of this contig,
           start and end of the contig can be thus found using:
           my $start = $contigobj->range()->start();
           my $end   = $contigobj->range()->end();
 Args    : none to get, OR Bio::Range to set

=cut

sub range {
    my ($self) = shift;
    return $self->{'_range'} = shift if @_;
    return $self->{'_range'};
}

=head2 position

 Title   : position
 Usage   : $ctgpos = $contigobj->position();
 Function: Get/set the position of the contig in the group
 Returns : scalar representing the position of the contig in the group
 Args    : none to get, OR string to set

=cut

sub position {
    my ($self) = shift;
    $self->{'_position'} = shift if @_;
    return $self->{'_position'} || 0;
}

=head2 anchor

 Title   : anchor
 Usage   : $ctganchor = $contig->anchor();
 Function: Get/set the anchor value for this Contig (True | False)
 Returns : scalar representing the anchor (1 | 0) for this contig
 Args    : none to get, OR string to set

=cut

sub anchor {
    my ($self) = shift;
    return $self->{'_anchor'} = shift if @_;
    return $self->{'_anchor'};
}

=head2 group

 Title   : group
 Usage   : $groupno = $contigobj->group();
 Function: Get/set the group number for this contig.
           This is a generic term, used for Linkage-Groups as well as for
           Chromosomes. 
 Returns : scalar representing the group number of this contig
 Args    : none

=cut

sub group {
    my ($self) = shift;
    $self->{'_group'} = shift if @_;
    return $self->{'_group'} || 0;
}

=head2 subgroup

 Title   : subgroup
 Usage   : $subgroup = $contig->subgroup();	
 Function: Get/set the subgroup for this contig. This is a generic term:
           subgroup here could represent subgroup of a Chromosome or of a
           Linkage Group. The user must take care of which subgroup he/she is
           querying for.	
 Returns : A scalar representing the subgroup of this contig
 Args    : none

=cut

sub subgroup {
    my ($self) = @_;
    return $self->{'_subgroup'} = shift if @_;
    return $self->{'_subgroup'} || 0;
}

=head2 each_cloneid

 Title   : each_cloneid
 Usage   : my @clones  = $map->each_cloneid();
 Function: retrieves all the clone ids in a map unordered
 Returns : list of strings (ids)
 Args    : none

 *** This only supplies the ids set with the set_clones method ***
 *** It has nothing to do with actual Bio::Map::MappableI objects ***

=cut

sub each_cloneid {
    my ($self) = @_;
    return $self->_each_element('clones');
}

=head2 each_markerid

 Title   : each_markerid
 Usage   : my @markers = $map->each_markerid();
 Function: retrieves all the marker ids in a map unordered
 Returns : list of strings (ids)
 Args    : none

 *** This only supplies the ids set with the set_markers method ***
 *** It has nothing to do with actual Bio::Map::MarkerI objects ***

=cut

sub each_markerid {
    my ($self) = @_;
    return $self->_each_element('markers');
}

sub _each_element {
    my ($self, $type) = @_;
    $type = 'clones' if (!defined($type));
    $type = lc("_$type");
    return keys %{$self->{$type} || {}};
}

=head2 set_clones

 Title   : set_clones
 Usage   : $marker->set_clones(\%clones)
 Function: Set the clones hashref
 Returns : None
 Args    : Hashref of clone ids

 *** This only sets a hash of ids ***
 *** It has nothing to do with actual Bio::Map::MappableI objects ***

=cut

sub set_clones {
   my ($self,$clones) = @_;
   if( defined $clones && ref($clones) =~ /HASH/ ) {
       $self->{'_clones'} = $clones;
   }
}

=head2 set_markers

 Title   : markers
 Usage   : $obj->set_markers($newval)
 Function: Set list of Markers (hashref)
 Returns : None
 Args    : Hashref of marker ids

 *** This only sets a hash of ids ***
 *** It has nothing to do with actual Bio::Map::MarkerI objects ***

=cut

sub set_markers {
    my ($self,$markers) = @_;
    if( defined $markers && ref($markers) =~ /HASH/ ) {
	$self->{'_markers'} = $markers;
    }
}

1;