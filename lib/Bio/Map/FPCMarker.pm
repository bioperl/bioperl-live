#
# BioPerl module for Bio::Map::fpcmarker
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Gaurav Gupta <gaurav@genome.arizona.edu>
#
# Copyright Gaurav Gupta
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::FPCMarker - An central map object representing a marker

=head1 SYNOPSIS

   # get the marker object of $marker from the Bio::Map::FPCMarker
   my $markerobj = $physical->get_markerobj($marker);

   # acquire all the clones that hit this marker
   foreach my $clone ($markerobj->each_cloneid()) {
       print "   +++$clone\n";
   }

   # find the position of this marker in $contig
   print "Position in contig $contig"," = ",$markerobj->position($contig),
         "\n";

   # find the group of the marker
   print "Group : ",$markerobj->group();


See L<Bio::Map::Position> and L<Bio::Map::PositionI> for more information.

=head1 DESCRIPTION

This object handles the notion of a marker.
This object is intended to be used by a map parser like fpc.pm.

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

package Bio::Map::FPCMarker;
use strict;
use Bio::Map::Position;
use Time::Local;

use base qw(Bio::Root::Root Bio::Map::MappableI);

=head2 new

 Title   : new
 Usage   : my $clone = Bio::Map::FPCMarker->new
                      (
		       -name    => $marker,
		       -type    => $type,
		       -global  => $global,
		       -frame   => $frame,
		       -group   => $group,
		       -subgroup=> $subgroup,
		       -anchor  => $anchor,
		       -clones  => \%clones,
		       -contigs => \%contigs,
		       -position => \%markerpos,
               -remark => $remark
		       );

 Function: Initialize a new Bio::Map::FPCMarker object
           Most people will not use this directly but get Markers
           through L<Bio::MapIO::fpc>
 Returns : L<Bio::Map::FPCMarker> object
 Args    : -name     => marker name string,
	       -type     => type string,
	       -global   => global position for marker,
	       -frame    => boolean if marker is framework or placement,
	       -group    => group number for marker,
	       -subgroup => subgroup number of marker,
	       -anchor   => boolean if marker is anchored,
	       -clones   => all the clone elements in map (hashref),
	       -contigs  => all the contig elements (hasref),
	       -position => mapping of marker names to map position (hasref),
           -remark   => remarks, separated by newlines

=cut

sub new {
   my ($class,@args) = @_;
   my $self= $class->SUPER::new(@args);

   my ($name,$type,$global,$frame,$group,
       $subgroup, $anchor, $clones,$contigs,
       $positions, $remark) = $self->_rearrange([qw(NAME TYPE GLOBAL FRAME
					   GROUP SUBGROUP ANCHOR
					   CLONES CONTIGS POSITIONS REMARK)],@args);

   $self->name($name)                  if defined $name;
   $self->type($type)                  if defined $type;
   $self->global($global)              if defined $global;
   $self->group($group)                if defined $group;
   $self->subgroup($group)             if defined $subgroup;
   $self->anchor($anchor)              if defined $anchor;
   $self->remark($remark)              if defined $remark;

   $self->set_clones($clones)          if defined $clones;
   $self->set_contigs($contigs)        if defined $contigs;
   $self->set_positions($positions)    if defined $positions;

   return $self;
}

=head1 Access Methods

These methods let you get and set the member variables

=head2 name

 Title   : name
 Usage   : my $name = $markerobj->name();
 Function: Get/set the name for this marker
 Returns : scalar representing the current name of this marker
 Args    : none to get, OR string to set

=cut

sub name {
    my ($self) = shift;
    return $self->{'_name'} = shift if @_;
    return $self->{'_name'};
}

=head2 type

 Title   : type
 Usage   : my $type = $markerobj->type();
 Function: Get/set the type for this marker
 Returns : scalar representing the current type of this marker
 Args    : none to get, OR string to set

=cut

sub type {
    my ($self) = shift;
    return $self->{'_type'} = shift if @_;
    return $self->{'_type'};
}

=head2 global

 Title   : global
 Usage   : my $type = $markerobj->global();
 Function: Get/set the global position for this marker
 Returns : scalar representing the current global position of this marker
 Args    : none to get, OR string to set

=cut

sub global {
    my ($self) = shift;
    return $self->{'_global'} = shift if @_;
    return $self->{'_global'};
}

=head2 anchor

 Title   : anchor
 Usage   : my $anchor = $markerobj->anchor();
 Function: indicate if the Marker is anchored or not (True | False)
 Returns : scalar representing the anchor (1 | 0) for this marker
 Args    : none to get, OR 1|0 to set

=cut

sub anchor {
    my ($self) = shift;
    return $self->{'_anchor'} = shift if @_;
    return $self->{'_anchor'};
}

=head2 framework

 Title   : framework
 Usage   : $frame = $markerobj->framework();
 Function: indicate if the Marker is framework or placement (1 | 0)
 Returns : scalar representing if the marker is framework
           (1 if framework, 0 if placement)
 Args    : none to get, OR 1|0 to set

=cut

sub framework {
    my ($self) = shift;
    return $self->{'_frame'} = shift if @_;
    return $self->{'_frame'};
}

=head2 group

 Title   : group
 Usage   : $grpno = $markerobj->group();
 Function: Get/set the group number for this marker. This is a generic term,
           used for Linkage-Groups as well as for Chromosomes.
 Returns : scalar representing the group number of this marker
 Args    : none to get, OR string to set

=cut

sub group {
    my ($self) = shift;
    $self->{'_group'} = shift if @_;
    return $self->{'_group'} || 0;
}

=head2 subgroup

 Title   : subgroup
 Usage   : $subgroup = $marker->subgroup();	
 Function: Get/set the subgroup for this marker. This is a generic term:
           subgroup here could represent subgroup of a Chromosome or of a
           Linkage Group. The user must take care of which subgroup he/she is
           querying for.	
 Returns : scalar representing the subgroup of this marker
 Args    : none to get, OR string to set

=cut

sub subgroup {
    my ($self) = shift;
    $self->{'_subgroup'} = shift if @_;
    return $self->{'_subgroup'} || 0;
}

=head2 position

 Title   : position
 Usage   : $markerpos = $markerobj->position($ctg);
 Function: get the position of the marker in the contig
 Returns : scalar representing the position of the markernumber of
           the contig
 Args    : $ctg is necessary to look for the position of the marker
           in that contig.

 *** This has nothing to do with an actual Bio::Map::PositionI object ***

=cut

sub position {
    my ($self,$ctg) = @_;
    return 0 unless defined $ctg;

    return 0 unless( defined $self->{'_position'} &&
		     defined $self->{'_position'}{$ctg});
    return $self->{'_position'}{$ctg};
}

=head2 remark

 Title   : remark
 Usage   : $markerremark = $markerobj->remark();
 Function: get the remarks for this marker
 Returns : scalar of newline-separated markers
 Args    : none

=cut

sub remark {
    my ($self) = shift;
    return $self->{'_remark'} = shift if @_;
    return $self->{'_remark'};
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

=head2 each_contigid

 Title   : each_contigid
 Usage   : my @contigs = $map->each_contigid();
 Function: retrieves all the contig ids in a map unordered
 Returns : list of strings (ids)
 Args    : none

 *** This only supplies the ids set with the set_contigs method ***
 *** It has nothing to do with actual Bio::Map::MapI objects ***

=cut

sub each_contigid {
    my ($self) = @_;
    return $self->_each_element('contigs');
}

sub _each_element{
    my ($self, $type) = @_;

    $type = 'clones' unless defined $type;
    $type = lc("_$type");

    return keys %{$self->{$type} || {}};
}

=head2 set_clones

 Title   : set_clones
 Usage   : $marker->set_clones(\%clones)
 Function: Set the clone ids hashref
 Returns : None
 Args    : Hashref of clone ids

 *** This only sets a hash of ids ***
 *** It has nothing to do with actual Bio::Map::MappableI objects ***

=cut

sub set_clones{
   my ($self,$clones) = @_;
   if( defined $clones && ref($clones) =~ /HASH/ ) {
       $self->{'_clones'} = $clones;
   }
}

=head2 set_contigs

 Title   : set_contigs
 Usage   : $marker->set_contigs(\%contigs)
 Function: Set the contig ids hashref
 Returns : None
 Args    : Hashref of contig ids

 *** This only sets a hash of ids ***
 *** It has nothing to do with actual Bio::Map::MapI objects ***

=cut

sub set_contigs{
   my ($self,$contigs) = @_;
   if( defined $contigs && ref($contigs) =~ /HASH/ ) {
       $self->{'_contigs'} = $contigs;
   }
}

=head2 set_positions

 Title   : set_positions
 Usage   : $marker->set_positions(\%markerpos)
 Function: Set the positions hashref
 Returns : None
 Args    : Hashref of marker positions

 *** This only sets a hash of numbers ***
 *** It has nothing to do with actual Bio::Map::PositionI objects ***

=cut

sub set_positions{
   my ($self,$pos) = @_;
   if( defined $pos && ref($pos) =~ /HASH/ ) {
       $self->{'_positions'} = $pos;
   }
}

1;
