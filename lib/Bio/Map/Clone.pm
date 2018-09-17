#
# BioPerl module for Bio::Map::clone
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

Bio::Map::Clone - An central map object representing a clone

=head1 SYNOPSIS

   # get the clone object of $clone from the Bio::Map::Clone
   my $cloneobj = $physical->get_cloneobj($clone);

   # acquire all the markers that hit this clone
   foreach my $marker ($cloneobj->each_markerid()) {
       print "   +++$marker\n";
   }

See L<Bio::Map::Position> and L<Bio::Map::PositionI> for more information.

=head1 DESCRIPTION

This object handles the notion of a clone. This clone will
have a name and a position in a map.

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

package Bio::Map::Clone;
use strict;
use Bio::Map::Position;

use base qw(Bio::Root::Root Bio::Map::MappableI);

=head2 new

 Title   : new
 Usage   : my $clone = Bio::Map::Clone->new
                      (
		       -name    => $clone,
		       -markers => \@markers,
		       -contig  => $contig,
		       -type    => $type,
		       -bands   => $bands,
		       -gel     => $gel,
		       -group   => $group,
		       -remark  => $remark,
		       -fpnumber=> $fp_number,
		       -sequencetype  => $seq_type,
		       -sequencestatus=> $seq_status,
		       -fpcremark => $fpc_remark,
		       -matche    => \@ematch,
		       -matcha    => \@amatch,
		       -matchp    => \@pmatch,
		       -range     => Bio::Range->new(-start => $startrange,
						     -end   => $endrange)
		       );
 Function: Initialize a new Bio::Map::Clone object
           Most people will not use this directly but get Clones 
           through L<Bio::MapIO::fpc>
 Returns : L<Bio::Map::Clone> object
 Args    :   -name => marker name string,
	     -markers => array ref of markers,
	     -contig  => contig name string,
	     -type    => type string,
	     -bands   => band string,
	     -gel     => gel string,
	     -group   => group name string,
	     -remark  => remark string,
	     -fpnumber=> FP number string,
	     -sequencetype  => seq type string,
	     -sequencestatus=> seq status string,
	     -fpcremark => FPC remark,
	     -matche    => array ref,
	     -matcha    => array ref,
	     -matchp    => array ref,
	     -range     => L<Bio::Range> object,

=cut

sub new {
   my ($class,@args) = @_;
   my $self= $class->SUPER::new(@args);
   
   my ($name,$markers,$contig,$type,$bands,$gel,$group,
       $remark,$fpnumber,$seqtype,$seqstatus,$fpcremark,
       $matche,$matcha,$matchp,
       $range) = $self->_rearrange([qw(NAME  MARKERS CONTIG TYPE
				       BANDS GEL GROUP REMARK FPNUMBER
				       SEQUENCETYPE SEQUENCESTATUS
				       FPCREMARK MATCHE MATCHA MATCHP
				       RANGE)],@args);

   $self->name($name)                  if defined $name;
   $self->markers($markers)            if defined $markers;
   $self->contigid($contig)            if defined $contig;
   $self->type($type)                  if defined $type;
   $self->bands($bands)                if defined $bands;
   $self->gel($gel)                    if defined $gel;
   $self->group($group)                if defined $group;
   $self->remark($remark)              if defined $remark;
   $self->fp_number($fpnumber)         if defined $fpnumber;
   $self->sequence_type($seqtype)     if defined $seqtype;
   $self->sequence_status($seqstatus) if defined $seqstatus;
   $self->fpc_remark($fpcremark)       if defined $fpcremark;
   $self->range($range)                if defined $range;

   $self->set_match('approx', $matcha) if defined $matcha;
   $self->set_match('pseudo', $matchp) if defined $matchp;
   $self->set_match('exact',  $matche) if defined $matche; 

   return $self;
}

=head1 Access Methods

These methods let you get and set the member variables

=head2 name

 Title   : name
 Usage   : my $name = $cloneobj->name();
 Function: Get/set the name for this Clone
 Returns : scalar representing the current name of this clone
 Args    : none to get, OR string to set

=cut

sub name {
    my ($self) = shift;    
    return $self->{'_name'} = shift if @_;
    return $self->{'_name'};
}

=head2 type

 Title   : type
 Usage   : my $type = $cloneobj->type();
 Function: Get/set the type for this clone
 Returns : scalar representing the current type of this clone
 Args    : none to get, OR string to set

=cut

sub type {
    my ($self) = shift;
    return $self->{'_type'} = shift if @_;
    return $self->{'_type'};
}

=head2 range

 Title   : range
 Usage   : my $range = $cloneobj->range();
 Function: Get/set the range of the contig that this clone covers
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

=head2 match

 Title   : match
 Usage   : @eclone = $cloneobj->match('exact');
           @aclone = $cloneobj->match('approximate');
           @pclone = $cloneobj->match('pseudo');
 Function: get all matching clones
 Returns : list 
 Args    : scalar representing the type of clone to be 
           queried.

=cut

sub match {
  my ($self,$type) = @_;

  $type = "_match" . lc(substr($type, 0, 1));
  return @{$self->{$type} || []};
}

=head2 each_match

 Title   : each_match
 Function: Synonym of the match() method.

=cut

*each_match = \&match;

=head2 set_match

 Title   : set_match
 Usage   : $clone->set_match($type,$values);
 Function: Set the Matches per type
 Returns : None
 Args    : type (one of 'exact' 'approx' 'pseudo')
           array ref of match values

=cut

sub set_match{
   my ($self,$type,$val) = @_;
   $type = "_match" . lc(substr($type, 0, 1));
   $self->{$type} = $val;
}

=head2 gel

 Title   : gel
 Usage   : $clonegel = $cloneobj->gel();
 Function: Get/set the gel number for this clone
 Returns : scalar representing the gel number of this clone
 Args    : none to get, OR string to set

=cut

sub gel {
    my ($self) = shift;
    return $self->{'_gel'} = shift if @_;
    return $self->{'_gel'};
}

=head2 remark

 Title   : remark
 Usage   : $cloneremark = $cloneobj->remark();
 Function: Get/set the remark for this clone
 Returns : scalar representing the current remark of this clone
 Args    : none to get, OR string to set

=cut

sub remark {
    my ($self) = shift;
    return $self->{'_remark'} = shift if @_;
    return $self->{'_remark'};
}

=head2 fp_number

 Title   : fp_number
 Usage   : $clonefpnumber = $cloneobj->fp_number();
 Function: Get/set the fp number for this clone
 Returns : scalar representing the fp number of this clone
 Args    : none to get, OR string to set

=cut

sub fp_number {
    my ($self) = shift;
    return $self->{'_fpnumber'} = shift if @_;
    return $self->{'_fpnumber'};
}

=head2 sequence_type

 Title   : sequence_type
 Usage   : $cloneseqtype = $cloneobj->sequence_type();
 Function: Get/set the sequence type for this clone
 Returns : scalar representing the sequence type of this clone
 Args    : none to get, OR string to set

=cut

sub sequence_type {
    my ($self) = shift;
    return $self->{'_sequencetype'} = shift if @_;
    return $self->{'_sequencetype'};
}

=head2 sequence_status

 Title   : sequence_status
 Usage   : $cloneseqstatus = $cloneobj->sequence_status();
 Function: Get/set the sequence status for this clone
 Returns : scalar representing the sequence status of this clone
 Args    : none to get, OR string to set

=cut

sub sequence_status {
    my ($self) = shift;
    return $self->{'_sequencestatus'} = shift if @_;
    return $self->{'_sequencestatus'};
}

=head2 fpc_remark

 Title   : fpc_remark
 Usage   : $clonefpcremark = $cloneobj->fpc_remark();
 Function: Get/set the fpc remark for this clone
 Returns : scalar representing the fpc remark of this clone
 Args    : none to get, OR string to set

=cut

sub fpc_remark {
    my ($self) = shift;
    return $self->{'_fpcremark'} = shift if @_;
    return $self->{'_fpcremark'};
}

=head2 bands

 Title   : bands
 Usage   : @clonebands = $cloneobj->bands();
 Function: Get/set the bands for this clone
 Returns : liat representing the band of this clone, if 
           readcor = 1 while creating the MapIO object and the
           .cor exists
 Args    : none to get, OR string to set

=cut

sub bands {
    my ($self) = shift; 
    return $self->{'_bands'} = shift if @_;
    return $self->{'_bands'};
}

=head2 group

 Title   : group
 Usage   : $cloneobj->group($chrno);
 Function: Get/set the group number for this clone.
           This is a generic term, used for Linkage-Groups as well as for
           Chromosomes.
 Returns : scalar representing the group number of this clone
 Args    : none to get, OR string to set

=cut

sub group {
    my ($self) = shift;
    return $self->{'_group'} = shift if @_;    
    return $self->{'_group'};
}

=head2 contigid

 Title   : contigid
 Usage   : my $ctg = $cloneobj->contigid();
 Function: Get/set the contig this clone belongs to
 Returns : scalar representing the contig
 Args    : none to get, OR string to set

=cut

sub contigid {
    my ($self) = shift;
    $self->{'_contig'} = shift if @_;
    return $self->{'_contig'} || 0;
}

=head2 each_markerid

 Title   : each_markerid
 Usage   : @markers = $cloneobj->each_markerid();
 Function: retrieves all the elements in a map unordered
 Returns : list of strings (ids)
 Args    : none

 *** This only supplies the ids set with the set_markers method ***
 *** It has nothing to do with actual Bio::Map::MarkerI objects ***

=cut

sub each_markerid {
  my ($self,$value) = @_;   
  return @{$self->{"_markers"}};
}

=head2 set_markers

 Title   : markers
 Usage   : $obj->set_markers($newval)
 Function: Set list of Marker ids (arrayref)
 Returns : None
 Args    : arrayref of strings (ids)

 *** This only sets a list of ids ***
 *** It has nothing to do with actual Bio::Map::MarkerI objects ***

=cut

sub set_markers {
    my ($self,$markers) = @_;
    if( defined $markers && ref($markers) =~ /ARRAY/ ) { 
	$self->{'_markers'} = $markers;
    }
}

1;
