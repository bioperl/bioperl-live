# $Id$
#
# bioperl module for Bio::Coordinate::Collection
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Coordinate::Collection - Continuous match between two coordinate sets

=head1 SYNOPSIS

  # create Bio::Coordinate::Pairs somehow
  $pair1; $pair2;

  # add them into a Collection
  $collection = Bio::Coordinate::Collection->new;
  $collection->add_mapper($pair1);
  $collection->add_mapper($pair2);

  # create a position and map it
  $pos = Bio::Location::Simple->new (-start => 5, -end => 9 );
  ok $res = $collection->map($pos);
  ok $res->match->start, 1;
  ok $res->match->end, 5;

  # if mapping is many to one (*>1) or many-to-many (*>*)
  # you have to give seq_id not get unrelevant entries
  $pos = Bio::Location::Simple->new
      (-start => 5, -end => 9 -seq_id=>'clone1');

=head1 DESCRIPTION

Generic, context neutral mapper to provide coordinate transforms
between two B<disjoint> coordinate systems. It brings into Bioperl the
functionality from Ewan Birney's Bio::EnsEMBL::Mapper ported into
current bioperl usage.

This class is aimed for representing mapping between whole chromosomes
and contigs, or between contigs and clones, or between sequencing
reads and assembly.

To map coordinates to other direction, you have to swap() the
collection. Keeping track of teh direction and when to id restricitons
are left to the calling code.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                         - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki@ebi.ac.uk
Address:

     EMBL Outstation, European Bioinformatics Institute
     Wellcome Trust Genome Campus, Hinxton
     Cambs. CB10 1SD, United Kingdom

=head1 CONTRIBUTORS

Ewan Birney, birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Coordinate::Collection;
use vars qw(@ISA );
use strict;

# Object preamble - inherits from Bio::Root::Root
use Bio::Root::Root;
use Bio::Coordinate::MapperI;
use Bio::Coordinate::Result;

@ISA = qw(Bio::Root::Root Bio::Coordinate::MapperI);


sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->{'_mappers'} = [];

    my($in, $out) =
	$self->_rearrange([qw(IN
                              OUT
			     )],
			 @args);

    $in  && $self->in($in);
    $out  && $self->out($out);
    return $self; # success - we hope!
}

=head2 add_mapper

 Title   : add_mapper
 Usage   : $obj->add_mapper($mapper)
 Function: Pushes one Bio::Coodinate::MapperI into the list of mappers.
           Sets _is_sorted() to false.
 Example : 
 Returns : 1 when succeeds, 0 for failure.
 Args    : mapper object

=cut

sub add_mapper {
  my ($self,$value) = @_;

  $self->throw("Is not a Bio::Coordinate::MapperI but a [$self]")
      unless defined $value && $value->isa('Bio::Coordinate::MapperI');

  # test pair range lengths
  $self->warn("Coodinates in pair [". $value . ":" .
	      $value->in->seq_id . "/". $value->in->seq_id .
	      "] are not right.")
		  unless $value->test;

  $self->_is_sorted(0);
  push(@{$self->{'_mappers'}},$value);
}


=head2 each_mapper

 Title   : each_mapper
 Usage   : $obj->each_mapper();
 Function: Returns a list of mappers.
 Example : 
 Returns : list of mappers
 Args    : none

=cut

sub each_mapper{
   my ($self,@args) = @_;
   return @{$self->{'_mappers'}}; 
}


=head2 swap

 Title   : swap
 Usage   : $obj->swap;
 Function: Swap the direction of mapping;input <-> output
 Example :
 Returns : 1
 Args    : 

=cut

sub swap {
   my ($self) = @_;

   #@{$self->{'mappers'}} = 
   map {$_->swap} @{$self->{'_mappers'}};
   $self->_sort;
}

=head2 test

 Title   : test
 Usage   : $obj->test;
 Function: test that both components of all pairs are of the same length.
           Ran automatically.
 Example :
 Returns : boolean
 Args    :

=cut

sub test {
   my ($self) = @_;

   my $res = 1;

   foreach my $mapper ($self->each_mapper) {
       $self->warn("Coodinates in pair [". $mapper . ":" .
		   $mapper->in->seq_id . "/". $mapper->in->seq_id .
		   "] are not right.") && ($res = 0)
	   unless $mapper->test;
   }
   $res;
}


=head2 map

 Title   : map
 Usage   : $newpos = $obj->map(5);
 Function: Map the location from the input coordinate system 
           to a new value in the output coordinate system.
 Example :
 Returns : new value in the output coordiante system
 Args    : integer

=cut

sub map {
   my ($self,$value) = @_;

   $self->throw("Need to pass me a value.")
       unless defined $value;
   $self->throw("I need a Bio::Location, not [$value]")
       unless $value->isa('Bio::LocationI');
   $self->throw("No coordinate mappers!")
       unless $self->each_mapper;

   $self->_sort unless $self->_is_sorted;

   my $result = new Bio::Coordinate::Result;

IDMATCH: {

       # bail out now we if are forcing the use of an ID 
       # and it is not in this collection
       last IDMATCH if defined $value->seq_id && 
	   ! $self->{'_in_ids'}->{$value->seq_id};

       foreach my $pair ($self->each_mapper) {

	   # if we are limiting input to a certain ID
	   next if defined $value->seq_id && $value->seq_id ne $pair->in->seq_id;

	   # if we haven't even reached the start, move on
	   next if $pair->in->end < $value->start;
	   # if we have over run, break
	   last if $pair->in->start > $value->end;

	   my $subres = $pair->map($value);
	   $result->add_result($subres);
       }
   }

   unless ($result->each_location) {
       #build one gap;
       my $gap = Bio::Location::Simple->new(-start => $value->start,
					    -end => $value->end,
					    -strand => $value->strand,
					    -location_type => $value->location_type
					   );
       $gap->seq_id($value->seq_id) if defined $value->seq_id;
       bless $gap, 'Bio::Coordinate::Result::Gap';
       $result->add_location($gap);
   }
   return $result;
}


=head2 _sort

 Title   : _sort
 Usage   : $obj->_sort;
 Function: Sort function so that all mappings are sorted by
           input coordinate start
 Example :
 Returns : 1
 Args    : 

=cut

sub _sort{
   my ($self) = @_;

   @{$self->{'_mappers'}} =
       sort { $a->in->start <=> $b->in->start } @{$self->{'_mappers'}};

   #create hashes for sequence ids
   $self->{'_in_ids'} = ();
   $self->{'_out_ids'} = ();
   foreach ($self->each_mapper) {
       $self->{'_in_ids'}->{$_->in->seq_id} = 1;
       $self->{'_out_ids'}->{$_->out->seq_id} = 1;
   }

   $self->_is_sorted(1);
}

=head2 _is_sorted

 Title   : _is_sorted
 Usage   : $newpos = $obj->_is_sorted;
 Function: toggle for whether the (internal) coodinate mapper data are sorted
 Example :
 Returns : boolean
 Args    : boolean

=cut

sub _is_sorted{
   my ($self,$value) = @_;

   $self->{'_is_sorted'} = 1 if defined $value && $value;
   return $self->{'_is_sorted'};
}

1;

