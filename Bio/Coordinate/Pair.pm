#
# bioperl module for Bio::Coordinate::Pair
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Coordinate::Pair - Continuous match between two coordinate sets

=head1 SYNOPSIS

  use Bio::Location::Simple;
  use Bio::Coordinate::Pair;

  my $match1 = Bio::Location::Simple->new 
      (-seq_id => 'propeptide', -start => 21, -end => 40, -strand=>1 );
  my $match2 = Bio::Location::Simple->new
      (-seq_id => 'peptide', -start => 1, -end => 20, -strand=>1 );
  my $pair = Bio::Coordinate::Pair->new(-in => $match1,
  					-out => $match2
                                        );
  # location to match
  $pos = Bio::Location::Simple->new 
      (-start => 25, -end => 25, -strand=> -1 );

  # results are in a Bio::Coordinate::Result
  # they can be Matches and Gaps; are  Bio::LocationIs
  $res = $pair->map($pos);
  $res->isa('Bio::Coordinate::Result');
  $res->each_match == 1;
  $res->each_gap == 0;
  $res->each_Location == 1;
  $res->match->start == 5;
  $res->match->end == 5;
  $res->match->strand == -1;
  $res->match->seq_id eq 'peptide';


=head1 DESCRIPTION

This class represents a one continuous match between two coordinate
systems represented by Bio::Location::Simple objects. The relationship
is directed and reversible. It implements methods to ensure internal
consistency, and map continuous and split locations from one
coordinate system to another.

The map() method returns Bio::Coordinate::Results with
Bio::Coordinate::Result::Gaps. The calling code have to deal (process
or ignore) them.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

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

report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Coordinate::Pair;
use strict;

# Object preamble - inherits from Bio::Root::Root
use Bio::Coordinate::Result;
use Bio::Coordinate::Result::Match;
use Bio::Coordinate::Result::Gap;

use base qw(Bio::Root::Root Bio::Coordinate::MapperI);


sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my($in, $out) =
	$self->_rearrange([qw(IN
                              OUT
			     )],
			 @args);

    $in  && $self->in($in);
    $out  && $self->out($out);
    return $self; # success - we hope!
}

=head2 in

 Title   : in
 Usage   : $obj->in('peptide');
 Function: Set and read the input coordinate system.
 Example :
 Returns : value of input system
 Args    : new value (optional), Bio::LocationI

=cut

sub in {
   my ($self,$value) = @_;
   if( defined $value) {
       $self->throw("Not a valid input Bio::Location [$value] ")
	   unless $value->isa('Bio::LocationI');
       $self->{'_in'} = $value;
   }
   return $self->{'_in'};
}


=head2 out

 Title   : out
 Usage   : $obj->out('peptide');
 Function: Set and read the output coordinate system.
 Example :
 Returns : value of output system
 Args    : new value (optional), Bio::LocationI

=cut

sub out {
   my ($self,$value) = @_;
   if( defined $value) {
       $self->throw("Not a valid output coordinate Bio::Location [$value] ")
	   unless $value->isa('Bio::LocationI');
       $self->{'_out'} = $value;
   }
   return $self->{'_out'};
}


=head2 swap

 Title   : swap
 Usage   : $obj->swap;
 Function: Swap the direction of mapping; input <-> output
 Example :
 Returns : 1
 Args    : 

=cut

sub swap {
   my ($self) = @_;
   ($self->{'_in'}, $self->{'_out'}) = ($self->{'_out'}, $self->{'_in'});
   return 1;
}

=head2 strand

 Title   : strand
 Usage   : $obj->strand;
 Function: Get strand value for the pair
 Example :
 Returns : ( 1 | 0 | -1 )
 Args    :

=cut

sub strand {
   my ($self) = @_;
   $self->warn("Outgoing coordinates are not defined")
       unless $self->out;
   $self->warn("Incoming coordinates are not defined")
       unless $self->in;

   return ($self->in->strand || 0) * ($self->out->strand || 0);
}

=head2 test

 Title   : test
 Usage   : $obj->test;
 Function: test that both components are of the same length
 Example :
 Returns : ( 1 | undef )
 Args    :

=cut

sub test {
   my ($self) = @_;
   $self->warn("Outgoing coordinates are not defined")
       unless $self->out;
   $self->warn("Incoming coordinates are not defined")
       unless $self->in;
   return ($self->in->end - $self->in->start) == ($self->out->end - $self->out->start);
}


=head2 map

 Title   : map
 Usage   : $newpos = $obj->map($pos);
 Function: Map the location from the input coordinate system
           to a new value in the output coordinate system.
 Example :
 Returns : new Bio::LocationI in the output coordinate system or undef
 Args    : Bio::LocationI object

=cut

sub map {
   my ($self,$value) = @_;

   $self->throw("Need to pass me a value.")
       unless defined $value;
   $self->throw("I need a Bio::Location, not [$value]")
       unless $value->isa('Bio::LocationI');
   $self->throw("Input coordinate system not set")
       unless $self->in;
   $self->throw("Output coordinate system not set")
       unless $self->out;


   if ($value->isa("Bio::Location::SplitLocationI")) {

       my $result = Bio::Coordinate::Result->new();
       foreach my $loc ( $value->sub_Location(1) ) {
           my $res = $self->_map($loc);
           map { $result->add_sub_Location($_) } $res->each_Location;
       }
       return $result;
   } else {
       return $self->_map($value);
   }
}


=head2 _map

 Title   : _map
 Usage   : $newpos = $obj->_map($simpleloc);
 Function: Internal method that does the actual mapping. Called
           multiple times by map() if the location to be mapped is a
           split location
 Example :
 Returns : new location in the output coordinate system or undef
 Args    : Bio::Location::Simple

=cut

sub _map {
   my ($self,$value) = @_;

   my $result = Bio::Coordinate::Result->new();

   my $offset = $self->in->start - $self->out->start;
   my $start  = $value->start - $offset;
   my $end    = $value->end - $offset;
   
   my $match = Bio::Location::Simple->new;
   $match->location_type($value->location_type);
   $match->strand($self->strand);

   #within
   #       |-------------------------|
   #            |-|
   if ($start >= $self->out->start and $end <= $self->out->end) {

       $match->seq_id($self->out->seq_id);
       $result->seq_id($self->out->seq_id);

       if ($self->strand >= 0) {
	   $match->start($start);
	   $match->end($end);
       } else {
	   $match->start($self->out->end - $end + $self->out->start);
	   $match->end($self->out->end - $start + $self->out->start);
       }
       if ($value->strand) {
	   $match->strand($match->strand * $value->strand);
	   $result->strand($match->strand);
       }
       bless $match, 'Bio::Coordinate::Result::Match';
       $result->add_sub_Location($match);
   }
   #out
   #       |-------------------------|
   #   |-|              or              |-|
   elsif ( ($end < $self->out->start or $start > $self->out->end ) or
	   #insertions just outside the range need special settings
	   ($value->location_type eq 'IN-BETWEEN' and 
	    ($end = $self->out->start or $start = $self->out->end)))  {

       $match->seq_id($self->in->seq_id);
       $result->seq_id($self->in->seq_id);
       $match->start($value->start);
       $match->end($value->end);
       $match->strand($value->strand);

       bless $match, 'Bio::Coordinate::Result::Gap';
       $result->add_sub_Location($match);
   }
   #partial I
   #       |-------------------------|
   #   |-----|
   elsif ($start < $self->out->start and $end <= $self->out->end ) {

       $result->seq_id($self->out->seq_id);
       if ($value->strand) {
	   $match->strand($match->strand * $value->strand);
	   $result->strand($match->strand);
       }
       my $gap = Bio::Location::Simple->new;
       $gap->start($value->start);
       $gap->end($self->in->start - 1);
       $gap->strand($value->strand);
       $gap->seq_id($self->in->seq_id);

       bless $gap, 'Bio::Coordinate::Result::Gap';
       $result->add_sub_Location($gap);

       # match
       $match->seq_id($self->out->seq_id);

       if ($self->strand >= 0) {
	   $match->start($self->out->start);
	   $match->end($end);
       } else {
	   $match->start($self->out->end - $end + $self->out->start);
	   $match->end($self->out->end);
       }
       bless $match, 'Bio::Coordinate::Result::Match';
       $result->add_sub_Location($match);
   }
   #partial II
   #       |-------------------------|
   #                             |------|
   elsif ($start >= $self->out->start and $end > $self->out->end ) {

       $match->seq_id($self->out->seq_id);
       $result->seq_id($self->out->seq_id);
       if ($value->strand) {
	   $match->strand($match->strand * $value->strand);
	   $result->strand($match->strand);
       }
       if ($self->strand >= 0) {
	   $match->start($start);
	   $match->end($self->out->end);
       } else {
	   $match->start($self->out->start);
	   $match->end($self->out->end - $start + $self->out->start);
       }
       bless $match, 'Bio::Coordinate::Result::Match';
       $result->add_sub_Location($match);

       my $gap = Bio::Location::Simple->new;
       $gap->start($self->in->end + 1);
       $gap->end($value->end);
       $gap->strand($value->strand);
       $gap->seq_id($self->in->seq_id);
       bless $gap, 'Bio::Coordinate::Result::Gap';
       $result->add_sub_Location($gap);

   }
   #enveloping
   #       |-------------------------|
   #   |---------------------------------|
   elsif ($start < $self->out->start and $end > $self->out->end ) {

       $result->seq_id($self->out->seq_id);
       if ($value->strand) {
	   $match->strand($match->strand * $value->strand);
	   $result->strand($match->strand);
       }
       # gap1
       my $gap1 = Bio::Location::Simple->new;
       $gap1->start($value->start);
       $gap1->end($self->in->start - 1);
       $gap1->strand($value->strand);
       $gap1->seq_id($self->in->seq_id);
       bless $gap1, 'Bio::Coordinate::Result::Gap';
       $result->add_sub_Location($gap1);

       # match
       $match->seq_id($self->out->seq_id);

       $match->start($self->out->start);
       $match->end($self->out->end);
       bless $match, 'Bio::Coordinate::Result::Match';
       $result->add_sub_Location($match);

       # gap2
       my $gap2 = Bio::Location::Simple->new;
       $gap2->start($self->in->end + 1);
       $gap2->end($value->end);
       $gap2->strand($value->strand);
       $gap2->seq_id($self->in->seq_id);
       bless $gap2, 'Bio::Coordinate::Result::Gap';
       $result->add_sub_Location($gap2);

   } else {
       $self->throw("Should not be here!");
   }
   return $result;
}


1;
