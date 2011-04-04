#
# bioperl module for Bio::Coordinate::ExtrapolatingPair
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

Bio::Coordinate::ExtrapolatingPair - Continuous match between two coordinate sets

=head1 SYNOPSIS


  use Bio::Location::Simple;
  use Bio::Coordinate::ExtrapolatingPair;


  $match1 = Bio::Location::Simple->new 
    (-seq_id => 'propeptide', -start => 21, -end => 40, -strand=>1 );
  $match2 = Bio::Location::Simple->new
    (-seq_id => 'peptide', -start => 1, -end => 20, -strand=>1 );

  $pair = Bio::Coordinate::ExtrapolatingPair->
    new(-in => $match1,
    	-out => $match2,
    	-strict => 1
       );

  $pos = Bio::Location::Simple->new 
      (-start => 40, -end => 60, -strand=> 1 );
  $res = $pair->map($pos);
  $res->start eq 20;
  $res->end eq 20;

=head1 DESCRIPTION

This class represents a one continuous match between two coordinate
systems represented by Bio::Location::Simple objects. The relationship
is directed and reversible. It implements methods to ensure internal
consistency, and map continuous and split locations from one
coordinate system to another.

This class is an elaboration of Bio::Coordinate::Pair. The map
function returns only matches which is the mode needed most of
tehtime. By default the matching regions between coordinate systems
are boundless, so that you can say e.g. that gene starts from here in
the chromosomal coordinate system and extends indefinetely in both
directions. If you want to define the matching regions exactly, you
can do that and set strict() to true.


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

package Bio::Coordinate::ExtrapolatingPair;
use strict;

# Object preamble - inherits from Bio::Root::Root
use Bio::Root::Root;
use Bio::LocationI;

use base qw(Bio::Coordinate::Pair);


sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my($strict) =
	$self->_rearrange([qw(STRICT
			     )],
			 @args);

    $strict  && $self->strict($strict);
    return $self;
}


=head2 strict

 Title   : strict
 Usage   : $obj->strict(1);
 Function: Set and read the strictness of the coordinate system.
 Example :
 Returns : value of input system
 Args    : boolean

=cut

sub strict {
   my ($self,$value) = @_;
   if( defined $value) {
       $self->{'_strict'} = 1 if $value;
   }
   return $self->{'_strict'};
}


=head2 map

 Title   : map
 Usage   : $newpos = $obj->map($loc);
 Function: Map the location from the input coordinate system
           to a new value in the output coordinate system.

           In extrapolating coodinate system there is no location zero.
           Locations are...
 Example :
 Returns : new location in the output coordinate system or undef
 Args    : Bio::Location::Simple

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

   my $match;

   if ($value->isa("Bio::Location::SplitLocationI")) {

       my $split = Bio::Coordinate::Result->new(-seq_id=>$self->out->seq_id);
       foreach my $loc ( sort { $a->start <=> $b->start }
                         $value->sub_Location ) {

           $match = $self->_map($loc);
           $split->add_sub_Location($match) if $match;

       }
       $split->each_Location ? (return $split) : return ;

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

   my ($offset, $start, $end);

   if ($self->strand == -1) {
       $offset = $self->in->end + $self->out->start;
       $start = $offset - $value->end;
       $end = $offset - $value->start ;
   } else { # undef, 0 or 1
       $offset = $self->in->start - $self->out->start;
       $start = $value->start - $offset;
       $end = $value->end - $offset;
   }

   # strict prevents matches outside stated range
   if ($self->strict) {
       return if $start < 0 and $end < 0;
       return if $start > $self->out->end;
       $start = 1 if $start < 0;
       $end = $self->out->end if $end > $self->out->end;
   }

   my $match = Bio::Location::Simple->
       new(-start => $start,
	   -end => $end,
	   -strand => $self->strand,
	   -seq_id => $self->out->seq_id,
	   -location_type => $value->location_type
	  );
   $match->strand($match->strand * $value->strand) if $value->strand;
   bless $match, 'Bio::Coordinate::Result::Match';

   return $match;
}

1;
