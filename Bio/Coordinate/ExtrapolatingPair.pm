# $Id$
#
# bioperl module for Bio::Coordinate::ExtrapolatingPair
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Coordinate::ExtrapolatingPair - Continuous match between two coordinate sets

=head1 SYNOPSIS

  # to use
  use Bio::Coordinate::ExtrapolatingPair;

  $a  = Bio::Coordinate::ExtrapolatingPair->new();
  $b  = Bio::Coordinate::ExtrapolatingPair -> new ( -id => 3 );

=head1 DESCRIPTION

Class

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
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki@ebi.ac.uk
Address:

     EMBL Outstation, European Bioinformatics Institute
     Wellcome Trust Genome Campus, Hinxton
     Cambs. CB10 1SD, United Kingdom

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Coordinate::ExtrapolatingPair;
use vars qw(@ISA );
use strict;

# Object preamble - inherits from Bio::Root::Root
use Bio::Root::Root;
use Bio::LocationI;
#use Bio::Coordinate::Result;
use Bio::Coordinate::Result::Match;
use Bio::Coordinate::Pair;

@ISA = qw(Bio::Coordinate::Pair);


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
 Function: Set and read the strictput coordinate system.
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
 Usage   : $newpos = $obj->map(5);
 Function: Map the location from the input coordinate system 
           to a new value in the output coordinate system.

           In extrapolating coodinate system there is no location zero.
           Locations are...
 Example :
 Returns : new value in the output coordinate system
 Args    : integer

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

#   my $result = new Bio::Coordinate::Result;

   my ($offset, $start, $end);
   if ($self->strand == -1) {
       $start = -1 * ($value->end - $self->in->end  - 1);
       $end = -1* ($value->start - $self->in->end  - 1);
   } else { # undef, 0 or 1
       $offset = $self->in->start - $self->out->start;
       $start = $value->start - $offset;
       $end = $value->end - $offset;
   }

   # strict prevents matches outside stated range
   if ($self->strict) {
       return undef if $start < 0 and $end < 0;
       return undef if $start > $self->out->end;
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
#   $result->add_Location($match);
   return $match;
}

1;
