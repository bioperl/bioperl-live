# $Id$
#
# bioperl module for Bio::Coordinate::Pair
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Coordinate::Pair - Continuous match between two coordinate sets

=head1 SYNOPSIS

  # to use
  use Bio::Coordinate::Pair;

  $a  = Bio::Coordinate::Pair->new();
  $b  = Bio::Coordinate::Pair -> new ( -id => 3 );

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
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki@ebi.ac.uk
Address:

     EMBL Outstation, European Bioinformatics Institute
     Wellcome Trust Genome Campus, Hinxton
     Cambs. CB10 1SD, United Kingdom

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Coordinate::Pair;
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
       $self->throw("Not a valid input coordinate system name [$value] ")
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
       $self->throw("Not a valid output coordinate system name [$value] ")
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

   my $tmp = $self->out;
   $self->out($self->in);
   $self->in($tmp);
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

   return $self->in->strand * $self->out->strand;
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

   1 if $self->in->end - $self->in->start == $self->out->end - $self->out->start;
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
   $self->throw("Input coordinate system not set")
       unless $self->in;
   $self->throw("Output coordinate system not set")
       unless $self->out;

   my $result = new Bio::Coordinate::Result;

   my $offset = $self->in->start - $self->out->start;
   my $start = $value->start - $offset;
   my $end = $value->end - $offset;

   #within
   #       |-------------------------|
   #            |-|
   if ($start >= $self->out->start and $end <= $self->out->end) {

       my $match = Bio::Location::Simple->new(-start => $start,
					      -end => $end,
					      -strand => $self->strand,
					      -seq_id => $self->out->seq_id,
					      -location_type => $value->location_type
					     );
       $match->strand($match->strand * $value->strand) if $value->strand;
       bless $match, 'Bio::Coordinate::Result::Match';
       $result->add_location($match);
   }
   #out
   #       |-------------------------|
   #   |-|              or              |-|
   elsif ( ($end < $self->out->start or $start > $self->out->end ) or
	   #insertions just outside the range need special settings
	   ($value->location_type eq 'IN-BETWEEN' and 
	    ($end = $self->out->start or $start = $self->out->end)))  {
       #print "$start ", $self->out->start,  " $end ", $self->out->end, "\n";
       my $gap = Bio::Location::Simple->new(-start => $value->start,
					    -end => $value->end,
					    -strand => $self->strand,
					    -seq_id => $self->in->seq_id,
					    -location_type => $value->location_type
					   );
       bless $gap, 'Bio::Coordinate::Result::Gap';
       $result->add_location($gap);
   }
   #partial I
   #       |-------------------------|
   #   |-----|
   elsif ($start < $self->out->start and $end <= $self->out->end ) {

       my $match = Bio::Location::Simple->new(-start => $self->out->start,
					      -end => $end,
					      -strand => $self->strand,
					      -seq_id => $self->out->seq_id,
					      -location_type => $value->location_type
					     );
       $match->strand($match->strand * $value->strand) if $value->strand;
       bless $match, 'Bio::Coordinate::Result::Match';
       $result->add_location($match);

       my $gap = Bio::Location::Simple->new(-start => $value->start,
					    -end => $self->in->start - 1,
					    -strand => $value->strand,
					    -seq_id => $self->in->seq_id
					   );
       bless $gap, 'Bio::Coordinate::Result::Gap';
       $result->add_location($gap);

   }
   #partial II
   #       |-------------------------|
   #                             |------|
   elsif ($start >= $self->out->start and $end > $self->out->end ) {

       my $match = Bio::Location::Simple->new(-start => $start,
					      -end => $self->out->end,
					      -strand => $self->strand,
					      -seq_id => $self->out->seq_id
					     );

       $match->strand($match->strand * $value->strand) if $value->strand;
       bless $match, 'Bio::Coordinate::Result::Match';
       $result->add_location($match);
       $result->match($match);

       my $gap = Bio::Location::Simple->new(-start => $self->in->end + 1,
					    -end => $value->end,
					    -strand => $value->strand,
					    -seq_id => $self->in->seq_id
					   );
       bless $gap, 'Bio::Coordinate::Result::Gap';
       $result->add_location($gap);

   }
   #enveloping
   #       |-------------------------|
   #   |---------------------------------|
   elsif ($start < $self->out->start and $end > $self->out->end ) {

       my $match = Bio::Location::Simple->new(-start => $self->out->start,
					      -end => $self->out->end,
					      -strand => $self->strand,
					      -seq_id => $self->out->seq_id
					     );
       $match->strand($match->strand * $value->strand) if $value->strand;
       bless $match, 'Bio::Coordinate::Result::Match';
       $result->add_location($match);

       my $gap1 = Bio::Location::Simple->new(-start => $value->start,
					     -end => $self->in->start - 1,
					     -strand => $value->strand,
					     -seq_id => $self->in->seq_id
					    );
       bless $gap1, 'Bio::Coordinate::Result::Gap';
       $result->add_location($gap1);

       my $gap2 = Bio::Location::Simple->new(-start => $self->in->end + 1,
					     -end => $value->end,
					     -strand => $value->strand,
					     -seq_id => $self->in->seq_id
					    );
       bless $gap2, 'Bio::Coordinate::Result::Gap';
       $result->add_location($gap2);


   } else {
       #use Data::Dumper; print $start, " ", $end, Dumper $value, $self;
       $self->throw("Should not be here!");
   }
   return $result;
}

1;

