# $Id$
#
# BioPerl module for Bio::Location::SplitLocationI
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::SplitLocationI - Abstract interface of a Location on a Sequence
which has multiple locations (start/end points)

=head1 SYNOPSIS

  # get a SplitLocationI somehow
    print $splitlocation->start, "..", $splitlocation->end, "\n";
    my @sublocs = $splitlocation->sub_Location();

    my $count = 1;
    # print the start/end points of the sub locations
    foreach my $location ( sort { $a->start <=> $b->start } 
			   @sublocs ) {
	printf "sub feature %d [%d..%d]\n", $location->start,$location->end;
        $count++;
    }

=head1 DESCRIPTION

This interface encapsulates the necessary methods for representing the
location of a sequence feature that has more that just a single
start/end pair.  Some examples of this are the annotated exons in a
gene or the annotated CDS in a sequence file.

=head1 FEEDBACK

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bio.perl.org/MailList.html - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@chg.mc.duke.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::Location::SplitLocationI;
use vars qw(@ISA);
use strict;

use Bio::LocationI;
use Carp;

@ISA = qw(Bio::LocationI);

# utility method Prints out a method like: 
# Abstract method stop defined in interface Bio::LocationI not
# implemented by package You::BadLocation

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  
  my $msg = "Abstract method '$caller' defined in interface Bio::ComplexLocationI but not implemented by package $package";
  if( $self->can('throw') ) {
      $self->throw($msg);
  } else {
      confess($msg);
  }
}

=head2 sub_Location

 Title   : sub_Location
 Usage   : @locations = $feat->sub_Location();
 Function: Returns an array of LocationI objects
 Returns : An array
 Args    : none

=cut

sub sub_Location {
    my ($self,@args) = @_;
    $self->_abstractDeath;
}

=head2 add_sub_Location

 Title   : add_sub_Location
 Usage   : $feat->add_sub_Location($locationIobj);
 Function: add an additional sublocation
 Returns : # of current sub locations
 Args    : LocationI object to add

=cut

sub add_sub_Location {
    my ($self,@args) = @_;
    $self->_abstractDeath;
}

=head2

  Title   : min_start
  Usage   : $min_start = $fuzzy->min_start();
  Function: get the minimum starting point
  Returns : the minimum starting point from the contained sublocations
  Args    : none

=cut

sub min_start {
    my ($self, $value) = @_;
    $self->_abstractDeath();
}

=head2

  Title   : max_end
  Usage   : $max_end = $fuzzy->max_end();
  Function: get the maximum ending point
  Returns : the maximum ending point from the contained sublocations
  Args    : none

=cut

sub max_end {
    my ($self, $value) = @_;
    $self->_abstractDeath();
}


# we'll need to override the RangeI methods since our locations will
# not be contiguous.

1;

