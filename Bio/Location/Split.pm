# $Id$
#
# BioPerl module for Bio::Location::SplitLocation
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Location::Split - Implementation of a Location on a Sequence
which has multiple locations (start/end points)

=head1 SYNOPSIS

    my $splitlocation = new Bio::Location::Split();
    $splitlocation->add_sub_Location(new Bio::Location::Simple(-start=>1,
							       -end=>30,
							       -strand=>1));
    $splitlocation->add_sub_Location(new Bio::Location::Simple(-start=>50,
							       -end=>61,
							       -strand=>1));   
    my @sublocs = $splitlocation->sub_Location();

    my $count = 1;
    # print the start/end points of the sub locations
    foreach my $location ( sort { $a->start <=> $b->start } 
			   @sublocs ) {
	printf "sub feature %d [%d..%d]\n", $location->start,$location->end;
        $count++;
    }

=head1 DESCRIPTION

This implementation handles locations which span more than one
start/end location.

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


package Bio::Location::Split;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::Location::SplitLocationI;
use Bio::Location::Simple;

@ISA = qw(Bio::Location::Simple Bio::Location::SplitLocationI );

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ( $locations ) = $self->_rearrange([qw(LOCATIONS)], @args);
    if( defined $locations && ref($locations) =~ /array/i ) {
	$self->{'_sublocations'} = $locations;
    } else { 
	$self->{'_sublocations'} = [];
    }
    return $self;
}

=head2 sub_Location

 Title   : sub_Location
 Usage   : @locations = $feat->sub_Location();
 Function: Returns an array of LocationI objects
 Returns : An array
 Args    : none

=cut

sub sub_Location {
    my ($self) = @_;
    return @{$self->{'_sublocations'}};
}

=head2 add_sub_Location

 Title   : add_sub_Location
 Usage   : $feat->add_sub_Location(@locationIobjs);
 Function: add an additional sublocation
 Returns : # of current sub locations
 Args    : list of LocationI object(s) to add

=cut

sub add_sub_Location {
    my ($self,@args) = @_;
    my @locs;
    foreach my $l ( @args ) {
	if( ref($l) && $l->isa('Bio::LocationI')) {
	    push @locs, $l;
	} else {
	    $self->warn("($l) was not a valid Bio::LocationI object to add as a subLocation");
	}
    }
    push @{$self->{'_sublocations'}}, @locs;
}


# we'll probably need to override the RangeI methods since our locations will
# not be contiguous.

1;

