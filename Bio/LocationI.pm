# $Id$
#
# BioPerl module for Bio::LocationI
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::LocationI - Abstract interface of a Location on a Sequence

=head1 SYNOPSIS

    # get a LocationI somehow
    printf( "start = %d, end = %d, strand = %s, seq_id = %s\n", 
	    $location->start, $location->end, $location->strand,
	    $location->seq_id);
    print "location str is ", $location->to_FTstring(), "\n"; 


=head1 DESCRIPTION

This Interface defines the methods for a Bio::LocationI, an object
which encapsulates a location on a biological sequence.  Locations
need not be attached to actual sequences as they are stand alone
objects.  LocationI objects are used by Bio::SeqFeatureI objects to
manage and represent locations for a Sequence Feature.

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

package Bio::LocationI;
use vars qw(@ISA $coord_policy);
use strict;

use Bio::RangeI;
use Bio::Location::WidestCoordPolicy;
use Carp;

@ISA = qw(Bio::RangeI);

BEGIN {
    $coord_policy = Bio::Location::WidestCoordPolicy->new();
}

# utility method Prints out a method like: 
# Abstract method stop defined in interface Bio::LocationI not
# implemented by package You::BadLocation

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  
  my $msg = "Abstract method '$caller' defined in interface Bio::LocationI but not implemented by package $package";
  if( $self->can('throw') ) {
      $self->throw($msg);
  } else {
      confess($msg);
  }
}

=head2 start

  Title   : start
  Usage   : $start = $location->start();
  Function: Get the start coordinate of this location as defined by the
            currently active coordinate computation policy. In simple cases,
            this will return the same number as min_start() and max_start(),
            in more ambiguous cases like fuzzy locations the number may be
            equal to one or neither of both.

            We override this here from RangeI in order to delegate 'get' to
            a Bio::Location::CoordinatePolicy implementing object. Implementing
            classes may also wish to provide 'set' functionality, in which
            case they *must* override this method. The implementation
            provided here will throw an exception if called with arguments.

  Returns : A positive integer value.
  Args    : none

=cut

sub start {
    my ($self,@args) = @_;

    $self->_abstractDeath() if @args;
    return $self->coordinate_policy()->start($self);
}

=head2 end

  Title   : end
  Usage   : $end = $location->end();
  Function: Get the end coordinate of this location as defined by the
            currently active coordinate computation policy. In simple cases,
            this will return the same number as min_end() and max_end(),
            in more ambiguous cases like fuzzy locations the number may be
            equal to one or neither of both.

            We override this here from RangeI in order to delegate 'get' to
            a Bio::Location::CoordinatePolicy implementing object. Implementing
            classes may also wish to provide 'set' functionality, in which
            case they *must* override this method. The implementation
            provided here will throw an exception if called with arguments.

  Returns : A positive integer value.
  Args    : none

=cut

sub end {
    my ($self,@args) = @_;

    $self->_abstractDeath() if @args;
    return $self->coordinate_policy()->end($self);
}

=head2 min_start

  Title   : min_start
  Usage   : my $minstart = $location->min_start();
  Function: Get minimum starting point of feature.

            Note that an implementation must not call start() in this method.

  Returns : integer or undef if no minimum starting point.
  Args    : none

=cut

sub min_start {
    my($self) = @_;
    $self->_abstractDeath();
}

=head2 max_start

  Title   : max_start
  Usage   : my $maxstart = $location->max_start();
  Function: Get maximum starting point of feature.

            Note that an implementation must not call start() in this method
            unless start() is overridden such as not to delegate to the
            coordinate computation policy object.

  Returns : integer or undef if no maximum starting point.
  Args    : none

=cut

sub max_start {
    my($self) = @_;
    $self->_abstractDeath();
}

=head2 start_pos_type

  Title   : start_pos_type
  Usage   : my $start_pos_type = $location->start_pos_type();
  Function: Get start position type encoded as text

            Known valid values are 'BEFORE' (<5..100), 'AFTER' (>5..100), 
            'EXACT' (5..100), 'WITHIN' ((5.10)..100), 'BETWEEN', (5^6), with
            their meaning best explained by their GenBank/EMBL location string
            encoding in brackets.

  Returns : string ('BEFORE', 'AFTER', 'EXACT','WITHIN', 'BETWEEN')
  Args    : none

=cut

sub start_pos_type {
    my($self) = @_;
    $self->_abstractDeath();
}

=head2 min_end

  Title   : min_end
  Usage   : my $minend = $location->min_end();
  Function: Get minimum ending point of feature. 

            Note that an implementation must not call end() in this method
            unless end() is overridden such as not to delegate to the
            coordinate computation policy object.

  Returns : integer or undef if no minimum ending point.
  Args    : none

=cut

sub min_end {
    my($self) = @_;
    $self->_abstractDeath();
}

=head2 max_end

  Title   : max_end
  Usage   : my $maxend = $location->max_end();
  Function: Get maximum ending point of feature.

            Note that an implementation must not call end() in this method
            unless end() is overridden such as not to delegate to the
            coordinate computation policy object.

  Returns : integer or undef if no maximum ending point.
  Args    : none

=cut

sub max_end {
    my($self) = @_;
    $self->_abstractDeath();
}

=head2 end_pos_type

  Title   : end_pos_type
  Usage   : my $end_pos_type = $location->end_pos_type();
  Function: Get end position encoded as text.

            Known valid values are 'BEFORE' (5..<100), 'AFTER' (5..>100), 
            'EXACT' (5..100), 'WITHIN' (5..(90.100)), 'BETWEEN', (5^6), with
            their meaning best explained by their GenBank/EMBL location string
            encoding in brackets.

  Returns : string ('BEFORE', 'AFTER', 'EXACT','WITHIN', 'BETWEEN')
  Args    : none

=cut

sub end_pos_type {
    my($self) = @_;
    $self->_abstractDeath();
}

=head2 seq_id

  Title   : seq_id
  Usage   : my $seqid = $location->seq_id();
  Function: Get/Set seq_id that location refers to
  Returns : seq_id (a string)
  Args    : [optional] seq_id value to set

=cut

sub seq_id {
    my ($self, $seqid) = @_;
    if( defined $seqid ) {
	$self->{'_seqid'} = $seqid;
    }
    return $self->{'_seqid'};
}

=head2 coordinate_policy

  Title   : coordinate_policy
  Usage   : $policy = $location->coordinate_policy();
            $location->coordinate_policy($mypolicy); # set may not be possible
  Function: Get the coordinate computing policy employed by this object.

            See Bio::Location::CoordinatePolicyI for documentation about
            the policy object and its use.

            The interface *does not* require implementing classes to accept
            setting of a different policy. The implementation provided here
            does, however, allow to do so.

            Implementors of this interface are expected to initialize every
            new instance with a CoordinatePolicyI object. The implementation
            provided here will return a default policy object if none has
            been set yet. To change this default policy object call this
            method as a class method with an appropriate argument. Note that
            in this case only subsequently created Location objects will be
            affected.
            
  Returns : A Bio::Location::CoordinatePolicyI implementing object.
  Args    : On set, a Bio::Location::CoordinatePolicyI implementing object.

=cut

sub coordinate_policy {
    my ($self, $policy) = @_;

    if(defined($policy)) {
	if(! $policy->isa('Bio::Location::CoordinatePolicyI')) {
	    $self->throw("Object of class ".ref($policy)." does not implement".
			 " Bio::Location::CoordinatePolicyI");
	}
	if(ref($self)) {
	    $self->{'_coordpolicy'} = $policy;
	} else {
	    # called as class method
	    $coord_policy = $policy;
	}
    }
    return (ref($self) && exists($self->{'_coordpolicy'}) ?
	    $self->{'_coordpolicy'} : $coord_policy);
}

=head2 to_FTstring

  Title   : to_FTstring
  Usage   : my $locstr = $location->to_FTstring()
  Function: returns the FeatureTable string of this location
  Returns : string
  Args    : none

=cut

sub to_FTstring { 
    my($self) = @_;
    $self->_abstractDeath();
}
1;

