# BioPerl module for Bio::Map::Microsatellite
#
# Cared for by Chad Matsalla <bioinformatics1@dieselwurks.com>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::Microsatellite - An object representing a Microsatellite marker.

=head1 SYNOPSIS

  $o_usat = new Bio::Map::Microsatellite
      (-name=>'Chad Super Marker 2',
       -sequence => 'gctgactgatcatatatatatatatatatatatatatatatcgcgatcgtga',
       -motif => 'at',
       -repeats => 15,
       -repeat_start_position => 11
       );

  $sequence_before_usat = $o_usat->get_leading_flank();
  $sequence_after_usat = $o_usat->get_trailing_flank();


=head1 DESCRIPTION

This object handles the notion of an Microsatellite. This microsatellite can
be placed on a (linear) Map or used on its own.  If this Microsatellites
will be used in a mapping context (it doesn't have to, you know) it can have
multiple positions in a map. For information about a Microsatellite's position
in a map one must query the associate PositionI object which is accessible
through the position() method.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Chad Matsalla

Email bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

Heikki Lehvaslaiho heikki@ebi.ac.uk
Lincoln Stein      lstein@cshl.org
Jason Stajich      jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::Microsatellite;
use vars qw(@ISA);
use strict;
use Bio::Root::Root;
use Bio::Map::Marker;

@ISA = qw(Bio::Map::Marker);

=head2 new

 Title   : new
 Usage   : $o_usat = 
 Function: Builds a new Bio::Map::Microsatellite object
 Returns : Bio::Map::Microsatellite
 Args    :
	-name    => name of this microsatellite (optional, string,
		default 'Unknown microsatellite')
        -positions => position(s) for this marker in maps[optional],
                An array reference of tuples (array refs themselves)
                Each tuple conatins a Bio::Map::MapI-inherited object and a 
		Bio::Map::PositionI-inherited obj, no default)
	-sequence => the sequence of this microsatellite (optional,
		 scalar, no default)
	-motif => the repeat motif of this microsatellite (optional,
		 scalar, no default)
	-repeats => the number of motif repeats for this microsatellite
		(optional, scalar, no default)
	-repeat_start_position => the starting position of the
		microsatellite in this sequence. The first base of the
		sequence is position "1". (optional, scalar, no default)

 Note    : Creating a Bio::Map::Microsatellite object with no position
	might be useful for microsatellite people wanting to embrace
	and extend this module. <raising hand> Me! Me! Me!
	- using repeat_start_position will trigger a mechinism to
	calculate a value for repeat_end_position. 


=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($map, $position, $sequence, $motif, $repeats, $start) = 
	$self->_rearrange([qw(MAP
			      POSITION
			      SEQUENCE 
			      MOTIF 
			      REPEATS 
			      REPEAT_START_POSITION
			      )], @args);
    if( ! $self->name ) { 
	$self->name('Unnamed microsatellite');
    }
    $map && $self->map($map);
    $position && $self->position($position);
    $sequence && $self->sequence($sequence);
    $self->motif(defined $motif ? $motif : 'Unknown motif'); 
    $repeats && $self->repeats($repeats);
    $start && $self->repeat_start_position($start);
    return $self;
}

=head2 Bio::Map::Marker methods

=cut

=head2 position

 Title   : position
 Usage   : my $position = $mappable->position($map); OR
           $mappable->position($map,$position); OR
 Function: Get/Set the Bio::Map::PositionI for a mappable element
           in a specific Map
 Returns : Bio::Map::PositionI
 Args    : $map =Bio::Map::MapI # Map we are talking about
           $position = Bio::Map::PositionI # Position we want to set


=head2 name($new_name)

 Title   : name($new_name)
 Usage   : $o_usat->name($new_name) _or_
	   my $name = $o_usat->name()
 Function: Get/Set the name for this Microsatellite
 Returns : A scalar representing the current name of this Microsatellite
 Args    : If provided, the current name of this Microsatellite
	   will be set to $new_name.

=head2 motif($new_motif)

 Title   : motif($new_motif)
 Usage   : my $motif = $o_usat->motif($new_motif) _or_
	my $motif = $o_usat->motif()
 Function: Get/Set the repeat motif for this Microsatellite
 Returns : A scalar representing the current repeat motif of this
	Microsatellite.
 Args    : If provided, the current repeat motif of this Microsatellite
	will be set to $new_motif.

=cut

sub motif {
	my ($self,$motif) = @_;
	if ($motif) {
		$self->{'_motif'} = $motif;
	}
	return $self->{'_motif'};	
}

=head2 sequence($new_sequence)

 Title   : sequence($new_sequence)
 Usage   : my $sequence = $o_usat->sequence($new_sequence) _or_
	my $sequence = $o_usat->sequence()
 Function: Get/Set the sequence for this Microsatellite
 Returns : A scalar representing the current sequence of this
	Microsatellite.
 Args    : If provided, the current sequence of this Microsatellite
	will be set to $new_sequence.

=cut

sub sequence {
	my ($self,$sequence) = @_;
	if ($sequence) {
		$self->{'_sequence'} = $sequence;
	}
	return $self->{'_sequence'};	
}

=head2 repeats($new_repeats)

 Title   : repeats($new_repeats)
 Usage   : my $repeats = $o_usat->repeats($new_repeats) _or_
	my $repeats = $o_usat->repeats()
 Function: Get/Set the repeat repeats for this Microsatellite
 Returns : A scalar representing the current number of repeats of this
	Microsatellite.
 Args    : If provided, the current number of repeats of this
	Microsatellite will be set to $new_repeats.

=cut

sub repeats {
	my ($self,$repeats) = @_;
	if ($repeats) {
		$self->{'_repeats'} = $repeats;
	}
	return $self->{'_repeats'};	
}

=head2 repeat_start_position($new_repeat_start_position)

 Title   : repeat_start_position($new_repeat_start_position)
 Usage   : my $repeat_start_position =
	$o_usat->repeat_start_position($new_repeat_start_position) _or_
	my $repeat_start_position = $o_usat->repeat_start_position()
 Function: Get/Set the repeat repeat_start_position for this
	Microsatellite
 Returns : A scalar representing the repeat start position for this 
	Microsatellite.
 Args    : If provided, the current repeat start position of this
	Microsatellite will be set to $new_repeat_start_position.
	This method will also try to set the repeat end position. This
	depends on having information for the motif and the number of
	repeats. If you want to use methods like get_trailing_flank or
	get_leading flank, be careful to include the right information.

=cut

sub repeat_start_position {
	my ($self,$repeat_start_position) = @_;
	if ($repeat_start_position) {
		$self->{'_repeat_start_position'} = $repeat_start_position;
		$self->repeat_end_position("set");
	}
	return $self->{'_repeat_start_position'};	
}

=head2 repeat_end_position($value)

 Title   : repeat_end_position($set)
 Usage   : $new_repeat_end_position =
		$o_usat->repeat_end_position("set"); _or_
	$new_repeat_end_position =
		$o_usat->repeat_end_position($value); _or_
	$current_repeat_end_position = $o_usat->repeat_end_position();
 Function: get/set the end position of the repeat in this sequence
 Returns : A scalar representing the base index of the end of the
	repeat in this Microsatellite. The first base in the sequence
	is base 1.
 Args    : A scalar representing a value, the string "set", or no
	argument (see Notes).
 Notes   : If you do not provide an argument to this method, the current
	end position of the repeat in this Microsatellite will be
	returned (a scalar).
	If you provide the string "set" to this method it will set the
	end position based on the start position, the length of the
	motif, and the nuimber of repeats.
	If you specify a value the current end position of the repeat
	will be set to that value. This is a really bad idea. Don't do
	it.

=cut

#'
sub repeat_end_position {
    my ($self,$caller) = @_;
    if( defined $caller ) { 
	if ($caller eq "set") {
	    $self->{'_repeat_end_position'} = 
		$self->{'_repeat_start_position'} + 
		    (length($self->motif()) * $self->repeats());
	}
	elsif ($caller) {
	    $self->{'_repeat_end_position'} = $caller;
	}
    }
    return $self->{'_repeat_end_position'};
}

=head2 equals

 Title   : equals
 Usage   : if( $mappable->equals($mapable2)) ...
 Function: Test if a position is equal to another position
 Returns : boolean
 Args    : Bio::Map::MappableI

=cut

sub equals {
	my ($self,@args) = @_;
	$self->warn("equals is not yet implemented in ".
		    ref($self)." yet. Check back real soon!");
}

=head2 less_than

 Title   : less_than
 Usage   : if( $mappable->less_than($m2) ) ...
 Function: Tests if a position is less than another position
 Returns : boolean
 Args    : Bio::Map::MappableI

=cut

sub less_than {
	my ($self,@args) = @_;
	$self->warn("less_then is not yet implemented in ".
		    ref($self)." yet. Check back real soon!");
}

=head2 greater_than

 Title   : greater_than
 Usage   : if( $mappable->greater_than($m2) ) ...
 Function: Tests if position is greater than another position
 Returns : boolean
 Args    : Bio::Map::MappableI

=cut

sub greater_than {
	my ($self,@args) = @_;
	$self->warn("greater_then is not yet implemented in ".
		    ref($self)." yet. Check back real soon!");
}

=head2 get_leading_flank()

 Title   : get_leading_flank()
 Usage   : $leading_sequence = $o_usat->get_leading_flank();
 Returns : A scalar representing the sequence before the repeats in this
	Microsatellite.
 Args    : None.

=cut

sub get_leading_flank {
	my $self = shift;
	return substr $self->sequence(),0,$self->repeat_start_position-1;

}

=head2 get_trailing_flank()

 Title   : get_trailing_flank()
 Usage   : $trailing_flank = $o_usat->get_trailing_flank();
 Returns : A scalar representing the sequence after the repeats in this
	Microsatellite.
 Args    : None.

=cut

sub get_trailing_flank {
	my $self = shift;
	return substr $self->sequence(),$self->repeat_end_position()-1;
}


1;


