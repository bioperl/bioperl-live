#--------------------------------------------------------------------
# PACKAGE : PersonManager.pm
# PURPOSE : Dummy object for testing Bio::Root::Vector.pm.
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 4 May 1997 (sac@genome.stanford.edu)
# REVISION: $Id$
#
# A PersonManager object may have one or more Person.pm objects.
# The _people data member holds the Bio::Root::Vector.pm object.
#--------------------------------------------------------------------

package PersonManager;

use Bio::Root::Global qw(:devel);
use Bio::Root::Object ();
use Person            ();

@ISA = qw( Bio::Root::Object );

my $ID = 'PersonManager';

##
## Using superclass constructor.
##

sub destroy { 
    my $self = shift;
    $DEBUG==2 && warn "DESTROYING $self ${\$self->name}";
    $self->{'_people'}->remove_all; 
    undef $self->{'_people'};
    $self->SUPER::destroy;
}

###################################
sub add_person {

    my($self,%param) = @_;

    my $p = eval{ $self->_make_person(%param); };

    if( not $@) {
	## Expand the linked list of Person objects if necessary.
	$self->{'_people'}->add($p) if defined $self->{'_people'};
    
	## Set $self->{'_people'} to the current Person object.
	$self->{'_people'} = $p;
    } else {
	my $err = $@;
	$self->throw("Can't add person: fatal error during Person construction", $err);
    }

}

###################################
sub _make_person {
    my($self,%param) = @_;

    my $p = new Person(%param);
}

###################################
sub insert_person {
    my($self,%param) = @_;
    my ($new, $before, $after) = $self->_rearrange([qw(NEW BEFORE AFTER)], %param);

    my $p = $self->_make_person(-NAME=>$new) or return;

    my($targetName,$target);
    if($targetName = $before) {
	$target = $self->get_person($targetName) || $self->get_person('last');
	$target->insert($p,'before');

    } elsif($targetName = $after) {
	$target = $self->get_person($targetName) || $self->get_person('last');
	$target->insert($p,'after');

    } else {
	$target = $self->get_person($targetName) || $self->get_person('last');
	$target->insert($p,'after');
    }
}

###################################
sub remove_person {

    my($self,$name) = @_;
    my %param = (-RET=>'last',-UPDATE=>1);

    if($name eq 'first') {
	$self->{'_people'}->first->remove(%param);
    } elsif($name eq 'last') {
	$self->{'_people'} = $self->{'_people'}->last->remove(%param);
    } else {
	my $target = $self->get_person($name) || return;

	### Remove the desired data and update $self->{'_people'}. 
	### remove() returns the new last VectorType object.

	$self->{'_people'} = $target->remove(%param);
    }
}

###################################
sub get_person {

    my($self,$name) = @_;
    my ($target);

    if($name eq 'first') {
	$target = $self->{'_people'}->first;
    } elsif($name eq 'last') {
	$target = $self->{'_people'}->last;
    } elsif(defined $name) {
	$target = $self->{'_people'}->get($name);
	defined $target || $self->throw("Can't locate person $name.");
    }
    $target;
}

###################################
sub sort_data {
    my( $self, $by ) = @_;
    $self->{'_people'} = $self->{'_people'}->sort($by);
}

###################################
sub _display_stats {

    my( $self, $OUT ) = @_;
    
    printf( $OUT "%12s: %-12s\n", "TOTAL", $self->{'_people'}->size);
    printf( $OUT "%12s: %-12s\n", "FIRST", $self->{'_people'}->first->name);
    printf( $OUT "%12s: %-12s\n", "LAST", $self->{'_people'}->last->name);
    printf( $OUT "%12s: %-12s\n", "RANKBY", $self->{'_people'}->rank_by);
    print $OUT "\n";
    $self->{'_people'}->_display_stats($OUT);
}

######################################
1;

