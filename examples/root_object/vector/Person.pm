#--------------------------------------------------------
# PACKAGE : Person.pm
# PURPOSE : Dummy object for testing Bio::Root::Vector.pm.
# AUTHOR  : Steve Chervitz (sac@bioperl.org)
# CREATED : 4 May 1997 (sac@bioperl.org)
# REVISION: $Id$
#--------------------------------------------------------

package Person;

use Bio::Root::Global qw(:devel);
use Bio::Root::Object ();
use Bio::Root::Vector ();

@ISA = qw( Bio::Root::Object Bio::Root::Vector);

my $ID = 'Person';


####################################
sub _initialize {
	
    my($self, %param) = @_;
    
    $DEBUG and do{ print ">>>> Initializing $ID (${\ref($self)}) ",$param{-NAME}||'anon';<STDIN>};
    
    $self->SUPER::_initialize(%param);
    $self->has_name() || $self->throw("Can't create anonymous person.");

    $DEBUG and print "---> Initialized $ID (${\ref($self)}) ",$self->name(),"\n";
}


sub _display_stats {
    # Prints most recent Person first.
    
    my( $self, $OUT ) = @_;
    my $object = $self->first();
    my $numerate = $object->size() > 1;

    printf ($OUT "  %-4s %-10s %-10s %-10s\n", 'Rank','Name','Prev','Next');
    printf ($OUT "  %-4s %-10s %-10s %-10s\n", '----','----','----','----');
    do{ 
	if($object->is_first()) {
	    if($object->is_last()) {
		printf $OUT ("  %-4d %-10s %-10s %-10s\n", $object->rank(), $object->name(), 'First', 'Last');
	    } else {
		printf $OUT ("  %-4d %-10s %-10s %-10s\n", 
			     $object->rank(), $object->name(), 'First', $object->next()->name());
	    }
	} elsif($object->is_last()) {
	    printf $OUT ("  %-4d %-10s %-10s %-10s\n", 
		    $object->rank(), $object->name(), $object->prev()->name(), 'Last');
	} else {
	    printf $OUT ("  %-4d %-10s %-10s %-10s\n", 
		    $object->rank(), $object->name(), $object->prev()->name(),$object->next()->name());
	}

    } while($object = $object->next());
}

######################################
1;

