package Bar;

use Bio::Root::Global  qw(:devel);
use Bio::Root::Object   ();

@ISA = qw( Bio::Root::Object );

my $ID = 'Bar';

use strict;

#----------------------------------------------------------------------
# PACKAGE  : Bar.pm
# PURPOSE  : Dummy object for testing Bio::Root::Object.pm and Bio::Root::Err.pm.
# AUTHOR   : Steve A. Chervitz (sac@genome.stanford.edu)
# REVISION : $Id$
#----------------------------------------------------------------------

####################################
sub _initialize {

    my( $self, %param) = @_;
    
    $DEBUG and do{ print ">>>> Initializing $ID (${\ref($self)}) ",$param{-NAME}||'anon';<STDIN>};
    
    #----------------------
    # DATA MEMBERS:
    # bardat
    # flavor
    #----------------------
    $self->SUPER::_initialize( %param );
    my ($bar, $flavor) = $self->_rearrange([qw(BAR FLAVOR)], %param);

    $self->_set_data( $bar );
    $self->_set_flavor( $flavor );
    $self->index();

    $DEBUG and print "---> Initialized $ID (${\ref($self)}) ",$self->name(),"\n";
}

###################################
sub _set_data {
    
    my ($self, $b ) = @_;	
    
    $b =~ /^\d+$/ or $self->throw("Bar data not defined or not an integer: $b");
    $self->{'bardat'} = $b;

}

###################################
sub _set_flavor {
    
    my ($self, $f ) = @_;	
    
    $DEBUG and do{ print "Bar::_set_flavor() $f"; <STDIN>; };

    if($self->err) {
	print "$ID: has an error:\n";
	$self->print_err; <STDIN>;
    }

    defined $f || $self->throw("Flavor not defined.",
			       "All $ID objects must have a flavor.");

    $self->{'flavor'} = $f;
}

###################################
sub flavor { my $self = shift; defined $self->{'flavor'} ? $self->{'flavor'} : 'unknown'; }
sub data { my $self = shift;  $self->{'bardat'} }

###################################
sub _display_stats {

    my ($self, $OUT ) = @_;

    $self->SUPER::_display_stats($OUT);

    printf( $OUT "%-15s: %d\n", 'DATA', $self->data());
    printf( $OUT "%-15s: %s\n", 'FLAVOR',$self->flavor());

    print $OUT "\n";
}

######################################
1;

