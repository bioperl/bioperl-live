package Foo;

use Bio::Root::Global qw(:devel);
use Bio::Root::Object ();
use Bar               ();

@ISA = qw( Bio::Root::Object );

my $ID            = 'Foo';
my $DEFAULT_COLOR = 'gray';
my @COLORS        = qw('red' 'green' 'blue' 'gray');

use strict;

#----------------------------------------------------------------------
# PACKAGE  : Foo.pm
# PURPOSE  : Dummy object for testing Bio::Root::Object.pm and Bio::Root::Err.pm.
# AUTHOR   : Steve Chervitz (sac@bioperl.org)
# REVISION : $Id$
#----------------------------------------------------------------------

####################################
sub _initialize {
	
    my($self, %param) = @_;
    
    $DEBUG && do{ printf ">>>> Initializing %s (%s) %s\n",$ID,ref($self),$param{-NAME}||'anon';
		  <STDIN>; };
    
    #----------------------
    # DATA MEMBERS:
    # foodat
    # color
    # bar
    #----------------------

    $self->SUPER::_initialize( %param );
    my ($foo, $color) = $self->_rearrange([qw(FOO COLOR)], %param);

    $self->set_data( $foo);    
    $self->set_color( $color);
    $self->set_bar( %param) if defined $param{-BAR};          
#    $self->index();

    $DEBUG and printf "---> Initialized %s (%s) %s\n",$ID, ref($self), $self->name;

}

###################################
sub destroy {
    my $self = shift;
    if(ref($self->{'Bar'})) {
	$self->{'Bar'}->destroy;
	undef $self->{'Bar'};
    } 
    $self->SUPER::destroy();
}

###################################
sub set_data {

    my ($self, $d) = @_;	

    if(not defined ($d) or $d =~ /\D/) {
	$self->throw("Undefined or invalid foo data: $d",
		     "All $ID objects must be initialized with numeric data.");
    }
    $self->{'Foodat'} = $d;
}

###################################
sub set_color {
    my ($self, $c) = @_;	

    $DEBUG and  print "$ID: setting color()"; 

    $c || do{ $self->warn("Color not defined.",
			  "All $ID objects must have a color.".
			  "\tSetting to $DEFAULT_COLOR");
	      $c = $DEFAULT_COLOR; };

    $self->{'Color'} = "\L$c\E";
}

###################################
sub change_color {
    my ($self,$c) = @_;
    
    if( grep /$c/i, @COLORS) {
	$self->set_color($c);
    } else {
	$self->throw("Invalid or unspecified color = $c", 
		     "Color not changed ($self->{'Color'}).\nAcceptable colors: @COLORS");
    }
}

###################################
sub set_bar {

    my ($self, %param) = @_;
    my $bar = undef;

    $DEBUG and print "Foo::set_bar()\n"; 
    
    $param{-PARENT} = $self;
    my $bar_name = 'BAR_'.$self->{'Foodat'};
    $param{-NAME} = $bar_name;

    eval { $bar = new Bar(%param); };

    if($@) {
	## If Bar throws an exception, assimilate the exception as a warning 
	## and add a note to it.
	$self->warn(-MSG  =>$@,
		    -NOTE =>"${\$self->name()} can't build Bar object \"$bar_name\": Invalid object.");
    } else {
	$self->{'Bar'} = $bar;
    }
}

###################################
sub color { my $self = shift; defined $self->{'Color'} ? $self->{'Color'} : 'unknown'; }
sub data { my $self = shift;  $self->{'Foodat'} }
sub bar { my $self = shift;  $self->{'Bar'} }

###################################
sub _display_stats {

    my ($self, $OUT ) = @_;

    $self->SUPER::_display_stats($OUT);

    printf( $OUT "%-15s: %d\n", 'DATA',$self->data());
    printf( $OUT "%-15s: %s\n", 'COLOR',$self->color());

    (defined $self->{'Bar'})
	? $self->{'Bar'}->display(-WHERE=>$OUT, -HEADER=>1) 
	: ( printf( $OUT "%-15s: %s\n", 'BAR','undefined'));

    print $OUT "\n";
}

######################################
1;

