package Outer;

use Bio::Root::Global qw(:devel);
use Bio::Root::Object ();
use Foo               ();

@ISA = qw( Bio::Root::Object );

my $ID = 'Outer';

use strict;

#-------------------------------------------------------------------
# PACKAGE  : Outer.pm
# PURPOSE  : Dummy object for testing Bio::Root::Object.pm and Bio::Root::Err.pm.
# AUTHOR   : Steve Chervitz (sac@bioperl.org)
# REVISION : $Id$
#-------------------------------------------------------------------

####################################
sub _initialize {
	
    my($self, %param) = @_;
    
    $DEBUG && do{ printf ">>>> Initializing %s (%s) %s\n",$ID,ref($self),$param{-NAME}||'anon';
		  <STDIN>; };
    
    #----------------------
    # DATA MEMBERS:
    # Outerdat
    # color
    # bar
    #----------------------

    $self->SUPER::_initialize( %param );
    my ($foodat) = $self->_rearrange([qw(FOODAT)], %param);
    $self->set_foo( @$foodat );          
    $self->index();

    $DEBUG and printf "---> Initialized %s (%s) %s\n",$ID, ref($self), $self->name;
    
}

###################################
sub destroy {
    my $self = shift;
    if(ref($self->{'Foo'})) {
	$self->{'Foo'}->destroy;
	undef $self->{'Foo'};
    } 
    $self->SUPER::destroy();
}

###################################
sub set_foo {

    my ($self, %param) = @_;
    my $foo = undef;

    $DEBUG and print "$ID: set_foo()\n"; 

    $param{-PARENT} = $self;
    
    eval { $foo = new Foo( %param ); };
    if($@) {
	## If Foo throws an exception, assimilate the exception as a warning 
	## and add a note to it.
	$self->warn(-MSG  =>$@,
		    -NOTE =>"${\$self->name()} can't build Foo object: Invalid object."
		    );
    } else {
	$self->{'Foo'} = $foo;
    }
}

###################################
sub foo { my $self = shift;  $self->{'Foo'} }

###################################
sub _display_stats {

    my ($self, $OUT ) = @_;

    $self->SUPER::_display_stats($OUT);

    (defined $self->{'Foo'})
	? $self->{'Foo'}->display(-WHERE=>$OUT, -HEADER=>1) 
	: ( printf( $OUT "%-15s: %s\n", 'Foo','undefined'));

    print $OUT "\n";
}

######################################
1;

