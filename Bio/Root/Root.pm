package Bio::Root::Root;
use strict;

# $Id$

=head1 NAME

Bio::Root::Root - Hash-based implementation of Bio::Root::Root

=head1 SYNOPSIS

  # any bioperl or bioperl compliant object is a RootI 
  # compliant object

  $obj->throw("This is an exception");

  eval {
      $obj->throw("This is catching an exception");
  };

  if( $@ ) {
      print "Caught exception";
  } else {
      print "no exception";
  }

=head1 DESCRIPTION

This is a hashref-based implementation of the Bio::Root::Root
interface.  Most bioperl objects should inherit from this.

See the documentation for Bio::Root::Root for most of the methods
implemented by this module.  Only overridden methods are described
here.

=head1 CONTACT

Functions originally from Steve Chervitz. Refactored by Ewan Birney.
Re-refactored by Lincoln Stein.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

use vars qw(@ISA $DEBUG $ID $Revision $VERSION $VERBOSITY);
use strict;
use Bio::Root::RootI;
use Carp 'confess';
@ISA = 'Bio::Root::RootI';

BEGIN { 
    $ID        = 'Bio::Root::Root';
    $VERSION   = 1.0;
    $Revision  = '$Id$ ';
    $DEBUG     = 0;
    $VERBOSITY = 0;
}

=head2 _create_object()

 Title   : _create_object()
 Usage   : $obj->create_object(@args)
 Function: Method which actually creates the blessed  hashref
 Returns : Blessed hashref
 Args    : Ignored

Override this method, the new() method, or _initialize() to make a
custom constructor.

=cut

sub _create_object {
  my ($class,@args) = @_;
  $class = ref($class) if ref($class);
  return bless {},$class;
}

sub verbose {
   my ($self,$value) = @_;

   if(ref($self) && (defined $value || ! defined $self->{'_root_verbose'}) ) {
       $value = 0 unless defined $value;
       $self->{'_root_verbose'} = $value;
   }
   return (ref($self) ? $self->{'_root_verbose'} : $VERBOSITY);
}

sub _register_for_cleanup {
  my ($self,$method) = @_;
  if($method) {
    if(! exists($self->{'_root_cleanup_methods'})) {
      $self->{'_root_cleanup_methods'} = [];
    }
    push(@{$self->{'_root_cleanup_methods'}},$method);
  }
}

sub _unregister_for_cleanup {
  my ($self,$method) = @_;
  my @methods = grep {$_ ne $method} $self->_cleanup_methods;
  $self->{'_root_cleanup_methods'} = \@methods;
}


sub _cleanup_methods {
  my $self = shift;
  my $methods = $self->{'_root_cleanup_methods'} or return;
  @$methods;
}

1;
