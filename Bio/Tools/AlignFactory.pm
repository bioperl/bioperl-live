
#
# BioPerl module for Bio::Tools::AlignFactory
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::AlignFactory - Base object for alignment factories

=head1 SYNOPSIS

You wont be using this as an object, but using a dervied class
like Bio::Tools::pSW

=head1 DESCRIPTION

Holds common Alignment Factory attributes in place

=head1 CONTACT

http://bio.perl.org/ or birney@sanger.ac.uk 

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::AlignFactory;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;

BEGIN {
    eval {
	require bp_sw;
    };
    if ( $@ ) {
	print STDERR ("\nThe C-compiled engine for Smith Waterman alignments (bp_sw) has not been installed.\n Please read the installation instructions for bioperl for using the compiled extensions\n\n");
	exit(1);
    }
}


@ISA = qw(Bio::Root::Object Exporter);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@p) = @_;
  my $make = $self->SUPER::_initialize(@p);
  
  # set up defaults

  $self->{'kbyte'} = 20000;
  $self->{'report'} = 0;
# set stuff in self from @args
 return $make; # success - we hope!
}


=head2 kbyte

 Title     : kbyte()
 Usage     : set/gets the amount of memory able to be used
 Function  : 
           : $factory->kbyte(200);
           :
 Returns   : 
 Argument  : memory in kilobytes

=cut

sub kbyte {
    my ($self,$value) = @_;
    
    if( defined $value ) {
	$self->{'kbyte'} = $value;
    } 

    return $self->{'kbyte'};
}


=head2 report

 Title     : report()
 Usage     : set/gets the report boolean to issue reports or not
 Function  : 
           : $factory->report(1); # reporting goes on
           :
 Returns   : n/a
 Argument  : 1 or 0

=cut

sub report {
    my ($self,$value) = @_;
    

    if( defined $value ) {
	if( $value != 1 && $value != 0 ) {
	    $self->throw("Attempting to modify AlignFactory Report with no boolean value!");
	}
	$self->{'report'} = $value;
    } 

    return $self->{'report'};
}

=head2 set_memory_and_report

 Title   : set_memory_and_report
 Usage   : Only used by subclasses.
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub set_memory_and_report{
   my ($self) = @_;

   if( $self->{'kbyte'} < 5 ) {
       $self->throw("You can suggest aligning things with less than 5kb");
   }

   &bp_sw::change_max_BaseMatrix_kbytes($self->{'kbyte'});

   if( $self->{'report'} == 0 ) {
       &bp_sw::error_off(16);
   } else {
       &bp_sw::error_on(16);
   }
}



