#
# BioPerl module for Bio::Tools::AlignFactory
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Tools::AlignFactory;
use strict;

use base qw(Bio::Root::Root);

BEGIN {
    eval {
	require Bio::Ext::Align;
    };
    if ( $@ ) {
	print STDERR ("\nThe C-compiled engine for Smith Waterman alignments (Bio::Ext::Align) has not been installed.\n Please install the bioperl-ext package\n\n");
	exit(1);
    }
}

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->_initialize(@args);
  # set up defaults
  
  $self->{'kbyte'} = 20000;
  $self->{'report'} = 0;  
  return $self;
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

   &Bio::Ext::Align::change_max_BaseMatrix_kbytes($self->{'kbyte'});

   if( $self->{'report'} == 0 ) {
       &Bio::Ext::Align::error_off(16);
   } else {
       &Bio::Ext::Align::error_on(16);
   }
}

1;
