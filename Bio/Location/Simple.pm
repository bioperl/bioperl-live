# $Id $
#
# BioPerl module for Bio::Location::Simple
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Location::Simple - Abstract interface of a Location on Sequence

=head1 SYNOPSIS

# get a Location::Simple somehow

=head1 DESCRIPTION

Descript to follow

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

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::Location::Simple;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::LocationI;


@ISA = qw(Bio::Root::RootI Bio::LocationI);

sub new { 
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($start,$end,$strand) = $self->_rearrange([qw(START 
						     END 
						     STRAND)],@args);
    $start && $self->start($start);
    $end && $self->end($end);
    $strand && $self->strand($strand);

    return $self;
}

=head2

  Title   : start
  Usage   : $start = $range->start();
  Function: get/set the start of this range
  Returns : the start of this range
  Args    : optionaly allows the start to be set
          : using $range->start($start)

=cut

sub start {
  my ($self, $value) = @_;
  if( defined $value || !defined $self->{'_start'} ) {
      $value = 0 unless defined ( $value );
      $self->{'_start'} = $value;
  }
  return $self->{'_start'};
}

=head2

  Title   : end
  Usage   : $end = $range->end();
  Function: get/set the end of this range
  Returns : the end of this range
  Args    : optionaly allows the end to be set
          : using $range->end($start)

=cut

sub end {
  my ($self, $value) = @_;
  if( defined $value || !defined $self->{'_end'}) {
      $value = 0 unless defined ( $value );
      $self->{'_end'} = $value;
  }
  return $self->{'_end'};
}

=head2

  Title   : strand
  Usage   : $strand = $range->strand();
  Function: get/set the strand of this range
  Returns : the strandidness (-1, 0, +1)
  Args    : optionaly allows the strand to be set
          : using $range->strand($strand)

=cut

sub strand {
  my ($self, $value) = @_;

  if ( defined $value ||  
       ! defined $self->{'_strand'} ) {
      # let's go ahead and force to '0' if
      # we are requesting the strand without it
      # having been set previously

       $value = 0 unless defined($value);

       if ( $value eq '+' ) { $value = 1; }
       elsif ( $value eq '-' ) { $value = -1; }
       elsif ( $value eq '.' ) { $value = 0; }
       elsif ( $value != -1 && $value != 1 && $value != 0 ) {
	   $self->throw("$value is not a valid strand info");
       }
       $self->{'_strand'} = $value
   }
   return $self->{'_strand'};
}

=head2 length

 Title   : length
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub length {
   my ($self) = @_;
   return $self->end() - $self->start() + 1;
}

1;

