# $Id$
#
# BioPerl module for Bio::LocationI
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::LocationI - Abstract interface of a Location on a Sequence

=head1 SYNOPSIS

    # get a LocationI somehow
    if( $location->strand == -1 ) {
	printf "complement(%d..%d)\n", $location->start, $location->end;
    } else {
	printf "%d..%d\n", $location->start, $location->end;
    }

=head1 DESCRIPTION


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

=head1 AUTHOR - Jason Stajich

Email jason@chg.mc.duke.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::LocationI;
use vars qw(@ISA);
use strict;

use Bio::RangeI;
use Carp;

@ISA = qw(Bio::RangeI);

# utility method Prints out a method like: 
# Abstract method stop defined in interface Bio::LocationI not
# implemented by package You::BadLocation

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  
  my $msg = "Abstract method '$caller' defined in interface Bio::LocationI but not implemented by package $package";
  if( $self->can('throw') ) {
      $self->throw($msg);
  } else {
      confess($msg);
  }
}

=head2 min_start

  Title   : min_start
  Usage   : my $minstart = $location->min_start();
  Function: Get minimum starting location of feature startpoint   
  Returns : integer or undef if no maximum starting point.
  Args    : none

=cut

sub min_start {
    my($self) = @_;
    $self->_abstractDeath();
}

=head2 max_start

  Title   : max_start
  Usage   : my $maxstart = $location->max_start();
  Function: Get maximum starting location of feature startpoint  
  Returns : integer or undef if no maximum starting point.
  Args    : none

=cut

sub max_start {
    my($self) = @_;
    $self->_abstractDeath();
}

=head2 start_pos_type

  Title   : start_pos_type
  Usage   : my $start_pos_type = $location->start_pos_type();
  Function: Get start position type encoded as text
  Returns : string ('BEFORE', 'AFTER', 'EXACT','WITHIN', 'BETWEEN')
  Args    : none

=cut

sub start_pos_type {
    my($self) = @_;
    $self->_abstractDeath();
}

=head2 min_end

  Title   : min_end
  Usage   : my $minend = $location->min_end();
  Function: Get minimum ending location of feature endpoint 
  Returns : integer or undef if no minimum ending point.
  Args    : none

=cut

sub min_end {
    my($self) = @_;
    $self->_abstractDeath();
}

=head2 max_end

  Title   : max_end
  Usage   : my $maxend = $location->max_end();
  Function: Get maximum ending location of feature endpoint 
  Returns : integer or undef if no maximum ending point.
  Args    : none

=cut

sub max_end {
    my($self) = @_;
    $self->_abstractDeath();
}

=head2 end_pos_type

  Title   : end_pos_type
  Usage   : my $end_pos_type = $location->end_pos_type();
  Function: Get end position encoded as text 
  Returns : string ('BEFORE', 'AFTER', 'EXACT','WITHIN', 'BETWEEN')
  Args    : none

=cut

sub end_pos_type {
    my($self) = @_;
    $self->_abstractDeath();
}

=head2 seq_id

  Title   : seq_id
  Usage   : my $seqid = $location->seq_id();
  Function: Get/Set seq_id that location refers to
  Returns : seq_id
  Args    : [optional] seq_id value to set

=cut

sub seq_id {
    my ($self, $seqid) = @_;
    if( defined $seqid ) {
	$self->{'_seqid'} = $seqid;
    }
    return $self->{'_seqid'};
}

=head2 to_FTstring

  Title   : to_FTstring
  Usage   : my $locstr = $location->to_FTstring()
  Function: returns the FeatureTable string of this location
  Returns : string
  Args    : none

=cut

sub to_FTstring { 
    my($self) = @_;
    $self->_abstractDeath();
}
1;

