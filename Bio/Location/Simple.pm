# $Id$
#
# BioPerl module for Bio::Location::Simple
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Location::Simple - Implementation of a Simple Location on a Sequence

=head1 SYNOPSIS

    use Bio::Location::Simple;

    my $location = new Bio::Location::Simple(-start => 1, -end => 100,
					     -strand => 1 );

    if( $location->strand == -1 ) {
	printf "complement(%d..%d)\n", $location->start, $location->end;
    } else {
	printf "%d..%d\n", $location->start, $location->end;
    }

=head1 DESCRIPTION

This is an implementation of Bio::LocationI to manage simple location
information on a Sequence.

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


package Bio::Location::Simple;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::LocationI;


@ISA = qw(Bio::Root::Root Bio::LocationI);

sub new { 
    my ($class, @args) = @_;
    my $self = {};

    bless $self,$class;

    my ($v,$start,$end,$strand,$seqid) = $self->_rearrange([qw(VERBOSE
							       START 
							    END 
							    STRAND
							    SEQID)],@args);
    defined $v && $self->verbose($v);
    defined $strand && $self->strand($strand);
    defined $start  && $self->start($start);
    defined $end    && $self->end($end);
    if( defined $self->start && defined $self->end &&
	$self->start > $self->end && $self->strand != -1 ) {
	$self->warn("When building a location start ($start) is expected to be less than end ($end), however it was not was not. Switching start and end and setting strand to -1");

	$self->strand(-1);
	my $e = $self->end;
	my $s = $self->start;
	$self->start($e);
	$self->end($s);
    }
    $seqid          && $self->seq_id($seqid);

    return $self;
}

=head2 start

  Title   : start
  Usage   : $start = $loc->start();
  Function: get/set the start of this range
  Returns : the start of this range
  Args    : optionaly allows the start to be set
          : using $loc->start($start)

=cut

sub start {
  my ($self, $value) = @_;
  $self->min_start($value) if( defined $value );
  return $self->SUPER::start();
}

=head2 end

  Title   : end
  Usage   : $end = $loc->end();
  Function: get/set the end of this range
  Returns : the end of this range
  Args    : optionaly allows the end to be set
          : using $loc->end($start)

=cut

sub end {
  my ($self, $value) = @_;

  $self->min_end($value) if( defined $value );
  return $self->SUPER::end();
}

=head2 strand

  Title   : strand
  Usage   : $strand = $loc->strand();
  Function: get/set the strand of this range
  Returns : the strandidness (-1, 0, +1)
  Args    : optionaly allows the strand to be set
          : using $loc->strand($strand)

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
 Usage   : $len = $loc->length();
 Function: get the length in the coordinate space this location spans
 Example :
 Returns : an integer
 Args    : none


=cut

sub length {
   my ($self) = @_;
   return abs($self->end() - $self->start()) + 1;
}

=head2 min_start

  Title   : min_start
  Usage   : my $minstart = $location->min_start();
  Function: Get minimum starting location of feature startpoint   
  Returns : integer or undef if no minimum starting point.
  Args    : none

=cut

sub min_start {
    my ($self,$value) = @_;

    if(defined($value)) {
	$self->{'_start'} = $value;
    }
    return $self->{'_start'};
}

=head2 max_start

  Title   : max_start
  Usage   : my $maxstart = $location->max_start();
  Function: Get maximum starting location of feature startpoint.

            In this implementation this is exactly the same as min_start().

  Returns : integer or undef if no maximum starting point.
  Args    : none

=cut

sub max_start {
    my ($self,@args) = @_;
    return $self->min_start(@args);
}

=head2 start_pos_type

  Title   : start_pos_type
  Usage   : my $start_pos_type = $location->start_pos_type();
  Function: Get start position type (ie <,>, ^).

            In this implementation this will always be 'EXACT'.

  Returns : type of position coded as text 
            ('BEFORE', 'AFTER', 'EXACT','WITHIN', 'BETWEEN')
  Args    : none

=cut

sub start_pos_type {
    my($self) = @_;
    return 'EXACT';
}

=head2 min_end

  Title   : min_end
  Usage   : my $minend = $location->min_end();
  Function: Get minimum ending location of feature endpoint 
  Returns : integer or undef if no minimum ending point.
  Args    : none

=cut

sub min_end {
    my($self,$value) = @_;

    if(defined($value)) {
	$self->{'_end'} = $value;
    }
    return $self->{'_end'};
}

=head2 max_end

  Title   : max_end
  Usage   : my $maxend = $location->max_end();
  Function: Get maximum ending location of feature endpoint 

            In this implementation this is exactly the same as min_end().

  Returns : integer or undef if no maximum ending point.
  Args    : none

=cut

sub max_end {
    my($self,@args) = @_;
    return $self->min_end(@args);
}

=head2 end_pos_type

  Title   : end_pos_type
  Usage   : my $end_pos_type = $location->end_pos_type();
  Function: Get end position type (ie <,>, ^) 

            In this implementation this will always be 'EXACT'.

  Returns : type of position coded as text 
            ('BEFORE', 'AFTER', 'EXACT','WITHIN', 'BETWEEN')
  Args    : none

=cut

sub end_pos_type {
    my($self) = @_;
    return 'EXACT';
}

=head2 location_type

  Title   : location_type
  Usage   : my $location_type = $location->location_type();
  Function: Get location type encoded as text
  Returns : string ('EXACT', 'WITHIN', 'BETWEEN')
  Args    : none

=cut

sub location_type {
    my ($self) = @_;
    return 'EXACT';
}

=head2 is_remote

 Title   : is_remote
 Usage   : $self->is_remote($newval)
 Function: Getset for is_remote value
 Returns : value of is_remote
 Args    : newvalue (optional)


=cut

sub is_remote {
   my $self = shift;
   if( @_ ) {
       my $value = shift;
       $self->{'is_remote'} = $value;
   }
   return $self->{'is_remote'};

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
    if( $self->start == $self->end ) {
	return $self->start;
    }
    my $str = $self->start . ".." . $self->end;
    if( $self->strand == -1 ) {
	$str = sprintf("complement(%s)", $str);
    }
    return $str;
}


sub trunc {
  my ($self,$start,$end,$relative_ori) = @_;

  
  my $newstart  = $self->start - $start+1;
  my $newend    = $self->end   - $start+1;
  my $newstrand = $relative_ori * $self->strand;

  my $out;
  if( $newstart < 1 || $newend > ($end-$start+1) ) {
    $out = Bio::Location::Simple->new();
    $out->start($self->start);
    $out->end($self->end);
    $out->strand($self->strand);
    $out->seq_id($self->seqid);
    $out->is_remote(1);
  } else {
    $out = Bio::Location::Simple->new();
    $out->start($newstart);
    $out->end($newend);
    $out->strand($newstrand);
    $out->seq_id();
  }

  return $out;
}

1;

