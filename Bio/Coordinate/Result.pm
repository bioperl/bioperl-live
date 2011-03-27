#
# bioperl module for Bio::Coordinate::Result
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Coordinate::Result - Results from coordinate transformation

=head1 SYNOPSIS

  use Bio::Coordinate::Result;

  #get results from a Bio::Coordinate::MapperI
  $matched = $result->each_match;

=head1 DESCRIPTION

The results from Bio::Coordinate::MapperI are kept in an object which
itself is a split location, See L<Bio::Location::Split>. The results
are either Matches or Gaps.  See L<Bio::Coordinate::Result::Match> and
L<Bio::Coordinate::Result::Gap>.

If only one Match is returned, there is a convenience method of
retrieving it or accessing its methods. Same holds true for a Gap.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Coordinate::Result;
use strict;


use base qw(Bio::Location::Split Bio::Coordinate::ResultI);


=head2 add_location

 Title   : add_sub_Location
 Usage   : $obj->add_sub_Location($variant)
 Function: 

           Pushes one Bio::LocationI into the list of variants.

 Example : 
 Returns : 1 when succeeds
 Args    : Location object

=cut

sub add_sub_Location {
  my ($self,$value) = @_;
  if( ! $value ) {
      $self->warn("provding an empty value for location\n");
      return;
  }
  $self->throw("Is not a Bio::LocationI but [$value]")
      unless $value->isa('Bio::LocationI');

  $self->{'_match'} = $value
      if $value->isa('Bio::Coordinate::Result::Match');

  $self->{'_gap'} = $value
      if $value->isa('Bio::Coordinate::Result::Gap');

  $self->SUPER::add_sub_Location($value);

}

=head2 add_result

 Title   : add_result
 Usage   : $obj->add_result($result)
 Function: Adds the contents of one Bio::Coordinate::Result
 Example : 
 Returns : 1 when succeeds
 Args    : Result object

=cut

sub add_result {
  my ($self,$value) = @_;

  $self->throw("Is not a Bio::Coordinate::Result but [$value]")
      unless $value->isa('Bio::Coordinate::Result');

  map { $self->add_sub_Location($_) } $value->each_Location;
}

=head2 seq_id

  Title   : seq_id
  Usage   : my $seqid = $location->seq_id();
  Function: Get/Set seq_id that location refers to

            We override this here in order to propagate to all sublocations
            which are not remote (provided this root is not remote either)

  Returns : seq_id
  Args    : [optional] seq_id value to set


=cut

sub seq_id {
    my ($self, $seqid) = @_;

    my @ls = $self->each_Location;
    if (@ls) {
	return $ls[0]->seq_id;
    } else {
	return;
    }
}


=head2 Convenience methods

These methods are shortcuts to Match and Gap locations.

=cut

=head2 each_gap

 Title   : each_gap
 Usage   : $obj->each_gap();
 Function: 

            Returns a list of Bio::Coordianate::Result::Gap objects.

 Returns : list of gaps
 Args    : none

=cut

sub each_gap {
   my ($self) = @_;

   my @gaps;
   foreach my $gap ($self->each_Location) {
       push @gaps, $gap if $gap->isa('Bio::Coordinate::Result::Gap');
   }
   return @gaps;

}


=head2 each_match

 Title   : each_match
 Usage   : $obj->each_match();
 Function: 

            Returns a list of Bio::Coordinate::Result::Match objects.

 Returns : list of Matchs
 Args    : none

=cut

sub each_match {
   my ($self) = @_;

   my @matches;
   foreach my $match ($self->each_Location) {
       push @matches, $match if $match->isa('Bio::Coordinate::Result::Match');
   }
   return @matches;
}

=head2 match

 Title   : match
 Usage   : $match_object = $obj->match(); #or
           $gstart = $obj->gap->start;
 Function: Read only method for retrieving or accessing the match object.
 Returns : one Bio::Coordinate::Result::Match
 Args    : 

=cut

sub match {
   my ($self) = @_;

   $self->warn("More than one match in results")
       if $self->each_match > 1 and $self->verbose > 0;
   unless (defined $self->{'_match'} ) {
       my @m = $self->each_match;
       $self->{'_match'} = $m[-1];
   }
   return $self->{'_match'};
}

=head2 gap

 Title   : gap
 Usage   : $gap_object = $obj->gap(); #or
           $gstart = $obj->gap->start;
 Function: Read only method for retrieving or accessing the gap object.
 Returns : one Bio::Coordinate::Result::Gap
 Args    : 

=cut

sub gap {
   my ($self) = @_;

   $self->warn("More than one gap in results")
       if $self->each_gap > 1 and $self->verbose > 0;
   unless (defined $self->{'_gap'} ) {
       my @m = $self->each_gap;
       $self->{'_gap'} = $m[-1];
   }
   return $self->{'_gap'};
}


=head2 purge_gaps

 Title   : purge_gaps
 Usage   : $gap_count = $obj->purge_gaps;
 Function: remove all gaps from the Result
 Returns : count of removed gaps
 Args    : 

=cut

sub purge_gaps {
    my ($self) = @_;
    my @matches;
    my $count = 0;

    foreach my $loc ($self->each_Location) {
        if ($loc->isa('Bio::Coordinate::Result::Match')) {
            push @matches, $loc;
        } else {
            $count++
        }
    }
    @{$self->{'_sublocations'}} = ();
    delete $self->{'_gap'} ;
    push @{$self->{'_sublocations'}}, @matches;
    return $count;
}


1;
