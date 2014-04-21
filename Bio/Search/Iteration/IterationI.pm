#-----------------------------------------------------------------
#
# BioPerl module Bio::Search::Iteration::IterationI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Iteration::IterationI - Abstract interface to an
iteration from an iterated search result, such as PSI-BLAST.

=head1 SYNOPSIS

    # Bio::Search::Iteration::IterationI objects cannot be 
    # instantiated since this module defines a pure interface.
    # Given an object that implements the 
    # Bio::Search::Iteration::IterationI interface, 
    # you can do the following things with it:

    # First, open up a SearchIO stream
    use Bio::SearchIO;
    my $file = shift or die "Usage: $0 <BLAST-report-file>\n";
    my $in = Bio::SearchIO->new(-format => 'blast',
                               -file => $file # comment out this line to read STDIN
                              );
    # Iterate over all results in the input stream
    while (my $result = $in->next_result) {

        printf "Result #%d: %s\n", $in->result_count, $result->to_string;
        printf "Total Iterations: %d\n", $result->num_iterations();

        # Iterate over all iterations and process old and new hits
        # separately.

        while( my $it = $result->next_iteration) { 
            printf "\nIteration %d\n", $it->number;
            printf "Converged: %d\n", $it->converged;

            # Print out the hits not found in previous iteration
            printf "New hits: %d\n", $it->num_hits_new;
            while( my $hit = $it->next_hit_new ) {
                printf "  %s, Expect=%g\n", $hit->name, $hit->expect; 
            }

            # Print out the hits found in previous iteration
            printf "Old hits: %d\n", $it->num_hits_old; 
            while( my $hit = $it->next_hit_old ) {
                printf "  %s, Expect=%g\n", $hit->name, $hit->expect; 
            }
        }
        printf "%s\n\n", '-' x 50;
    }

    printf "Total Reports processed: %d: %s\n", $in->result_count;

    __END__

    # NOTE: The following functionality is just proposed
    # (does not yet exist but might, given sufficient hew and cry):

    # Zero-in on the new hits found in last iteration.
    # By default, iteration() returns the last one.

    my $last_iteration = $result->iteration();
    while( my $hit = $last_iteration->next_hit) {
        # Do something with new hit...
    }

    # Get the first iteration

    my $first_iteration = $result->iteration(1);


=head1 DESCRIPTION

Bio::Search::Result::ResultI objects are data structures containing
the results from the execution of a search algorithm.  As such, it may
contain various algorithm specific information as well as details of
the execution, but will contain a few fundamental elements, including
the ability to return Bio::Search::Hit::HitI objects.

=head2 Classification of Hits

Within a given iteration, the hits can be classified into a number of
useful subsets based on whether or not the hit appeard in a previous
iteration and whether or not the hit is below the threshold E-value
for inclusion in the score matrix model.

                           All hits
                             (A)
               _______________|_________________
               |                               |
            New hits                        Old hits
              (B)                             (C)
      _________|________                _______|_________
      |                |                |               |
    Below            Above             Below          Above
  threshold        threshold         threshold      threshold
     (D)              (E)              (F)             (G)
                               _________|___________
                               |                   |
                         Occurred in a         Occurred in a
                         previous iteration    previous iteration
                         below threshold       above threshold
                              (H)                  (I)

Notes: The term I<threshold> in the diagram and descriptions below
refer to this inclusion threshold. I<Below threshold> actually means
I<at or below threshold>.

The IterationI interface defines a number of methods for extracting
these subsets of hits.

=over 4

=item * newhits_below_threshold() [subset D]

Hits that did not appear in a previous iteration and are below
threshold in the current iteration.

=item * newhits_not_below_threshold() [subset E]

Hits that did not appear in a previous iteration and are not below
threshold in the current iteration.

=item * newhits() [subset B]

All newly found hits, below and above the inclusion threshold.  This
is the union of newhits_below_threshold() + newhits_not_below_threshold()
[subset D + subset E].

=item * oldhits_below_threshold() [subset H]

Hits that appeared in a previous iteration below threshold and are
still below threshold in the current iteration.

=item * oldhits_newly_below_threshold() [subset I]

Hits that appeared in a previous iteration above threshold but are
below threshold in the current iteration. (Not applicable to the first
iteration.)

=item * oldhits_not_below_threshold() [subset G]

Hits that appeared in a previous iteration not below threshold and
are still not below threshold in the current iteration.

=item * oldhits()  [subset C]

All hits that occured in a previous iteration, whether below or above
threshold in the current iteration. Union of oldhits_below_threshold()
+ oldhits_newly_below_threshold() + oldhits_not_below_threshold()
[subset H + subset I + subset G]. (Not applicable to the first
iteration.)

=item * hits_below_threshold() [subset D + subset F]

All hits, old and new, that are below the inclusion threshold in this
iteration. This is the union of newhits_below_threshold() +
oldhits_below_threshold() + oldhits_newly_below_threshold()
[subset D + subset H + subset I].

=item * hits() [subset A]

The union of newhits() and oldhits() [subset B + subset C].

=back

For the first iteration, the methods L<oldhits>, L<oldhits_below_threshold>,
L<oldhits_newly_below_threshold>, and oldhits_not_below_threshold()
will return empty lists.

Iterator and numbers-of-hit methods are provided for subsets A, B, and C:

=over 4

=item * next_hit_new(), num_hits_new() [subset B]

=item * next_hit_old(), num_hits_old() [subset C]

=item * next_hit(), num_hits() [subset A]

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR 

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports and comments.

=head1 COPYRIGHT

Copyright (c) 2003 Steve Chervitz. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...


package Bio::Search::Iteration::IterationI;

use strict;


use base qw(Bio::Root::RootI);

=head2 number

 Title   : number
 Usage   : $it_number = $iteration->number();
 Purpose : returns the number of the iteration (a.k.a "round") 
           within the Result.
 Returns : integer
 Args    : [optional] integer to set the number of the iteration

=cut

sub number {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

=head2 converged

 Title   : converged
 Usage   : $it_converged = $iteration->converged();
 Purpose : Indicates whether or not the iteration has converged 
 Returns : boolean 
 Args    : [optional] boolean value to set the converged of the iteration

=cut

sub converged {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

=head2 next_hit

 Title   : next_hit
 Usage   : while( $hit = $iteration->next_hit( [$found_again]) ) { ... }
 Purpose : Iterates through all of the HitI objects
           including new hits and old hits found in a previous iteration
           and both below and above the inclusion threshold.
           Corresponds to subset A in the "Classification of Hits"
           documentation section of this module.
 Returns : A Bio::Search::Hit::HitI object or undef if there are no more.
           Hits will be returned in the order in which they occur in the report
           unless otherwise specified.
 Args    : none

See Also: L<hits>, L<Classification of Hits>

next_hit() iterates through all hits, including the new ones
for this iteration and those found in previous iterations.
You can interrogate each hit using L<Bio::Search::Hit::HitI::found_again>
to determine whether it is new or old.

To get just the new hits, use L<next_hit_new>.
To get just the old hits, use L<next_hit_old>.

=cut

sub next_hit {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

=head2 next_hit_new

 Title   : next_hit_new
 Usage   : while( $hit = $iteration->next_hit_new() ) { ... }
 Purpose : Iterates through all newly found hits (did not occur in a
           previous iteration) and are either below or above the inclusion threshold.
           Corresponds to subset B in the "Classification of Hits"
           documentation section of this module.
 Returns : A Bio::Search::Hit::HitI object or undef if there are no more.
           Hits will be returned in the order in which they occur in the report
           unless otherwise specified.
 Args    : none

See Also: L<next_hit>, L<next_hit_old>, L<newhits>, L<Classification of Hits>

=cut

sub next_hit_new {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

=head2 next_hit_old

 Title   : next_hit_old
 Usage   : while( $hit = $iteration->next_hit_old() ) { ... }
 Purpose : Iterates through the Hit objects representing just the
           hits that have been found in a previous iteration, whether
           below or above the inclusion threshold.
           Corresponds to subset C in the "Classification of Hits"
           documentation section of this module.
 Returns : A Bio::Search::Hit::HitI object or undef if there are no more.
           Hits will be returned in the order in which they occur in the report
           unless otherwise specified.
 Args    : none

See Also: L<next_hit>, L<next_hit_old>, L<oldhits>, L<Classification of Hits>

=cut

sub next_hit_old {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

=head2 num_hits

 Title   : num_hits
 Usage   : my $hitcount_total = $iteration->num_hits
 Purpose : Returns the total number of hits for this query result, including new and old
           below and above inclusion threshold.
 Returns : integer
 Args    : none

See Also: L<num_hits_new>, L<num_hits_old>, L<Classification of Hits>

=cut

sub num_hits {
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 num_hits_new

 Title   : num_hits_new
 Usage   : my $hitcount_new = $result->num_hits_new;
         : my $hitcount_new_below_thresh = $result->num_hits_new( 1 );
 Purpose : Returns the number of new hits in this iteration that were not
           found in a previous iteration and are either below or above the
           the inclusion threshold.
           Corresponds to subset B in the "Classification of Hits"
           documentation section of this module.
 Returns : integer
 Args    : (optional) boolean, true if you want to get a count of just the new hits
           that are below the inclusion threshold.


See Also: L<num_hits>, L<num_hits_old>, L<Classification of Hits>

=cut

sub num_hits_new {
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 num_hits_old

 Title   : num_hits_old
 Usage   : my $hitcount_old = $result->num_hits_old;
         : my $hitcount_old_below_thresh = $result->num_hits_old( 1 );
 Purpose : Returns the number of new hits in this iteration that were
           found in a previous iteration and are either below or above the
           the inclusion threshold.
           Corresponds to subset C in the "Classification of Hits"
           documentation section of this module.
 Returns : integer
 Args    : (optional) boolean, true if you want to get a count of just the old hits
           that are below the inclusion threshold.

See Also: L<num_hits>, L<num_hits_new>, L<Classification of Hits>

=cut

sub num_hits_old {
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 hits

 Title    : hits
 Usage    : foreach( $obj->hits() ) { ... };
 Purpose  : Provides access to all hits, both new and old, and either
            below or above the inclusion threshold.
            Corresponds to subset A in the "Classification of Hits"
            documentation section of this module.
 Returns  : An array containing all HitI objects.
            Hits will be ordered according to their occurrence in the report
            unless otherwise specified.
 Args     : none

See Also: L<newhits>, L<oldhits>, L<Classification of Hits>

=cut

sub hits  { shift->throw_not_implemented(); }

=head2 newhits

 Title    : newhits
 Usage    : foreach( $obj->newhits() ) { ... };
 Purpose  : Provides access to hits that were not found in a previous iteration
            and may be either below or above the inclusion threshold.
            Corresponds to subset B in the "Classification of Hits"
            documentation section of this module.
 Returns  : An array containing Bio::Search::Hit::HitI objects.
            Hits will be ordered according to their occurrence in the report
            unless otherwise specified.
 Args     : none

See Also: L<hits>, L<oldhits>, L<newhits_below_threshold> + L<newhits_not_below_threshold>, L<Classification of Hits>

=cut

sub newhits  { shift->throw_not_implemented(); }

=head2 oldhits

 Title    : oldhits
 Usage    : foreach( $obj->oldhits() ) { ... };
 Purpose  : Provides access to hits that were found in a previous iteration
            and are either below or above the inclusion threshold in the current iteration.
            Corresponds to subset C in the "Classification of Hits"
            documentation section of this module.
 Returns  : An array containing Bio::Search::Hit::HitI objects.
            Hits will be ordered according to their occurrence in the report
            unless otherwise specified.
 Args     : none

See Also: L<hits>, L<newhits>, L<oldhits_below_threshold>, L<oldhits_newly_below_threshold>, L<oldhits_not_below_threshold>, L<Classification of Hits>

=cut

sub oldhits  { shift->throw_not_implemented(); }

=head2 newhits_below_threshold

 Title   : newhits_below_threshold
 Usage   : foreach( $obj->newhits_below_threshold() ) { ... };
 Purpose : Provides access to hits that did not appear in a 
           previous iteration and are below threshold.
           Corresponds to subset D in the "Classification of Hits"
           documentation section of this module.
 Returns : An array containing Bio::Search::Hit::HitI objects.
           Hits will be returned in the order in which they occur in the report
           unless otherwise specified.
 Args    : none

See Also: L<newhits_not_below_threshold>, L<oldhits_newly_below_threshold>, L<newhits>, L<Classification of Hits>

=cut

sub newhits_below_threshold  { shift->throw_not_implemented(); }

=head2 oldhits_below_threshold

 Title   : oldhits_below_threshold
 Usage   : foreach( $obj->oldhits_below_threshold() ) { ... };
 Purpose : Provides access to hits that appeared in a 
           previous iteration below inclusion threshold and are still below threshold.
           Corresponds to subset H in the "Classification of Hits"
           documentation section of this module.
 Returns : An array containing Bio::Search::Hit::HitI objects.
           Hits will be returned in the order in which they occur in the report
           unless otherwise specified.
 Args    : none

See Also: L<oldhits_not_below_threshold>, L<oldhits_newly_below_threshold>, L<oldhits>, L<Classification of Hits>

=cut

sub oldhits_below_threshold  { shift->throw_not_implemented(); }

=head2 oldhits_newly_below_threshold

 Title   : oldhits_newly_below_threshold
 Usage   : foreach( $obj->oldhits_newly_below_threshold() ) { ... };
 Purpose : Provides access to hits that appeared in a previous
           iteration above threshold but are below threshold in the 
           current iteration. Not applicable to the first iteration.
           Corresponds to subset I in the "Classification of Hits"
           documentation section of this module.
 Returns : An array containing Bio::Search::Hit::HitI objects.
           Hits will be returned in the order in which they occur in the report
           unless otherwise specified.
 Args    : none

See Also: L<newhits_below_threshold>, L<oldhits>, L<Classification of Hits>

=cut

sub oldhits_newly_below_threshold  { shift->throw_not_implemented(); }

=head2 oldhits_not_below_threshold

 Title   : oldhits_not_below_threshold
 Usage   : foreach( $obj->oldhits_not_below_threshold() ) { ... };
 Purpose : Provides access to hits that appeared in a previous iteration
           not below threshold and are still not below threshold.
           Corresponds to subset G in the "Classification of Hits"
           documentation section of this module.
 Returns : An array containing Bio::Search::Hit::HitI objects.
           Hits will be returned in the order in which they occur in the report
           unless otherwise specified.
 Args    : none

See Also: L<oldhits_below_threshold>, L<oldhits>, L<Classification of Hits>

=cut

sub oldhits_not_below_threshold  { shift->throw_not_implemented(); }

=head2 newhits_not_below_threshold

 Title   : newhits_not_below_threshold
 Usage   : foreach( $obj->newhits_not_below_threshold() ) { ... };
 Purpose : Provides access to hits that did not appear in a 
           previous iteration and are not below threshold 
           in the current iteration.
           Corresponds to subset E in the "Classification of Hits"
           documentation section of this module.
 Returns : An array containing Bio::Search::Hit::HitI objects.
           Hits will be returned in the order in which they occur in the report
           unless otherwise specified.
 Args    : none

See Also: L<newhits_below_threshold>, L<newhits>, L<Classification of Hits>

=cut

sub newhits_not_below_threshold  { shift->throw_not_implemented(); }

=head2 hits_below_threshold

 Title   : hits_below_threshold
 Usage   : foreach( $obj->hits_below_threshold() ) { ... };
 Purpose : Provides access to all hits, old and new, that are below the inclusion threshold.
           Corresponds to the union of subset D and subset F in the 
           "Classification of Hits" documentation section of this module.
 Returns : An array containing Bio::Search::Hit::HitI objects.
           Hits will be returned in the order in which they occur in the report
           unless otherwise specified.
 Args    : none

See Also: L<newhits_below_threshold>, L<oldhits_newly_below_threshold>, L<oldhits_below_threshold>, L<Classification of Hits>

=cut

sub hits_below_threshold  { shift->throw_not_implemented(); }


=head2 add_hit

 Title   : add_hit
 Usage   : $report->add_hit(-hit             =>$hit_obj,
                            -old             =>$boolean,
                            -below_threshold =>$boolean,
                            -newly_below     =>$boolean )
 Purpose : Adds a HitI to the stored list of hits
 Returns : Number of HitI currently stored for the class of the added hit.
 Args    : Tagged values, the only required one is -hit. All others are used
           only for PSI-BLAST reports.
           -hit => Bio::Search::Hit::HitI object
           -old => boolean, true indicates that the hit was found 
                   in a previous iteration. Default=false.
           -below_threshold => boolean, true indicates that the hit is below
                   the inclusion threshold.
           -newly_below => boolean, true indicates that the hit is below
                   the inclusion threshold in this iteration but was above
                   the inclusion threshold in a previous iteration. 
                   Only appropriate for old hits. Default=false.
 Throws  : Bio::Root::BadParameter if the hit is not a
           Bio::Search::Hit::HitI.
           Bio::Root::BadParameter if -old=>false and -newly_below=>true.

=cut

sub add_hit { shift->throw_not_implemented }



=head2 get_hit

 Title   : get_hit
 Usage   : $hit = $report->get_hit( $hit_name )
 Purpose : Gets a HitI object given its name 
           if a hit with this name exists within this Iteration.
 Returns : Bio::Search::Hit::HitI object or undef if there is no such hit.
 Args    : $hit_name = string containing name of the hit
 Throws  : n/a

The name string must be the same as that returned by
Bio::Search::Hit::HitI::name().

The lookup should be case-insensitive.

=cut

sub get_hit { shift->throw_not_implemented }


1;


