#
# BioPerl module for Bio::Search::Result::BlastResult
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# Copyright Steve Chervitz
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Result::BlastResult - Blast-specific subclass of Bio::Search::Result::GenericResult

=head1 SYNOPSIS

    # Working with iterations (PSI-BLAST results)

    $result->next_iteration();
    $result->num_iterations();
    $result->iteration();
    $result->iterations();

# See Bio::Search::Result::GenericResult for information about working with Results.

# See L<Bio::Search::Iteration::IterationI|Bio::Search::Iteration::IterationI>
# for details about working with iterations.

# TODO:
#     * Show how to configure a SearchIO stream so that it generates
#       BlastResult objects.


=head1 DESCRIPTION

This object is a subclass of Bio::Search::Result::GenericResult
and provides some operations that facilitate working with BLAST
and PSI-BLAST results.

For general information about working with Results, see 
Bio::Search::Result::GenericResult.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Steve Chervitz

Email sac@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::Result::BlastResult;
use strict;

use Bio::Search::BlastStatistics;

use base qw(Bio::Search::Result::GenericResult);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::Result::BlastResult->new();
 Function: Builds a new Bio::Search::Result::BlastResult object
 Returns : Bio::Search::Result::BlastResult
 Args    : See Bio::Search::Result::GenericResult();
           The following parameters are specific to BlastResult:
             -iterations  => array ref of Bio::Search::Iteration::IterationI objects
             -inclusion_threshold => e-value threshold for inclusion in the
                                     PSI-BLAST score matrix model (blastpgp)

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  $self->{'_iterations'} = [];
  $self->{'_iteration_index'} = 0;
  $self->{'_iteration_count'} = 0;

  my( $iters, $ithresh ) = $self->_rearrange([qw(ITERATIONS
                                                 INCLUSION_THRESHOLD)],@args);

  $self->{'_inclusion_threshold'} = $ithresh;  # This is a read-only variable

  if( defined $iters  ) {
      $self->throw("Must define arrayref of Iterations when initializing a $class\n") unless ref($iters) =~ /array/i;

      foreach my $i ( @{$iters} ) {
          $self->add_iteration($i);
      }
  } 
  else {
      # This shouldn't get called with the new SearchIO::blast.
      print STDERR "BlastResult::new(): Not adding iterations.\n";
      $self->{'_no_iterations'} = 1;
  }

  return $self;
}


=head2 hits

 Title   : hits
 Usage   : my @hits = $result->hits
 Function: Returns the available hits for this Result
 Returns : Array of L<Bio::Search::Hit::HitI> objects
 Args    : none
 Note    : This method overrides L<Bio::Search::Result::GenericResult::hits> to
           take into account the possibility of multiple iterations, as occurs
           in PSI-BLAST reports.
           If there are multiple iterations, all 'new' hits for all iterations
           are returned. These are the hits that did not occur in a previous
           iteration.

           See Also: L<Bio::Search::Result::GenericResult::hits>

=cut

sub hits {
   my ($self) = shift;
   if ($self->{'_no_iterations'}) {
       return $self->SUPER::hits;
   }
   my @hits = ();
   foreach my $it ($self->iterations) {
       push @hits, $it->hits;
   }
   return @hits;
}

=head2 next_hit

 Title   : next_hit
 Usage   : while( $hit = $result->next_hit()) { ... }
 Function: Returns the next available Hit object, representing potential
           matches between the query and various entities from the database.
 Returns : a Bio::Search::Hit::HitI object or undef if there are no more.
 Args    : none
 Note    : This method overrides L<Bio::Search::Result::GenericResult::next_hit>
           to take into account the possibility of multiple iterations, as
           occurs in PSI-BLAST reports.

           If there are multiple iterations, calling next_hit() traverses the
           all of the hits, old and new, for each iteration, calling next_hit()
           on each iteration.

           See Also: L<Bio::Search::Iteration::GenericIteration::next_hit>

=cut

sub next_hit {
    my ($self,@args) = @_;
    if ($self->{'_no_iterations'}) {
        return $self->SUPER::next_hit(@args);
    }

    my $iter_index;
    if (not defined $self->{'_last_hit'}) {
        $iter_index = $self->{'_iter_index'} = $self->_next_iteration_index;
    } else {
        $iter_index = $self->{'_iter_index'};
    }

    return if $iter_index >= scalar @{$self->{'_iterations'}};

    my $it = $self->{'_iterations'}->[$iter_index];
    my $hit = $self->{'_last_hit'} = $it->next_hit;

    return defined($hit) ? $hit : $self->next_hit;
}


=head2 num_hits

 Title   : num_hits
 Usage   : my $hitcount= $result->num_hits
 Function: returns the number of hits for this query result
 Returns : integer
 Args    : none
 Note    : This method overrides L<Bio::Search::Result::GenericResult::num_hits>
           to take into account the possibility of multiple iterations, as
           occurs in PSI-BLAST reports.

           If there are multiple iterations, calling num_hits() returns the
           number of 'new' hits for each iteration. These are the hits that did
           not occur in a previous iteration.

           See Also: L<Bio::Search::Result::GenericResult::num_hits>

=cut

sub num_hits{
   my ($self) = shift;
   if ($self->{'_no_iterations'}) {
       return $self->SUPER::num_hits;
   }
   if (not defined $self->{'_iterations'}) {
       $self->throw("Can't get Hits: data not collected.");
    }
    return scalar( $self->hits );
}

=head2 add_hit

 Title   : add_hit
 Usage   : $report->add_hit($hit)
 Function: Adds a HitI to the stored list of hits
 Returns : Number of HitI currently stored
 Args    : Bio::Search::Hit::HitI

=cut

sub add_hit {
    my ($self,$hit) = @_;
    my $iter = $self->iteration;
    if( $hit->isa('Bio::Search::Hit::HitI') ) { 
	return $iter->add_hit(-hit => $hit);
    } else { 
        $self->throw("Passed in a " .ref($hit). 
                     " as a Iteration which is not a Bio::Search::Hit::HitI.");
    }
    return $iter->num_hits;
}

=head2 add_iteration

 Title   : add_iteration
 Usage   : $report->add_iteration($iteration)
 Function: Adds a IterationI to the stored list of iterations
 Returns : Number of IterationI currently stored
 Args    : Bio::Search::Iteration::IterationI

=cut

sub add_iteration {
    my ($self,$i) = @_;
    if( $i->isa('Bio::Search::Iteration::IterationI') ) { 
        push @{$self->{'_iterations'}}, $i;
        $self->{'_iteration_count'}++;
    } else { 
        $self->throw("Passed in a " .ref($i). 
                     " as a Iteration which is not a Bio::Search::Iteration::IterationI.");
    }
    return scalar @{$self->{'_iterations'}};
}


=head2 next_iteration

 Title   : next_iteration
 Usage   : while( $it = $result->next_iteration()) { ... }
 Function: Returns the next Iteration object, representing all hits
           found within a given PSI-Blast iteration.
 Returns : a Bio::Search::Iteration::IterationI object or undef if there are no more.
 Args    : none

=cut

sub next_iteration {
    my ($self) = @_;

   unless($self->{'_iter_queue_started'}) {
       $self->{'_iter_queue'} = [$self->iterations()];
       $self->{'_iter_queue_started'} = 1;
   }
   return shift @{$self->{'_iter_queue'}};
}

=head2 iteration

 Usage     : $iteration = $blast->iteration( $number );
 Purpose   : Get an IterationI object for the specified iteration
             in the search result (PSI-BLAST).
 Returns   : Bio::Search::Iteration::IterationI object
 Throws    : Bio::Root::NoSuchThing exception if $number is not within 
             range of the number of iterations in this report.
 Argument  : integer (optional, if not specified get the last iteration)
             First iteration = 1

=cut

sub iteration {
    my ($self,$num) = @_;
    $num = scalar @{$self->{'_iterations'}} unless defined $num;
    unless ($num >= 1 and $num <= scalar $self->{'_iteration_count'}) {
        $self->throw(-class=>'Bio::Root::NoSuchThing',
                     -text=>"No such iteration number: $num. Valid range=1-$self->{'_iteration_count'}",
                     -value=>$num);
    }
    return $self->{'_iterations'}->[$num-1];
}

=head2 num_iterations

 Usage     : $num_iterations = $blast->num_iterations; 
 Purpose   : Get the number of iterations in the search result (PSI-BLAST).
 Returns   : Total number of iterations in the report
 Argument  : none (read-only)

=cut

sub num_iterations { shift->{'_iteration_count'} }

# Methods provided for consistency with BPpsilite.pm (now deprecated);
# these are now merely synonyms

=head2 number_of_iterations

 Usage     : $num_iterations = $blast->number_of_iterations; 
 Purpose   : Get the number of iterations in the search result (PSI-BLAST).
 Returns   : Total number of iterations in the report
 Argument  : none (read-only)
 Note      : Alias of L<num_iterations>.

=cut

sub number_of_iterations { shift->num_iterations }

=head2 round

 Usage     : $round = $blast->round( $number );
 Purpose   : Get an IterationI object for the specified iteration
             in the search result (PSI-BLAST).
 Returns   : Bio::Search::Iteration::IterationI object
 Throws    : Bio::Root::NoSuchThing exception if $number is not within 
             range of the number of iterations in this report.
 Argument  : integer (optional, if not specified get the last iteration)
             First iteration = 1
 Note      : Alias of L<iteration>.

=cut

sub round { shift->iteration(@_) }


=head2 iterations

 Title   : iterations
 Usage   : my @iterations = $result->iterations
 Function: Returns the IterationI objects contained within this Result
 Returns : Array of L<Bio::Search::Iteration::IterationI> objects
 Args    : none

=cut

sub iterations { 
    my $self = shift;
    my @its = ();
    if( ref($self->{'_iterations'}) =~ /ARRAY/i ) {
       @its = @{$self->{'_iterations'}};
    }
    return @its;
}

=head2 psiblast

 Usage     : if( $blast->psiblast ) { ... }
 Purpose   : Set/get a boolean indicator whether or not the report 
             is a PSI-BLAST report.
 Returns   : 1 if PSI-BLAST, undef if not.
 Argument  : 1 (when setting)

=cut

#----------------
sub psiblast {
#----------------
    my ($self, $val ) = @_;
    if( $val ) {
        $self->{'_psiblast'} = 1;
    }
    return $self->{'_psiblast'};
}


=head2 no_hits_found

 Usage     : $nohits = $blast->no_hits_found( $iteration_number );
 Purpose   : Get boolean indicator indicating whether or not any hits
             were present in the report.

             This is NOT the same as determining the number of hits via
             the hits() method, which will return zero hits if there were no
             hits in the report or if all hits were filtered out during the parse.

             Thus, this method can be used to distinguish these possibilities
             for hitless reports generated when filtering.

 Returns   : Boolean
 Argument  : (optional) integer indicating the iteration number (PSI-BLAST)
             If iteration number is not specified and this is a PSI-BLAST result,
             then this method will return true only if all iterations had
             no hits found.

=cut

sub no_hits_found {
    my ($self, $round) = @_;

    my $result = 0;   # final return value of this method.
    # Watch the double negative! 
    # result = 0 means "yes hits were found"
    # result = 1 means "no hits were found" (for the indicated iteration or all iterations)

    # If a iteration was not specified and there were multiple iterations,
    # this method should return true only if all iterations had no hits found.
    if( not defined $round ) {
        if( $self->{'_iterations'} > 1) {
            $result = 1;
            foreach my $i( 1..$self->{'_iterations'} ) {
                if( not defined $self->{"_iteration_$i"}->{'_no_hits_found'} ) {
                    $result = 0;
                    last;
                }
            }
        }
        else {
            $result = $self->{"_iteration_1"}->{'_no_hits_found'};
        }
    }
    else {
        $result = $self->{"_iteration_$round"}->{'_no_hits_found'};
    }

    return $result;
}


=head2 set_no_hits_found

 Usage     : $blast->set_no_hits_found( $iteration_number ); 
 Purpose   : Set boolean indicator indicating whether or not any hits
             were present in the report.
 Returns   : n/a
 Argument  : (optional) integer indicating the iteration number (PSI-BLAST)

=cut

sub set_no_hits_found {
    my ($self, $round) = @_;
    $round ||= 1;
    $self->{"_iteration_$round"}->{'_no_hits_found'} = 1;
}

=head2 _next_iteration_index

 Title   : _next_iteration_index
 Usage   : private

=cut

sub _next_iteration_index{
   my ($self,@args) = @_;
   return $self->{'_iteration_index'}++;
}


=head2 rewind

 Title   : rewind
 Usage   : $result->rewind;
 Function: Allow one to reset the Iteration iterator to the beginning
           Since this is an in-memory implementation
 Returns : none
 Args    : none

=cut

sub rewind {
   my $self = shift;
   $self->SUPER::rewind(@_);
   $self->{'_iteration_index'} = 0;
   foreach ($self->iterations) {
       $_->rewind;
   }
}


=head2 inclusion_threshold

 Title   : inclusion_threshold
 Usage   : my $incl_thresh = $result->inclusion_threshold; (read-only)
 Function: Gets the e-value threshold for inclusion in the PSI-BLAST 
           score matrix model (blastpgp) that was used for generating the report
           being parsed.
 Returns : number (real) or undef if not a PSI-BLAST report.
 Args    : none

=cut

sub inclusion_threshold {
    my $self = shift;
    return $self->{'_inclusion_threshold'};
}

1;
