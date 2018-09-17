#
# BioPerl module for Bio::Search::Iteration::GenericIteration
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# Copyright Steve Chervitz
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

# TODO: Consider calling this BlastIteration (strongly) and maybe simplifying IterationI.

=head1 NAME

Bio::Search::Iteration::GenericIteration - A generic implementation of the Bio::Search::Iteration::IterationI interface.

=head1 SYNOPSIS

    use Bio::Search::Iteration::GenericIteration;
    my $it = Bio::Search::GenericIteration->new(
                              -number => 1,
                              -converged => 0,
                              -newhits_unclassified => [@newhits_unclass],
                              -newhits_below => [@newhits_below_threshold],
                              -newhits_not_below => [@newhits_not_below_threshold],
                              -oldhits_below => [@oldhits_below_threshold],
                              -oldhits_newly_below => [@oldhits_newly_below_threshold],
                              -oldhits_not_below => [@oldhits_not_below_threshold],
                                        );

# TODO: Describe how to configure a SearchIO stream so that it generates
#       GenericIteration objects.


=head1 DESCRIPTION

This module acts as a container for Bio::Search::Hit::HitI objects,
allowing a Search::Result::ResultI object to partition its hits based
on which iteration the hit occurred in (e.g., a PSI-BLAST round).

Unless you're writing a parser, you won't ever need to create a
GenericIteration or any other IterationI-implementing object. If you use
the SearchIO system, IterationI objects are created automatically from
a SearchIO stream which returns Bio::Search::Result::ResultI objects
and you get the IterationI objects via the ResultI API.

For documentation on what you can do with GenericIteration (and other IterationI
objects), please see the API documentation in
L<Bio::Search::Iteration::IterationI|Bio::Search::Iteration::IterationI>.

Bio::Search::Iteration::GenericIteration is similar in spirit to the deprecated
Bio::Tools::BPlite::Iteration modules in bioperl releases prior to 1.6, except
that Bio::Search::Iteration::GenericIteration is a pure container, without any
parsing functionality as is in Bio::Tools::BPlite::Iteration.

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


package Bio::Search::Iteration::GenericIteration;
use strict;


use base qw(Bio::Root::Root Bio::Search::Iteration::IterationI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::Iteration->new(%args);
 Function: Builds a new Bio::Search::Iteration object 
 Returns : Bio::Search::Iteration::GenericIteration object
 Args    : -number => integer for the number of this iteration (required)
           -converged => boolean value whether or not the iteration converged
           -newhits_unclassified => array reference to hits that were not found
                       in a previous iteration for the iteration and have not been 
                       classified with regard to the inclusion threshold

           # The following are only used for PSI-BLAST reports:

           -newhits_below => array reference to hits were not found in a 
                        previous iteration and are below the inclusion threshold.
           -newhits_not_below => array reference to hits that were not found in a 
                        previous iteration below threshold that and are not below 
                        the inclusion threshold threshold.
           -oldhits_below => array reference to hits that were found
                        in a previous iteration below inclusion threshold and are
                        still below threshold in the current iteration.
           -oldhits_newly_below => array reference to hits that were found
                        in a previous iteration above threshold but are below
                        threshold in the current iteration.
           -oldhits_not_below => array reference to hits that were found in a
                        previous iteration above threshold that and are still above
                        the inclusion threshold threshold.

           -hit_factory => Bio::Factory::ObjectFactoryI capable of making
                        Bio::Search::Hit::HitI objects

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($number, $newhits_unclassified, $newhits_below, $newhits_not_below,
      $oldhits_below, $oldhits_newly_below, $oldhits_not_below, $converged,
      $h_f) =
      $self->_rearrange([qw(NUMBER
                            NEWHITS_UNCLASSIFIED
                            NEWHITS_BELOW
                            NEWHITS_NOT_BELOW
                            OLDHITS_BELOW
                            OLDHITS_NEWLY_BELOW
                            OLDHITS_NOT_BELOW
                            CONVERGED
                            HIT_FACTORY
                           )], @args);

  if( ! defined $number ) { 
      $self->throw(-class=>'Bio::Root::BadParameter',
                   -text=>"Iteration number not specified.");
  } else { 
      $self->number($number);
  }

  defined $converged && $self->converged($converged);

  # TODO: Performance optimization test calling add_hit() vs. simple assignment:
  #       push @{$self->{'_hits_new'}}, @{$newhits};
  #             vs.
  #       foreach(@{$newhits_below}) {$self->add_hit(-hit=>$_, -old=>0, -below=>1);}

  if(defined $newhits_unclassified ) {
    if( ref($newhits_unclassified) =~ /ARRAY/i) {
         push @{$self->{'_newhits_unclassified'}}, @{$newhits_unclassified};
    } else {
      $self->throw(-class=>'Bio::Root::BadParameter',
                   -text=>"Parameter NEWHITS is not an array ref: $newhits_unclassified");
    }
  } else {
      $self->{'_newhits_unclassified'} = [];
  }

  if(defined $newhits_below ) {
    if( ref($newhits_below) =~ /ARRAY/i) {
        push @{$self->{'_newhits_below_threshold'}}, @{$newhits_below};
    } else {
      $self->throw(-class=>'Bio::Root::BadParameter',
                   -text=>"Parameter NEWHITS_BELOW is not an array ref: $newhits_below");
    }
  } else {
      $self->{'_newhits_below_threshold'} = [];
  }

  if(defined $newhits_not_below ) {
    if( ref($newhits_not_below) =~ /ARRAY/i) {
         push @{$self->{'_newhits_not_below_threshold'}}, @{$newhits_not_below};
    } else {
      $self->throw(-class=>'Bio::Root::BadParameter',
                   -text=>"Parameter NEWHITS_NOT_BELOW is not an array ref: $newhits_not_below");
    }
  } else {
      $self->{'_newhits_not_below_threshold'} = [];
  }

  if(defined $oldhits_below ) {
    if( ref($oldhits_below) =~ /ARRAY/i) {
         push @{$self->{'_oldhits_below_threshold'}}, @{$oldhits_below};
    } else {
      $self->throw(-class=>'Bio::Root::BadParameter',
                   -text=>"Parameter OLDHITS_BELOW is not an array ref: $oldhits_below");
    }
  } else {
      $self->{'_oldhits_below_threshold'} = [];
  }

  if(defined $oldhits_newly_below ) {
    if( ref($oldhits_newly_below) =~ /ARRAY/i) {
         push @{$self->{'_oldhits_newly_below_threshold'}}, @{$oldhits_newly_below};
    } else {
      $self->throw(-class=>'Bio::Root::BadParameter',
                   -text=>"Parameter OLDHITS_NEWLY_BELOW is not an array ref: $oldhits_newly_below");
    }
  } else {
      $self->{'_oldhits_newly_below_threshold'} = [];
  }

  if(defined $oldhits_not_below ) {
    if( ref($oldhits_not_below) =~ /ARRAY/i) {
         push @{$self->{'_oldhits_not_below_threshold'}}, @{$oldhits_not_below};
    } else {
      $self->throw(-class=>'Bio::Root::BadParameter',
                   -text=>"Parameter OLDHITS_NOT_BELOW is not an array ref: $oldhits_not_below");
    }
  } else {
      $self->{'_oldhits_not_below_threshold'} = [];
  }
  
  $self->hit_factory($h_f) if $h_f;
  
  return $self;
}


=head2 number

See documentation in Bio::Search::Iteration::IterationI.

=cut

sub number {
    my ($self,$value) = @_;
    my $previous = $self->{'_number'};
    if( defined $value || ! defined $previous ) {
        $value = $previous = '' unless defined $value;
        $self->{'_number'} = $value;
    } 
    return $previous;
}

=head2 converged

See documentation in Bio::Search::Iteration::IterationI.

=cut

sub converged {
    my ($self,$value) = @_;
    my $previous = $self->{'_converged'};
    if( defined $value || ! defined $previous ) {
        $value = $previous = '' unless defined $value;
        $self->{'_converged'} = $value;
    } 
    return $previous;
}


=head2 hit_factory

 Title   : hit_factory
 Usage   : $hit->hit_factory($hit_factory)
 Function: Get/set the factory used to build HitI objects if necessary.
 Returns : Bio::Factory::ObjectFactoryI
 Args    : Bio::Factory::ObjectFactoryI

=cut

sub hit_factory {
    my $self = shift;
    if (@_) { $self->{_hit_factory} = shift }
    return $self->{_hit_factory} || return;
}

=head2 next_hit

This iterates through all old hits as returned by L<oldhits> 
followed by all new hits as returned by L<newhits>.

For more documentation see L<Bio::Search::Iteration::IterationI::next_hit()|Bio::Search::Iteration::IterationI>.

=cut

sub next_hit {
   my ($self) = @_;

   unless($self->{'_hit_queue_started'}) {
       $self->{'_hit_queue'} = ( [$self->oldhits(), $self->newhits()] );
       $self->{'_hit_queue_started'} = 1;
   }
   return shift @{$self->{'_hit_queue'}};
}

=head2 next_hit_new

See documentation in L<Bio::Search::Iteration::IterationI::next_hit_new()|Bio::Search::Iteration::IterationI>.

=cut

sub next_hit_new {
   my ($self) = @_;

   unless($self->{'_hit_queue_new_started'}) {
       $self->{'_hit_queue_new'} = [$self->newhits()];
       $self->{'_hit_queue_new_started'} = 1;
   }
   return shift @{$self->{'_hit_queue_new'}};
}

=head2 next_hit_old

See documentation in L<Bio::Search::Iteration::IterationI::next_hit_old()|Bio::Search::Iteration::IterationI>.

=cut

sub next_hit_old {
   my ($self,$found_again) = @_;

   unless($self->{'_hit_queue_old_started'}) {
       $self->{'_hit_queue_old'} = [$self->oldhits()];
       $self->{'_hit_queue_old_started'} = 1;
   }
   return shift @{$self->{'_hit_queue_old'}};
}

=head2 rewind

 Title   : rewind
 Usage   : $iteration->rewind;
 Function: Allow one to reset the Hit iterators to the beginning
           Since this is an in-memory implementation
 Returns : none
 Args    : none

=cut

sub rewind {
   my $self = shift;
   $self->{'_hit_queue_started'} = 0;
   $self->{'_hit_queue_new_started'} = 0;
   $self->{'_hit_queue_old_started'} = 0;
   foreach ($self->hits) {
      $_->rewind;
   }
}


=head2 num_hits

See documentation in L<Bio::Search::Iteration::IterationI::num_hits()|Bio::Search::Iteration::IterationI>.

=cut

sub num_hits {
   my $self = shift;

   return $self->num_hits_old + $self->num_hits_new;
}

=head2 num_hits_new

See documentation in L<Bio::Search::Iteration::IterationI::num_hits_new()|Bio::Search::Iteration::IterationI>.

=cut

sub num_hits_new {
   my $self = shift;

    return scalar $self->newhits();
}

=head2 num_hits_old

See documentation in L<Bio::Search::Iteration::IterationI::num_hits_old()|Bio::Search::Iteration::IterationI>.

=cut

sub num_hits_old {
   my ($self,$found_again) = @_;

   return scalar $self->oldhits();
}

=head2 add_hit

See documentation in L<Bio::Search::Iteration::IterationI::add_hit()|Bio::Search::Iteration::IterationI>.

=cut

sub add_hit { 
    my ($self,@args) = @_;
    my( $hit, $old, $below, $newly_below ) = 
        $self->_rearrange([qw(HIT
                              OLD
                              BELOW_THRESHOLD
                              NEWLY_BELOW
                             )], @args);
    my $count = 0;

    unless( ref($hit) eq 'HASH' || $hit->isa('Bio::Search::Hit::HitI') ) { 
        $self->throw(-class=>'Bio::Root::BadParameter',
                     -text=>"Passed in " .ref($hit). 
                    " as a Hit which is not a Bio::Search::Hit::HitI.");
    }

    if($old) {
        if ($newly_below) {
            push @{$self->{'_oldhits_newly_below_threshold'}}, $hit;
            $count = scalar @{$self->{'_oldhits_newly_below_threshold'}};
        } elsif ($below) {
            push @{$self->{'_oldhits_below_threshold'}}, $hit;
            $count = scalar @{$self->{'_oldhits_below_threshold'}};
        } else {
            push @{$self->{'_oldhits_not_below_threshold'}}, $hit;
            $count = scalar @{$self->{'_oldhits_not_below_threshold'}};
        }
    } elsif (defined $old) {
        # -old is defined but false, so this is a new PSI-BLAST hit
        if ($below) {
            push @{$self->{'_newhits_below_threshold'}}, $hit;
            $count = scalar @{$self->{'_newhits_below_threshold'}};
        } elsif (defined $below) {
            push @{$self->{'_newhits_not_below_threshold'}}, $hit;
            $count = scalar @{$self->{'_newhits_not_below_threshold'}};
        } else {
            # -below not defined, PSI-BLAST threshold may not be known
            push @{$self->{'_newhits_unclassified'}}, $hit;
            $count = scalar @{$self->{'_newhits_unclassified'}};
        }
    } else {
        # -old not defined, so it's non-PSI-BLAST
        push @{$self->{'_newhits_unclassified'}}, $hit;
        $count = scalar @{$self->{'_newhits_unclassified'}};
    }
    return $count;
}

=head2 hits

See Documentation in InterfaceI.

=cut

sub hits  { 
    my $self = shift;
#    print STDERR "Called GenericIteration::hits()\n";
    my @new = $self->newhits;
    my @old = $self->oldhits;
    return ( @new, @old );
}

=head2 newhits

Returns a list containing all newhits in this order:

newhits_below_threshold
newhits_not_below_threshold
newhits_unclassified

See more documentation in InterfaceI.

=cut

sub newhits  { 
    my $self = shift;
    my @hits = $self->newhits_below_threshold;
    push @hits, $self->newhits_not_below_threshold;
    push @hits, $self->newhits_unclassified;
    return @hits;
}

=head2 newhits_below_threshold

See documentation in L<Bio::Search::Iteration::IterationI::newhits_below_threshold()|Bio::Search::Iteration::IterationI>.

=cut

sub newhits_below_threshold  { 
    my $self = shift;
    if (ref $self->{'_newhits_below_threshold'} ) {
        my $factory = $self->hit_factory || return @{$self->{'_newhits_below_threshold'}};
        for (0..$#{$self->{'_newhits_below_threshold'}}) {
            ref(${$self->{'_newhits_below_threshold'}}[$_]) eq 'HASH' || next;
            ${$self->{'_newhits_below_threshold'}}[$_] = $factory->create_object(%{${$self->{'_newhits_below_threshold'}}[$_]});
        }
        return @{$self->{'_newhits_below_threshold'}};
    }
    return;
}

=head2 newhits_not_below_threshold

See documentation in L<Bio::Search::Iteration::IterationI::newhits_not_below_threshold()|Bio::Search::Iteration::IterationI>.

=cut

sub newhits_not_below_threshold  { 
    my $self = shift;
    if (ref $self->{'_newhits_not_below_threshold'} ) {
        my $factory = $self->hit_factory || return @{$self->{'_newhits_not_below_threshold'}};
        for (0..$#{$self->{'_newhits_not_below_threshold'}}) {
            ref(${$self->{'_newhits_not_below_threshold'}}[$_]) eq 'HASH' || next;
            ${$self->{'_newhits_not_below_threshold'}}[$_] = $factory->create_object(%{${$self->{'_newhits_not_below_threshold'}}[$_]});
        }
        return @{$self->{'_newhits_not_below_threshold'}};
    }
    return;
}

=head2 newhits_unclassified

 Title   : newhits_unclassified
 Usage   : foreach( $iteration->hits_unclassified ) {...}
 Function: Gets all newhits that have not been partitioned into
           sets relative to the inclusion threshold.
 Returns : Array of Bio::Search::Hit::HitI objects.
 Args    : none

=cut

sub newhits_unclassified  { 
    my $self = shift;
    if (ref $self->{'_newhits_unclassified'} ) {
        my $factory = $self->hit_factory || return @{$self->{'_newhits_unclassified'}};
        for (0..$#{$self->{'_newhits_unclassified'}}) {
            ref(${$self->{'_newhits_unclassified'}}[$_]) eq 'HASH' || next;
            ${$self->{'_newhits_unclassified'}}[$_] = $factory->create_object(%{${$self->{'_newhits_unclassified'}}[$_]});
        }
        return @{$self->{'_newhits_unclassified'}};
    }
    return;
}

=head2 oldhits

Returns a list containing all oldhits in this order:

oldhits_below_threshold
oldhits_newly_below_threshold
oldhits_not_below_threshold

See more documentation in InterfaceI.

=cut

sub oldhits  { 
    my $self = shift;
    my @hits = $self->oldhits_below_threshold;
    push @hits, $self->oldhits_newly_below_threshold;
    push @hits, $self->oldhits_not_below_threshold;
    return @hits;
}

=head2 oldhits_below_threshold

See documentation in L<Bio::Search::Iteration::IterationI::oldhits_below_threshold()|Bio::Search::Iteration::IterationI>.

=cut

sub oldhits_below_threshold  { 
    my $self = shift;
    if (ref $self->{'_oldhits_below_threshold'} ) {
        my $factory = $self->hit_factory || return @{$self->{'_oldhits_below_threshold'}};
        for (0..$#{$self->{'_oldhits_below_threshold'}}) {
            ref(${$self->{'_oldhits_below_threshold'}}[$_]) eq 'HASH' || next;
            ${$self->{'_oldhits_below_threshold'}}[$_] = $factory->create_object(%{${$self->{'_oldhits_below_threshold'}}[$_]});
        }
        return @{$self->{'_oldhits_below_threshold'}};
    }
    return;
}

=head2 oldhits_newly_below_threshold

See documentation in L<Bio::Search::Iteration::IterationI::oldhits_newly_below_threshold()|Bio::Search::Iteration::IterationI>.

=cut

sub oldhits_newly_below_threshold  { 
    my $self = shift;
    if (ref $self->{'_oldhits_newly_below_threshold'} ) {
        my $factory = $self->hit_factory || return @{$self->{'_oldhits_newly_below_threshold'}};
        for (0..$#{$self->{'_oldhits_newly_below_threshold'}}) {
            ref(${$self->{'_oldhits_newly_below_threshold'}}[$_]) eq 'HASH' || next;
            ${$self->{'_oldhits_newly_below_threshold'}}[$_] = $factory->create_object(%{${$self->{'_oldhits_newly_below_threshold'}}[$_]});
        }
        return @{$self->{'_oldhits_newly_below_threshold'}};
    }
    return;
}

=head2 oldhits_not_below_threshold

See documentation in L<Bio::Search::Iteration::IterationI::oldhits_not_below_threshold()|Bio::Search::Iteration::IterationI>.

=cut

sub oldhits_not_below_threshold  { 
    my $self = shift;
    if (ref $self->{'_oldhits_not_below_threshold'} ) {
        my $factory = $self->hit_factory || return @{$self->{'_oldhits_not_below_threshold'}};
        for (0..$#{$self->{'_oldhits_not_below_threshold'}}) {
            ref(${$self->{'_oldhits_not_below_threshold'}}[$_]) eq 'HASH' || next;
            ${$self->{'_oldhits_not_below_threshold'}}[$_] = $factory->create_object(%{${$self->{'_oldhits_not_below_threshold'}}[$_]});
        }
        return @{$self->{'_oldhits_not_below_threshold'}};
    }
    return;
}

=head2 hits_below_threshold

See documentation in L<Bio::Search::Iteration::IterationI::hits_below_threshold()|Bio::Search::Iteration::IterationI>.

=cut

sub hits_below_threshold  {
    my $self = shift;
    my @hits = $self->newhits_below_threshold;
    push @hits, $self->oldhits_newly_below_threshold;
    return @hits;
}

=head2 get_hit

See documentation in L<Bio::Search::Iteration::IterationI::get_hit()|Bio::Search::Iteration::IterationI>.

To free up the memory used by the get_hit() functionality, call free_hit_lookup().

This functionality might be useful at the Result level, too.
BlastResult::get_hit() would return a list of HitI objects for hits 
that occur in multiple iterations.

=cut

sub get_hit {
    my ($self,$name) = @_;
    $self->_create_hit_lookup() unless defined $self->{'_hit_lookup'};

    return $self->{'_hit_lookup'}->{"\U$name"};
}

# Internal method.
sub _create_hit_lookup {
    my $self = shift;
    foreach ($self->hits) {
        my $hname = $_->name;
        $self->{'_hit_lookup'}->{"\U$hname"} = $_;
    }
}

=head2 free_hit_lookup

 Purpose : Frees up the memory used by the get_hit() functionality.
           For the memory-conscious.

=cut

sub free_hit_lookup {
    my $self = shift;
    undef $self->{'_hit_lookup'};
}

1;
