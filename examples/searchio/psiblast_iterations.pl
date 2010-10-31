#!/usr/bin/perl 

# Demonstrates the use of a SearchIO parser for processing
# the iterations within a PSI-BLAST report.
#
# Usage:
#   STDIN:  none; supply filename of PSI-BLAST report on command-line
#   STDOUT: information parsed from the input data.
#   STDERR: errors.
#
# For more documentation about working with Iteration objects,
# see docs for:
#   Bio::Search::Iteration::IterationI
#
# Author: Steve Chervitz <sac@bioperl.org>

use strict;
use lib '../../';

use Bio::SearchIO;

my $file = shift or die "Usage: $0 <BLAST-report-file>\n";
my $in = new Bio::SearchIO(-format => 'blast',
                           -file => $file, #comment this out to read STDIN
                           #-fh => \*ARGV,  #uncomment this to read STDIN
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

