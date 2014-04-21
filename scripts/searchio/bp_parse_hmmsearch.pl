#!/usr/bin/perl

use strict;
use warnings;

=head1 NAME

bp_parse_hmmsearch - parse single/multiple HMMSEARCH results file(s) with
                  different output options

=head1 SYNOPSIS

bp_parse_hmmsearch [--po] [--ps] -s hmmsearch_file

bp_parse_hmmsearch [--po] [--ps] -m index_file

=head1 DESCRIPTION

=head2 Mandatory Options:

  -s  HMMSEARCH file to parse.
  -m  INDEX file that contains a list of HMMSEARCH files for multiple
      parsing.

=head2 Special Options:

  --po    Print only the hits that have positive scores.
  --ps    Print the total of positive scores found.
  --help  Show this documentation.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list. Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

  Mauricio Herrera Cuadra <mauricio at open-bio.org>

=cut

# Modules, pragmas and variables to use
use Bio::SearchIO;
use Getopt::Long;
use vars qw($opt_s $opt_m $opt_po $opt_ps $opt_help);

# Gets options from the command line
GetOptions qw(-s:s -m:s --po --ps --help);

# Print documentation if help switch was given
exec('perldoc', $0) and exit() if $opt_help;

# If no mandatory options are given prints an error and exits
if (!$opt_s && !$opt_m) {
    print "ERROR: No HMMSEARCH or INDEX file has been specified.\n       Use
'--help' switch for documentation.\n" and exit();
} elsif ($opt_s && $opt_m) {
    print "ERROR: You must select only one option (-s or -m) for input.\n      
Use '--help' switch for documentation.\n" and exit();
}

# Initializes a counter for the domain positive scores if the option
# was given
my $pos_scores = 0 if $opt_ps;

# If single file mode was selected
if ($opt_s) {
    parse_hmmsearch($opt_s);

    # Prints the total domain positive scores if the option was given
    if ($opt_ps) {
        print "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
- - - -\n";
        print "Total domain positive scores: $pos_scores\n";
    }

# If multiple files mode was selected
} elsif ($opt_m) {

    # Opens the INDEX file sent as input
    open my $FH, '<', $opt_m or die "Could not read INDEX file '$opt_m': $!\n";

    # Cycle that extracts one line for every loop until finding the
    # end of file
    while (my $line = <$FH>) {

        # Deletes the new line characters from the line
        chomp $line;

        # Parses the result file in turn
        parse_hmmsearch($line);
        print "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
= = = =\n";
    }

    # Prints the total domain positive scores if the option was given
    print "Total domain positive scores: $pos_scores\n" if $opt_ps;

    # Closes INDEX files
    close $FH;
}

# Exits the program
exit();

# Subroutine that parses a HMMSEARCH results file
sub parse_hmmsearch {

    # Gets the parameters sent to the function
    my ($file) = @_;

    # Creates a new Bio::SearchIO object
    my $in = new Bio::SearchIO(
        -format => 'hmmer',
        -file   => $file,
    );

    # Loops through the results file
    while (my $result = $in->next_result()) {

        # Prints program name and version (these are values from
        # Bio::Search::Result::GenericResult methods)
        print $result->algorithm(), " ", $result->algorithm_version(), "\n";
        print "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
- - - -\n";

        # Prints HMM file and sequence database (these are values from
        # Bio::Search::Result::HMMERResult methods)
        print "HMM file:\t\t\t", $result->hmm_name(), "\n";
        print "Sequence database:\t\t", $result->sequence_file(), "\n";
        print "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
-\n";

        # Prints some values from Bio::Search::Result::GenericResult
        # methods
        print "Query HMM:\t\t\t", $result->query_name(), "\n";
        print "Accession:\t\t\t", $result->query_accession(), "\n";
        print "Description:\t\t\t", $result->query_description(), "\n";
        print "Total hits:\t\t\t", $result->num_hits(), "\n";

        # Loops through the sequence in turn
        while (my $hit = $result->next_hit()) {

            # If only positive scores option was given and the score
            # in turn is greater than zero
            if ($opt_po) {
                printHits($hit) if ($hit->score() >= 0);

            # Prints all hits otherwise
            } else {
                printHits($hit);
            }
        }
    }
}

# Subroutine that prints the values from a Bio::Search::Hit::HitI
# object
sub printHits {

    # Gets the parameters sent to the function
    my ($hit) = @_;

    # Prints some values from Bio::Search::Hit::HitI methods
    print "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
    print "Hit ", $hit->rank(), "\n";
    print "Sequence:\t\t\t", $hit->name(), "\n";
    print "Description:\t\t\t", $hit->description(), "\n";
    print "Score:\t\t\t\t", $hit->score(), "\n";
    print "E-value:\t\t\t", $hit->significance(), "\n";
    print "Number of domains:\t\t", $hit->num_hsps(), "\n";

    # Loops through the domain in turn
    while (my $hsp = $hit->next_hsp()) {

        # Prints some values from Bio::Search::HSP::HSPI methods
        print "   - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
        print "   Domain:\t\t\t", $hsp->rank(), " of ", $hit->num_hsps(), "\n";
        print "   seq-f:\t\t\t", $hsp->start('hit'), "\n";
        print "   seq-t:\t\t\t", $hsp->end('hit'), "\n";
        print "   hmm-f:\t\t\t", $hsp->start(), "\n";
        print "   hmm-t:\t\t\t", $hsp->end(), "\n";
        print "   score:\t\t\t", $hsp->score(), "\n";
        $pos_scores++ if ($hsp->score() >= 0) && $opt_ps;
        print "   E-value:\t\t\t", $hsp->evalue(), "\n";
        my $hmm_string = $hsp->query_string();
        $hmm_string =~ s/<-\*$//;
        print "   hmm string:\t\t\t", $hmm_string, "\n";
        print "   homology string:\t\t", $hsp->homology_string(), "\n";
        print "   hit string:\t\t\t", $hsp->hit_string(), "\n";
    }
}
