#!/usr/bin/perl
use strict;
use warnings;

=head1 NAME

bp_find-blast-matches.pl - extract DNA sequences based on BLAST hits

=head1 SYNOPSIS

bp_find-blast-matches.pl [-h -e -p -5 -n -o -3 -header] -blast <BLAST_FILE> -fasta <FASTA_FILE>

=head1 OPTIONS

=head2 Mandatory:

=over

=item B<-blast>

BLAST output file to read from. The alignment should use the file specified by
'-fasta' option ideally

=item B<-fasta>

FASTA file to read from. This is where the sequence will be extracted from

=back

=head2 Optional:

=over

=item B<-h>

Displays this help message

=item B<-e>

Maximum e-value for matches (0.01 by default)

=item B<-p>

Number of base pairs of 5' region to be included with the sequence

=item B<-5>

Number of base pairs of 5' region only, excluding the regular sequence

=item B<-3>

Number of base pairs of 3' region only, excluding the regular sequence

=item B<-n>

Number of top hits to display, starting with the best hit

=item B<-o>

Exact match to display (this option can't be used in conjunction with '-n'

=item B<-header>

The FASTA header to display instead of the default

=back

=head1 DESCRIPTION

This script takes a BLAST output file and a FASTA file as arguments, 
given after the '-blast' and '-fasta' options respectively. The BLAST output 
file should have been generated with your sequence of interest and the 
FASTA file supplied as an argument.
Example: find-blast-matches.pl -blast BLAST_FILE -fasta FASTA_FILE

It parses through the BLAST file to check for high quality matches, 
which are then searched for in the FASTA file.  The sequence may vary 
from you candidate sequence, hence the BLAST search prior. 

The sequence from the FASTA file is then displayed to STDOUT.
Optional arguments can be used, such as to extract the 5' or 3' region.

=head1 AUTHOR

Gabriel Abud - E<lt>gabriel.jabud-at-gmail.comE<gt>

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via 
email or the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 EDIT HISTORY

2014-08-04 - Gabriel Abud
    First features added

=head1 DEPENDANCIES

Getopt::long, Pod::Usage, Bio::SearchIO, Bio::Seq, Bio::SeqIO,
File::Basename

=cut


# Modules
use Bio::SearchIO qw(new);
use Bio::Seq qw(new);
use Bio::SeqIO qw(new);
use File::Basename qw(basename);
use Getopt::Long qw(GetOptions);
use Pod::Usage;

# Variables
my $baseProg = basename($0);    # Program name
my $line;
my @scaffolds;
my $inputFile;
my $blastFile;
my $scaffoldFind;
my $baseList;
my @start;
my @end;
my @strand;
my @arrayBases;
my $query_name;
my @query_names;
my @accessions;
my $accession;
my $baseFile;
my $total_size;
my $title;
my $default_header;
my $in;
my $out;

# Command line options
my $exact_match;    # Undef by default
my $e_value = 0.01; # Default e-value
my $matches = 1;    # Default number of matches
my $promoter;       # Undef by default
my $three_prime;    # Undef by default
my $promoter_only;  # Undef by default
my $opt_help;       # Undef by default
my $header;         # Undef by default

# Functions for proper command line usage
sub syntax {
    print STDERR
      "Usage: $baseProg -blast 'BLAST_file' -fasta 'FASTA_file' [OPTIONS]\n";
    print STDERR "Try '$baseProg --help' for more information.\n";
    exit;
}

sub help {
    print STDERR
      "\nNAME:\n",
      "\t$baseProg - extract a DNA sequence based on BLAST hits\n\n",
      "SYNTAX:\n",
      "\t$baseProg -blast 'BLAST_file' -fasta 'FASTA__file'  [OPTIONS]\n\n",
      "OPTIONS:\n",
      "\t(All options require an additional number argument [ie: -e 0.01])\n\n",
      "\t-e, maximum e-value for matches (0.01 by default)\n\n",
      "\t-p, number of base pairs of 5' region to be included with the\n",
      "\t sequence of interest\n\n",
      "\t-5, number of base pairs of 5' region, excluding the sequence\n",
      "\t of interest (unlike '-p')\n\n",
      "\t-n, number of top hits to display, starting with the highest hit\n",
      "\t(1 by default)\n\n",
      "\t-o, exact match to display (this option can't be used in conjuction\n",
      "\t with '-n')\n\n",
      "\t-3, number of base pairs of 3' region to display\n\n",
      "\t-header, the fasta header to display instead of the default\n\n";
    exit;
}

# Get command line options
GetOptions(
    'e=f'      => \$e_value,
    'p=i'      => \$promoter,
    'n=i'      => \$matches,
    'o=i'      => \$exact_match,
    '5=i'      => \$promoter_only,
    '3=i'      => \$three_prime,
    'help|h'     => \$opt_help,
    'header|head=s' => \$header,
    'blast|b=s'  => \$blastFile,
    'fasta|f=s'  => \$inputFile
) or pod2usage(-exitstatus => 0, -verbose => 1);

# Help screen
pod2usage(-exitstatus => 0, -verbose => 2) if $opt_help;
#help() if defined $opt_help;

# Checks for required arguments
pod2usage(-exitstatus => 0, -verbose => 1, 
          -msg => "You must specify the FASTA and BLAST files to read from!\n")
    if (!defined $blastFile || !defined $inputFile);
#syntax() if ( !defined $blastFile || !defined $inputFile );

# Checks for any negative numbers
#syntax()
pod2usage(-exitstatus => 0, -verbose => 1,
          -msg => "You must use positive numbers as values to options!")
  if ( (defined $e_value && $e_value < 0)
    || (defined $promoter && $promoter < 0)
    || (defined $matches && $matches < 0)
    || (defined $exact_match && $exact_match < 0)
    || (defined $promoter_only && $promoter_only < 0)
    || (defined $three_prime && $three_prime < 0 ) );

if ( $matches > 1 && defined $exact_match ) {
    print STDERR "Cannot use both options '-n' and '-o' at the same time\n";
    print STDERR "(Type '$baseProg --help' for more information\n)";
    exit;
}

if ( defined $promoter && defined $promoter_only ) {
    print STDERR "Cannot use both options '-p' and '-5' at the same time\n";
    print STDERR "(Type '$baseProg --help' for more information)\n";
    exit;
}

if ( defined $three_prime && ( defined $promoter || defined $promoter_only ) ) {
    print STDERR "Cannot use both '-3' with '-p' or '-5' at the same time\n";
    print STDERR "(Type $baseProg --help' for more information)\n";
    exit;
}

# Class used to search through the blast file
eval { $in = Bio::SearchIO->new( -file => $blastFile, -format => "blast" ); };
if ($@) {
    die "'$blastFile' does not appear to be a BLAST output file! Exiting...\n";
}
$out = Bio::SeqIO->new( -fh => \*STDOUT, -format => "fasta" );

# Creates arrays of all the scaffold names and the coordinates of where those scaffolds are found
my $n = 0;
OUTERLOOP: while ( my $result = $in->next_result ) {
  LOOP: while ( my $hit = $result->next_hit ) {
        while ( my $hsp = $hit->next_hsp ) {

            # Finds all matches with an evalue <= $e_value (or 0.01 by default)
            if ( $hsp->evalue <= $e_value ) {
                ( $start[$n], $end[$n] ) = $hsp->range('hit');
                $scaffolds[$n] = $hit->name, $strand[$n] = $hsp->strand('hit');
                $n += 1;
            }

            # If an exact_match option is given
            if ( defined($exact_match) && $exact_match == $n ) {
                ( $start[0], $end[0] ) = $hsp->range('hit');
                $scaffolds[0] = $hit->name, $strand[0] = $hsp->strand('hit');
                last OUTERLOOP;
            }
            elsif ( !defined($exact_match) && $n >= $matches )
            { # Exits after the correct amount of matches have been found (1 by default)
                last LOOP;
            }
        }
    }
}

$baseFile = basename $inputFile;

open INFILE, "<$inputFile"
  or die "Couldn't open the input file '$inputFile'! Exiting...\n";

$scaffoldFind = 0;
my %scaffoldList;
my %baseList;

# Extracts only the scaffolds of interest, avoiding duplicates
if ( defined $exact_match ) {
    $scaffoldList{ $scaffolds[0] } = 1;
}
else {
    foreach my $num ( 0 .. $#scaffolds ) {
        if ( !defined $scaffoldList{ $scaffolds[$num] } ) {
            $scaffoldList{ $scaffolds[$num] } = 1;
        }
    }
}
my $remaining_scaffolds = my $unique_scaffolds = keys(%scaffoldList);

local $/ = "\n>";

# Reads the FASTA file here, storing FASTA headers where the sequence is found
while ( $line = <INFILE> ) {
    chomp($line);
    next if ( $line =~ m/^\s*?$/ );

    foreach my $scaffold ( sort keys %scaffoldList ) {
        if ( $line =~ /^[ \t]{0,5}$scaffold/ )
        {    # True if FASTA segment contains the scaffold
            $line =~ s/^.*?\n//;
            $line =~ s/\s//g;
            $baseList{$scaffold} = $line;
            $scaffoldFind++;
            $remaining_scaffolds--;
        }
    }
    last unless ($remaining_scaffolds);
}

if ( $scaffoldFind != $unique_scaffolds ) {
    print STDERR "The scaffold specified in the BLAST file was not found.\n";
    print STDERR "Make sure you are using the correct FASTA file.\n";
    exit;
}

for my $m ( 0 .. $#scaffolds )
{    # Runs a loop as many times as there are scaffolds
    $baseList = $baseList{ $scaffolds[$m] };

    # Print title line for each scaffold
    $accession      = $baseFile;
    $default_header = "(BLAST hit:$scaffolds[$m]|";

    my $real_start = $start[$m];
    my $real_end   = $end[$m];

    # For "normal", positive strands (+/+)
    if ( $strand[$m] == 1 ) {

        # If the -p flag was used
        if ( defined $promoter ) {    # 5' region specified with -p flag
             # If 5' region is too large for sequence, don't include the 5' region
            if ( $promoter >= $start[$m] ) {
                print STDERR "ERROR: 5' region is too big!! (max promoter = ",
                  $start[$m] - 1, " )\n",
                  "Showing sequence without 5' region...\n";
                $baseList = substr $baseList, $start[$m], $end[$m] - $start[$m];
            }
            else {
                my $real_start = $start[$m] - $promoter;
                $baseList = substr $baseList, $start[$m] - $promoter - 1,
                  $promoter + $end[$m] - $start[$m] + 1;
            }
        }

        # If the -3 flag was used
        elsif ( defined $three_prime ) {
            $total_size = length($baseList);

            if ( ( $three_prime + $end[$m] ) > $total_size ) {
                die "ERROR: 3' region is too big!! ",
                  "(max 3' region = ", $total_size - $end[$m], ")\n",
                  ;
            }
            else {
                $real_start = $end[$m] + 1;
                $real_end   = $end[$m] + $three_prime;
                $baseList   = substr $baseList, $end[$m], $three_prime;
            }
        }

        # If the -5 flag was used
        elsif ( defined $promoter_only )
        {    # 5' region was specified with -5 flag
             # If 5' region is too large for sequence, don't include the 5' region
            if ( $promoter_only >= $start[$m] ) {
                die "ERROR: 5' region is too big!! (max promoter = ",
                  $start[$m] - 1, " )\n",
                  ;
            }
            else {
                $real_start = $start[$m] - $promoter_only;
                $real_end   = $start[$m] - 1;
                $baseList   = substr $baseList, $start[$m] - $promoter_only - 1,
                  $promoter_only;
            }
        }
        else {    # Default: just the BLAST hit (no 5' or 3')
            $baseList = substr $baseList, $start[$m] - 1,
              $end[$m] - $start[$m] + 1;
        }
    }

    # Checks to see if the hit sequence is the reverse compliment
    elsif ( $strand[$m] == -1 ) {

        # If -p flag was used:
        if ( defined $promoter ) {
            $total_size = length($baseList);

            if ( ( $promoter + $end[$m] ) > $total_size ) {
                print STDERR "ERROR: 5' region is too big!! ",
                  "(max 5' region = ", $total_size - $end[$m], ")\n",
                  "Showing sequence without 5' region...\n";
                $baseList = substr $baseList, $start[$m], $end[$m] - $start[$m];
            }
            else {
                $real_start = $start[$m] + 1;
                $real_end   = $end[$m] + $promoter;
                $baseList   = substr $baseList, $start[$m] - 1,
                  $end[$m] - $start[$m] + $promoter + 1;
            }

        }

        # If -3 flag was used
        elsif ( defined $three_prime ) {
            if ( $three_prime >= $start[$m] ) {
                die "ERROR: 3' region is too big!! (max = ", $start[$m] - 1,
                  " )\n";
            }
            else {
                $real_start = $start[$m] - $three_prime;
                $real_end   = $start[$m] - 1;
                $baseList   = substr $baseList, $start[$m] - $three_prime - 1,
                  $three_prime;
            }
        }

        # If -5 flag was used
        elsif ( defined $promoter_only ) {
            $total_size = length($baseList);

            if ( ( $promoter_only + $end[$m] ) > $total_size ) {
                die "ERROR: 5' region is too big!! ",
                  "(max promoter = ", $total_size - $end[$m], ")\n",
                  ;
            }
            else {
                $real_start = $end[$m] + 1;
                $real_end   = $end[$m] + $promoter_only;
                $baseList   = substr $baseList, $end[$m], $promoter_only;
            }

        }
        else {    # If 5' region wasn't specified at all
            $baseList = substr $baseList, $start[$m] - 1,
              $end[$m] - $start[$m] + 1;
        }
    }
    $default_header .= "$real_start-$real_end)";
    $default_header .= " showing 5' region ($promoter_only bp) only"
      if defined($promoter_only);
    $default_header .= " showing 3' region ($three_prime bp)"
      if defined($three_prime);
    $default_header .= " with 5' region ($promoter bp)" if defined($promoter);

    my $seq_obj = Bio::Seq->new( -seq => "$baseList", -alphabet => 'dna' );
    if ($strand[$m] == -1) {
        $seq_obj = $seq_obj->revcom();
    }

    if ( defined $header ) {
        $header =~ s/^\s*//;
        $header =~ s/\s*$//;
        $header =~ s/^>//;
        $header =~ s/^([^\s]+)//;
        my $default_name = $1;
        $seq_obj->display_id($default_name);
        $seq_obj->desc($header);
    }
    else {
        $seq_obj->desc($default_header);
        $seq_obj->display_id($accession);
    }

    # Prints the sequence to STDOUT
    $out->write_seq($seq_obj);
    print "\n";

   # If 'exact_match' option is specified, exit after first (and only) iteration
    if ( defined($exact_match) ) {
        close INFILE;
        exit;
    }

    close INFILE;
}    # End of major for loop

__END__
