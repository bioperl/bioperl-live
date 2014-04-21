#!perl
use strict;
use warnings;

=head1 NAME

bp_split_seq - splits a sequence into equal sized chunks with an optional
            overlapping range

=head1 SYNOPSIS

bp_split_seq -c 10000 [-o 1000] [-i] -f seq.in

=head1 DESCRIPTION 

The script will split sequences into chunks

Mandatory Options:

  -c  Desired length of the resulting sequences.
  -f  Input file (must be FASTA format).

Special Options:

  -o  Overlapping range between the resulting sequences.
  -i  Create an index file with the resulting sequence files. This is
      useful if you want to pass this list as input arguments into
      another programs (i.e. CLUSTAL, HMMER, etc.).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHORS

  Ewan Birney E<lt>birney-at-ebi.ac.ukE<gt>
  Mauricio Herrera Cuadra E<lt>mauricio at open-bio.orgE<gt>
  (some enhancements)

=cut

# Modules, pragmas and variables to use
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use vars qw($opt_c $opt_o $opt_i $opt_f $index_file);

# Gets options from the command line
GetOptions qw(-c=i -o:i -i -f=s);

# If no mandatory options are given prints an error and exits
if (!$opt_c) {
    print "ERROR: No chunk size has been specified.\n" and exit();
} elsif (!$opt_f) {
    print "ERROR: No FASTA file has been specified.\n" and exit();
}

# Declares offset size
my $offset = $opt_o ? $opt_o : "0";

# Opens the FASTA file
my $in = Bio::SeqIO->new(
    -file   => "$opt_f",
    -format => "Fasta",
);
print "==> Opening FASTA file:\t\t\t\t$opt_f\n";

# Reads the next sequence object
while (my $seq = $in->next_seq()) {

    # Reads the ID for the sequence and prints it
    my $id = $seq->id();
    print "--> The ID for this sequence is:\t\t$id\n";

    # Reads the description for the sequence and prints it
    my $desc = $seq->desc();
    print "--> The description for this sequence is:\t$desc\n";

    # Gets sequence length and prints it
    my $seq_length = $seq->length();
    print "--> The length of this sequence is:\t\t$seq_length\n";

    # If the chunk size is bigger than the sequence length prints the error and exits
    (print "ERROR: Specified chunk size is bigger than sequence length.\n" and exit()) if ($opt_c > $seq_length);

    # Creates a directory for writing the resulting files
    mkdir("split", 0755) unless -e "split" and -d "split";

    # Creates the INDEX file if the option was given
    my $FH;
    if ($opt_i) {
        $index_file = "$id.c$opt_c.o$offset.INDEX";
        open $FH, '>', $index_file or die "Could not write file '$index_file': $!\n";
    }

    # Loops through the sequence
    for (my $i = 1; $i < $seq_length; $i += $opt_c) {
        my $end = (($i + $opt_c) > $seq_length) ? ($seq_length + 1) : ($i + $opt_c);
        my $seq_range = (($i + $opt_c) > $seq_length) ? "$i-".($end - 1) : "$i-$end";
        my $id = $seq->id();
        $id .= "_$seq_range";

        # Stores chunk into its corresponding FASTA file
        my $out = Bio::SeqIO->new(
            -file   => ">split/$id.faa",
            -format => "Fasta",
        );
        my $trunc_seq = $seq->trunc($i, $end - 1);
        $trunc_seq->id($id);
        $out->write_seq($trunc_seq);
        print "==> Sequence chunk:\t$seq_range\tstored in file:\tsplit/$id.faa\n";

        # Prints the current file name into the INDEX file if the option was given
        print $FH "split/$id.faa\n" if $opt_i;

        # Decreases the $i value with the offset value
        $i -= $offset;
    }

    # Closes the INDEX file if the option was given
    if ($opt_i) {
        print "==> INDEX stored in file:\t\t\t$index_file\n";
        close $FH;
    }
}

__END__
