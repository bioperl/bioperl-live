#!perl
# -*-Perl-*-
#
# Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
# Finding flanking sequences for a variant.
#
#
#   v. 1     16 Mar 2001
#   v. 1.1    9 Aug 2001 interface change, more info in fasta header
#   v. 2.0   23 Nov 2001 new code from the flanks CGI program
#                         support for EMBL-like positions
#   v. 3.0   21 Feb 2003 new command line interface


use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::DB::EMBL;
use Bio::DB::GenBank;
use Getopt::Long;
use strict;
use warnings;

use constant VERSION => '3.0';

my $help = '';
my $flank = 100;              # flank length on both sides of the region
my $in_format = 'EMBL';       # format of the file to read in
my @pos;                      # position(s) in the sequence


GetOptions ("help" => \$help, "flanklength=i" => \$flank,
            "position=s" => \@pos );

@pos = split(/,/,join(',',@pos));

system("perldoc $0") if $help;
system("perldoc $0") unless @ARGV;
print STDERR "\nYou need to provide --position option\n" and system("perldoc $0") 
    unless @pos;

my $file = shift;
$file || system("perldoc $0");

my $seq = get_seq($file);
exit 0 unless $seq;

&extract($seq, \@pos, $flank);

#
## end main
#

sub get_seq {
    my ($file) = @_;
    my $IN_FORMAT = 'EMBL';	# format of the local file on disk

    if (-e $file ) {		# local file
	my $in  = Bio::SeqIO->new('-file' => $file,
				  '-format' => $IN_FORMAT);
	$seq = $in->next_seq();
    }
    elsif ($file =~ /\./) {	# sequence version from GenBank
	eval {
	    my $gb = new Bio::DB::GenBank;
	    $seq = $gb->get_Seq_by_version($file);
	};
    } else {			# plain accession mumber from more reliable EMBL
	eval {
	    my $gb = new Bio::DB::EMBL;
	    $seq = $gb->get_Seq_by_acc($file);
	};
	
    }
    print STDERR "Could not find sequence [$file]" && return unless $seq;
    return $seq;
}

sub extract {
    my ($seq, $pos, $flank) = @_;
    my ($out_seq);
    my $OUT_FORMAT = 'FASTA';	# output format, going into STDOUT
    my $strand = 1;		# default for the forward strand

    my $out = Bio::SeqIO->new('-format' => $OUT_FORMAT);

    my $count = 1;
    foreach my $idpos (@$pos) {

	my ($id, $pos_range, $start, $end, $allele_len);
	my $inbetween = 0;	# handle 23^24 notation as well as plain integer (24)
                          	# but set flag and make corrections when needed

	if ($idpos =~ /:/ ) {	# id and position separator
	    ($id, $pos_range) = split /:/, $idpos;
	} else {		# no id
	    $id = $count;
	    $count++;
	    $pos_range = $idpos;
	}
	$strand = -1 if $pos_range =~ /-$/; # opposite strand
	$pos_range = $1 if $pos_range =~ /(.+)-/; # remove trailing '-'

	if ($pos_range =~ /\^/) { # notation 23^24 used
	    ($start, $end) = split /\^/, $pos_range;
	    print STDERR $id, ": Give adjacent nucleotides flanking '^' character, not [",
	        $start, "] and [", $end, "]\n" and next
		unless $end == $start + 1;
	    $end = $start;
	    $inbetween = 1;
	} else {		#  notation 23..24 used
	    ($start, $end) = split /\.\./, $pos_range;
	}
	$end ||= $start;	# notation 23 used
	print STDERR $id, ": Start can not be larger than end. Not [",
	        $start, "] and [", $end, "]\n" and next
		if $start > $end;
	$allele_len = $end - $start;

	# sanity checks
	next unless defined $start && $start =~ /\d+/ && $start != 0;
	print STDERR "Position '$start' not in sequence '$file'\n",  and next
	    if $start < 1 or $start > $seq->length;
	print STDERR "End position '$end' not in sequence '$file'\n",  and next
	    if $end < 1 or $end > $seq->length;
	
	# determine nucleotide positions
	# left edge
	my $five_start = $start - $flank;
	$five_start = 1 if $five_start < 1; # not outside the sequence
	# right edge
	my $three_end = $start + $allele_len + $flank;
	$three_end = $seq->length if $start + $allele_len + $flank > $seq->length;
	$three_end-- if $inbetween;

	# extract the sequences
	my $five_prime = lc $seq->subseq($five_start , $start - 1); # left flank
	my $snp = uc $seq->subseq($start, $end); # allele (always > 0 length)
	$snp = lc $snp if $inbetween;

	my $three_prime;	# right flank
	if ($end < $seq->length) { # make sure we are not beyond reference sequece
	    $three_prime = lc $seq->subseq($end + 1, $three_end);
	} else {
	    $three_prime = '';
	}

	# allele positions in local, extracted coordinates
	my $locpos = length($five_prime) + 1;
	my $loc_end;
	if ($allele_len) {
	    $loc_end = "..". ($locpos+$allele_len);
	} else {
	    $loc_end = '';
	    $loc_end = '^'. ($locpos+1) if $inbetween;
	}
	# build FASTA id and description line
	my $fastaid = uc($id). "_". uc($file).
	    " oripos=$pos_range strand=$strand allelepos=$locpos$loc_end";

	#build BioPerl sequence objects
	if ($strand == -1) {
	    my $five_prime_seq = new Bio::PrimarySeq(-alphabet=>'dna',-seq=>$five_prime);
	    my $snp_seq = new Bio::PrimarySeq(-alphabet=>'dna',-seq=>$snp);
	    my $three_prime_seq = new Bio::PrimarySeq(-alphabet=>'dna',-seq=>$three_prime);

	    my $str = $three_prime_seq->revcom->seq. " ".
		$snp_seq->revcom->seq. " ". $five_prime_seq->revcom->seq;
	    $str =~ s/ //g;
	    $out_seq = new Bio::PrimarySeq (-id => $fastaid,
					    -alphabet=>'dna',
					    -seq => $str );
	} else {
	    my $str = $five_prime. " ". $snp. " ". $three_prime;
	    $str =~ s/ //g;
	    $out_seq = new Bio::PrimarySeq (-id => $fastaid,
					    -alphabet=>'dna',
					    -seq => $str );
	}
	$out->write_seq($out_seq); # print sequence out
    }
}



=head1 NAME

bp_flanks - finding flanking sequences for a variant in a sequence position

=head1 SYNOPSIS

  bp_flanks --position POS [-p POS ...]  [--flanklen INT]
         accession | filename

=head1 DESCRIPTION

This script allows you to extract a subsequence around a region of
interest from an existing sequence. The output if fasta formatted
sequence entry where the header line contains additional information
about the location.

=head1 OPTIONS

The script takes one unnamed argument which be either a file name in
the local file system or a nucleotide sequence accession number.


  -p         Position uses simple nucleotide sequence feature table
  --position notation to define the region of interest, typically a
             SNP or microsatellite repeat around which the flanks are
             defined.

             There can be more than one position option or you can
             give a comma separated list to one position option.

             The format of a position is:

                 [id:] int | range | in-between [-]

             The optional id is the name you want to call the new
             sequence. If it not given in joins running number to the
             entry name with an underscore.

             The position is either a point (e.g. 234), a range (e.g
             250..300) or insertion point between nucleotides
             (e.g. 234^235)

             If the position is not completely within the source
             sequence the output sequence will be truncated and it
             will print a warning.

             The optional hyphen [-] at the end of the position
             indicates that that you want the retrieved sequence to be
             in the opposite strand.


  -f         Defaults to 100. This is the length of the nucleotides
  --flanklen sequence retrieved on both sides of the given position.

             If the source file does not contain 

=head1 OUTPUT FORMAT

The output is a fasta formatted entry where the description file
contains tag=value pairs for information about where in the original
sequence the subsequence was taken.

The ID of the fasta entry is the name given at the command line joined
by hyphen to the filename or accesion of the source sequence. If no id
is given a series of consequtive integers is used.

The tag=value pairs are:

=over 3

=item oripos=int

position in the source file

=item strand=1|-1

strand of the sequence compared to the source sequence

=item allelepos=int

position of the region of interest in the current entry.
The tag is the same as used by dbSNP@NCBI

=back

The sequence highlights the allele variant position by showing it in
upper case and rest of the sequence in lower case characters.

=head1 EXAMPLE

  % bp_flanks ~/seq/ar.embl

  >1_/HOME/HEIKKI/SEQ/AR.EMBL oripos=100 strand=1 allelepos=100
  taataactcagttcttatttgcacctacttcagtggacactgaatttggaaggtggagga
  ttttgtttttttcttttaagatctgggcatcttttgaatCtacccttcaagtattaagag
  acagactgtgagcctagcagggcagatcttgtccaccgtgtgtcttcttctgcacgagac
  tttgaggctgtcagagcgct


=head1 TODO

The input files are assumed to be in EMBL format and the sequences are
retrieved only from the EMB database. Make this more generic and use
the registry.


head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  E<lt>heikki-at-bioperl-dot-orgE<gt>

=cut
