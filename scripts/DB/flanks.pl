#!/usr/bin/perl -- # -*-Perl-*-
#
# $Id$
#
# Heikki Lehvaslaiho <heikki@ebi.ac.uk>
# Finding flanking sequences for a variant.
# 
# Requires BioPerl >  0.7 
#
#   v. 1     16 Mar 2001
#   v. 1.1    9 Aug 2001 interface change, more info in fasta header 
#   v. 2.0   23 Nov 2001 new code from the flanks CGI program
#                         support for EMBL-like positions


use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::DB::EMBL;
use Bio::DB::GenBank;
use strict;

use constant VERSION => '2.0';

my $flank;                    # Will default to 100
my ($file, $pos, $strand);
my $in_format = 'EMBL';       # format of the local file on disk
my $out_format = 'FASTA';     # output format, going into STDOUT

if (@ARGV != 2 and @ARGV != 3 and @ARGV != 4) {
    print "Usage: flanks (accession or filename) position [length]\n";
    exit;
}

$file = shift;
$pos = shift;
$flank = shift;

$flank ||= 100;          # defaults to 100

&extract($file, $pos, $flank);

#
## end main
#

sub extract {
    my ($file, $ids, $flank) = @_;
    my ($seq, $out_seq);
    my $IN_FORMAT = 'EMBL';	# format of the local file on disk
    my $OUT_FORMAT = 'FASTA';	# output format, going into STDOUT
    my $strand = 1;		# default for the forward strand


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

    my $out = Bio::SeqIO->new('-format' => $OUT_FORMAT);
    
    my $count = 1;
    foreach my $idpos (split /\s+/, $ids) {

	my ($id, $pos_range, $start, $end, $allele_len);
	my $inbetween = 0;	# handle 23^24 notation like plain integer (24) 
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

	# extract the sequnces
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
