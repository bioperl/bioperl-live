=encoding utf8

=head1 NAME

Bio::BPWrapper::AlnManipulations - Functions for L<bioaln>

=head1 SYNOPSIS

    use Bio::BPWrapper::AlnManipulations;
    # Set options hash ...
    initialize(\%opts);
    write_out(\%opts);

=cut

package Bio::BPWrapper::AlnManipulations;

use strict;
use warnings;
use 5.010;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Data::Dumper;
use List::Util qw(shuffle);
use Bio::Align::Utilities qw(:all);
use Exporter ();
use Bio::SearchIO;
#use Bio::Tools::GuessSeqFormat;

if ($ENV{'DEBUG'}) { use Data::Dumper }

use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA         = qw(Exporter);

# FIXME: some of these have too generic names like
# "upper_case" or "concat".
@EXPORT      = qw(initialize can_handle handle_opt write_out write_out_paml
phylip_non_interleaved split_cdhit trim_ends gap_states
gap_states_matrix print_avp_id bootstrap draw_codon_view del_seqs
remove_gaps print_length print_match print_num_seq pick_seq
change_ref aln_slice get_uniq binary_informative variables_sites
avg_id_by_win concat conserve_blocks get_consensus dns_to_protein
remove_gapped_cols_in_one_seq colnum_from_residue_pos
list_ids premute_states protein_to_dna sample_seqs
shuffle_sites random_slice select_third_sites remove_third_sites
upper_case );

use Bio::BPWrapper;
# Package global variables
my ($in, $out, $aln, %opts, $file, $in_format, $out_format, @alns, $binary);

my $VERSION = $Bio::BPWrapper::VERSION;

## For new options, just add an entry into this table with the same
## key as in the GetOpts function in the main program. Make the key be
## a reference to the handler subroutine (defined below), and test
## that it works.
my %opt_dispatch = (
    "avg-pid" => \&print_avp_id,
    "boot" => \&bootstrap,
    "codon-view" => \&draw_codon_view,
    "delete" => \&del_seqs,
    "gap-char" => \&gap_char,
    "no-gaps" => \&remove_gaps,
    "length" => \&print_length,
    "match" => \&print_match,
    "num-seq" => \&print_num_seq,
    "pick" => \&pick_seqs,
    "ref-seq" => \&change_ref,
    "slice" => \&aln_slice,
    "split-cdhit" => \&split_cdhit,
    "uniq" => \&get_unique,
    "var-sites" => \&variable_sites,
    "window" => \&avg_id_by_win,
    "concat" => \&concat,
    "con-blocks" => \&conserved_blocks,
    "consensus" => \&get_consensus,
    "dna2pep" => \&dna_to_protein,
    "rm-col" => \&remove_gapped_cols_in_one_seq,
    "aln-index" => \&colnum_from_residue_pos,
    "list-ids" => \&list_ids,
    "pair-diff" => \&pair_diff,
    "pair-diff-ref" => \&pair_diff_ref,
    "permute-states" => \&permute_states,
    "pep2dna" => \&protein_to_dna,
    "resample" => \&sample_seqs,
    "shuffle-sites" => \&shuffle_sites,
    "select-third" => \&select_third_sites,
    "remove-third" => \&remove_third_sites,
    "random-slice" => \&random_slice,
    "upper" => \&upper_case,
    "gap-states" => \&gap_states,
    "gap-states2" => \&gap_states_matrix,
    "trim-ends" => \&trim_ends,
    "bin-inform" => \&binary_informative,
    "bin-ref" => \&binary_ref,
    "phy-nonint" => \&phylip_non_interleaved
   );


##################### initializer & option handlers ###################

## TODO Formal testing!

=head1 SUBROUTINES

=head2 initialize()

Sets up most of the actions to be performed on an alignment.

Call this right after setting up an options hash.

Sets package variables: C<$in_format>, C<$binary>, C<$out_format>, and C<$out>.

=cut

sub initialize {
    my $opts_ref = shift;
    Bio::BPWrapper::common_opts($opts_ref);
    %opts = %{$opts_ref};

    # This is the format that aln-manipulations expects by default
    my $default_format = "clustalw";

    # assume we're getting input from standard input

#    my $in_format = $opts{"input"} || $default_format;
#    my $in_format;
#    use IO::Scalar;
#    my $s;
#    my ($guesser);
#    if ($file eq "STDIN") {
#	my $line_ct = 0; 
#	my $lines;
#	while(<>) { $lines .= $_; $line_ct++; last if $line_ct >= 100 } # read the first 100 lines
#	$guesser = Bio::Tools::GuessSeqFormat->new( -text => $lines );
#   } else {
#	open $ifh, "<", $file or die $!;
#	$guesser = Bio::Tools::GuessSeqFormat->new( -file => $file );
#    }
#    $in_format  = $guesser->guess();
#    die "unknown file format. Try specify with -i flag.\n" unless $in_format;
#    seek (STDIN, 0, 0);
#    warn "$in_format\n";

    my $in_format = $opts{'input'} || 'clustalw';
    if ($opts{"concat"}) {
#	foreach my $file (glob @ARGV) {
	while ($file = shift @ARGV) {
#	    warn "reading $file\n";
#	       $guesser = Bio::Tools::GuessSeqFormat->new( -file => $file);
#	       $in_format  = $guesser->guess;
	       $in = Bio::AlignIO->new(-file => $file, -format => $in_format);
	       while ($aln=$in->next_aln()) { push @alns, $aln }
	}
    } else {
	$file = shift @ARGV || "STDIN";    # If no more arguments were given on the command line
	if ($in_format && $in_format =~ /blast/) { # guess blastoutput as "phylip", so -i 'blast' is needed
#	if ($opts{"input"} && $opts{"input"} =~ /blast/) { # "blastxml" (-outfmt 5 ) preferred
	    my $searchio = Bio::SearchIO->new( -format => 'blast', ($file eq "STDIN")? (-fh => \*STDIN) : (-file => $file)); # works for regular blast output
#	    my $searchio = Bio::SearchIO->new( -format => 'blast', -fh => $ifh);
	    while ( my $result = $searchio->next_result() ) {
		while( my $hit = $result->next_hit ) {
 		    my $hsp = $hit->next_hsp; # get first hit; others ignored
		    $aln = $hsp->get_aln();
		}
	    }
	} else { # would throw error if format guessed wrong
#	    $in = Bio::AlignIO->new(-format => $in_format, ($file eq "STDIN")? (-fh => \*STDIN) : (-file => $file));
#	    $in = Bio::AlignIO->new(-format => $in_format, -fh => $ifh);
	    $in = Bio::AlignIO->new(-format=>$in_format, ($file eq "STDIN")? (-fh => \*STDIN) : (-file => $file) );
	    $aln = $in->next_aln()
	}
    }
    
    $binary = $opts{"binary"} ? 1 : 0;
    
    #### Options which *require an output FH* go *after* this ####
    $out_format = $opts{"output"} || $default_format;
    $out = Bio::AlignIO->new(-format => $out_format, -fh => \*STDOUT) unless $out_format eq 'paml'
}

sub can_handle {
    my $option = shift;
    return defined($opt_dispatch{$option})
}

sub handle_opt {
    my $option = shift;
    # This passes option name to all functions
    $opt_dispatch{$option}->($option)
}

=head2 write_out_paml()

Writes output in PAML format.

=cut

sub write_out_paml() {
    my @seq;
    my $ct=0;

    foreach my $seq ($aln->each_seq()) {
        my $id = $seq->display_id();
        if ($seq->seq() =~ /^-+$/) { print STDERR "all gaps: $file\t$id\n"; next }
        $ct++;
        push @seq, $seq
    }

    die "No computable sequences: less than 2 seq.\n" unless $ct >= 2;
    print $ct, "\t", $aln->length(), "\n";
    foreach (@seq) {
	   print $_->display_id(), "\n";
	   print $_->seq(), "\n"
    }
}

=head2 write_out()

Performs the bulk of the alignment actions actions set via
L<C<initialize(\%opts)>|/initialize> and calls
L<C<$AlignIO-E<gt>write_aln()>|https://metacpan.org/pod/Bio::AlignIO#write_aln>
or L<C<write_out_paml()>|/write_out_paml>.

Call this after calling C<#initialize(\%opts)>.

=cut

sub write_out($) {

    my $opts = shift;
    for my $option (keys %{$opts}) {
	next if ($option eq 'input') || ($option eq 'output') || ($option eq 'noflatname') || ($option eq 'binary'); # Don't process these options: they are for AlignIO

	if (can_handle($option)) { handle_opt($option) } # If there is a function to handle the current option, execute it
	else { warn "Missing handler for: $option\n" }
    }


    $aln->set_displayname_flat() unless $opts{"no-flat"};
    if ($out_format eq 'paml') { &write_out_paml($aln) }
    else { $out->write_aln($aln) }
}

sub phylip_non_interleaved {
    my @seq;
    my $ct=0;

    foreach my $seq ($aln->each_seq()) {
        my $id = $seq->display_id();
        if ($seq->seq() =~ /^-+$/) { print STDERR "all gaps: $file\t$id\n"; next }
        $ct++;
        push @seq, $seq
    }

    die "No computable sequences: less than 2 seq.\n" unless $ct >= 2;
    print "\t", $ct, "\t", $aln->length(), "\n";
    foreach (@seq) {
	   printf "%-50s", $_->display_id();
	   print $_->seq(), "\n"
    }
    exit;
}

###################### subroutine ######################

sub gap_char {
    my $char = $opts{'gap-char'};
    die "gap-char takes a single character\n" unless length($char) == 1;
    foreach my $seq ($aln->each_seq()) { 
	my $seq_str = $seq->seq();
	$seq_str =~ s/[\.-]/$char/g;
	$seq->seq($seq_str); 
    }
}

sub pair_diff_ref {
    my $refId = $opts{'pair-diff-ref'};
    my (@seqs, $refSeq);
    foreach my $seq ($aln->each_seq()) { 
	if ($seq->id eq $refId) {
	    $refSeq = $seq;
	} else {
	    push @seqs, $seq;
	}
    }

    die "ref seq $refId not found\n" unless $refSeq;

    @seqs = sort { $a->id() cmp $b->id() } @seqs;
    for (my $i=0; $i < $#seqs; $i++) {
	my $idB = $seqs[$i]->id();
	my $seqB = $seqs[$i];
	my $pair = new Bio::SimpleAlign;
	$pair->add_seq($refSeq);
	$pair->add_seq($seqB);
#	$pair = $pair->remove_gaps(); # not relialbe, e.g., n/N for DNA seqs
	my $ct_diff = 0;
	my $ct_valid = 0;
#	my $matchLine = $pair->match_line();
#	my @match_symbols = split //, $matchLine;
	for (my $j = 1; $j <= $pair->length; $j++) {
	    my $mA = $refSeq->subseq($j,$j);
	    my $mB = $seqB->subseq($j,$j);
	    next if $refSeq->alphabet eq 'dna' && $seqB->alphabet eq 'dna' && ($mA !~ /^[ATCG]$/i || $mB !~ /^[ATCG]$/i);
#		next if $match_symbols[$i] eq '*';
	    next if $mA eq '-' || $mB eq '-'; 
	    $ct_valid++;
	    $ct_diff++ unless $mA eq $mB;
#	    next if $match_symbols[$i] eq '*'; 
	}
	my $pairdiff = $pair->percentage_identity();
	print join "\t", ($refId, $idB, $ct_diff, $ct_valid, $pair->length());
	printf "\t%.4f\t%.4f\n", $pairdiff, 1-$pairdiff/100;
    }
    exit;
}

sub pair_diff {
    my $alnBack = $aln;
    my $lenTotal = $alnBack->length();
#    $alnBack = $alnBack->remove_gaps();
    my $matchLineFull = $alnBack->match_line();
    my @match_symbols_full = split //, $matchLineFull;
    my $num_var = 0; # contain gaps
#    my $num_var = 0; # de-gapped variable sites
    for (my $i = 0; $i < $alnBack->length; $i++) {
	next if $match_symbols_full[$i] eq '*'; 
	$num_var++;
    }

    my (@seqs);
    print join "\t", qw(seq_1 seq_2 len_gapless len_variable diff_gapless diff_variable len_aln identitty fac_diff frac_gapless frac_variable);
    print "\n";
    foreach my $seq ($aln->each_seq()) { push @seqs, $seq }
    @seqs = sort { $a->id() cmp $b->id() } @seqs;
    for (my $i=0; $i < $#seqs; $i++) {
	my $idA = $seqs[$i]->id();
	my $seqA = $seqs[$i];
	for (my $j=$i+1; $j <= $#seqs; $j++) {
	    my $idB = $seqs[$j]->id();
	    my $seqB = $seqs[$j];
	    my $pair = new Bio::SimpleAlign;
	    $pair->add_seq($seqA);
	    $pair->add_seq($seqB);
#	    $pair = $pair->remove_gaps(); # unreliable: e.g., "n", "N" in DNA seqs
#	    my $mask = $seqA->seq ^ $seqB->seq; #  (exclusive or) operator: returns "\0" if same
	    my $ct_diff = 0;
	    my $ct_valid = 0;
	    my $gap_included_diff = 0;
#	    my $matchLine = $pair->match_line();
#	    my @match_symbols = split //, $matchLine;
	    for (my $j = 1; $j <= $pair->length; $j++) {
		my $mA = $seqA->subseq($j,$j);
		my $mB = $seqB->subseq($j,$j);
		next if $seqA->alphabet eq 'dna' && $seqB->alphabet eq 'dna' && ($mA !~ /^[ATCG]$/i || $mB !~ /^[ATCG]$/i);
		$gap_included_diff++ unless $mA eq $mB;
#		next if $match_symbols[$i] eq '*'; 
		next if $mA eq '-' || $mB eq '-';
		$ct_valid++;
		$ct_diff++ unless $mA eq $mB;
	    }
#	    while ($mask =~ /[^\0]/g) { $ct_diff++ }
	    my $pairdiff = $pair->percentage_identity();
	    print join "\t", ($idA, $idB, $ct_valid, $num_var, $ct_diff, $gap_included_diff, $pair->length());
	    printf "\t%.4f\t%.4f\t%.4f\t%.4f\n", $pairdiff, 1-$pairdiff/100, $ct_diff/$ct_valid, $gap_included_diff/$num_var;
	}
    }
    exit;
}


sub split_cdhit {
    my $cls_file = $opts{'split-cdhit'};
    open IN, "<" . $cls_file || die "cdhit clstr file not found: $cls_file\n";
    my %clusters;
    my $cl_id;
    my @mem;
    while (<IN>) {
	my $line = $_;
	chomp $line;
	if ($line =~ /^>(\S+)\s+(\d+)/) {
	    $cl_id = $1 . "_" . $2;
	    my @mem = ();
	    $clusters{$cl_id} = \@mem;
	} else {
	    my $ref = $clusters{$cl_id};
	    my @mems = @$ref;
	    my @els = split /\s+/, $line;
	    my $seq_id = $els[2];
	    $seq_id =~ s/>//;
	    $seq_id =~ s/\.\.\.$//;
	    push @mems, $seq_id;
	    $clusters{$cl_id} = \@mems;
	}
    }
#    print Dumper(\%clusters);
    my %seqs;
    foreach my $seq ($aln->each_seq()) {
        my $id = $seq->display_id();
	$seqs{$id} = $seq;
    }

    foreach my $id (keys %clusters) {
	my $out = Bio::AlignIO->new( -file => ">" . $file . "-". $id . ".aln", -format => 'clustalw');
	my $new_aln = Bio::SimpleAlign->new();
	my @seqids = @{ $clusters{$id} };
	foreach my $seq_id (@seqids) {
	    $new_aln->add_seq($seqs{$seq_id});
	}
	$new_aln->set_displayname_flat();
#	$new_aln = &_remove_common_gaps($new_aln);
	$out->write_aln($new_aln);
    }
    exit;
}

#sub _remove_common_gaps {
#}

sub trim_ends {
    my (@seqs, @gaps);
    foreach my $seq ($aln->each_seq()) {
        my $id = $seq->display_id();
        my @nts = split //, $seq->seq();
	my $gap_start = 0;
	my $new;
	for (my $i=0; $i< $aln->length(); $i++) {
	    if ($nts[$i] eq '-') {
		if ($gap_start) { # gap -> gap
		    $new->{end}++;
		} else { # nt -> gap
		    $gap_start = 1;
		    $new = { 'start' => $i+1, 'end' => $i+1, 'seq_name' => $id }
		}
	    } else {
		if ($gap_start) { # gap -> nt
		    $gap_start = 0;
		    push @gaps, $new;
		} else { # nt -> nt
		    next;
		}
	    }
	}
	push @gaps, $new if $gap_start;
    }

    my (@three_end_gaps, @five_end_gaps);

    foreach my $gap (@gaps) {
	$gap->{length} = $gap->{end} - $gap->{start} + 1;
	push @three_end_gaps, $gap if $gap->{start} == 1;
	push @five_end_gaps, $gap if $gap->{end} == $aln->length;
    }

    return unless @three_end_gaps or @five_end_gaps;

    my $longest_three_end = 0;
    my $longest_five_start = 0;
    my $longest_three_length = 0;
    my $longest_five_length = 0;

    foreach my $gap (@three_end_gaps) {
	if ($gap->{length} > $longest_three_length) {
	    $longest_three_end = $gap->{end};
	    $longest_three_length = $gap->{length};
	}
    }

    foreach my $gap (@five_end_gaps) {
	if ($gap->{length} > $longest_five_length) {
	    $longest_five_start = $gap->{start};
	    $longest_five_length = $gap->{length};
	}
    }

#    print STDERR $longest_three, "\t", $longest_five, "\n";
    if (@three_end_gaps) {
	print STDERR Dumper(\@three_end_gaps);
	print STDERR $longest_three_end, "\n";
	$aln = $aln->slice($longest_three_end + 1, $aln->length);
    }

    if (@five_end_gaps) {
	print STDERR Dumper(\@five_end_gaps);
	print STDERR $longest_five_start, "\n";
	$aln = $aln->slice(1, $longest_five_start - $longest_three_end - 1);
    }
}

sub gap_states {
    my (@seqs, @gaps);
    foreach my $seq ($aln->each_seq()) {
        my $id = $seq->display_id();
        my @nts = split //, $seq->seq();
	my $gap_start = 0;
	my $new;
	for (my $i=0; $i< $aln->length(); $i++) {
	    if ($nts[$i] eq '-') {
		if ($gap_start) { # gap -> gap
		    $new->{end}++;
		} else { # nt -> gap
		    $gap_start = 1;
		    $new = { 'start' => $i+1, 'end' => $i+1, 'seq_name' => $id }
		}
	    } else {
		if ($gap_start) { # gap -> nt
		    $gap_start = 0;
		    push @gaps, $new;
		} else { # nt -> nt
		    next;
		}
	    }
	}
	push @gaps, $new if $gap_start;
    }
    my (%gap_freqs, @uniq_gaps);
    foreach my $gap (@gaps) {
	my $id = $gap->{start} . "-" . $gap->{end};
	$gap->{id} = $id;
	$gap_freqs{$id}++;
    }

    foreach my $id (keys %gap_freqs) {
	my ($start, $end) = split /-/, $id;
	push @uniq_gaps, {
	    'start' => $start,
	    'end' => $end,
	    'is_edge' => ($start == 1 || $end == $aln->length) ? 1 : 0,
	    'in_frame' => ($end - $start + 1) % 3 ? 0 : 1,
	    'counts' => $gap_freqs{$id},
	};
    }

    foreach my $gap (@uniq_gaps) { say join "\t", ($file, $gap->{start}, $gap->{end}, $gap->{is_edge}, $gap->{in_frame}, $gap->{counts}, $aln->length()) }
#    print Dumper(\@uniq_gaps);
    exit;
}

sub gap_states_matrix {
    my (@seq_ids, @gaps);
    foreach my $seq ($aln->each_seq()) {
        my $id = $seq->display_id();
	push @seq_ids, $id;
        my @nts = split //, $seq->seq();
	my $gap_start = 0;
	my $new;
	for (my $i=0; $i< $aln->length(); $i++) {
	    if ($nts[$i] eq '-') {
		if ($gap_start) { # gap -> gap
		    $new->{end}++;
		} else { # nt -> gap
		    $gap_start = 1;
		    $new = { 'start' => $i+1, 'end' => $i+1, 'seq_name' => $id }
		}
	    } else {
		if ($gap_start) { # gap -> nt
		    $gap_start = 0;
		    push @gaps, $new;
		} else { # nt -> nt
		    next;
		}
	    }
	}
	push @gaps, $new if $gap_start;
    }
    my (%gap_freqs, @uniq_gaps, %gap_presence);
    foreach my $gap (@gaps) {
	my $id = $gap->{start} . "-" . $gap->{end};
	$gap->{id} = $id;
	$gap_freqs{$id}++;
	$gap_presence{$id}->{$gap->{seq_name}} = 1;
    }

    foreach my $id (keys %gap_freqs) {
	my ($start, $end) = split /-/, $id;
	push @uniq_gaps, {
	    'start' => $start,
	    'end' => $end,
	    'is_edge' => ($start == 1 || $end == $aln->length) ? 1 : 0,
	    'in_frame' => ($end - $start + 1) % 3 ? 0 : 1,
	    'counts' => $gap_freqs{$id},
	    'id' => $id
	};
    }

    my @gaps_sorted = sort {$a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}} @uniq_gaps;
    foreach (@gaps_sorted) { print "\t", $_->{id}}
    print "\n";
    foreach my $sid (sort @seq_ids) {
	print $sid;
	foreach my $u_gap (@gaps_sorted) {
	    print "\t", $gap_presence{$u_gap->{id}}->{$sid} || 0;
	}
	print "\n";
    }

#    foreach my $gap (@uniq_gaps) { say join "\t", ($file, $gap->{start}, $gap->{end}, $gap->{is_edge}, $gap->{in_frame}, $gap->{counts}, $aln->length()) }
#    print Dumper(\@uniq_gaps);
    exit;
}

=head2 print_avp_id

Print the average percent identity of an alignment.

Wraps
L<Bio::SimpleAlign-E<gt>average_percentage_identity()|https://metacpan.org/pod/Bio::SimpleAlign#average_percentage_identity>.


=cut

sub print_avp_id {
    printf "%.4f\n", $aln->average_percentage_identity();
    exit
}

=head2 boostrap()

Produce a bootstrapped
alignment. L<Bio::Align::Utilities-E<gt>bootstrap()|https://metacpan.org/pod/Bio::Align::Utilities#bootstrap_replicates>.

=cut

sub bootstrap {
    my $replicates = bootstrap_replicates($aln,1);
    $aln = shift @$replicates
}

=head2 draw_codon_view()

Print a CLUSTALW-like alignment, but separated by codons. Intended for
use with DNA sequences. Block-final position numbers are printed at
the end of every alignment block at the point of wrapping, and
block-initial counts appear over first nucleotide in a block.

=cut

sub draw_codon_view {
#    my $aln = shift;
    # Is 20 by default. Blocks are measured in CODONS, so mult by 3
    my $block_length = 3 * $opts{"codon-view"};
    my $aln_length   = $aln->length();
    my $num_seqs     = $aln->num_sequences();
    my $min_pad = 4;    # Minimum padding between sequence and ID
    my $seq_matrix;
    my @seqs = ($aln->each_seq);
    my @display_ids;

    # Find longest id length, add id/sequence padding
    my $max_id_len = _find_max_id_len(\@seqs);

    # id length includes padding
    $max_id_len += $min_pad;

    # Extract display_ids and sequences from AlignIO object.
    foreach my $seq (@seqs) {
        my @seq_str = split '', $seq->seq();
        push @$seq_matrix, \@seq_str;
        push @display_ids, $seq->display_id;

       # Pad display ids so that space between them and sequence is consistent
        $display_ids[-1] = _pad_display_id($display_ids[-1], $max_id_len)
    }

    my $nuc_count = 0;

    # Loop over each sequence.
    for (my $i = 0; $i < $num_seqs; $i++) {

        # Print count at end of block when we are starting out a new block
        _print_positions($nuc_count, $aln_length, $max_id_len) if $i == 0;

        # Loop over nucleotides
        for (my $j = $nuc_count; $j < $aln_length; $j++) {

            # When we're starting, or starting a new block, print the display id's.
            print $display_ids[$i] if $j % $block_length == 0;

            print "$$seq_matrix[$i]->[$j]";
            print " " if ((($j + 1) % 3) == 0);

            # When we've reached the end of the alignment or a block
            if ($j + 1 == $aln_length || (($j + 1) % $block_length) == 0) {
                if ($i + 1 == $num_seqs) { $nuc_count = $j + 1 }  # If this is the last sequence, save the ending (next) position.
                else { print "\n" } # Otherwise, start on the next line.
                last  # In either case, need to exit this loop.
            }
        }    # END for LOOP OVER NUCLEOTIDES

        # Finish if we've reached the end of the alignment, and the last sequence
        if (($i + 1 == $num_seqs) && ($nuc_count == $aln_length)) { print "\n"; last }

      # If we haven't reached the end of the alignment, but we've run through
      # all sequences, print final block position and start at first sequence.
        elsif (($i + 1 == $num_seqs) && ($nuc_count < $aln_length)) {
            $i = -1;  # Always increments after a loop; next increment sets to 0.
            print "\n\n"
        }
    }    # END for LOOP OVER SEQUENCES

  # Can't let script terminate normally: produces traditional alignment output
    exit 0
}

=head2 del_seqs()

Delete sequences based on their id. Option takes a comma-separated list of ids.
The list of sequences to delete is in C<$opts{"delete"}> which is set via
L<C<#initilize(\%opts)>|/initialize>

=cut

sub del_seqs {
    _del_or_pick($opts{"delete"}, "remove_seq", 0)
}

=head2 remove_gaps()

Remove gaps (and returns an de-gapped alignment). Wraps
L<Bio::SimpleAlign-E<gt>remove_gaps()|https://metacpan.org/pod/Bio::SimpleAlign#remove_gaps>.

=cut

sub remove_gaps {
    $aln = $aln->remove_gaps()
}

=head2 print_length()

Print alignment length. Wraps L<Bio::SimpleAlign-E<gt>length()|https://metacpan.org/pod/Bio::SimpleAlign#length>.

=cut

sub print_length {
    say $aln->length();
    exit
}

=head2 print_match()

Go through all columns and change residues identical to the reference
sequence to be the match character, '.' Wraps
L<Bio::SimpleAlign-E<gt>match()|https://metacpan.org/pod/Bio::SimpleAlign#match>.

=cut

sub print_match {
    $aln->match()
}

=head2 print_num_seq()

Print number of sequences in alignment.

=cut

sub print_num_seq {
    say $aln->num_sequences();
    exit
}

=head2 pick_seqs()

Pick sequences based on their id. Option takes a comma-separated list
of ids.  The sequences to pick set in C<$opts{"pick"}> which is set via
L<C<#initilize(\%opts)>|/initialize>.

=cut

sub pick_seqs {
    _del_or_pick($opts{"pick"}, "add_seq", 1)
}

=head2 change_ref()

Change the reference sequence to be what is in C<$opts{"refseq"}>
which is set via L<C<#initilize(\%opts)>|/initialize>. Wraps
L<Bio::SimpleAlign-E<gt>set_new_reference()|https://metacpan.org/pod/Bio::SimpleAlign#set_new_reference>.

=cut

sub change_ref {
    my @newAlns;
    if ($opts{'concat'}) {
	foreach (@alns) {
	    push @newAlns, $_->set_new_reference($opts{"ref-seq"})
	}
	@alns = @newAlns;
    } else {
	$aln = $aln->set_new_reference($opts{"ref-seq"})
    }
}


=head2 aln_slice()

Get a slice of the alignment.  The slice is specified
C<$opts{"slice"}> which is set via L<C<#initilize(\%opts)>|/initialize>.

Wraps
L<Bio::SimpleAlign-E<gt>slice()|https://metacpan.org/pod/Bio::SimpleAlign#slice>
with improvements.

=cut


sub aln_slice {    # get alignment slice
    my ($begin, $end) = split(/\s*,\s*/, $opts{"slice"});

    # Allow for one parameter to be omitted. Default $begin to the
    # beginning of the alignment, and $end to the end.
    $begin = 1            if $begin eq "-";
    $end   = $aln->length if $end   eq "-";
    $aln = $aln->slice($begin, $end)
}

=head2 get_unique()

Extract the alignment of unique sequences. Wraps
L<Bio::SimpleAlign-E<gt>uniq_seq()|https://metacpan.org/pod/Bio::SimpleAlign#uniq_seq>.

=cut


sub get_unique() {
    $aln->verbose(1);
    $aln = $aln->uniq_seq();
}

sub _has_gap {
    my $ref = shift;
    foreach (@$ref) {
	return 1 if $_ eq '-';
    }
    return 0;
}

sub _has_singleton {
    my $ref = shift;
    foreach my $key (keys %$ref) {
	return 1 if $ref->{$key} == 1;
    }
    return 0;
}

=head2 binary_informative

extract binary and informative sites (for clique): discard constant,
3/4-states, non-informative

=cut

sub binary_informative {
    my $new_aln = Bio::SimpleAlign->new();
    my $len=$aln->length();
    my (@seq_ids, @inf_sites, %bin_chars);

    # Go through each column and save variable sites
    my $ref_bases = &_get_a_site_v2(); #print Dumper($ref_bases); exit;
    foreach (sort keys %$ref_bases) { push @seq_ids, $_ }
    for (my $i=1; $i<=$len; $i++) {
	my (%seen, @bases);
	foreach my $id (@seq_ids) { push @bases, $ref_bases->{$id}->{$i}; }
	%seen = %{&_seen_bases(\@bases)};
	next if &_has_gap( [ values %seen ] );
	next if keys %seen != 2;
	next if &_has_singleton(\%seen);
	my ($base1, $base2) = sort keys %seen;
	$bin_chars{$i}{$base1} = 0;
	$bin_chars{$i}{$base2} = 1;
	push @inf_sites, $i;
    }

    die "informative sites not found\n" unless @inf_sites;
    foreach (@inf_sites) { warn $_, "\n" }

    # Recreate the object for output
    foreach my $id (@seq_ids) {
        my $seq_str;
        foreach my $i (@inf_sites) {
            $seq_str .= $binary ? $bin_chars{$i}->{$ref_bases->{$id}->{$i}} : $ref_bases->{$id}->{$i};
        }
        my $loc_seq = Bio::LocatableSeq->new(-seq => $seq_str, -id => $id, -start => 1);
        my $end = $loc_seq->end;
        $loc_seq->end($end);
        $new_aln->add_seq($loc_seq)
    }

    $aln = $new_aln
}

sub binary_ref {
    my $new_aln = Bio::SimpleAlign->new();
    my $len=$aln->length();
    my (@seq_ids, @inf_sites, %bin_chars, @ref_states);
    my $refId = $opts{'bin-ref'} || die "need ref id as an argument\n";
    # Go through each column and save variable sites
    my $ref_bases = &_get_a_site_v2(); #print Dumper($ref_bases); exit;
    my $seenRef = 0;
    foreach (sort keys %$ref_bases) { 
	push @seq_ids, $_;
	next unless $_ eq $refId;
	$seenRef++;
    }
    die "ref seq not found: $refId\n" unless $seenRef;

    for (my $i=1; $i<=$len; $i++) {
	my (%seen, @bases);
	foreach my $id (@seq_ids) { push @bases, $ref_bases->{$id}->{$i}; }
	%seen = %{&_seen_bases(\@bases)};
	next if &_has_gap( [ values %seen ] ); # skip gaps
	next if keys %seen != 2; # skip multi-states or constant sites
	my ($base1, $base2) = sort keys %seen;
	if ($base1 eq $ref_bases->{$refId}->{$i}) { # base1 is ref
	    $bin_chars{$i}{$base1} = 1;
	    $bin_chars{$i}{$base2} = 0;
	} else { # base2 is ref
	    $bin_chars{$i}{$base1} = 0;
	    $bin_chars{$i}{$base2} = 1;
	}
	push @inf_sites, $i;
    }

    die "no binary sites\n" unless @inf_sites;
    foreach (@inf_sites) { warn $_, "\n" }

    # Recreate the object for output
    foreach my $id (@seq_ids) {
        my $seq_str;
        foreach my $i (@inf_sites) {
            $seq_str .= $bin_chars{$i}->{$ref_bases->{$id}->{$i}};
        }
        my $loc_seq = Bio::LocatableSeq->new(-seq => $seq_str, -id => $id, -start => 1);
        my $end = $loc_seq->end;
        $loc_seq->end($end);
        $new_aln->add_seq($loc_seq)
    }
    $aln = $new_aln
}


=head2 variable_sites()

Extracts variable sites.

=cut
sub variable_sites {
    $aln = $aln->remove_gaps();
    my $new_aln = Bio::SimpleAlign->new();
    my $len=$aln->length();
    my (%seq_ids, @sites, @var_sites);

    # Go through each column and save variable sites
    for (my $i=1; $i<=$len; $i++) {
        my ($ref_bases, $ref_ids) = &_get_a_site($i);
        %seq_ids = %{$ref_ids};
        my $is_constant = &_is_constant(&_paste_nt($ref_bases));
        if ($is_constant < 1) { push @sites, $ref_bases; push @var_sites, $i }
    }

    foreach (@var_sites) { warn $_, "\n" }

    # Recreate the object for output
    foreach my $id (sort keys %seq_ids) {
        my $seq_str;
        foreach my $aln_site (@sites) {
            foreach (@$aln_site) { $seq_str .= $_->{nt} if $_->{id} eq $id }
        }

        my $loc_seq = Bio::LocatableSeq->new(-seq => $seq_str, -id => $id, -start => 1);
        my $end = $loc_seq->end;
        $loc_seq->end($end);
        $new_aln->add_seq($loc_seq)
    }

    $aln = $new_aln
}

=head2 avg_id_by_win()

Calculate pairwise average sequence difference by windows (overlapping
windows with fixed step of 1). The window size is set in
C<$opts{"window"}> which is set via L<C<#initilize(\%opts)>|/initialize>.

=cut


sub avg_id_by_win {
    my $window_sz = $opts{"window"};
    for my $i (1 .. ($aln->length() - $window_sz + 1)) {
        my $slice = $aln->slice($i, $i + $window_sz - 1);
        my $pi = (100 - $slice->average_percentage_identity()) / 100;
        printf "%d\t%d\t%.4f\n", $i, $i + $window_sz - 1, $pi
    }
    exit
}

=head2 concat()

Concatenate multiple alignments sharing the same set of unique
IDs. This is normally used for concatenating individual gene
alignments of the same set of samples to a single one for making a
"supertree". Wraps
L<Bio::Align::UtilitiesE<gt>cat()|https://metacpan.org/pod/Bio::Align::Utilities#cat>.

=cut

##################################################################
# 2/26/2021: add position map (for rate4site applications)
##############################################################
sub concat {
    $aln = cat(@alns);
    warn "Alignment concated. Getting position maps...\n";
    my $refSeq = $opts{"ref-seq"} ? $aln->get_seq_by_id($opts{"ref-seq"}) : $aln->get_seq_by_pos(1);
#    my $refSeq = $aln->get_seq_by_pos(1);
    my @refGenesInOrder = map { $_ -> get_seq_by_id($refSeq->id) } @alns;
    # remap start & end for individual gene alignments, sync pos with concatenated aln:
    my $pos = 0;
    my %geneRange;
    for (my $i = 0; $i <= $#refGenesInOrder; $i++) {
	my $start = $pos + 1;
	$pos += $refGenesInOrder[$i]->length();
	my $end = $pos;
	$geneRange{$i+1} = {'start' => $start, 
			    'end' => $end, 
	};
    }

#    warn Dumper(\%geneRange);

    my @locTable;
    for (my $i = 1; $i <= $refSeq->length(); $i++) {
	next if $refSeq->subseq($i, $i) eq '-';
	my ($inGene, $posGene) = &__gene_order($i, \%geneRange);
	push @locTable, {
	    'pos_concat' => $i,
	    'pos_unaligned' => $refSeq->location_from_column($i)->start(),
	    'gene_order' => $inGene,
	    'pos_gene_aln' => $posGene,
	    'pos_gene_unaligned' => $refGenesInOrder[$inGene-1]->location_from_column($posGene)->start() 
	};
    }
    open LOG, ">concat.log";
    print LOG join "\t", ("seq_id", "pos_concat", "pos_residue", "gene", "pos_gene_aligned", "pos_gene");
    print LOG "\n"; 
    foreach (@locTable) {
	print LOG join "\t", (
	    $refSeq->id(),
	    $_->{pos_concat}, 
	    $_->{pos_unaligned}, 
	    $_->{gene_order},
	    $_->{pos_gene_aln},
	    $_->{pos_gene_unaligned}
	);
	print LOG "\n";
    }
    close LOG;
    warn "Position map of reference seq is saved in file concat.log\n";
}

sub __gene_order {
    my $pos = shift;
    my $ref = shift;
    my %range = %$ref;
    foreach my $gene (keys %range) {
	next unless $pos >= $range{$gene}->{start} && $pos <= $range{$gene}->{end};
	return ($gene, $pos - $range{$gene}->{'start'} + 1) 
    }
    die "Not found in any genes: position $pos\n";
}

#######################################

sub conserved_blocks {
    my $len=$aln->length();
    my $nseq = $aln->num_sequences();
    my $min_block_size = $opts{"con-blocks"};
    my %seq_ids;

    die "Alignment contains only one sequence: $file\n" if $nseq < 2;

    my (@blocks, $block);
    my $in_block=0;
    for (my $i=1; $i<=$len; $i++) {
        my ($ref_bases, $ref_ids) = &_get_a_site($i);
        %seq_ids = %{$ref_ids};
        my $is_constant = &_is_constant(&_paste_nt($ref_bases));
        if ($in_block) { # previous site is a contant one
            if ($is_constant) {
                $block->{length} ++;
                my @sites = @{$block->{sites}};
                push @sites, $ref_bases;
                $block->{sites} = \@sites;
                if ($i == $len) {
                    warn "Leaving a constant block at the end of alignment: $i\n";
                    push @blocks, $block if $block->{length} >= $min_block_size
                }
            } else {
                $in_block = 0;
                push @blocks, $block if $block->{length} >= $min_block_size;
                warn "Leaving a constant block at $i\n"
            }
        } else { # previous site not a constant one
            if ($is_constant) { # entering a block
                warn "Entering a constant block at site $i ...\n";
                $in_block=1;
                $block = {start => $i, length => 1, num_seq => $nseq, sites => [($ref_bases)]}  # start a new block
            }
        }
    }

    foreach my $bl (@blocks) {
        my $out = Bio::AlignIO->new(-file=> ">$file" . ".slice-". $bl->{start} . ".aln" , -format=>'clustalw');
        my $block_aln = Bio::SimpleAlign->new();
        foreach my $id (sort keys %seq_ids) {
            my ($seq_str, $ungapped_start, $ungapped_end);
            my @sites = @{$bl->{sites}};
            for (my $i = 0; $i <= $#sites; $i++) {
                my $ref_chars = $sites[$i];
                foreach (@$ref_chars) {
                    next unless $_->{id} eq $id;
                    $ungapped_start = $_->{ungapped_pos} if $i == 0;
                    $ungapped_end = $_->{ungapped_pos} if $i == $#sites;
                    $seq_str .= $_->{nt}
                }
            }

            my $loc_seq = Bio::LocatableSeq->new(-seq => $seq_str, -id => $id, -start => $ungapped_start, -end => $ungapped_end);
            $block_aln->add_seq($loc_seq)
        }
        $out->write_aln($block_aln)
    }
    exit
}

sub get_consensus {
    my $percent_threshold = $opts{"consensus"};
    my $consense = Bio::LocatableSeq->new(
        -seq   => $aln->consensus_string($percent_threshold),
        -id    => "Consensus_$percent_threshold",
        -start => 1,
        -end   => $aln->length()
   );
    $aln->add_seq($consense)
}

=head2 dna_to_protein()

Align CDS sequences according to their corresponding protein
alignment. Wraps
L<Bio::Align::Utilities-E<gt>aa_to_dna_aln()|https://metacpan.org/pod/Bio::Align::Utilities#aa_to_dna_aln>.

=cut

sub dna_to_protein {
    $aln = dna_to_aa_aln($aln)
}

sub remove_gapped_cols_in_one_seq {
    my $id = $opts{"rm-col"};
    my $nmatch=0;
    my $ref_seq;
    foreach ($aln->each_seq) {
        if ($_->id() =~ /$id/) { $nmatch++; $ref_seq = $_ }
    }
    die "Quit. No ref seq found or more than one ref seq!\n" if !$nmatch || $nmatch > 1;
    my ($ct_gap, $ref) = &_get_gaps($ref_seq);
    warn "Original length: " . $aln->length() . "\n";
    if ($ct_gap) {
        my @args;
        push @args, [$_, $_] foreach @$ref;
        $aln = $aln->remove_columns(@args);
        warn "New length: " . $aln->length() . "\n"
    } else {
        warn "No gap: " . $aln->length() . "\n"
    }
}

sub colnum_from_residue_pos {
    my ($id, $pos) = split /\s*,\s*/, $opts{"aln-index"};
    print $aln->column_from_residue_number($id, $pos), "\n";
    exit
}

=head2 list_ids()

List all sequence ids.

=cut

sub list_ids {
    my @ids;
    foreach ($aln->each_seq) { push @ids, $_->display_id() }
    say join "\n", @ids;
    exit
}

sub permute_states {
    my $new_aln = Bio::SimpleAlign->new();
    my $len=$aln->length();
    my $nseq = $aln->num_sequences();
    my @seq_ids;

    die "Alignment contains only one sequence: $file\n" if $nseq < 2;

    my @sites;
    my $ref_bases = &_get_a_site_v2();
    foreach (sort keys %$ref_bases) { push @seq_ids, $_ }
    for (my $i=1; $i<=$len; $i++) {
        my @bases;
        foreach (keys %$ref_bases) { push @bases, $ref_bases->{$_}->{$i} }
        @bases = shuffle(@bases);
        for (my $j=0; $j<$nseq; $j++) { $ref_bases->{$seq_ids[$j]}->{$i} = $bases[$j] }
    }

    foreach my $id (@seq_ids) {
        my $seq_str;
        for (my $i=1; $i<=$len; $i++) { $seq_str .= $ref_bases->{$id}->{$i} }

        my $loc_seq = Bio::LocatableSeq->new(-seq => $seq_str, -id => $id, -start => 1);
        my $end = $loc_seq->end;
        $loc_seq->end($end);
        $new_aln->add_seq($loc_seq)
    }
    $aln = $new_aln
}

=head2 protein_to_dna()

Align CDS sequences according to their corresponding protein
alignment. Wraps
L<Bio::Align::Utilities-E<gt>aa_to_dna_aln()https://metacpan.org/pod/Bio::Align::Utilities#aa_to_dna_aln>.

=cut

sub protein_to_dna {
    use Bio::SeqIO;
    my $cds_in = Bio::SeqIO->new(-file=>$opts{pep2dna}, -format=>'fasta');
    my %CDSs;
    while (my $seq = $cds_in->next_seq()) { $CDSs{$seq->display_id()} = $seq }
    $aln = aa_to_dna_aln($aln, \%CDSs);
}

=head2 sample_seqs()

Picks I<n> random sequences from input alignment and produces a new
alignment consisting of those sequences.

If n is not given, default is the number of sequences in alignment
divided by 2, rounded down.

This functionality uses an implementation of Reservoir Sampling, based
on the algorithm found here:
http://blogs.msdn.com/b/spt/archive/2008/02/05/reservoir-sampling.aspx

=cut

sub sample_seqs {
    # If option was given with no number, take the integer part of num_sequences/2
    # Its OK to use int() here (especially since we want to round towards 0)
    my $num_seqs = $aln->num_sequences;
    my $sample_size = ($opts{"resample"} == 0) ? int($num_seqs / 2) : $opts{"resample"};

    die "Error: sample size ($sample_size) exceeds number of sequences in alignment: ($num_seqs)" if $sample_size > $num_seqs;

    # Use Reservoir Sampling to pick random sequences.
    my @sampled = (1 .. $sample_size);
    for my $j ($sample_size + 1 .. $num_seqs) {
        $sampled[ rand(@sampled) ] = $j if rand() <= ($sample_size / $j)
    }

    warn "Sampled the following sequences: @sampled\n\n";
    my $tmp_aln = $aln->select_noncont(@sampled);
    $aln = $tmp_aln
}

=head2 shuffle_sites()

Make a shuffled (not bootstrapped) alignment. This operation randomizes
alignment columns. It is used for testing the significance of long-runs
of conserved sites in an alignment (e.g., conserved intergenic spacers
[IGSs]).

=cut

sub shuffle_sites {
    my $new_aln = Bio::SimpleAlign->new();
    my $len = $aln->length();
    my $nseq = $aln->num_sequences();
    my %seq_ids;

    die "Alignment contains only one sequence: $file\n" if $nseq < 2;

    my @sites;
    for (my $i=1; $i<=$len; $i++) {
        my ($ref_bases, $ref_ids) = &_get_a_site($i);
        %seq_ids = %{$ref_ids};
        push @sites, $ref_bases
    }

    @sites = shuffle(@sites);

    my @order;
    push @order, $_->[0]->{pos} foreach @sites;
    print STDERR "Shuffled site order:\t", join(",", @order);
    print STDERR "\n";

    foreach my $id (sort keys %seq_ids) {
        my $seq_str;
        foreach my $aln_site (@sites) {
            foreach (@$aln_site) { $seq_str .= $_->{nt} if $_->{id} eq $id }
        }

        my $loc_seq = Bio::LocatableSeq->new(-seq => $seq_str, -id => $id, -start => 1);
        my $end = $loc_seq->end;
        $loc_seq->end($end);
        $new_aln->add_seq($loc_seq);
    }
    $aln = $new_aln;
}

sub random_slice {
    my $slice_length = $opts{'random-slice'};
    my $len=$aln->length();
    my $start = int(rand($len - $slice_length+1));
    my $end = $start + $slice_length - 1;
    $aln = $aln->slice($start, $end);
}

sub select_third_sites {
    my $new_aln = Bio::SimpleAlign->new();
    my $len=$aln->length();
    my $nseq = $aln->num_sequences();
    my @seq_ids;

    die "Alignment contains only one sequence: $file\n" if $nseq < 2;

    my $ref_bases = &_get_a_site_v2();
    foreach (sort keys %$ref_bases) { push @seq_ids, $_ }

    my @sites;
    for (my $i=3; $i<=$len; $i+=3) { push @sites, $i }

    foreach my $id (sort @seq_ids) {
        my $seq_str;
        $seq_str .= $ref_bases->{$id}->{$_} foreach @sites;

        my $loc_seq = Bio::LocatableSeq->new(-seq => $seq_str, -id => $id, -start => 1);
        my $end = $loc_seq->end;
        $loc_seq->end($end);

        $new_aln->add_seq($loc_seq)
    }
    $aln = $new_aln
}

sub remove_third_sites {
    my $new_aln = Bio::SimpleAlign->new();
    my $len=$aln->length();
    my $nseq = $aln->num_sequences();
    my @seq_ids;

    die "Alignment contains only one sequence: $file\n" if $nseq < 2;

    my $ref_bases = &_get_a_site_v2();
    foreach (sort keys %$ref_bases) { push @seq_ids, $_ }

    my @sites;
    for (my $i=1; $i<=$len; $i++) { push @sites, $i if $i % 3 }

    foreach my $id (sort @seq_ids) {
        my $seq_str;
        $seq_str .= $ref_bases->{$id}->{$_} foreach @sites;

        my $loc_seq = Bio::LocatableSeq->new(-seq => $seq_str, -id => $id, -start => 1);
        my $end = $loc_seq->end;
        $loc_seq->end($end);

        $new_aln->add_seq($loc_seq)
    }
    $aln = $new_aln
}

sub upper_case {
    $aln->uppercase()
}

########################## internal subroutine #######################

# For use in draw_codon_view
# Pad display ids with a minimum of 4 spaces using the longest display id as a reference point for length. Pass-by-reference, so don't return array.
# Return length of longest id plus padding.
sub _pad_display_id {
    my $display_id = shift;
    my $max_len    = shift;
    my $padding = ($max_len - length($display_id));
    $display_id .= " " x $padding;
    return $display_id
}

# Used by draw_codon_view. Calculates position of final position in alinged block, prints the current position there.
sub _print_positions {
    my $nuc_count    = shift;
    my $aln_length   = shift;
    my $block_length = 3 * $opts{"codon-view"};
    my $max_id_len   = shift;
    my $num_spaces   = 0;

    my $start_pos = $nuc_count + 1;
    my $last_pos  = 0;
    my $offset    = 0;
    if (($nuc_count + $block_length) >= $aln_length) {
        $last_pos = $aln_length;
        my $diff = $aln_length - $nuc_count;
        $offset = $diff + ($diff) / 3 + ($diff % 3) - 2 # $diff % 3 gives the number of extra non-codon nucleotides
    } else {
        $last_pos = $nuc_count + $block_length;
        $offset = $block_length + ($block_length) / 3 - 2
    }

    # -1 since we are also printing the starting position.
    $num_spaces += $offset - 1;

 # $last_pos_len = length of last_pos treated as a string (ie length(335) = 3)
    my $last_pos_len = length($last_pos);

    # Pad $start_pos with $num_blanks blanks if it is shorter than $last_pos
    my $num_blanks = $last_pos_len - length($start_pos);
    $start_pos = " " x $num_blanks . "$start_pos" if length($start_pos) < $last_pos_len;

    for (my $i = 0; $i < $last_pos_len; $i++) {
        print " " x $max_id_len . substr($start_pos, $i, 1) . " " x ($num_spaces) . substr($last_pos, $i, 1) . "\n"
    }
}

# Function: _del_or_pick
# Desc: Internal function. Generic code for either picking or deleting a sequence from an alignment. Used by del_seqs and pick_seqs.
# Input:
#   $id_list, a user-supplied string consisting of comma-separated seq id values
#   $method, the name of the Bio::SimpleAlign method to use (remove_seq or add_seq)
#   $need_new, a flag indicating whether a new Bio::SimpleAlign object is needed
# Returns: Nothing; uses the $aln global variable

sub _del_or_pick {
    my ($id_list, $method, $need_new) = @_;
    my $new_aln = ($need_new) ? Bio::SimpleAlign->new() : $aln;

    my @selected = split(/\s*,\s*/, $id_list);
    foreach my $seq ($aln->each_seq) {
        my $seqid = $seq->display_id();
        foreach (@selected) {
            next unless $seqid eq $_;
            $new_aln->$method($seq)
        }
    }
    $aln = $new_aln if $need_new == 1
}

sub _get_gaps {
    my $seq = shift;
    my $seq_str = $seq->seq();
    my @chars = split //, $seq_str;
    my $cts = 0;
    my @pos=();
    for (my $i=0; $i<=$#chars; $i++) {
	   if ($chars[$i] eq '-') { push @pos, $i; $cts++ }
    }
    warn "Found " . scalar(@pos) ." gaps at (@pos) on " . $seq->id() . "\n";
    return ($cts, \@pos)
}

sub _paste_nt {
    my $ref = shift;
    my @nts;
    push @nts, $_->{nt} foreach @$ref;
    return \@nts
}

sub _get_a_site {
    my $pos = shift;
    my (@chars, %seq_ids);

    foreach my $seq ($aln->each_seq) {
        my $ungapped = 0;
        $seq_ids{$seq->id()}++;
        my $state;
        for (my $i = 1; $i <= $pos; $i++) {
            $state = $seq->subseq($i, $i);
            $ungapped++ unless $state eq '-'
        }

        push @chars, {
	       nt => $seq->subseq($pos, $pos),
	       ungapped_pos => ($state eq '-') ? "gap" : $ungapped++,
	       id => $seq->id(),
	       pos => $pos,
	   }
    }
    return (\@chars, \%seq_ids)
}

sub _seen_bases {
    my %count;
    my $ref   = shift;
    my @array = @$ref;
    my $constant = 1;

    $count{$_}++ foreach @array;
    return \%count;
}

sub _is_constant {
    my %count;
    my $ref   = shift;
    my @array = @$ref;
    my $constant = 1;

    $count{$_}++ foreach @array;
    my @keys = keys %count;
    $constant = 0 if @keys > 1;
    return $constant
}

sub _column_status {
    my %count;
    my $ref   = shift;
    my @array = @$ref;
    my $st    = { gap => 0, informative => 1, constant => 1 };

    foreach (@array) {
        $count{$_}++;
        $st->{gap} = 1 if $_ =~ /[\-\?]/
    }

    my @keys = keys %count;
    foreach (values %count) {
        if ($_ < 2) { $st->{informative} = 0; last }
    }
    $st->{constant} = 0 if @keys > 1;
    return $st
}

sub _get_a_site_v2 {
    my %seq_ids;
    my $len = $aln->length();
    foreach my $seq ($aln->each_seq) {
        my $id = $seq->id();
        for (my $i = 1; $i <= $len; $i++) { $seq_ids{$id}{$i} = $seq->subseq($i, $i) }
    }
    return \%seq_ids
}

sub _find_max_id_len {
    my $seqs = shift;
    my @sorted_by_length = sort {length $a->display_id <=> length $b->display_id} @$seqs;
    return length $sorted_by_length[-1]->display_id
}

1;
__END__

=head1 EXTENDING THIS MODULE

We encourage BioPerl developers to add command-line interface to their BioPerl methods here.

Here is how to extend.  We'll use option C<--avpid> as an example.

=over 4

=item *
Create a new method like one of the above in the previous section.

=item *
Document your method in pod using C<=head2>. For example:

    =head2 print_avpid

    Print the average percent identity of an alignment.

    Wraps
    L<Bio::SimpleAlign-E<gt>average_percentage_identity()|https://metacpan.org/pod/Bio::SimpleAlign#average_percentage_identity>.

    =cut

See L<C<print_avpid()>|/print_avpid> for how this gets rendered.


=item *
Add the method to C<@EXPORT> list in C<AlnManipulations.pm>.

=item *
Add option to C<%opt_displatch> which maps the option used in C<bioaln> to the subroutine that
gets called here. For example:

    "avpid" => \&print_avp_id,

=item *
Add option in to C<bioaln> script. See the code that starts:

    GetOptions(
    ...
    "avpid|a",
    ...

This option has a short option name C<a> and takes no additional argument

=item *
Write a test for the option. See the file C<t/10test-bioaln.t> and L<Testing|https://github.com/bioperl/p5-bpwrapper/wiki/Testing>.

=item *
Share back. Create a pull request to the github repository and contact
Weigang Qiu, City University of New York, Hunter College (L<mailto:weigang@genectr.hunter.cuny.edu>)

=back

=head1 SEE ALSO

=over 4

=item *

L<bioaln>: command-line tool for using this

=item *

L<Qiu Lab wiki page|http://diverge.hunter.cuny.edu/labwiki/Bioutils>

=item *

L<Github project wiki page|https://github.com/bioperl/p5-bpwrapper>

=back

=head1 CONTRIBUTORS

=over 4

=item *
William McCaig <wmccaig at gmail dot com>

=item *
Girish Ramrattan <gramratt at gmail dot com>

=item  *
Che Martin <che dot l dot martin at gmail dot com>

=item  *
Yözen Hernández yzhernand at gmail dot com

=item *
Levy Vargas <levy dot vargas at gmail dot com>

=item  *
L<Weigang Qiu|mailto:weigang@genectr.hunter.cuny.edu> (Maintainer)

=item *
Rocky Bernstein

=back

=cut
