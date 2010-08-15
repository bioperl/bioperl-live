# BioPerl module for SimpleAlign
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code
#
#  History:
#	11/3/00 Added threshold feature to consensus and consensus_aa  - PS
#	May 2001 major rewrite - Heikki Lehvaslaiho
#	August 2010 refactoring - Jun Yin

=head1 NAME

Bio::SimpleAlign - Multiple alignments held as a set of sequences

=head1 SYNOPSIS

  # Use Bio::AlignIO to read in the alignment
  $str = Bio::AlignIO->new(-file => 't/data/testaln.pfam');
  $aln = $str->next_aln();

  # Describe
  print $aln->length;
  print $aln->num_residues;
  print $aln->is_flush;
  print $aln->num_sequences;
  print $aln->score;
  print $aln->percentage_identity;
  print $aln->consensus_string(50);

  # Find the position in the alignment for a sequence location
  $pos = $aln->column_from_residue_number('1433_LYCES', 14); # = 6;

  # Extract sequences and check values for the alignment column $pos
  foreach $seq ($aln->next_Seq) {
      $res = $seq->subseq($pos, $pos);
      $count{$res}++;
  }
  foreach $res (keys %count) {
      printf "Res: %s  Count: %2d\n", $res, $count{$res};
  }

  # Manipulate
  $aln->remove_LocatableSeq($seq);
  $mini_aln = $aln->select_columns([20..30,45..48])  # get a block of columns
  #more parameters
  $mini_aln = $aln->select_columns(-selection=>[20..30],-toggle=>0,-keepgaponly=>1)
  $mini_aln = $aln->select_Seqs([1,5..10,15]) # select certain sequences
  $aln->remove_columns([20..30]); # remove by position
  $aln->remove_columns(['mismatch']); # remove by property

  # Analyze
  $str = $aln->consensus_string($threshold_percent);
  $str = $aln->match_line();
  $str = $aln->cigar_line();
  $id = $aln->percentage_identity;

  # See the module documentation for details and more methods.

=head1 DESCRIPTION

SimpleAlign is an object that handles a multiple sequence alignment
(MSA). It is very permissive of types (it does not insist on sequences
being all same length, for example). Think of it as a set of sequences
with a whole series of built-in manipulations and methods for reading and
writing alignments.

SimpleAlign uses L<Bio::LocatableSeq>, a subclass of L<Bio::PrimarySeq>,
to store its sequences. These are subsequences with a start and end
positions in the parent reference sequence. Each sequence in the
SimpleAlign object is a Bio::LocatableSeq.

SimpleAlign expects the combination of name, start, and end for a
given sequence to be unique in the alignment, and this is the key for the
internal hashes (name, start, end are abbreviated C<nse> in the code).
However, in some cases people do not want the name/start-end to be displayed:
either multiple names in an alignment or names specific to the alignment
(ROA1_HUMAN_1, ROA1_HUMAN_2 etc). These names are called
C<displayname>, and generally is what is used to print out the
alignment. They default to name/start-end.

The SimpleAlign Module is derived from the Align module by Ewan Birney.
This package is refactored by Jun Yin as part of the Google Summer of 
Code project in 2010.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
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
the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

Ewan Birney, birney@ebi.ac.uk

=head1 CONTRIBUTORS

Allen Day, allenday-at-ucla.edu,
Richard Adams, Richard.Adams-at-ed.ac.uk,
David J. Evans, David.Evans-at-vir.gla.ac.uk,
Heikki Lehvaslaiho, heikki-at-bioperl-dot-org,
Allen Smith, allens-at-cpan.org,
Jason Stajich, jason-at-bioperl.org,
Anthony Underwood, aunderwood-at-phls.org.uk,
Xintao Wei & Giri Narasimhan, giri-at-cs.fiu.edu
Brian Osborne, bosborne at alum.mit.edu
Weigang Qiu, Weigang at GENECTR-HUNTER-CUNY-EDU
Hongyu Zhang, forward at hongyu.org
Jay Hannah, jay at jays.net
Alexandr Bezginov, albezg at gmail.com
Jun Yin, jun dot yin at ucd dot ie

=head1 SEE ALSO

L<Bio::LocatableSeq>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# 'Let the code begin...

package Bio::SimpleAlign;
use vars qw(%CONSERVATION_GROUPS);
use strict;

use Bio::LocatableSeq;  # uses Seq's as list

use Bio::Seq;
use Bio::SeqFeature::Generic;

BEGIN {
    # This data should probably be in a more centralized module...
    # it is taken from Clustalw documentation.
    # These are all the positively scoring groups that occur in the
    # Gonnet Pam250 matrix. The strong and weak groups are
    # defined as strong score >0.5 and weak score =<0.5 respectively.

    %CONSERVATION_GROUPS = (
            'strong' => [ qw(
						 STA
						 NEQK
						 NHQK
						 NDEQ
						 QHRK
						 MILV
						 MILF
						 HY
						 FYW )],
				'weak' => [ qw(
                      CSA
					       ATV
					       SAG
					       STNK
					       STPA
					       SGND
					       SNDEQK
					       NDEQHK
					       NEQHRK
					       FVLIM
					       HFY )],);
}

use base qw(Bio::Root::Root Bio::Align::AlignI Bio::AnnotatableI 
	    Bio::FeatureHolderI);

=head2 new

 Title     : new
 Usage     : my $aln = Bio::SimpleAlign->new();
 Function  : Creates a new simple align object
 Returns   : Bio::SimpleAlign
 Args      : -source     => string representing the source program
                            where this alignment came from
             -annotation => Bio::AnnotationCollectionI
             -seq_annotation => Bio::AnnotationCollectionI for sequences
                                 (requires -annotation also be set)
             -seqs       => array ref containing Bio::LocatableSeq 
                                                 or Bio::Seq::Meta
             -consensus  => consensus string
             -consensus_meta  => Bio::Seq::Meta object containing 
                                 consensus met information (kludge)

=cut


sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($src, $score, $id, $acc, $desc, $seqs, $feats, $coll, $sa, $con, $cmeta) = $self->_rearrange([qw(
                                            SOURCE
                                            SCORE
                                            ID
                                            ACCESSION
                                            DESCRIPTION
                                            SEQS
                                            FEATURES
                                            ANNOTATION
                                            SEQ_ANNOTATION
                                            CONSENSUS
                                            CONSENSUS_META
                                            )], @args);
  $src && $self->source($src);
  defined $score && $self->score($score);
  # we need to set up internal hashs first!

  $self->{'_seq'} = {};
  $self->{'_order'} = {};
  $self->{'_start_end_lists'} = {};
  $self->{'_dis_name'} = {};
  $self->{'_id'} = 'NoName';
  # maybe we should automatically read in from args. Hmmm...
  $id  && $self->id($id);
  $acc && $self->accession($acc);
  $desc && $self->description($desc);
  $coll && $self->annotation($coll);
  # sequence annotation is layered into a provided annotation collection (or dies)
  if ($sa) {
    $self->throw("Must supply an alignment-based annotation collection (-annotation) ".
                 "with a sequence annotation collection")
        if !$coll;
    $coll->add_Annotation('seq_annotation', $sa);
  }
  if ($feats && ref $feats eq 'ARRAY') {
    for my $feat (@$feats) {
        $self->add_SeqFeature($feat);
    }
  }
  $con && $self->consensus($con);
  $cmeta && $self->consensus_meta($cmeta);
  # assumes these are in correct alignment order
  if ($seqs && ref($seqs) eq 'ARRAY') {
    for my $seq (@$seqs) {
        $self->add_Seq($seq);
    }
  }

  return $self; # success - we hope!
}

=head1 Alignment modifier methods

These methods modify the original MSA by adding, removing or shuffling 
complete sequences.

=head2 add_Seq

 Title     : add_Seq
 Usage     : $aln->add_Seq($newseq);
             $aln->add_Seq(-SEQ=>$newseq, -ORDER=>5);
 Function  : Adds another sequence to the alignment. *Does not* align
             it - just adds it to the hashes.
             If -ORDER is specified, the sequence is inserted at the
             the position spec'd by -ORDER, and existing sequences
             are pushed down the storage array.
 Returns   : 1
 Args      : A Bio::LocatableSeq object
             Positive integer for the sequence position (optional)

See L<Bio::LocatableSeq> for more information

=cut

sub addSeq {
    my $self = shift;
    $self->deprecated("addSeq - deprecated method. Use add_Seq() instead.");
    $self->add_Seq(@_);
}

sub add_seq {
    my $self = shift;
    $self->deprecated("add_seq - deprecated method. Use add_Seq() instead.");
    $self->add_Seq(@_);
}

sub add_Seq {
    my $self = shift;
    my @args = @_;
    my ($seq, $order) = $self->_rearrange([qw(SEQ ORDER)], @args);
    my ($name,$id,$start,$end);

    unless ($seq) {
	$self->throw("LocatableSeq argument required");
    }
    if( ! ref $seq || ! $seq->isa('Bio::LocatableSeq') ) {
	$self->throw("Unable to process non locatable sequences [". ref($seq). "]");
    }
    !defined($order) and $order = 1 + keys %{$self->{'_seq'}}; # default 
    $order--; # jay's patch (user-specified order is 1-origin)
    
    if ($order < 0) {
	$self->throw("User-specified value for ORDER must be >= 1");
    }

    $id = $seq->id() ||$seq->display_id || $seq->primary_id;

    # build the symbol list for this sequence,
    # will prune out the gap and missing/match chars
    # when actually asked for the symbol list in the
    # symbol_chars
    # map { $self->{'_symbols'}->{$_} = 1; } split(//,$seq->seq) if $seq->seq;

    $name = $seq->get_nse;

    if( $self->{'_seq'}->{$name} ) {
	$self->warn("Replacing one sequence [$name]\n") unless $self->verbose < 0;
    }
    else {
	$self->debug( "Assigning $name to $order\n");

    my $ordh = $self->{'_order'};
    if ($ordh->{$order}) {
        # make space to insert
        # $c->() returns (in reverse order) the first subsequence 
        # of consecutive integers; i.e., $c->(1,2,3,5,6,7) returns
        # (3,2,1), and $c->(2,4,5) returns (2).
        my $c;
        $c = sub { return (($_[1]-$_[0] == 1) ? ($c->(@_[1..$#_]),$_[0]) : $_[0]); };
        map { 
     $ordh->{$_+1} = $ordh->{$_}
        } $c->(sort {$a <=> $b} grep {$_ >= $order} keys %{$ordh});

    }
    $ordh->{$order} = $name;

	unless( exists( $self->{'_start_end_lists'}->{$id})) {
	    $self->{'_start_end_lists'}->{$id} = [];
	}
	push @{$self->{'_start_end_lists'}->{$id}}, $seq;
    }

    $self->{'_seq'}->{$name} = $seq;
    
    1;

}


=head2 remove_LocatableSeq

 Title     : remove_LocatableSeq
 Usage     : $aln->remove_LocatableSeq($seq);
 Function  : Removes a single sequence from an alignment
 Returns   : 1
 Argument  : a Bio::LocatableSeq object

=cut

sub removeSeq {
    my $self = shift;
    $self->deprecated("removeSeq - deprecated method. Use remove_seq() instead.");
    $self->remove_LocatableSeq(@_);
}

sub remove_seq {
    my $self = shift;
    $self->deprecated("remove_seq - deprecated method. Use remove_LocatableSeq() instead.");
    $self->remove_LocatableSeq(@_);	
}


sub remove_LocatableSeq {
    my $self = shift;
    my $seq = shift;
    my ($name,$id,$start,$end);

    $self->throw("Need Bio::Locatable seq argument ")
	unless ref $seq && $seq->isa( 'Bio::LocatableSeq');

    $id = $seq->id();
    $start = $seq->start();
    $end  = $seq->end();
    $name = sprintf("%s/%d-%d",$id,$start,$end);
    

    if( !exists $self->{'_seq'}->{$name} ) {
	$self->throw("Sequence $name does not exist in the alignment to remove!");
    }

    delete $self->{'_seq'}->{$name};

    # we need to remove this seq from the start_end_lists hash

    if (exists $self->{'_start_end_lists'}->{$id}) {
	# we need to find the sequence in the array.

	my ($i, $found);;
	for ($i=0; $i < @{$self->{'_start_end_lists'}->{$id}}; $i++) {
	    if (${$self->{'_start_end_lists'}->{$id}}[$i] eq $seq) {
		$found = 1;
		last;
	    }
	}
	if ($found) {
	    splice @{$self->{'_start_end_lists'}->{$id}}, $i, 1;
	}
	else {
	    $self->throw("Could not find the sequence to remoce from the start-end list");
	}
    }
    else {
	$self->throw("There is no seq list for the name $id");
    }
    # we need to shift order hash
    my %rev_order = reverse %{$self->{'_order'}};
    my $no = $rev_order{$name};
    my $num_sequences = $self->num_sequences;
    for (; $no < $num_sequences; $no++) {
       $self->{'_order'}->{$no} = $self->{'_order'}->{$no+1};
    }
    delete $self->{'_order'}->{$no};
    return 1;
}

=head2 remove_Seqs

 Title     : remove_Seqs
 Usage     : $aln->remove_Seqs([1,3,5..7]);
 Function  : Removes specified sequences from the alignment
 Returns   : 1 
 Argument  : An reference list of positive integers for the selected 
             sequences. An optional parameter can be defined to toggle 
             the coordinate selection.

=cut

sub remove_Seqs {
	my $self=shift;
	my ($sel, $toggle) = $self->_rearrange([qw(SELECTION TOGGLE)], @_);
	@{$sel}=sort {$a<=>$b} @{$sel};

   $self->throw("Select start has to be a positive integer, not [".$sel->[0]."]")
	unless $sel->[0] =~ /^\d+$/ and $sel->[0] > 0;
    $self->throw("Select end has to be a positive integer, not [".$sel->[$#$sel]."]")
	unless $sel->[$#$sel] =~ /^\d+$/ and $sel->[$#$sel] > 0;
	
	#use select_Seqs to select removed Seq objects, then use remove_LocatableSeq to remove the selected Seqs
	my $removed_aln;
	if($toggle) {
		$removed_aln=select_Seqs(_toggle_selection($sel));
	}
	else {
		$removed_aln=select_Seqs($sel);
	}
	
	foreach my $seq ($removed_aln->next_Seq) {
		$self->remove_LocatableSeq($seq);
	}

	return 1;
}

=head2 remove_redundant_Seqs

 Title   : remove_redundant_Seqs
 Usage   : $aln->remove_redundant_Seqs(0.7);
 Function: Removes sequences above given sequence similarity
           This function will grind on large alignments. Beware!
 Example :
 Returns : An array of the removed sequences
 Args    : float, threshold for similarity

=cut

sub purge {
    my ($self,@args) = @_;
    $self->deprecated("purge - deprecated method. Use remove_redundant_Seqs() instead.");
    $self->remove_redundant_Seq(@args);
}

sub remove_redundant_Seqs {
	my ($self,$perc) = @_;
	my (%duplicate, @dups);

	my @seqs = $self->next_Seq();

	for (my $i=0;$i< @seqs - 1;$i++ ) { #for each seq in alignment
		my $seq = $seqs[$i];

		#skip if already in duplicate hash
		next if exists $duplicate{$seq->display_id} ;
		my $one = $seq->seq();

		my @one = split '', $one;	#split to get 1aa per array element

		for (my $j=$i+1;$j < @seqs;$j++) {
			my $seq2 = $seqs[$j];

			#skip if already in duplicate hash
			next if exists $duplicate{$seq2->display_id} ;

			my $two = $seq2->seq();
			my @two = split '', $two;

			my $count = 0;
			my $res = 0;
			for (my $k=0;$k<@one;$k++) {
				if ( $one[$k] ne '.' && $one[$k] ne '-' && defined($two[$k]) &&
					  $one[$k] eq $two[$k]) {
					$count++;
				}
				if ( $one[$k] ne '.' && $one[$k] ne '-' && defined($two[$k]) &&
					  $two[$k] ne '.' && $two[$k] ne '-' ) {
					$res++;
				}
			}

			my $ratio = 0;
			$ratio = $count/$res unless $res == 0;

			# if above threshold put in duplicate hash and push onto
			# duplicate array for returning to get_unique
			if ( $ratio > $perc ) {
				$self->warn("duplicate: ", $seq2->display_id) if $self->verbose > 0;
				$duplicate{$seq2->display_id} = 1;
				push @dups, $seq2;
			}
		}
	}
	foreach my $seq (@dups) {
		$self->remove_LocatableSeq($seq);
	}
	return @dups;
}

=head2 uniq_Seq

 Title     : uniq_Seq
 Usage     : $aln->uniq_Seq():  Remove identical sequences in
             in the alignment.  Ambiguous base ("N", "n") and
             leading and ending gaps ("-") are NOT counted as
             differences.
             $aln->uniq_Seq is a variation of $aln->remove_redundant_Seqs(1)
 Function  : Make a new alignment of unique sequence types (STs)
 Returns   : 1a. if called in a scalar context, 
                a new Bio::SimpleAlign object (all sequences renamed as "ST")
             1b. if called in an array context, 
                a new Bio::SimpleAlign object, and a hashref whose keys
                are sequence types, and whose values are arrayrefs to 
                lists of sequence ids within the corresponding sequence type
             2. if $aln->verbose > 0, ST of each sequence is sent to 
                STDERR (in a tabular format)
 Argument  : None

=cut

sub uniq_seq {
	my ($self,@args) = @_;
	$self->deprecated("uniq_seq - deprecated method. Use uniq_Seq() instead.");
	$self->uniq_Seq(@args);	
}

sub uniq_Seq {
    my ($self, $seqid) = @_;
    my $aln = $self->new;
    my (%member, %order, @seq, @uniq_str, $st);
    my $order=0;
    my $len = $self->length();
    $st = {};
    foreach my $seq ( $self->next_Seq() ) {
	my $str = $seq->seq();

# it's necessary to ignore "n", "N", leading gaps and ending gaps in
# comparing two sequence strings

    # 1st, convert "n", "N" to "?" (for DNA sequence only):
	$str =~ s/n/\?/gi if $str =~ /^[atcgn-]+$/i;
    # 2nd, convert leading and ending gaps to "?":
	$str = &_convert_leading_ending_gaps($str, '-', '?');
    # Note that '?' also can mean unknown residue.
    # I don't like making global class member changes like this, too
    # prone to errors... -- cjfields 08-11-18
    local $Bio::LocatableSeq::GAP_SYMBOLS = '-\?';
	my $new = Bio::LocatableSeq->new(
                     -id      => $seq->id(),
					 -alphabet=> $seq->alphabet,
					 -seq     => $str,
					 -start   => $seq->start,
					 -end     => $seq->end
					 );
	push @seq, $new;
    }

    foreach my $seq (@seq) {
	my $str = $seq->seq();
	my ($seen, $key) = &_check_uniq($str, \@uniq_str, $len);
	if ($seen) { # seen before
	    my @memb = @{$member{$key}};
	    push @memb, $seq;
	    $member{$key} = \@memb;
	} else {  # not seen
	    push @uniq_str, $key;
	    $order++;
	    $member{$key} = [ ($seq) ];
	    $order{$key} = $order;
	}
    }

    foreach my $str (sort {$order{$a} <=> $order{$b}} keys %order) { # sort by input order
# convert leading/ending "?" back into "-" ("?" throws errors by SimpleAlign):
	my $str2 = &_convert_leading_ending_gaps($str, '?', '-');
# convert middle "?" back into "N" ("?" throws errors by SimpleAlign):
	$str2 =~ s/\?/N/g if $str2 =~ /^[atcg\-\?]+$/i;
	my $gap='-';
	my $end= CORE::length($str2);
	$end -= CORE::length($1) while $str2 =~ m/($gap+)/g;
	my $new = Bio::LocatableSeq->new(-id   =>"ST".$order{$str},
					 -seq  =>$str2,
					 -start=>1,
					 -end  =>$end
					 );
	$aln->add_Seq($new);
	foreach (@{$member{$str}}) {
	    push @{$$st{$order{$str}}}, $_->id(); # per Tristan's patch/Bug #2805
        $self->debug($_->id(), "\t", "ST", $order{$str}, "\n");
        }
    }
    return wantarray ? ($aln, $st) : $aln;
}

sub _check_uniq {  # check if same seq exists in the alignment
    my ($str1, $ref, $length) = @_;
    my @char1=split //, $str1;
    my @array=@$ref;

    return (0, $str1) if @array==0; # not seen (1st sequence)

    foreach my $str2 (@array) {
	my $diff=0;
	my @char2=split //, $str2;
	for (my $i=0; $i<=$length-1; $i++) {
	    next if $char1[$i] eq '?';
	    next if $char2[$i] eq '?';
	    $diff++ if $char1[$i] ne $char2[$i];
	}
	return (1, $str2) if $diff == 0;  # seen before
    }

    return (0, $str1); # not seen
}

sub _convert_leading_ending_gaps {
    my $s=shift;
    my $sym1=shift;
    my $sym2=shift;
    my @array=split //, $s;
# convert leading char:
    for (my $i=0; $i<=$#array; $i++) {
	($array[$i] eq $sym1) ? ($array[$i] = $sym2):(last);
    }
# convert ending char:
    for (my $i = $#array; $i>= 0; $i--) {
	($array[$i] eq $sym1) ? ($array[$i] = $sym2):(last);
    }
    my $s_new=join '', @array;
    return $s_new;
}
=head2 remove_columns

 Title     : remove_columns
 Usage     : $aln->remove_columns(['mismatch','weak']) or
             $aln->remove_columns([3,6..8]) or
             $aln->remove_columns(-pos=>[3,6..8],-toggle=>T)
 Function  : Modify the aligment with columns removed corresponding to
             the specified type or by specifying the columns by number.
             If you dont want to modify the alignment, 
             you can use select_columns instead
 Returns   : 1
 Args      : Array ref of types ('match'|'weak'|'strong'|'mismatch') 
             or array ref where the referenced array
             contains a pair of integers that specify a range.
             use remove_gaps to remove columns containing gaps
             The first column is 1

=cut

sub remove_columns {
	my ($self,@args) = @_;
   @args || $self->throw("Must supply column ranges or column types");
   my ($aln,$cont_loc);
	if(@args) {
		if(ref($args[0]) && $args[0][0] =~ /^[a-z_]+$/i) {
			$self->_remove_columns_by_type($args[0]);
		}
		else {
			my ($sel, $toggle) = $self->_rearrange([qw(SELECTION TOGGLE)], @args);
			if($toggle) {
				$cont_loc=_cont_coords(_toggle_selection($sel,$self->length));
			}
			else {
				$cont_loc=_cont_coords($sel);
			}
			$self->_remove_columns_by_num($cont_loc);
		}
	}
	else {
   	$self->throw("You must pass array references to remove_columns(), not @args");
	}
	# fix for meta, sf, ann
	1;
}

sub _remove_col {
    my ($self,$remove) = @_;
    my @new;
    
    my $gap = $self->gap_char;
    
    # splice out the segments and create new seq
    my ($firststart,$firstend)=($remove->[0][0]-1,$remove->[0][1]-1);
    foreach my $seq($self->next_Seq){
		my $sequence = $seq->seq;
		my $orig = $sequence;
		#calculate the new start
		if ($firststart == 0) {
		  my $start_adjust = () = substr($orig, 0, $firstend + 1) =~ /$gap/g;
		  $seq->start($seq->start() + $firstend + 1 - $start_adjust);
		}
		else {
		  my $start_adjust = $orig =~ /^$gap+/;
		  if ($start_adjust) {
		      $start_adjust = $+[0] == $firststart;
		  }
		  $seq->start($seq->start + $start_adjust);
		}
        
        foreach my $pair(@{$remove}){
            my $start = $pair->[0]-1;
            my $end   = $pair->[1]-1;
            $sequence = $seq->seq unless $sequence;
            $orig = $sequence;
            my $head =  $start > 0 ? substr($sequence, 0, $start) : '';
            my $tail = ($end + 1) >= CORE::length($sequence) ? '' : substr($sequence, $end + 1);
            $sequence = $head.$tail;
				$seq->seq($sequence) if $sequence;
				
            #calculate the new end
            #if (($end + 1) >= CORE::length($orig)) {
            #	print STDERR "SA2:$end\t,",CORE::length($orig),"\n";
             #   my $end_adjust = () = substr($orig, $start) =~ /$gap/g;
             #   $seq->end($seq->end - (CORE::length($orig) - $start) + $end_adjust-1);
            #}
            #else {
            #	print STDERR "SA:",$seq->display_id(),"\t",$seq->_ungapped_len,"\t",$seq->start,"\n";
            	$seq->end($seq->_ungapped_len+$seq->start-1);
            #}
        }
        
        if ($seq->end < $seq->start) {
            # we removed all columns except for gaps: set to 0 to indicate no
            # sequence
            $seq->start(0);
            $seq->end(0);
        }
        
    }
    
    # fix for meta, sf, ann    
    return 1;
}

sub _remove_columns_by_type {
	my ($self,$type) = @_;
	my @remove;

	my $gap = $self->gap_char if (grep { $_ eq 'gaps'} @{$type});
	my $all_gaps_columns = $self->gap_char if (grep /all_gaps_columns/,@{$type});
	my %matchchars = ( 'match'           => '\*',
                       'weak'             => '\.',
                       'strong'           => ':',
                       'mismatch'         => ' ',
                     );
	# get the characters to delete against
	my $del_char;
	foreach my $type (@{$type}){
		$del_char.= $matchchars{$type};
	}

	my $length = 0;
	my $match_line = $self->match_line;
	# do the matching to get the segments to remove
	if($del_char){
		while($match_line =~ m/[$del_char]/g ){
			my $start = pos($match_line);
			$match_line=~/\G[$del_char]+/gc;
			my $end = pos($match_line);

			#have to offset the start and end for subsequent removes
			$start-=$length;
			$end  -=$length;
			$length += ($end-$start+1);
			push @remove, [$start,$end];
		}
	}

	# remove the segments
	$self->_remove_col(\@remove) if $#remove >= 0;

    # fix for meta, sf, ann    
	1;
}


sub _remove_columns_by_num {
	my ($self,$positions) = @_;
	
	my @remove;
	my $length = 0;
	for(my $num=0;$num<@{$positions};) {
		my ($start, $end) = ($positions->[$num],$positions->[$num+1]);
        
		#have to offset the start and end for subsequent removes
		$start-=$length;
		$end  -=$length;
		$length += ($end-$start+1);
		push @remove, [$start,$end];
		$num+=2;
    }

	#remove the segments
	$self->_remove_col(\@remove) if $#remove >= 0;
	# fix for meta, sf, ann    
	1;
}


sub _cont_coords {
	#This function is used to merge the coordinates from select and remove functions in order to reduce the number of calculations in select and #remove. For exmaple, if the input of remove_columns is remove_columns([2,5,7..10]), this function will transform ([2,5,7..10]) to 
 	# ([2,2,5,5,7,10]).
	
	my ($old_coords)=@_;
	@{$old_coords}=sort {$a<=>$b} @{$old_coords};
	
	my $cont_coords;
	
	push @{$cont_coords},$old_coords->[0];
	for(my $num=0;$num<@{$old_coords};) {
		if($num==@{$old_coords}-1) {
			#for single selection
			push @{$cont_coords},$old_coords->[$num];
			last;
		}
		if($old_coords->[$num+1]-$old_coords->[$num]>1) {
			if($num+2==@{$old_coords}) {
				push @{$cont_coords},$old_coords->[$num],$old_coords->[$num+1],$old_coords->[$num+1];
				last;
			}
			else {
				push @{$cont_coords},$old_coords->[$num],$old_coords->[$num+1];
			}
		}
		else {
			if ($num+2==@{$old_coords}) {
				push @{$cont_coords},$old_coords->[$num+1];
				last;
			}
		}
		$num++;
	}
	return $cont_coords;
}


sub _toggle_selection {
	#This function is used to toggle the selection of sequences or columns
	my ($old_coords,$length)=@_;
	my %hash=map {$_=>1} @{$old_coords};
	my $new_coords;
	for(my $num=1;$num<=$length;$num++) {
		unless(defined($hash{$num})) {
			push @{$new_coords},$num;
		}
	}
	return $new_coords;
}


=head2 sort_alphabetically

 Title     : sort_alphabetically
 Usage     : $ali->sort_alphabetically
 Function  : Changes the order of the alignment to alphabetical on name
             followed by numerical by number.
 Returns   : 1
 Argument  : None

=cut

sub sort_alphabetically {
    my $self = shift;
    my ($seq,$nse,@arr,%hash,$count);

    foreach $seq ( $self->next_Seq() ) {
	$nse = $seq->get_nse;
	$hash{$nse} = $seq;
    }

    $count = 0;

    %{$self->{'_order'}} = (); # reset the hash;

    foreach $nse ( sort _alpha_startend keys %hash) {
	$self->{'_order'}->{$count} = $nse;

	$count++;
    }
    1;
}

=head2 sort_by_list

 Title     : sort_by_list
 Usage     : $aln_ordered=$aln->sort_by_list($list_file)
 Function  : Arbitrarily order sequences in an alignment
 Returns   : A new Bio::SimpleAlign object
 Argument  : a file listing sequence names in intended order 
             (one name per line)

=cut

sub sort_by_list {
    my ($self, $list) = @_;
    my (@seq, @ids, %order);

    foreach my $seq ( $self->next_Seq() ) {
        push @seq, $seq;
        push @ids, $seq->display_id;
    }

    my $ct=1;
    open(my $listfh, '<', $list) || $self->throw("can't open file for reading: $list");
    while (<$listfh>) {
      chomp;
      my $name=$_;
      $self->throw("Not found in alignment: $name") unless &_in_aln($name, \@ids);
      $order{$name}=$ct++;
    }
    close($listfh);
    
    # use the map-sort-map idiom:
    my @sorted= map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [$order{$_->id()}, $_] } @seq;
    my $aln = $self->new;
    foreach (@sorted) { $aln->add_Seq($_) }
    return $aln;
}


=head2 sort_by_pairwise_identity

 Title     : sort_by_pairwise_identity()
 Usage     : $ali->sort_by_pairwise_identity(2)
 Function  : Changes the order of the alignment by the pairwise percentage 
             identity of the reference sequence
 Returns   : 1
 Argument  : Optional, the position or id of reference sequences to be compared 
             with. Default is the first sequence

=cut

sub sort_by_pairwise_identity {
	my ($self,$seqid) = @_;
	my ($seq,$nse,@arr,%hash,$count);
	$seqid=defined($seqid)?$seqid:1;
	my @pairwise_iden=$self->pairwise_percentage_identity($seqid);
	
	my $seqcount = 0;
	foreach $seq ( $self->next_Seq() ) {
		$nse = $seq->get_nse;
		$hash{$nse} = $pairwise_iden[$seqcount++];
	}

	$count=0;

	%{$self->{'_order'}} = (); # reset the hash;

	foreach $nse ( sort {$hash{$b}<=>$hash{$a}} keys %hash) {
		$self->{'_order'}->{$count} = $nse;
		$count++;
	}
	1;	
}

=head2 sort_by_length

 Title     : sort_by_length()
 Usage     : $ali->sort_by_length()
 Function  : Changes the order of the alignment by the ungapped length of 
              the sequences
 Returns   : 1
 Argument  : None

=cut

sub sort_by_length {
	my $self = shift;
	my %seq2length;
	
	foreach my $order (keys %{$self->{'_order'}}) {
		my $nse=$self->{'_order'}->{$order};
		my $seq=$self->{'_seq'}->{$nse};
		$seq2length{$seq->get_nse}=$seq->_ungapped_len;
	}
	
	%{$self->{'_order'}} = (); # reset the hash;
	
	my $count=0;
	foreach my $nse (sort {$seq2length{$a}<=>$seq2length{$b}} keys %seq2length) {
		$self->{'_order'}->{$count} = $nse;
		$count++;
	}
	
	1;
}

=head2 sort_by_start

 Title     : sort_by_start
 Usage     : $ali->sort_by_start
 Function  : Changes the order of the alignment to the start position of 
             each subalignment    
 Returns   : 1
 Argument  : None

=cut

sub sort_by_start {
    my $self = shift;
    my ($seq,$nse,@arr,%hash,$count);
    foreach $seq ( $self->next_Seq() ) {
        $nse = $seq->get_nse;
        $hash{$nse} = $seq;
    }
    $count = 0;
    %{$self->{'_order'}} = (); # reset the hash;
    foreach $nse ( sort _startend keys %hash) {
        $self->{'_order'}->{$count} = $nse;
        $count++;
    }
    1;
}

sub _startend
{
    my ($aname,$arange) = split (/[\/]/,$a);
    my ($bname,$brange) = split (/[\/]/,$b);
    my ($astart,$aend) = split(/\-/,$arange);
    my ($bstart,$bend) = split(/\-/,$brange);
    return $astart <=> $bstart;
}
=head2 set_new_reference

 Title     : set_new_reference
 Usage     : $aln->set_new_reference(3 or 'B31'):  Select the 3rd sequence, 
             or the sequence whoes name is "B31" (full, exact, and 
             case-sensitive), as the reference (1st) sequence
 Function  : Change/Set a new reference (i.e., the first) sequence
 Returns   : a new Bio::SimpleAlign object.
             Throws an exception if designated sequence not found
 Argument  : a positive integer of sequence order, or a sequence name
             in the original alignment

=cut

sub set_new_reference {
    my ($self, $seqid) = @_;
    my $aln = $self->new;
    my (@seq, @ids, @new_seq);
    my $is_num=0;
    foreach my $seq ( $self->next_Seq() ) {
	push @seq, $seq;
	push @ids, $seq->display_id;
    }

    if ($seqid =~ /^\d+$/) { # argument is seq position
	$is_num=1;
	$self->throw("The new reference sequence number has to be a positive integer >1 and <= num_sequences ") if ($seqid <= 1 || $seqid > $self->num_sequences);
    } else { # argument is a seq name
	$self->throw("The new reference sequence not in alignment ") unless &_in_aln($seqid, \@ids);
    }

    for (my $i=0; $i<=$#seq; $i++) {
	my $pos=$i+1;
        if ( ($is_num && $pos == $seqid) || ($seqid eq $seq[$i]->display_id) ) {
	    unshift @new_seq, $seq[$i];
	} else {
	    push @new_seq, $seq[$i];
	}
    }
    foreach (@new_seq) { $aln->add_Seq($_);  }
    return $aln;
}

sub _in_aln {  # check if input name exists in the alignment
    my ($str, $ref) = @_;
    foreach (@$ref) {
	return 1 if $str eq $_;
    }
    return 0;
}



=head1 Alignment selection methods

These methods are used to select the sequences or horizontal/vertical 
subsets of the current MSA.

=head2 next_Seq

 Title     : next_Seq
 Usage     : foreach $seq ( $aln->next_Seq() )
 Function  : Gets a Seq object from the alignment
 Returns   : Seq object
 Argument  :

=cut

sub eachSeq {
    my $self = shift;
    $self->deprecated("eachSeq - deprecated method. Use next_Seq() instead.");
    $self->next_Seq();
}

sub each_seq {
    my $self = shift;
    $self->deprecated("each_seq - deprecated method. Use next_Seq() instead.");
    $self->next_Seq();
}

sub next_Seq {
	my $self = shift;
	my (@arr,$order);

	foreach $order ( sort { $a <=> $b } keys %{$self->{'_order'}} ) {
		if( exists $self->{'_seq'}->{$self->{'_order'}->{$order}} ) {
			push(@arr,$self->{'_seq'}->{$self->{'_order'}->{$order}});
		}
	}
	return @arr;
}


=head2 next_alphabetically

 Title     : next_alphabetically
 Usage     : foreach $seq ( $aln->next_alphabetically() )
 Function  : Returns a sequence object, but the objects are returned
             in alphabetically sorted order.
             Does not change the order of the alignment.
 Returns   : Seq object
 Argument  : None

=cut

sub each_alphabetically {
    my $self = shift;
    $self->deprecated("each_alphabetically - deprecated method. Use next_alphabetically() instead.");
    $self->next_Seq();	
}


sub next_alphabetically {
	my $self = shift;
	my ($seq,$nse,@arr,%hash,$count);

	foreach $seq ( $self->next_Seq() ) {
		$nse = $seq->get_nse;
		$hash{$nse} = $seq;
	}

	foreach $nse ( sort _alpha_startend keys %hash) {
		push(@arr,$hash{$nse});
	}
	return @arr;
}

sub _alpha_startend {
    my ($aname,$astart,$bname,$bstart);
    ($aname,$astart) = split (/-/,$a);
    ($bname,$bstart) = split (/-/,$b);

    if( $aname eq $bname ) {
	return $astart <=> $bstart;
    }
    else {
	return $aname cmp $bname;
    }
}

=head2 next_Seq_with_id

 Title     : next_Seq_with_id
 Usage     : foreach $seq ( $aln->next_Seq_with_id() )
 Function  : Gets a Seq objects from the alignment, the contents
             being those sequences with the given name (there may be
             more than one)
 Returns   : Seq object
 Argument  : a seq name

=cut

sub eachSeqWithId {
    my $self = shift;
    $self->deprecated("eachSeqWithId - deprecated method. Use next_Seq_with_id() instead.");
    $self->next_Seq_with_id(@_);
}

sub each_seq_with_id {
    my $self = shift;
    $self->deprecated("each_seq_with_id - deprecated method. Use next_Seq_with_id() instead.");
    $self->next_Seq_with_id(@_);
}

sub next_Seq_with_id {
    my $self = shift;
    my $id = shift;

    $self->throw("Method next_Seq_with_id needs a sequence name argument")
	unless defined $id;

    my (@arr, $seq);

    if (exists($self->{'_start_end_lists'}->{$id})) {
	@arr = @{$self->{'_start_end_lists'}->{$id}};
    }
    return @arr;
}

=head2 get_Seq_by_pos

 Title     : get_Seq_by_pos
 Usage     : $seq = $aln->get_Seq_by_pos(3) # third sequence from the alignment
 Function  : Gets a sequence based on its position in the alignment.
             Numbering starts from 1.  Sequence positions larger than
             num_sequences() will thow an error.
 Returns   : a Bio::LocatableSeq object
 Args      : positive integer for the sequence position

=cut

sub get_seq_by_pos {
	my $self = shift;
	$self->deprecated("get_seq_by_pos - deprecated method. Use get_Seq_by_pos() instead.");
	$self->get_Seq_by_pos(@_);	
}

sub get_Seq_by_pos {

    my $self = shift;
    my ($pos) = @_;

    $self->throw("Sequence position has to be a positive integer, not [$pos]")
	unless $pos =~ /^\d+$/ and $pos > 0;
    $self->throw("No sequence at position [$pos]")
	unless $pos <= $self->num_sequences ;

    my $nse = $self->{'_order'}->{--$pos};
    return $self->{'_seq'}->{$nse};
}

=head2 get_Seq_by_id

 Title     : get_Seq_by_id
 Usage     : $seq = $aln->get_Seq_by_id($name) # seq named $name
 Function  : Gets a sequence based on its name.
             Sequences that do not exist will warn and return undef
 Returns   : a Bio::LocatableSeq object
 Args      : string for sequence name

=cut

sub get_seq_by_id {
	my $self = shift;
	$self->deprecated("get_seq_by_id - deprecated method. Use get_Seq_by_id() instead.");
	$self->get_Seq_by_id(@_);	
}


sub get_Seq_by_id {
    my ($self,$name) = @_;
    unless( defined $name ) {
      $self->warn("Must provide a sequence name");
      return;
    }
    for my $seq ( values %{$self->{'_seq'}} ) {
      if ( $seq->id eq $name) {
	return $seq;
      }
    }
    return;
}


=head2 select_Seqs

 Title     : select_Seqs
 Usage     : $aln2 = $aln->select_Seqs([1,5..10,15]) # three first sequences
 Function  : Creates a new alignment from a subset of
             sequences.  Numbering starts from 1.  Sequence positions
             larger than num_sequences() will thow an error.
 Returns   : a Bio::SimpleAlign object
 Args      : An reference list of positive integers for the selected 
             sequences. An optional parameter can be defined to toggle the 
             coordinate selection.
             

=cut

sub select_Seqs {
	my $self=shift;
	my ($sel, $toggle) = $self->_rearrange([qw(SELECTION TOGGLE)], @_);
	@{$sel}=sort {$a<=>$b} @{$sel};

   $self->throw("Select start has to be a positive integer, not [".$sel->[0]."]")
	unless $sel->[0] =~ /^\d+$/ and $sel->[0] > 0;
    $self->throw("Select end has to be a positive integer, not [".$sel->[$#$sel]."]")
	unless $sel->[$#$sel] =~ /^\d+$/ and $sel->[$#$sel] > 0;
	
	my @positions;
	if($toggle) {
		@positions=@{_toggle_selection($sel,$self->num_sequences)};
	}
	else {
		@positions=@{$sel};
	}
	
	my $aln = $self->new;
	foreach my $pos (@positions) {
		$aln->add_Seq($self->get_Seq_by_pos($pos));
	}
	$aln->id($self->id);
	# fix for meta, sf, ann    
	return $aln;	
}


sub select {
    my $self = shift;
    $self->deprecated("select - deprecated method. Use select_Seqs() instead.");
    $self->select_Seqs([$_[0]..$_[1]]);
}

sub select_noncont {
	my $self = shift;
	$self->deprecated("select_noncont - deprecated method. Use select_Seqs() instead.");
	if($_[0] eq 'nosort') {
    	shift;
   }
   $self->select_Seqs([@_]);
}

=head2 select_columns

 Title     : select_columns
 Usage     : $newaln = $aln->select_columns([20..30,45..48]);
             $newaln = $aln->select_columns(['mismatch']);
 Function  : Creates a slice from the alignment from the selected columns. 
             The first column in the alignment is denoted 1.
             Sequences with no residues in the slice are excluded from the
             new alignment and a warning is printed. Slice beyond the length
             of the sequence does not do padding.
 Returns   : A Bio::SimpleAlign object
 Args      : Positive integers for the selected colums, or specified type
             ('match'|'weak'|'strong'|'mismatch')
             First optional boolean can be defined to toggle the coordinate 
             selection. Second optional boolean which if true will keep gap-only 
             columns in the newly created slice. Example:

             $aln2 = $aln->select_columns([20..30],0,1)
             or $aln2 = $aln->select_columns(-selection=>[20..30],-toggle=>0,-keepgaponly=>1)

=cut

sub slice {
	my $self = shift;
	$self->deprecated("slice - deprecated method. Use select_columns() instead.");
	$self->select_columns([$_[0]..$_[1]],0,defined($_[2])?$_[2]:0);	
}


sub select_columns {
	my $self = shift;
	#my ($start, $end, $keep_gap_only) = @_;
	my ($sel, $toggle,$keep_gap_only) = $self->_rearrange([qw(SELECTION TOGGLE KEEPGAPONLY)], @_);
	
	if($sel->[0]=~/^[a-z]/i) {
		$sel=$self->_select_columns_by_type($sel);
	}
	
	@{$sel}=sort {$a<=$b} @{$sel};

	#warnings
	$self->throw("select_columns start has to be a positive integer, not [".$sel->[0]."]")
	  unless $sel->[0] =~ /^\d+$/ and $sel->[0] > 0;
	$self->throw("select_columns end has to be a positive integer, not [".$sel->[$#$sel]."]")
	  unless $sel->[$#$sel] =~ /^\d+$/ and $sel->[$#$sel] > 0;
	$self->throw("This alignment has only ". $self->length . " residues. select_columns start " .
					 "[".$sel->[0]."] is too big.") if $sel->[0] > $self->length;

	my $newcoords;
	
	if($toggle) {
		$newcoords=_cont_coords(_toggle_selection($sel,$self->length));
	}
	else {
		$newcoords=_cont_coords($sel);
	}

    my $cons_meta = $self->consensus_meta;
	my $aln = $self->new;
	$aln->id($self->id);
	my ($start,$end)=($newcoords->[0],$newcoords->[$#$newcoords]);
	
	foreach my $seq ( $self->next_Seq() ) {
	    my $new_seq = $seq->isa('Bio::Seq::MetaI') ?
            Bio::Seq::Meta->new
        (-id      => $seq->id,
		 -alphabet => $seq->alphabet,
		 -strand  => $seq->strand,
		 -verbose => $self->verbose) :
            Bio::LocatableSeq->new
        (-id      => $seq->id,
		 -alphabet => $seq->alphabet,
		 -strand  => $seq->strand,
		 -verbose => $self->verbose);
        
	    # seq
	    
	    my $seq_end = $end;
	    $seq_end = $seq->length if( $end > $seq->length );
	    # the actual selected seq

	    my $slice_seq;
      for(my $num=0;$num<$#$newcoords;) {
	    	$slice_seq.= $seq->subseq($newcoords->[$num], $newcoords->[$num+1]);
	    	$num+=2;
	   }
	    
	    $new_seq->seq( $slice_seq );

		#How to calculate the new start and new end
		#new_start=old_start+pre_start-1
		#new_end=old_end-after_end+1
		#in order to pass the $seq->end test, we should assign a fake start first
		
		#all the non residue chars
		my $nonres = join("",$self->gap_char, $self->match_char,$self->missing_char);
	   
	   my $slice_length=($slice_seq =~s/([^\Q$nonres\E])/$1/g);
	   my ($newstart, $newend, $fakestart, $oldend);
	   my ($pre_start_seq,$pre_start_length,$after_end_seq,$after_end_length);
	   
	   
	   if($end>=CORE::length($seq->seq())) {
	   	$after_end_length=0;
	   }
	   else {
	   	$after_end_seq = $seq->subseq($end+1,CORE::length($seq->seq()));	
	   	$after_end_length=($after_end_seq=~s/([^\Q$nonres\E])/$1/g); #need to check
	   }
		
	   if($start>1) {
	   	
		   my $pre_start_seq = $seq->subseq(1, $start - 1);
			my $pre_start_length=($pre_start_seq=~s/([^\Q$nonres\E])/$1/g); #need to check
			if (defined($seq->strand) && $seq->strand < 0){
				$newend=$seq->end-$pre_start_length;
				$fakestart=$newend-$slice_length+1;
				$newstart=$seq->start+$after_end_length;				
			}
			else {
				$newend=$seq->end-$after_end_length;
				$fakestart=$newend-$slice_length+1;
				$newstart=$seq->start+$pre_start_length;
			}
		}
		else {
			if (defined($seq->strand) && $seq->strand < 0){
				$newend=$seq->end;
				$fakestart=$newend-$slice_length+1;
				$newstart=$seq->start+$after_end_length;
			}
			else {
				$newend=$seq->end-$after_end_length;
				$fakestart=$newend-$slice_length+1;
				$newstart=$seq->start;
			}			
		}
		$new_seq->start($fakestart);
		$new_seq->end($newend);
		$new_seq->start($newstart);

        if ($new_seq->isa('Bio::Seq::MetaI')) {
            for my $meta_name ($seq->meta_names) {
                $new_seq->named_meta($meta_name, $seq->named_submeta($meta_name, $start, $end));
            }
        }

	    if ($new_seq->start and $new_seq->end >= $new_seq->start) {
            $aln->add_Seq($new_seq);
	    } else {
            if( $keep_gap_only ) {
                $aln->add_Seq($new_seq);
            } else {
                my $nse = $seq->get_nse();
                $self->warn("Slice [$start-$end] of sequence [$nse] contains no residues.".
                    " Sequence excluded from the new alignment.");
            }
	    }
	}
    if ($cons_meta) {
        my $new = Bio::Seq::Meta->new();
        for my $meta_name ($cons_meta->meta_names) {
            $new->named_meta($meta_name, $cons_meta->named_submeta($meta_name, $start, $end));
        }
        $aln->consensus_meta($new);
    }
    $aln->annotation($self->annotation);
    # fix for meta, sf, ann
	return $aln;
}


sub _select_columns_by_type {
	#code copied from _remove_columns_by_type to enable selecting columns by type
	my ($self,$type) = @_;
	my @selects;

	my $gap = $self->gap_char if (grep { $_ eq 'gaps'} @{$type});
	my $all_gaps_columns = $self->gap_char if (grep /all_gaps_columns/,@{$type});
	my %matchchars = ( 'match'           => '\*',
                       'weak'             => '\.',
                       'strong'           => ':',
                       'mismatch'         => ' ',
                     );
	# get the characters to delete against
	my $del_char;
	foreach my $type (@{$type}){
		$del_char.= $matchchars{$type};
	}

	my $length = 0;
	my $match_line = $self->match_line;
	# do the matching to get the segments to remove
	if($del_char){
		while($match_line =~ m/[$del_char]/g ){
			my $start = pos($match_line);
			$match_line=~/\G[$del_char]+/gc;
			my $end = pos($match_line);

			#have to offset the start and end for subsequent removes
			push @selects, $start..$end;
		}
	}

	# remove the segments
	return \@selects;
}



=head2 remove_gaps

 Title     : remove_gaps
 Usage     : $aln2 = $aln->remove_gaps(-reference=>5)
 Function  : Creates an aligment with gaps removed
 Returns   : a Bio::SimpleAlign object
 Args      : -GAPCHAR a gap character(optional) if none specified taken
                from $self->gap_char,
             -ALLGAPCOL $all_gaps_columns flag (1 or 0, default is 0)
                 indicates that only all-gaps columns should be deleted
             -REFERENCE splices all aligned sequences where the specified 
                 sequence has gaps.

=cut

sub remove_gaps {
	my $self=shift @_;
	my ($gapchar,$all_gaps_columns,$reference) = $self->_rearrange([qw(GAPCHAR ALLGAPCOL REFERENCE)], @_);
	my $gap_line;
   if($reference) {
   	my $refseq=$self->get_Seq_by_pos($reference);
    	$gap_line=$refseq->seq();
   }
   else {
		if ($all_gaps_columns) {
		  $gap_line = $self->all_gap_line($gapchar);
		} else {
		  $gap_line = $self->gap_line($gapchar);
		}
	}
    
    my $aln = $self->new;

    my $length = 0;
    my $del_char = $gapchar || $self->gap_char;
    # Do the matching to get the segments to remove
    my $removed_cols;
    while ($gap_line =~ m/[$del_char]/g) {
        push @{$removed_cols}, pos($gap_line);
    }
    #remove the segments
    $aln = $#$removed_cols >= 0 ? $self->select_columns($removed_cols,1) : $self;
    # fix for meta, sf, ann        
    return $aln;
}

sub splice_by_seq_pos{
    my $self = shift;
    $self->deprecated("splice_by_seq_pos - deprecated method. Use remove_gaps() instead.");
    $self->remove_gaps(-reference=>$_[0]);
}

=head2 mask_columns

 Title     : mask_columns
 Usage     : $aln2 = $aln->mask_columns([2,5,7..10],1)
 Function  : Masks slices of the alignment inclusive of the defined
             columns, and the first column in the alignment is denoted 1.
             Mask beyond the length of the sequence does not do padding.
 Returns   : A Bio::SimpleAlign object
 Args      : Positive integers should be used to defined the column numbers
             The mask character should be defined by $aln->mask_char() or "?" as 
             default. An optional parameter can be defined to toggle the 
             coordinate selection.
 Note      : 

=cut

sub mask_columns {
    #based on slice(), but did not include the Bio::Seq::Meta sections as I was not sure what it is doing
	my $self=shift;
	my ($sel, $toggle) = $self->_rearrange([qw(SELECTION TOGGLE)], @_);
	@{$sel}=sort {$a<=>$b} @{$sel};
	
	my $nonres = join("",$self->gap_char, $self->match_char,$self->missing_char);
	my $mask_char=$self->mask_char;
   
	#Exceptions
	$self->throw("Mask start has to be a positive integer and less than ".
                 "alignment length, not [".$sel->[0]."]")
	unless $sel->[0] =~ /^\d+$/ && $sel->[0] > 0 && $sel->[0] <= $self->length;	

   $self->throw("Mask end has to be a positive integer and less than ".
                 "alignment length, not [".$sel->[$#$sel]."]")
	unless $sel->[$#$sel] =~ /^\d+$/ && $sel->[$#$sel] > 0 && $sel->[$#$sel] <= $self->length;                 
	
	#calculate the coords, and make it in a continous way
	my @newcoords;
	if($toggle) {
		@newcoords=@{_cont_coords(_toggle_selection($sel,$self->length))};
	}
	else {
		@newcoords=@{_cont_coords($sel)};
	}
	
	my $aln = $self->new;
	$aln->id($self->id);
	

	foreach my $seq ( $self->next_Seq() ) {
		my $new_seq = Bio::LocatableSeq->new(-id => $seq->id,
				-alphabet => $seq->alphabet,
				-strand  => $seq->strand,
				-verbose => $self->verbose);
		my $dna_string=$seq->seq();
		my $masked_string;
		for(my $num=0;$num<@newcoords;) {
			my ($start,$end)=@newcoords[$num,$num+1];
			$masked_string = substr($dna_string, $start - 1, $end - $start + 1);
         $masked_string =~ s{[^$nonres]}{$mask_char}g;
         $dna_string = substr($dna_string,0,$start-1) . $masked_string . substr($dna_string,$end);    
			$num+=2;
		}
		$new_seq->seq($dna_string);
		$aln->add_Seq($new_seq);
	}
	return $aln;
}
=head2 seq_with_features

 Title   : seq_with_features
 Usage   : $seq = $aln->seq_with_features(-pos => 1,
                                          -consensus => 60
                                          -mask =>
           sub { my $consensus = shift;

                 for my $i (1..5){
                    my $n = 'N' x $i;
                    my $q = '\?' x $i;
                    while($consensus =~ /[^?]$q[^?]/){
                       $consensus =~ s/([^?])$q([^?])/$1$n$2/;
                    }
                  }
                 return $consensus;
               }
                                         );
 Function: produces a Bio::Seq object by first splicing gaps from -pos
           (by means of a remove_gaps(-reference=>1) call), then creating
           features using non-? chars (by means of a consensus_string()
           call with stringency -consensus).
 Returns : a Bio::Seq object
 Args    : -pos : required. sequence from which to build the Bio::Seq
             object
           -consensus : optional, defaults to consensus_string()'s
             default cutoff value
           -mask : optional, a coderef to apply to consensus_string()'s
             output before building features.  this may be useful for
             closing gaps of 1 bp by masking over them with N, for
             instance

=cut

sub seq_with_features{
   my ($self,%arg) = @_;

   #first do the preparatory splice
   $self->throw("must provide a -pos argument") unless $arg{-pos};
   $self->remove_gaps(-reference=>$arg{-pos});

   my $consensus_string = $self->consensus_string($arg{-consensus});
   $consensus_string = $arg{-mask}->($consensus_string)
	 if defined($arg{-mask});

   my(@bs,@es);

   push @bs, 1 if $consensus_string =~ /^[^?]/;

   while($consensus_string =~ /\?[^?]/g){
	 push @bs, pos($consensus_string);
   }
   while($consensus_string =~ /[^?]\?/g){
	 push @es, pos($consensus_string);
   }

   push @es, CORE::length($consensus_string) if $consensus_string =~ /[^?]$/;

   my $seq = Bio::Seq->new();

#   my $rootfeature = Bio::SeqFeature::Generic->new(
#                -source_tag => 'location',
#                -start      => $self->get_seq_by_pos($arg{-pos})->start,
#                -end        => $self->get_seq_by_pos($arg{-pos})->end,
#                                                  );
#   $seq->add_SeqFeature($rootfeature);

   while(my $b = shift @bs){
	 my $e = shift @es;
	 $seq->add_SeqFeature(
       Bio::SeqFeature::Generic->new(
         -start => $b - 1 + $self->get_Seq_by_pos($arg{-pos})->start,
         -end   => $e - 1 + $self->get_Seq_by_pos($arg{-pos})->start,
         -source_tag => $self->source || 'MSA',
       )
     );
   }

   return $seq;
}


=head1 Change sequences within the MSA

These methods affect characters in all sequences without changing the
alignment.


=head2 map_chars

 Title     : map_chars
 Usage     : $ali->map_chars('\.','-')
 Function  : Does a s/$arg1/$arg2/ on the sequences. Useful for gap
             characters

             Notice that the from (arg1) is interpretted as a regex,
             so be careful about quoting meta characters (eg
             $ali->map_chars('.','-') wont do what you want)
 Returns   :
 Argument  : 'from' rexexp
             'to' string

=cut

sub map_chars {
    my $self = shift;
    my $from = shift;
    my $to   = shift;
    my ($seq,$temp);

    $self->throw("Need exactly two arguments")
	unless defined $from and defined $to;

    foreach $seq ( $self->next_Seq() ) {
	$temp = $seq->seq();
	$temp =~ s/$from/$to/g;
	$seq->seq($temp);
    }
    return 1;
}


=head2 uppercase

 Title     : uppercase()
 Usage     : $ali->uppercase()
 Function  : Sets all the sequences to uppercase
 Returns   : 1
 Argument  : None

=cut

sub uppercase {
    my $self = shift;
    my $seq;
    my $temp;

    foreach $seq ( $self->next_Seq() ) {
      $temp = $seq->seq();
      $temp =~ tr/[a-z]/[A-Z]/;

      $seq->seq($temp);
    }
    return 1;
}

=head2 lowercase

 Title     : lowercase()
 Usage     : $ali->lowercase()
 Function  : Sets all the sequences to lowercase
 Returns   : 1
 Argument  : None

=cut

sub lowercase {
    my $self = shift;
    my $seq;
    my $temp;

    foreach $seq ( $self->next_Seq() ) {
      $temp = $seq->seq();
      $temp =~ tr/[A-Z]/[a-z]/;

      $seq->seq($temp);
    }
    return 1;
}


=head2 togglecase

 Title     : togglecase()
 Usage     : $ali->togglecase()
 Function  : Sets all the sequences to opposite case
 Returns   : 1
 Argument  : None

=cut

sub togglecase {
    my $self = shift;
    my $seq;
    my $temp;

    foreach $seq ( $self->next_Seq() ) {
      $temp = $seq->seq();
      $temp =~ tr/[A-Za-z]/[a-zA-Z]/;

      $seq->seq($temp);
    }
    return 1;
}


=head2 match

 Title     : match()
 Usage     : $aln->match()
 Function  : Goes through all columns and changes residues that are
             identical to residue in first sequence to match '.'
             character. Sets match_char.

             USE WITH CARE: Most MSA formats do not support match
             characters in sequences, so this is mostly for output
             only. NEXUS format (Bio::AlignIO::nexus) can handle
             it.
 Returns   : 1
 Argument  : a match character, optional, defaults to '.'
             If the character is defined, it will reset $aln->match_char 

=cut

sub match {
    my ($self, $match) = @_;

   if(defined $match) {
		$self->match_char($match);
	}
	else {
		$match=$self->match_char();
	}
    
    my ($matching_char) = $match;
    $matching_char = "\\$match" if $match =~ /[\^.$|()\[\]]/ ;  #';
    $self->map_chars($matching_char, '-');

    my @seqs = $self->next_Seq();
    return 1 unless scalar @seqs > 1;

    my $refseq = shift @seqs ;
    my @refseq = split //, $refseq->seq;
    my $gapchar = $self->gap_char;

    foreach my $seq ( @seqs ) {
	my @varseq = split //, $seq->seq();
	for ( my $i=0; $i < scalar @varseq; $i++) {
	    $varseq[$i] = $match if defined $refseq[$i] &&
		( $refseq[$i] =~ /[A-Za-z\*]/ ||
		  $refseq[$i] =~ /$gapchar/ )
		      && $refseq[$i] eq $varseq[$i];
	}
	$seq->seq(join '', @varseq);
    }
    
    return 1;
}


=head2 unmatch

 Title     : unmatch()
 Usage     : $ali->unmatch()
 Function  : Undoes the effect of method match. Unsets match_char.
 Returns   : 1
 Argument  : a match character, optional, defaults to '.'

See L<match> and L<match_char>

=cut

sub unmatch {
    my ($self, $match) = @_;

   if(defined $match) {
		$self->match_char($match);
	}
	else {
		$match=$self->match_char();
	}

    my @seqs = $self->next_Seq();
    return 1 unless scalar @seqs > 1;

    my $refseq = shift @seqs ;
    my @refseq = split //, $refseq->seq;
    my $gapchar = $self->gap_char;
    foreach my $seq ( @seqs ) {
	my @varseq = split //, $seq->seq();
	for ( my $i=0; $i < scalar @varseq; $i++) {
	    $varseq[$i] = $refseq[$i] if defined $refseq[$i] &&
		( $refseq[$i] =~ /[A-Za-z\*]/ ||
		  $refseq[$i] =~ /$gapchar/ ) &&
		      $varseq[$i] eq $match;
	}
	$seq->seq(join '', @varseq);
    }
    $self->match_char('');
    return 1;
}


=head1 Consensus sequences

Methods to calculate consensus sequences for the MSA

=head2 consensus_string

 Title     : consensus_string
 Usage     : $str = $ali->consensus_string($threshold_percent)
 Function  : Makes a strict consensus
 Returns   : Consensus string
 Argument  : Optional treshold ranging from 0 to 100.
             The consensus residue has to appear at least threshold %
             of the sequences at a given location, otherwise a '?'
             character will be placed at that location.
             (Default value = 0%)

=cut

sub consensus_string {
    my $self = shift;
    my $threshold = shift;

    my $out = "";
    my $len = $self->length - 1;

    foreach ( 0 .. $len ) {
	$out .= $self->_consensus_aa($_,$threshold);
    }
    return $out;
}

sub _consensus_aa {
    my $self = shift;
    my $point = shift;
    my $threshold_percent = shift || -1 ;
    my ($seq,%hash,$count,$letter,$key);
    my $gapchar = $self->gap_char;
    foreach $seq ( $self->next_Seq() ) {
	$letter = substr($seq->seq,$point,1);
	$self->throw("--$point-----------") if $letter eq '';
	($letter eq $gapchar || $letter =~ /\./) && next;
	# print "Looking at $letter\n";
	$hash{$letter}++;
    }
    my $number_of_sequences = $self->num_sequences();
    my $threshold = $number_of_sequences * $threshold_percent / 100. ;
    $count = -1;
    $letter = '?';

    foreach $key ( sort keys %hash ) {
	# print "Now at $key $hash{$key}\n";
	if( $hash{$key} > $count && $hash{$key} >= $threshold) {
	    $letter = $key;
	    $count = $hash{$key};
	}
    }
    return $letter;
}


=head2 consensus_iupac

 Title     : consensus_iupac
 Usage     : $str = $ali->consensus_iupac()
 Function  : Makes a consensus using IUPAC ambiguity codes from DNA
             and RNA. The output is in upper case except when gaps in
             a column force output to be in lower case.

             Note that if your alignment sequences contain a lot of
             IUPAC ambiquity codes you often have to manually set
             alphabet.  Bio::PrimarySeq::_guess_type thinks they
             indicate a protein sequence.
 Returns   : consensus string
 Argument  : none
 Throws    : on protein sequences

=cut

sub consensus_iupac {
    my $self = shift;
    my $out = "";
    my $len = $self->length-1;

    # only DNA and RNA sequences are valid
    foreach my $seq ( $self->next_Seq() ) {
	$self->throw("Seq [". $seq->get_nse. "] is a protein")
	    if $seq->alphabet eq 'protein';
    }
    # loop over the alignment columns
    foreach my $count ( 0 .. $len ) {
	$out .= $self->_consensus_iupac($count);
    }
    return $out;
}

sub _consensus_iupac {
    my ($self, $column) = @_;
    my ($string, $char, $rna);

    #determine all residues in a column
    foreach my $seq ( $self->next_Seq() ) {
	$string .= substr($seq->seq, $column, 1);
    }
    $string = uc $string;

    # quick exit if there's an N in the string
    if ($string =~ /N/) {
	$string =~ /\W/ ? return 'n' : return 'N';
    }
    # ... or if there are only gap characters
    return '-' if $string =~ /^\W+$/;

    # treat RNA as DNA in regexps
    if ($string =~ /U/) {
	$string =~ s/U/T/;
	$rna = 1;
    }

    # the following s///'s only need to be done to the _first_ ambiguity code
    # as we only need to see the _range_ of characters in $string

    if ($string =~ /[VDHB]/) {
	$string =~ s/V/AGC/;
	$string =~ s/D/AGT/;
	$string =~ s/H/ACT/;
	$string =~ s/B/CTG/;
    }

    if ($string =~ /[SKYRWM]/) {
	$string =~ s/S/GC/;
	$string =~ s/K/GT/;
	$string =~ s/Y/CT/;
	$string =~ s/R/AG/;
	$string =~ s/W/AT/;
	$string =~ s/M/AC/;
    }

    # and now the guts of the thing

    if ($string =~ /A/) {
        $char = 'A';                     # A                      A
        if ($string =~ /G/) {
            $char = 'R';                 # A and G (purines)      R
            if ($string =~ /C/) {
                $char = 'V';             # A and G and C          V
                if ($string =~ /T/) {
                    $char = 'N';         # A and G and C and T    N
                }
            } elsif ($string =~ /T/) {
                $char = 'D';             # A and G and T          D
            }
        } elsif ($string =~ /C/) {
            $char = 'M';                 # A and C                M
            if ($string =~ /T/) {
                $char = 'H';             # A and C and T          H
            }
        } elsif ($string =~ /T/) {
            $char = 'W';                 # A and T                W
        }
    } elsif ($string =~ /C/) {
        $char = 'C';                     # C                      C
        if ($string =~ /T/) {
            $char = 'Y';                 # C and T (pyrimidines)  Y
            if ($string =~ /G/) {
                $char = 'B';             # C and T and G          B
            }
        } elsif ($string =~ /G/) {
            $char = 'S';                 # C and G                S
        }
    } elsif ($string =~ /G/) {
        $char = 'G';                     # G                      G
        if ($string =~ /C/) {
            $char = 'S';                 # G and C                S
        } elsif ($string =~ /T/) {
            $char = 'K';                 # G and T                K
        }
    } elsif ($string =~ /T/) {
        $char = 'T';                     # T                      T
    }

    $char = 'U' if $rna and $char eq 'T';
    $char = lc $char if $string =~ /\W/;

    return $char;
}


=head2 consensus_meta

 Title     : consensus_meta
 Usage     : $seqmeta = $ali->consensus_meta()
 Function  : Returns a Bio::Seq::Meta object containing the consensus
             strings derived from meta data analysis.
 Returns   : Bio::Seq::Meta 
 Argument  : Bio::Seq::Meta 
 Throws    : non-MetaI object

=cut

sub consensus_meta {
    my ($self, $meta) = @_;
    if ($meta && (!ref $meta || !$meta->isa('Bio::Seq::MetaI'))) {
        $self->throw('Not a Bio::Seq::MetaI object');
    }
    return $self->{'_aln_meta'} = $meta if $meta;
    return $self->{'_aln_meta'} 
}


=head2 bracket_string

 Title     : bracket_string
 Usage     : my @params = (-refseq     => 'testseq',
                           -allele1    => 'allele1',
                           -allele2    => 'allele2',
                           -delimiters => '{}',
                           -separator  => '/');
             $str = $aln->bracket_string(@params)

 Function :  When supplied with a list of parameters (see below), returns a
             string in BIC format. This is used for allelic comparisons.
             Briefly, if either allele contains a base change when compared to
             the refseq, the base or gap for each allele is represented in
             brackets in the order present in the 'alleles' parameter.

             For the following data:

             >testseq
             GGATCCATTGCTACT
             >allele1
             GGATCCATTCCTACT
             >allele2
             GGAT--ATTCCTCCT

             the returned string with parameters 'refseq => testseq' and
             'alleles => [qw(allele1 allele2)]' would be:

             GGAT[C/-][C/-]ATT[C/C]CT[A/C]CT
 Returns   : BIC-formatted string
 Argument  : Required args
                refseq    : string (ID) of the reference sequence used
                            as basis for comparison
                allele1   : string (ID) of the first allele
                allele2   : string (ID) of the second allele
             Optional args
                delimiters: two symbol string of left and right delimiters.
                            Only the first two symbols are used
                            default = '[]'
                separator : string used as a separator.  Only the first
                            symbol is used
                            default = '/'
 Throws    : On no refseq/alleles, or invalid refseq/alleles.

=cut

sub bracket_string {
    my ($self, @args) = @_;
    my ($ref, $a1, $a2, $delim, $sep) =
        $self->_rearrange([qw(refseq allele1 allele2 delimiters separator)], @args);
    $self->throw('Missing refseq/allele1/allele2') if (!$a1 || !$a2 || !$ref);
    my ($ld, $rd);
    ($ld, $rd) = split('', $delim, 2) if $delim;
    $ld ||= '[';
    $rd ||= ']';
    $sep ||= '/';
    my ($refseq, $allele1, $allele2) =
        map {( $self->next_Seq_with_id($_) )} ($ref, $a1, $a2);
    if (!$refseq || !$allele1 || !$allele2) {
        $self->throw("One of your refseq/allele IDs is invalid!");
    }
    my $len = $self->length-1;
    my $bic = '';
    # loop over the alignment columns
    for my $column ( 0 .. $len ) {
        my $string;
        my ($compres, $res1, $res2) =
            map{substr($_->seq, $column, 1)} ($refseq, $allele1, $allele2);
        # are any of the allele symbols different from the refseq?
        $string = ($compres eq $res1 && $compres eq $res2) ? $compres :
                $ld.$res1.$sep.$res2.$rd;
        $bic .= $string;
    }
    return $bic;
}


=head2 cigar_line

 Title    : cigar_line()
 Usage    : %cigars = $align->cigar_line()
 Function : Generates a "cigar" (Compact Idiosyncratic Gapped Alignment
            Report) line for each sequence in the alignment. Examples are
            "1,60" or "5,10:12,58", where the numbers refer to conserved
            positions within the alignment. The keys of the hash are the
            NSEs (name/start/end) assigned to each sequence.
 Args     : threshold (optional, defaults to 100)
 Returns  : Hash of strings (cigar lines)

=cut

sub cigar_line {
	my $self = shift;
	my $thr=shift||100;
	my %cigars;

	my @consensus = split "",($self->consensus_string($thr));
	my $len = $self->length;
	my $gapchar = $self->gap_char;

	# create a precursor, something like (1,4,5,6,7,33,45),
	# where each number corresponds to a conserved position
	foreach my $seq ( $self->next_Seq ) {
		my @seq = split "", uc ($seq->seq);
		my $pos = 1;
		for (my $x = 0 ; $x < $len ; $x++ ) {
			if ($seq[$x] eq $consensus[$x]) {
				push @{$cigars{$seq->get_nse}},$pos;
				$pos++;
			} elsif ($seq[$x] ne $gapchar) {
				$pos++;
			}
		}
	}
	# duplicate numbers - (1,4,5,6,7,33,45) becomes (1,1,4,5,6,7,33,33,45,45)
	for my $name (keys %cigars) {
		splice @{$cigars{$name}}, 1, 0, ${$cigars{$name}}[0] if
		  ( ${$cigars{$name}}[0] + 1 < ${$cigars{$name}}[1] );
      push @{$cigars{$name}}, ${$cigars{$name}}[$#{$cigars{$name}}] if
           ( ${$cigars{$name}}[($#{$cigars{$name}} - 1)] + 1 <
		          ${$cigars{$name}}[$#{$cigars{$name}}] );
		for ( my $x = 1 ; $x < $#{$cigars{$name}} - 1 ; $x++) {
			if (${$cigars{$name}}[$x - 1] + 1 < ${$cigars{$name}}[$x]  &&
		       ${$cigars{$name}}[$x + 1]  > ${$cigars{$name}}[$x] + 1) {
	         splice @{$cigars{$name}}, $x, 0, ${$cigars{$name}}[$x];
			}
      }
	}
  # collapse series - (1,1,4,5,6,7,33,33,45,45) becomes (1,1,4,7,33,33,45,45)
  for my $name (keys %cigars) {
	  my @remove;
	  for ( my $x = 0 ; $x < $#{$cigars{$name}} ; $x++) {
		   if ( ${$cigars{$name}}[$x] == ${$cigars{$name}}[($x - 1)] + 1 &&
			     ${$cigars{$name}}[$x] == ${$cigars{$name}}[($x + 1)] - 1 ) {
		      unshift @remove,$x;
	      }
	   }
      for my $pos (@remove) {
		  	splice @{$cigars{$name}}, $pos, 1;
	   }
   }
   # join and punctuate
   for my $name (keys %cigars) {
 	  my ($start,$end,$str) = "";
 	  while ( ($start,$end) = splice @{$cigars{$name}}, 0, 2 ) {
 		  $str .= ($start . "," . $end . ":");
 	  }
 	  $str =~ s/:$//;
      $cigars{$name} = $str;
   }
   %cigars;
}


=head2 match_line

 Title    : match_line()
 Usage    : $line = $aln->match_line()
 Function : Generates a match line - much like consensus string
            except that a line indicating the '*' for a match.
 Args     : (optional) Match line characters ('*' by default)
            (optional) Strong match char (':' by default)
            (optional) Weak match char ('.' by default)
 Returns  : String

=cut

sub match_line {
	my ($self,$matchlinechar, $strong, $weak) = @_;
	my %matchchars = ('match'    => $matchlinechar || '*',
							  'weak'     => $weak          || '.',
							  'strong'   => $strong        || ':',
							  'mismatch' => ' ',
						  );

	my @seqchars;
	my $alphabet;
	foreach my $seq ( $self->next_Seq ) {
		push @seqchars, [ split(//, uc ($seq->seq)) ];
		$alphabet = $seq->alphabet unless defined $alphabet;
	}
	my $refseq = shift @seqchars;
	# let's just march down the columns
	my $matchline;
 POS:
	foreach my $pos ( 0..$self->length ) {
		my $refchar = $refseq->[$pos];
		my $char = $matchchars{'mismatch'};
		unless( defined $refchar ) {
			last if $pos == $self->length; # short circuit on last residue
			# this in place to handle jason's soon-to-be-committed
			# intron mapping code
			goto bottom;
		}
		my %col = ($refchar => 1);
		my $dash = ($refchar eq '-' || $refchar eq '.' || $refchar eq ' ');
		foreach my $seq ( @seqchars ) {
			next if $pos >= scalar @$seq;
			$dash = 1 if( $seq->[$pos] eq '-' || $seq->[$pos] eq '.' ||
							  $seq->[$pos] eq ' ' );
			$col{$seq->[$pos]}++ if defined $seq->[$pos];
		}
		my @colresidues = sort keys %col;

		# if all the values are the same
		if( $dash ) { $char =  $matchchars{'mismatch'} }
		elsif( @colresidues == 1 ) { $char = $matchchars{'match'} }
		elsif( $alphabet eq 'protein' ) { # only try to do weak/strong
			# matches for protein seqs
	    TYPE:
			foreach my $type ( qw(strong weak) ) {
				# iterate through categories
				my %groups;
				# iterate through each of the aa in the col
				# look to see which groups it is in
				foreach my $c ( @colresidues ) {
					foreach my $f ( grep { index($_,$c) >= 0 } @{$CONSERVATION_GROUPS{$type}} ) {
						push @{$groups{$f}},$c;
					}
				}
			 GRP:
				foreach my $cols ( values %groups ) {
					@$cols = sort @$cols;
					# now we are just testing to see if two arrays
					# are identical w/o changing either one
					# have to be same len
					next if( scalar @$cols != scalar @colresidues );
					# walk down the length and check each slot
					for($_=0;$_ < (scalar @$cols);$_++ ) {
						next GRP if( $cols->[$_] ne $colresidues[$_] );
					}
					$char = $matchchars{$type};
					last TYPE;
				}
			}
		}
	 bottom:
		$matchline .= $char;
	}
	return $matchline;
}


=head2 gap_line

 Title    : gap_line()
 Usage    : $line = $aln->gap_line()
 Function : Generates a gap line - much like consensus string
            except that a line where '-' represents gap
 Args     : (optional) gap line characters ('-' by default)
 Returns  : string

=cut

sub gap_line {
    my ($self,$gapchar) = @_;
    $gapchar = $gapchar || $self->gap_char;
    my %gap_hsh; # column gaps vector
    foreach my $seq ( $self->next_Seq ) {
		my $i = 0;
    	map {$gap_hsh{$_->[0]} = undef} grep {$_->[1] eq $gapchar}
		  map {[$i++, $_]} split(//, uc ($seq->seq));
    }
    my $gap_line;
    foreach my $pos ( 0..$self->length-1 ) {
	  $gap_line .= (exists $gap_hsh{$pos}) ? $gapchar:'.';
    }
    return $gap_line;
}

=head2 all_gap_line

 Title    : all_gap_line()
 Usage    : $line = $aln->all_gap_line()
 Function : Generates a gap line - much like consensus string
            except that a line where '-' represents all-gap column
 Args     : (optional) gap line characters ('-' by default)
 Returns  : string

=cut

sub all_gap_line {
    my ($self,$gapchar) = @_;
    $gapchar = $gapchar || $self->gap_char;
    my %gap_hsh;		# column gaps counter hash
    my @seqs = $self->next_Seq;
    foreach my $seq ( @seqs ) {
	my $i = 0;
    	map {$gap_hsh{$_->[0]}++} grep {$_->[1] eq $gapchar}
	map {[$i++, $_]} split(//, uc ($seq->seq));
    }
    my $gap_line;
    foreach my $pos ( 0..$self->length-1 ) {
	if (exists $gap_hsh{$pos} && $gap_hsh{$pos} == scalar @seqs) {
            # gaps column
	    $gap_line .= $gapchar;
	} else {
	    $gap_line .= '.';
	}
    }
    return $gap_line;
}

=head2 gap_col_matrix

 Title    : gap_col_matrix()
 Usage    : my $cols = $aln->gap_col_matrix()
 Function : Generates an array of hashes where
            each entry in the array is a hash reference
            with keys of all the sequence names and
            and value of 1 or 0 if the sequence has a gap at that column
 Args     : (optional) gap line characters ($aln->gap_char or '-' by default)

=cut

sub gap_col_matrix {
    my ($self,$gapchar) = @_;
    $gapchar = $gapchar || $self->gap_char;
    my %gap_hsh; # column gaps vector
    my @cols;
    foreach my $seq ( $self->next_Seq ) {
	my $i = 0;
	my $str = $seq->seq;
	my $len = $seq->length;
	my $ch;
	my $id = $seq->display_id;
	while( $i < $len ) {
	    $ch = substr($str, $i, 1);
	    $cols[$i++]->{$id} = ($ch eq $gapchar);
	}
    }
    return \@cols;
}


=head1 MSA attributes

Methods for setting and reading the MSA attributes.

Note that the methods defining character semantics depend on the user
to set them sensibly.  They are needed only by certain input/output
methods. Unset them by setting to an empty string ('').

=head2 id

 Title     : id
 Usage     : $myalign->id("Ig")
 Function  : Gets/sets the id field of the alignment
 Returns   : An id string
 Argument  : An id string (optional)

=cut

sub id {
    my ($self, $name) = @_;

    if (defined( $name )) {
	$self->{'_id'} = $name;
    }

    return $self->{'_id'};
}

=head2 accession

 Title     : accession
 Usage     : $myalign->accession("PF00244")
 Function  : Gets/sets the accession field of the alignment
 Returns   : An acc string
 Argument  : An acc string (optional)

=cut

sub accession {
    my ($self, $acc) = @_;

    if (defined( $acc )) {
	$self->{'_accession'} = $acc;
    }

    return $self->{'_accession'};
}

=head2 description

 Title     : description
 Usage     : $myalign->description("14-3-3 proteins")
 Function  : Gets/sets the description field of the alignment
 Returns   : An description string
 Argument  : An description string (optional)

=cut

sub description {
    my ($self, $name) = @_;

    if (defined( $name )) {
	$self->{'_description'} = $name;
    }

    return $self->{'_description'};
}

=head2 source

 Title   : source
 Usage   : $aln->source($newval)
 Function: sets the Alignment source program
 Example :
 Returns : value of source
 Args    : newvalue (optional)


=cut

sub source{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_source'} = $value;
    }
    return $self->{'_source'};
}

=head2 missing_char

 Title     : missing_char
 Usage     : $myalign->missing_char("?")
 Function  : Gets/sets the missing_char attribute of the alignment
             It is generally recommended to set it to 'n' or 'N'
             for nucleotides and to 'X' for protein.
 Returns   : An missing_char string,
 Argument  : An missing_char string (optional), default as '?'

=cut

sub missing_char {
    my ($self, $char) = @_;
	
	if (defined $char || ! defined $self->{'_missing_char'} ) {
		$char= '?' unless defined $char;
		$self->throw("Single gap character, not [$char]!") if CORE::length($char) > 1;
		$self->{'_missing_char'} = $char;
	}

    return $self->{'_missing_char'};
}

=head2 match_char

 Title     : match_char
 Usage     : $myalign->match_char('.')
 Function  : Gets/sets the match_char attribute of the alignment
 Returns   : An match_char string,
 Argument  : An match_char string (optional), default as '.'

=cut

sub match_char {
    my ($self, $char) = @_;

	if (defined $char || ! defined $self->{'_match_char'} ) {
		$char= '.' unless defined $char;
		$self->throw("Single gap character, not [$char]!") if CORE::length($char) > 1;
		$self->{'_match_char'} = $char;
	}

    return $self->{'_match_char'};
}

=head2 gap_char

 Title     : gap_char
 Usage     : $myalign->gap_char('-')
 Function  : Gets/sets the gap_char attribute of the alignment
 Returns   : An gap_char string, defaults to '-'
 Argument  : An gap_char string (optional), default as '-'

=cut

sub gap_char {
    my ($self, $char) = @_;

	if (defined $char || ! defined $self->{'_gap_char'} ) {
		$char= '-' unless defined $char;
		$self->throw("Single gap character, not [$char]!") if CORE::length($char) > 1;
		$self->{'_gap_char'} = $char;
	}
	return $self->{'_gap_char'};
}

=head2 mask_char

 Title     : mask_char
 Usage     : $aln->mask_char('?')
 Function  : Gets/sets the mask_char attribute of the alignment
 Returns   : An mask_char string,
 Argument  : An mask_char string (optional), default as '?'

=cut

sub mask_char {
	my ($self, $char) = @_;
	
	#may need to check whether $char is the same with _gap_char, _match_char and _missing_char
	if (defined $char || ! defined $self->{'_mask_char'} ) {
		$char= '?' unless defined $char;
		$self->throw("Single mask character, not [$char]!") if CORE::length($char) > 1;
		$self->{'_mask_char'} = $char;
	}
	return $self->{'_mask_char'};
}

=head2 symbol_chars

 Title   : symbol_chars
 Usage   : my @symbolchars = $aln->symbol_chars;
 Function: Returns all the seen symbols (other than gaps)
 Returns : array of characters that are the seen symbols
 Args    : boolean to include the gap/missing/match characters

=cut

sub symbol_chars{
   my ($self,$includeextra) = @_;

   unless ($self->{'_symbols'}) {
       foreach my $seq ($self->next_Seq) {
           map { $self->{'_symbols'}->{$_} = 1; } split(//,$seq->seq);
       }
   }
   my %copy = %{$self->{'_symbols'}};
   if( ! $includeextra ) {
       foreach my $char ( $self->gap_char, $self->match_char,
			  $self->missing_char, $self->mask_char) {
	   delete $copy{$char} if( defined $char );
       }
   }
   return keys %copy;
}



=head1 Alignment descriptors

These read only methods describe the MSA in various ways.


=head2 score

 Title     : score
 Usage     : $str = $ali->score()
 Function  : get/set a score of the alignment
 Returns   : a score for the alignment
 Argument  : an optional score to set

=cut

sub score {
  my $self = shift;
  $self->{score} = shift if @_;
  return $self->{score};
}



=head2 is_flush

 Title     : is_flush
 Usage     : if ( $ali->is_flush() )
 Function  : Tells you whether the alignment
           : is flush, i.e. all of the same length
 Returns   : 1 or 0
 Argument  : None

=cut

sub is_flush {
    my ($self,$report) = @_;
    my $seq;
    my $length = (-1);
    my $temp;

    foreach $seq ( $self->next_Seq() ) {
	if( $length == (-1) ) {
	    $length = CORE::length($seq->seq());
	    next;
	}

	$temp = CORE::length($seq->seq());
	if( $temp != $length ) {
	    $self->warn("expecting $length not $temp from ".
			$seq->display_id) if( $report );
	    $self->debug("expecting $length not $temp from ".
			 $seq->display_id);
	    $self->debug($seq->seq(). "\n");
	    return 0;
	}
    }

    return 1;
}


=head2 length

 Title     : length()
 Usage     : $len = $ali->length()
 Function  : Returns the maximum length of the alignment.
             To be sure the alignment is a block, use is_flush
 Returns   : Integer
 Argument  : None

=cut

sub length_aln {
    my $self = shift;
    $self->deprecated("length_aln - deprecated method. Use length() instead.");
    $self->length(@_);
}

sub length {
    my $self = shift;
    my $seq;
    my $length = -1;
    my $temp;
    
    foreach $seq ( $self->next_Seq() ) {
        $temp = $seq->length();
        if( $temp > $length ) {
            $length = $temp;
        }
    }

    return $length;
}


=head2 maxdisplayname_length

 Title     : maxdisplayname_length
 Usage     : $ali->maxdisplayname_length()
 Function  : Gets the maximum length of the displayname in the
             alignment. Used in writing out various MSA formats.
 Returns   : integer
 Argument  : None

=cut

sub maxname_length {
    my $self = shift;
    $self->deprecated("maxname_length - deprecated method.".
		      " Use maxdisplayname_length() instead.");
    $self->maxdisplayname_length();
}

sub maxnse_length {
    my $self = shift;
    $self->deprecated("maxnse_length - deprecated method.".
		      " Use maxnse_length() instead.");
    $self->maxdisplayname_length();
}

sub maxdisplayname_length {
    my $self = shift;
    my $maxname = (-1);
    my ($seq,$len);

    foreach $seq ( $self->next_Seq() ) {
	$len = CORE::length $self->displayname($seq->get_nse());

	if( $len > $maxname ) {
	    $maxname = $len;
	}
    }

    return $maxname;
}

=head2 max_metaname_length

 Title     : max_metaname_length
 Usage     : $ali->max_metaname_length()
 Function  : Gets the maximum length of the meta name tags in the
             alignment for the sequences and for the alignment.
             Used in writing out various MSA formats.
 Returns   : integer
 Argument  : None

=cut

sub max_metaname_length {
    my $self = shift;
    my $maxname = (-1);
    my ($seq,$len);
    
    # check seq meta first
    for $seq ( $self->next_Seq() ) {
        next if !$seq->isa('Bio::Seq::MetaI' || !$seq->meta_names);
        for my $mtag ($seq->meta_names) {
            $len = CORE::length $mtag;
            if( $len > $maxname ) {
                $maxname = $len;
            }
        }
    }
    
    # alignment meta
    for my $meta ($self->consensus_meta) {
        next unless $meta;
        for my $name ($meta->meta_names) {
            $len = CORE::length $name;
            if( $len > $maxname ) {
                $maxname = $len;
            }
        }
    }

    return $maxname;
}

=head2 num_residues

 Title     : num_residues
 Usage     : $no = $ali->num_residues
 Function  : number of residues in total in the alignment
 Returns   : integer
 Argument  : None
 Note      : replaces no_residues() 

=cut

sub num_residues {
    my $self = shift;
    my $count = 0;

    foreach my $seq ($self->next_Seq) {
	my $str = $seq->seq();

	$count += ($str =~ s/[A-Za-z]//g);
    }

    return $count;
}

=head2 num_sequences

 Title     : num_sequences
 Usage     : $depth = $ali->num_sequences
 Function  : number of sequence in the sequence alignment
 Returns   : integer
 Argument  : none
 Note      : replaces no_sequences()

=cut

sub num_sequences {
    my $self = shift;
    return scalar($self->next_Seq);
}


=head2 average_percentage_identity

 Title   : average_percentage_identity
 Usage   : $id = $align->average_percentage_identity
 Function: The function uses a fast method to calculate the average
           percentage identity of the alignment
 Returns : The average percentage identity of the alignment
 Args    : None
 Notes   : This method implemented by Kevin Howe calculates a figure that
           is designed to be similar to the average pairwise identity of 
           the alignment (identical in the absence of gaps), without having
           to explicitly calculate pairwise identities proposed by Richard 
           Durbin. Validated by Ewan Birney ad Alex Bateman.

=cut

sub average_percentage_identity{
   my ($self,@args) = @_;

	my %alphabet = map {$_ => undef} qw (A C G T U B D E F H I J K L M N O P Q R S V W X Y Z);
   
   my ($len, $total, $subtotal, $divisor, $subdivisor, @seqs, @countHashes);

   if (! $self->is_flush()) {
       $self->throw("All sequences in the alignment must be the same length");
   }

   @seqs = $self->next_Seq();
   $len = $self->length();

   # load the each hash with correct keys for existence checks

   foreach my $seq (@seqs)  {
       my @seqChars = split //, $seq->seq();
       for( my $column=0; $column < @seqChars; $column++ ) {
	   my $char = uc($seqChars[$column]);
	   if (exists $alphabet{$char}) {
	       $countHashes[$column]->{$char}++;
	   }
       }
   }

   $total = 0;
   $divisor = 0;
   for(my $column =0; $column < $len; $column++) {
       my %hash = %{$countHashes[$column]};
       $subdivisor = 0;
       foreach my $res (keys %hash) {
	   $total += $hash{$res}*($hash{$res} - 1);
	   $subdivisor += $hash{$res};
       }
       $divisor += $subdivisor * ($subdivisor - 1);
   }
   return $divisor > 0 ? ($total / $divisor )*100.0 : 0;
}

=head2 percentage_identity

 Title   : percentage_identity
 Usage   : $id = $align->percentage_identity
 Function: The function calculates the average percentage identity
           (aliased to average_percentage_identity)
 Returns : The average percentage identity
 Args    : None

=cut

sub percentage_identity {
    my $self = shift;
    return $self->average_percentage_identity();
}

=head2 overall_percentage_identity

 Title   : overall_percentage_identity
 Usage   : $id = $align->overall_percentage_identity
           $id = $align->overall_percentage_identity('short')
 Function: The function calculates the percentage identity of
           the conserved columns
 Returns : The percentage identity of the conserved columns
 Args    : length value to use, optional defaults to alignment length
                 possible values: 'align', 'short', 'long'

The argument values 'short' and 'long' refer to shortest and longest
sequence in the alignment. Method modification code by Hongyu Zhang.

=cut

sub overall_percentage_identity{
   my ($self, $length_measure) = @_;

   my %alphabet = map {$_ => undef} qw (A C G T U B D E F H I J K L M N O P Q R S V W X Y Z);

   my %enum = map {$_ => undef} qw (align short long);

   $self->throw("Unknown argument [$length_measure]") 
       if $length_measure and not exists $enum{$length_measure};
   $length_measure ||= 'align';

   if (! $self->is_flush()) {
       $self->throw("All sequences in the alignment must be the same length");
   }

   # Count the residues seen at each position
   my $len;
   my $total = 0; # number of positions with identical residues
   my @countHashes;
   my @seqs = $self->next_Seq;
   my $nof_seqs = scalar @seqs;
   my $aln_len = $self->length();
   for my $seq (@seqs)  {
       my $seqstr = $seq->seq;

       # Count residues for given sequence
       for my $column (0 .. $aln_len-1) {
           my $char = uc( substr($seqstr, $column, 1) );
           if ( exists $alphabet{$char} ) {

               # This is a valid char
               if ( defined $countHashes[$column]->{$char} ) {
                 $countHashes[$column]->{$char}++;
               } else {
                 $countHashes[$column]->{$char} = 1;
               }

               if ( $countHashes[$column]->{$char} == $nof_seqs ) {
                   # All sequences have this same residue
                   $total++;
               }

           }
       }

       # Sequence length
       if ($length_measure eq 'short' || $length_measure eq 'long') {
           my $seq_len = $seqstr =~ tr/[A-Za-z]//;
           if ($length_measure eq 'short') {
               if ( (not defined $len) || ($seq_len < $len) ) {
                   $len = $seq_len;
               }
           } elsif ($length_measure eq 'long') {
               if ( (not defined $len) || ($seq_len > $len) ) {
                   $len = $seq_len;
               }
           }
       }

   }

   if ($length_measure eq 'align') {
       $len = $aln_len;
   }

   return ($total / $len ) * 100.0;
}

=head2 pairwise_percentage_identity

 Title     : pairwise_percentage_identity()
 Usage     : @pairwiseiden=$ali->pairwise_percentage_identity(3)
 Function  : Returns pairwise percentage identity of each sequence to the 
             reference sequence(first sequence as default), or selected sequence
 				 See set_new_reference for information on reference sequence
 Returns   : A list of percentage identity to the reference sequence
 Argument  : A number for the position of the reference sequence or the 
             sequence name of the reference sequence

=cut

sub pairwise_percentage_identity {
	my ($self,$seqid)=@_;

	my %alphabet = map {$_ => undef} qw (A C G T U B D E F H I J K L M N O P Q R S V W X Y Z);

   my ($refseq, @seqs,@ids,@pairwise_iden);

   if (! $self->is_flush()) {
       $self->throw("All sequences in the alignment must be the same length");
   }

	foreach my $seq ( $self->next_Seq() ) {
		push @seqs, $seq;
		push @ids, $seq->display_id;
	}

   # load the each hash with correct keys for existence checks
   if(defined($seqid)) {
   	my $is_num=0;
	   if ($seqid =~ /^\d+$/) { # argument is seq position
			$is_num=1;
			$self->throw("The new reference sequence number has to be a positive integer >=1 and <= num_sequences ") if ($seqid < 1 || $seqid > $self->num_sequences);
		}
		else { # argument is a seq name
			$self->throw("The new reference sequence not in alignment ") unless &_in_aln($seqid, \@ids);
	   }
		#return the requested seq as the reference seq in the pairwise comparison		
		for (my $i=0; $i<=$#seqs; $i++) {
			my $pos=$i+1;
			if ( ($is_num && $pos == $seqid) || ($seqid eq $seqs[$i]->display_id) ) {
				$refseq= $seqs[$i];
				last;
			}
		}
	}
	else {
		$refseq=$seqs[0];
	}
	
	#calculate the length of the reference sequence
   my @refseqchars=split //, uc($refseq->seq());
   my $reflength=$refseq->seq()=~tr/A-Za-z//;;

   if($reflength==0) {
   	$self->throw("The reference sequence should be non-zero length ");
   }
   
   #calculate the pairwise identity
   foreach my $seq (@seqs)  {
   	my $idcount=0;
		my @seqChars = split //, $seq->seq();
		for( my $column=0; $column < @seqChars; $column++ ) {
			my $char = uc($seqChars[$column]);
			if(exists $alphabet{$char} && exists $alphabet{$refseqchars[$column]} && $char eq $refseqchars[$column]) {
				$idcount++;
			}
		}
		push @pairwise_iden,$idcount/$reflength;
	}
	
	return @pairwise_iden;
}


=head2 column_from_residue_number

 Title   : column_from_residue_number
 Usage   : $col = $ali->column_from_residue_number( $seqname, $resnumber)
 Function: This function gives the position in the alignment
           (i.e. column number) of the given residue number in the
           sequence with the given name. For example, for the
           alignment

    	     Seq1/91-97 AC..DEF.GH.
   	     Seq2/24-30 ACGG.RTY...
  	        Seq3/43-51 AC.DDEF.GHI

           column_from_residue_number( "Seq1", 94 ) returns 6.
           column_from_residue_number( "Seq2", 25 ) returns 2.
           column_from_residue_number( "Seq3", 50 ) returns 10.

           An exception is thrown if the residue number would lie
           outside the length of the aligment
           (e.g. column_from_residue_number( "Seq2", 22 )

      	  Note: If the the parent sequence is represented by more than
	        one alignment sequence and the residue number is present in
	        them, this method finds only the first one.

 Returns : A column number for the position in the alignment of the
           given residue in the given sequence (1 = first column)
 Args    : A sequence id/name (not a name/start-end)
           A residue number in the whole sequence (not just that
           segment of it in the alignment)

=cut

sub column_from_residue_number {
    my ($self, $name, $resnumber) = @_;

    $self->throw("No sequence with name [$name]") unless $self->{'_start_end_lists'}->{$name};
    $self->throw("Second argument residue number missing") unless $resnumber;

    foreach my $seq ($self->next_Seq_with_id($name)) {
	my $col;
	eval {
	    $col = $seq->column_from_residue_number($resnumber);
	};
	next if $@;
	return $col;
    }

    $self->throw("Could not find a sequence segment in $name ".
		 "containing residue number $resnumber");

}


=head1 Sequence names

Methods to manipulate the display name. The default name based on the
sequence id and subsequence positions can be overridden in various
ways.

=head2 displayname

 Title     : displayname
 Usage     : $myalign->displayname("Ig", "IgA")
 Function  : Gets/sets the display name of a sequence in the alignment
 Returns   : A display name string
 Argument  : name of the sequence
             displayname of the sequence (optional)

=cut

sub get_displayname {
    my $self = shift;
    $self->deprecated("get_displayname - deprecated method. Use displayname() instead.");
    $self->displayname(@_);
}

sub set_displayname {
    my $self = shift;
    $self->deprecated("set_displayname - deprecated method. Use displayname() instead.");
    $self->displayname(@_);
}

sub displayname {
    my ($self, $name, $disname) = @_;

    $self->throw("No sequence with name [$name]")
        unless defined $self->{'_seq'}->{$name};

    if(  $disname and  $name) {
	$self->{'_dis_name'}->{$name} = $disname;
	return $disname;
    }
    elsif( defined $self->{'_dis_name'}->{$name} ) {
	return  $self->{'_dis_name'}->{$name};
    } else {
	return $name;
    }
}

=head2 set_displayname_count

 Title     : set_displayname_count
 Usage     : $ali->set_displayname_count
 Function  : Sets the names to be name_# where # is the number of
             times this name has been used.
 Returns   : 1, on success
 Argument  : None

=cut

sub set_displayname_count {
    my $self= shift;
    my (@arr,$name,$seq,$count,$temp,$nse);

    foreach $seq ( $self->each_alphabetically() ) {
	$nse = $seq->get_nse();

	#name will be set when this is the second
	#time (or greater) is has been seen

	if( defined $name and $name eq ($seq->id()) ) {
	    $temp = sprintf("%s_%s",$name,$count);
	    $self->displayname($nse,$temp);
	    $count++;
	} else {
	    $count = 1;
	    $name = $seq->id();
	    $temp = sprintf("%s_%s",$name,$count);
	    $self->displayname($nse,$temp);
	    $count++;
	}
    }
    return 1;
}

=head2 set_displayname_flat

 Title     : set_displayname_flat
 Usage     : $ali->set_displayname_flat()
 Function  : Makes all the sequences be displayed as just their name,
             not name/start-end
 Returns   : 1
 Argument  : None

=cut

sub set_displayname_flat {
    my $self = shift;
    my ($nse,$seq);

    foreach $seq ( $self->next_Seq() ) {
	$nse = $seq->get_nse();
	$self->displayname($nse,$seq->id());
    }
    return 1;
}

=head2 set_displayname_normal

 Title     : set_displayname_normal
 Usage     : $ali->set_displayname_normal()
 Function  : Makes all the sequences be displayed as name/start-end
 Returns   : 1, on success
 Argument  : None

=cut

sub set_displayname_normal {
    my $self = shift;
    my ($nse,$seq);

    foreach $seq ( $self->next_Seq() ) {
	$nse = $seq->get_nse();
	$self->displayname($nse,$nse);
    }
    return 1;
}

=head2 set_displayname_safe

 Title     : set_displayname_safe
 Usage     : ($new_aln, $ref_name)=$ali->set_displayname_safe(4)
 Function  : Assign machine-generated serial names to sequences in input order.
             Designed to protect names during PHYLIP runs. Assign 10-char string
             in the form of "S000000001" to "S999999999". Restore the original
             names using "restore_displayname".
 Returns   : 1. a new $aln with system names;
             2. a hash ref for restoring names
 Argument  : Number for id length (default 10)

=cut

sub set_displayname_safe {
    my $self = shift;
    my $idlength = shift || 10;
    my ($seq, %phylip_name);
    my $ct=0;
    my $new=Bio::SimpleAlign->new();
    foreach $seq ( $self->next_Seq() ) {
	$ct++;
	my $pname="S". sprintf "%0" . ($idlength-1) . "s", $ct;
	$phylip_name{$pname}=$seq->id();
	my $new_seq= Bio::LocatableSeq->new(-id       => $pname,
					    -seq      => $seq->seq(),
					    -alphabet => $seq->alphabet,
					    -start    => $seq->{_start},
					    -end      => $seq->{_end}
					    );
	$new->add_Seq($new_seq);
    }

    $self->debug("$ct seq names changed. Restore names by using restore_displayname.");
    return ($new, \%phylip_name);
}

=head2 restore_displayname

 Title     : restore_displayname
 Usage     : $aln_name_restored=$ali->restore_displayname($hash_ref)
 Function  : Restore original sequence names (after running
             $ali->set_displayname_safe)
 Returns   : a new $aln with names restored.
 Argument  : a hash reference of names from "set_displayname_safe".

=cut

sub restore_displayname {
    my $self = shift;
    my $ref=shift;
    my %name=%$ref;
    my $new=Bio::SimpleAlign->new();
    foreach my $seq ( $self->next_Seq() ) {
      $self->throw("No sequence with name") unless defined $name{$seq->id()};
      my $new_seq= Bio::LocatableSeq->new(-id       => $name{$seq->id()},
					  -seq      => $seq->seq(),
					  -alphabet => $seq->alphabet,
					  -start    => $seq->{_start},
					  -end      => $seq->{_end}
					  );
      $new->add_Seq($new_seq);
    }
    return $new;
}


=head2 methods implementing Bio::FeatureHolderI

FeatureHolderI implementation to support labeled character sets like one
would get from NEXUS represented data.

=head2 get_SeqFeatures

 Usage   : @features = $aln->get_SeqFeatures
 Function: Get the feature objects held by this feature holder.
 Example :
 Returns : an array of Bio::SeqFeatureI implementing objects
 Args    : optional filter coderef, taking a Bio::SeqFeatureI 
         : as argument, returning TRUE if wanted, FALSE if 
         : unwanted

=cut

sub get_SeqFeatures {
    my $self = shift;
    my $filter_cb = shift;
    $self->throw("Arg (filter callback) must be a coderef") unless 
	!defined($filter_cb) or ref($filter_cb) eq 'CODE';
    if( !defined $self->{'_as_feat'} ) {
	$self->{'_as_feat'} = [];
    }
    if ($filter_cb) {
	return grep { $filter_cb->($_) } @{$self->{'_as_feat'}};
    }
    return @{$self->{'_as_feat'}};
}

=head2 add_SeqFeature

 Usage   : $aln->add_SeqFeature($subfeat);
 Function: adds a SeqFeature into the SeqFeature array.
 Example :
 Returns : true on success
 Args    : a Bio::SeqFeatureI object
 Note    : This implementation is not compliant
           with Bio::FeatureHolderI

=cut

sub add_SeqFeature {
   my ($self,@feat) = @_;

   $self->{'_as_feat'} = [] unless $self->{'_as_feat'};

   foreach my $feat ( @feat ) {
       if( !$feat->isa("Bio::SeqFeatureI") ) {
           $self->throw("$feat is not a SeqFeatureI and that's what we expect...");
       }

       push(@{$self->{'_as_feat'}},$feat);
   }
   return 1;
}


=head2 remove_SeqFeatures

 Usage   : $obj->remove_SeqFeatures
 Function: Removes all SeqFeatures.  If you want to remove only a subset,
           remove that subset from the returned array, and add back the rest.
 Returns : The array of Bio::SeqFeatureI features that was
           deleted from this alignment.
 Args    : none

=cut

sub remove_SeqFeatures {
    my $self = shift;

    return () unless $self->{'_as_feat'};
    my @feats = @{$self->{'_as_feat'}};
    $self->{'_as_feat'} = [];
    return @feats;
}

=head2 feature_count

 Title   : feature_count
 Usage   : $obj->feature_count()
 Function: Return the number of SeqFeatures attached to the alignment
 Returns : integer representing the number of SeqFeatures
 Args    : None

=cut

sub feature_count {
    my ($self) = @_;

    if (defined($self->{'_as_feat'})) {
        return ($#{$self->{'_as_feat'}} + 1);
    } else {
        return 0;
    }
}

=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : 
 Function: Get all SeqFeatures.
 Example :
 Returns : an array of Bio::SeqFeatureI implementing objects
 Args    : none
 Note    : Falls through to Bio::FeatureHolderI implementation.

=cut

=head2 methods for Bio::AnnotatableI

AnnotatableI implementation to support sequence alignments which
contain annotation (NEXUS, Stockholm).

=head2 annotation

 Title   : annotation
 Usage   : $ann = $aln->annotation or 
           $aln->annotation($ann)
 Function: Gets or sets the annotation
 Returns : Bio::AnnotationCollectionI object
 Args    : None or Bio::AnnotationCollectionI object

See L<Bio::AnnotationCollectionI> and L<Bio::Annotation::Collection>
for more information

=cut

sub annotation {
    my ($obj,$value) = @_;
    if( defined $value ) {
        $obj->throw("object of class ".ref($value)." does not implement ".
                "Bio::AnnotationCollectionI. Too bad.")
            unless $value->isa("Bio::AnnotationCollectionI");
        $obj->{'_annotation'} = $value;
    } elsif( ! defined $obj->{'_annotation'}) {
        $obj->{'_annotation'} = Bio::Annotation::Collection->new();
    }
    return $obj->{'_annotation'};
}


sub no_residues {
	my $self = shift;
	$self->deprecated(-warn_version => 1.0069,
					  -throw_version => 1.0075,
                      -message => 'Use of method no_residues() is deprecated, use num_residues() instead');
    $self->num_residues(@_);
}

sub no_sequences {
	my $self = shift;
	$self->deprecated(-warn_version => 1.0069,
					  -throw_version => 1.0075,
                      -message => 'Use of method no_sequences() is deprecated, use num_sequences() instead');
    $self->num_sequences(@_);
}



1;
