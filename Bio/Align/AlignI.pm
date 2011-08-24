#
# BioPerl module for Bio::Align::AlignI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
#   Copyright Jason Stajich
#	August 2010 refactoring - Jun Yin
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Align::AlignI - An interface for describing sequence alignments.

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

  https://redmine.open-bio.org/projects/bioperl/

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
Jun Yin, jun.yin at ucd.ie

=head1 SEE ALSO

L<Bio::LocatableSeq>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::Align::AlignI;
use strict;


use base qw(Bio::Root::RootI);


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

sub add_Seq {
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 remove_LocatableSeq

 Title     : remove_LocatableSeq
 Usage     : $aln->remove_LocatableSeq($seq);
 Function  : Removes a single sequence from an alignment
 Returns   : 1
 Argument  : a Bio::LocatableSeq object

=cut

sub remove_LocatableSeq {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 remove_Seqs

 Title     : remove_Seqs
 Usage     : $aln->remove_Seqs([1,3,5..7]);
             $aln->remove_Seqs(-selection=>[1],-toggle=>1);
 Function  : Removes specified sequences from the alignment
 Returns   : 1 
 Argument  : An reference list of positive integers for the selected 
             sequences. An optional parameter can be defined to toggle 
             the coordinate selection.
             See also select_Seqs

=cut

sub remove_Seqs {
    my ($self) = @_;
    $self->throw_not_implemented();
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

sub remove_redundant_Seqs {
    my ($self) = @_;
    $self->throw_not_implemented();
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


sub uniq_Seq {
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
     my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 sort_by_pairwise_identity

 Title     : sort_by_pairwise_identity()
 Usage     : $ali->sort_by_pairwise_identity(2)
 Function  : Changes the order of the alignment by the pairwise percentage 
             identity of the reference sequence
 Returns   : 1
 Argument  : Optional, the position or id of reference sequences to be compared 
             with. Default is the first sequence
             See also set_new_reference

=cut

sub sort_by_pairwise_identity {
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
     my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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


sub next_Seq {
    my ($self) = @_;
    $self->throw_not_implemented();
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


sub next_alphabetically {
    my ($self) = @_;
    $self->throw_not_implemented();
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

sub next_Seq_with_id {
    my ($self) = @_;
    $self->throw_not_implemented();
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


sub get_Seq_by_pos {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 get_Seq_by_id

 Title     : get_Seq_by_id
 Usage     : $seq = $aln->get_Seq_by_id($name) # seq named $name
 Function  : Gets a sequence based on its name.
             Sequences that do not exist will warn and return undef
 Returns   : a Bio::LocatableSeq object
 Args      : string for sequence name

=cut

sub get_Seq_by_id {
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 select_Seqs

 Title     : select_Seqs
 Usage     : $aln2 = $aln->select_Seqs([1,5..10,15]) # three first sequences
             $aln2 = $aln->select_Seqs(-selection=>[2..4,11..14],-toggle=>1) # toggle selection
 Function  : Creates a new alignment from a subset of
             sequences.  Numbering starts from 1.  Sequence positions
             larger than num_sequences() will thow an error.
 Returns   : a Bio::SimpleAlign object
 Args      : An reference list of positive integers for the selected 
             sequences. An optional parameter can be defined to toggle the 
             coordinate selection.
             

=cut

sub select_Seqs {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 select_columns

 Title     : select_columns
 Usage     : $newaln = $aln->select_columns([20..30,45..48]);
             $newaln = $aln->select_columns(['mismatch']);
             $newaln = $aln->select_columns(-selection=>['mismatch'],-toggle=>1);
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

sub select_columns {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 remove_columns

 Title     : remove_columns
 Usage     : my $newaln=$aln->remove_columns(['mismatch','weak']) or
             my $newaln=$aln->remove_columns([3,6..8]) or
             my $newaln=$aln->remove_columns(-selection=>[3,6..8],-toggle=>1,-keepgaponly=>0)
 Function  : Modify the aligment with columns removed corresponding to
             the specified type or by specifying the columns by number.
             remove_columns is a variance of select_columns(-toogle=>1)
 Returns   : 1
 Args      : Array ref of types ('match'|'weak'|'strong'|'mismatch') 
             or array ref where the referenced array
             contains a pair of integers that specify a range.
             use remove_gaps to remove columns containing gaps
             The first column is 1
             See also select_columns

=cut

sub remove_columns {
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 mask_columns

 Title     : mask_columns
 Usage     : $aln2 = $aln->mask_columns([3..8,11..22])
             $aln2 = $aln->mask_columns(-selection=>[2,5,7..10],-toggle=>1)
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 uppercase

 Title     : uppercase()
 Usage     : $ali->uppercase()
 Function  : Sets all the sequences to uppercase
 Returns   : 1
 Argument  : None

=cut

sub uppercase {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 lowercase

 Title     : lowercase()
 Usage     : $ali->lowercase()
 Function  : Sets all the sequences to lowercase
 Returns   : 1
 Argument  : None

=cut

sub lowercase {
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 togglecase

 Title     : togglecase()
 Usage     : $ali->togglecase()
 Function  : Sets all the sequences to opposite case
 Returns   : 1
 Argument  : None

=cut

sub togglecase {
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 consensus_conservation

 Title     : consensus_conservation
 Usage     : @conservation = $ali->consensus_conservation();
 Function  : Conservation (as a percent) of each position of alignment
 Returns   : Array of percentages [0-100]. Gap columns are 0% conserved.
 Argument  : 
 
=cut

sub consensus_conservation {
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
     my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 accession

 Title     : accession
 Usage     : $myalign->accession("PF00244")
 Function  : Gets/sets the accession field of the alignment
 Returns   : An acc string
 Argument  : An acc string (optional)

=cut

sub accession {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 description

 Title     : description
 Usage     : $myalign->description("14-3-3 proteins")
 Function  : Gets/sets the description field of the alignment
 Returns   : An description string
 Argument  : An description string (optional)

=cut

sub description {
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 missing_char

 Title     : missing_char
 Usage     : $myalign->missing_char("&")
 Function  : Gets/sets the missing_char attribute of the alignment
             It is generally recommended to set it to 'n' or 'N'
             for nucleotides and to 'X' for protein.
 Returns   : An missing_char string,
 Argument  : An missing_char string (optional), default as '&'

=cut

sub missing_char {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 match_char

 Title     : match_char
 Usage     : $myalign->match_char('.')
 Function  : Gets/sets the match_char attribute of the alignment
 Returns   : An match_char string,
 Argument  : An match_char string (optional), default as '.'

=cut

sub match_char {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 gap_char

 Title     : gap_char
 Usage     : $myalign->gap_char('-')
 Function  : Gets/sets the gap_char attribute of the alignment
 Returns   : An gap_char string, defaults to '-'
 Argument  : An gap_char string (optional), default as '-'

=cut

sub gap_char {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 mask_char

 Title     : mask_char
 Usage     : $aln->mask_char('?')
 Function  : Gets/sets the mask_char attribute of the alignment
 Returns   : An mask_char string,
 Argument  : An mask_char string (optional), default as '?'

=cut

sub mask_char {
    my ($self) = @_;
    $self->throw_not_implemented();
}		
		
		



=head2 symbol_chars
 Title   : symbol_chars
 Usage   : my @symbolchars = $aln->symbol_chars;
 Function: Returns all the seen symbols (other than gaps)
 Returns : array of characters that are the seen symbols
 Args    : boolean to include the gap/missing/match characters

=cut

sub symbol_chars{
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
}

sub length {
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 maxdisplayname_length

 Title     : maxdisplayname_length
 Usage     : $ali->maxdisplayname_length()
 Function  : Gets the maximum length of the displayname in the
             alignment. Used in writing out various MSA formats.
 Returns   : integer
 Argument  : None

=cut

sub maxdisplayname_length {
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();

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


sub displayname {
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
     my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 set_displayname_normal

 Title     : set_displayname_normal
 Usage     : $ali->set_displayname_normal()
 Function  : Makes all the sequences be displayed as name/start-end
 Returns   : 1, on success
 Argument  : None

=cut

sub set_displayname_normal {
     my ($self) = @_;
    $self->throw_not_implemented();
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
     my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    $self->throw_not_implemented();
}


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
    my ($self) = @_;
    $self->throw_not_implemented();
}

1;

