#
# BioPerl module for Bio::Align::AlignI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Align::AlignI - An interface for describing sequence alignments.

=head1 SYNOPSIS

  # get a Bio::Align::AlignI somehow - typically using Bio::AlignIO system
  # some descriptors
  print $aln->length, "\n";
  print $aln->num_residues, "\n";
  print $aln->is_flush, "\n";
  print $aln->num_sequences, "\n";
  print $aln->percentage_identity, "\n";
  print $aln->consensus_string(50), "\n";

  # find the position in the alignment for a sequence location
  $pos = $aln->column_from_residue_number('1433_LYCES', 14); # = 6;

  # extract sequences and check values for the alignment column $pos
  foreach $seq ($aln->each_seq) {
      $res = $seq->subseq($pos, $pos);
      $count{$res}++;
  }
  foreach $res (keys %count) {
      printf "Res: %s  Count: %2d\n", $res, $count{$res};
  }

=head1 DESCRIPTION

This interface describes the basis for alignment objects.

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

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 CONTRIBUTORS

Ewan Birney, birney@ebi.ac.uk
Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Align::AlignI;
use strict;


use base qw(Bio::Root::RootI);

=head1 Modifier methods

These methods modify the MSE by adding, removing or shuffling complete
sequences.

=head2 add_seq

 Title     : add_seq
 Usage     : $myalign->add_seq($newseq);
 Function  : Adds another sequence to the alignment. *Does not* align
             it - just adds it to the hashes.
 Returns   : None
 Argument  : a Bio::LocatableSeq object
             order (optional)

See L<Bio::LocatableSeq> for more information.

=cut

sub add_seq {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 remove_seq

 Title     : remove_seq
 Usage     : $aln->remove_seq($seq);
 Function  : Removes a single sequence from an alignment
 Returns   :
 Argument  : a Bio::LocatableSeq object

=cut

sub remove_seq {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 purge

 Title   : purge
 Usage   : $aln->purge(0.7);
 Function:

           Removes sequences above whatever %id.

           This function will grind on large alignments. Beware!
           (perhaps not ideally implemented)

 Example :
 Returns : An array of the removed sequences
 Argument:


=cut

sub purge {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 sort_alphabetically

 Title     : sort_alphabetically
 Usage     : $ali->sort_alphabetically
 Function  : 

             Changes the order of the alignment to alphabetical on name 
             followed by numerical by number.

 Returns   : an array
 Argument  : 

=cut

sub sort_alphabetically {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head1 Sequence selection methods

Methods returning one or more sequences objects.

=head2 each_seq

 Title     : each_seq
 Usage     : foreach $seq ( $align->each_seq() ) 
 Function  : Gets an array of Seq objects from the alignment
 Returns   : an array
 Argument  : 

=cut

sub each_seq {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 each_alphabetically

 Title     : each_alphabetically
 Usage     : foreach $seq ( $ali->each_alphabetically() )
 Function  :

             Returns an array of sequence object sorted alphabetically 
             by name and then by start point.
             Does not change the order of the alignment

 Returns   : 
 Argument  : 

=cut

sub each_alphabetically {
    my($self) = @_;
    $self->throw_not_implemented();
}

=head2 each_seq_with_id

 Title     : each_seq_with_id
 Usage     : foreach $seq ( $align->each_seq_with_id() ) 
 Function  : 

             Gets an array of Seq objects from the
             alignment, the contents being those sequences
             with the given name (there may be more than one)

 Returns   : an array
 Argument  : a seq name

=cut

sub each_seq_with_id {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 get_seq_by_pos

 Title     : get_seq_by_pos
 Usage     : $seq = $aln->get_seq_by_pos(3) # third sequence from the alignment
 Function  : 

             Gets a sequence based on its position in the alignment.
             Numbering starts from 1.  Sequence positions larger than
             num_sequences() will throw an error.

 Returns   : a Bio::LocatableSeq object
 Argument  : positive integer for the sequence position

=cut

sub get_seq_by_pos {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head1 Create new alignments

The result of these methods are horizontal or vertical subsets of the
current MSE.

=head2 select

 Title     : select
 Usage     : $aln2 = $aln->select(1, 3) # three first sequences
 Function  : 

             Creates a new alignment from a continuous subset of
             sequences.  Numbering starts from 1.  Sequence positions
             larger than num_sequences() will throw an error.

 Returns   : a Bio::SimpleAlign object
 Argument  : positive integer for the first sequence
             positive integer for the last sequence to include (optional)

=cut

sub select {
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 select_noncont

 Title     : select_noncont
 Usage     : $aln2 = $aln->select_noncont(1, 3) # first and 3rd sequences
 Function  : 

             Creates a new alignment from a subset of
             sequences.  Numbering starts from 1.  Sequence positions
             larger than num_sequences() will throw an error.

 Returns   : a Bio::SimpleAlign object
 Args      : array of integers for the sequences

=cut

sub select_noncont {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 slice

 Title     : slice
 Usage     : $aln2 = $aln->slice(20, 30)
 Function  : 

             Creates a slice from the alignment inclusive of start and
             end columns.  Sequences with no residues in the slice are
             excluded from the new alignment and a warning is printed.
             Slice beyond the length of the sequence does not do
             padding.

 Returns   : a Bio::SimpleAlign object
 Argument  : positive integer for start column 
             positive integer for end column 

=cut

sub slice {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head1 Change sequences within the MSE

These methods affect characters in all sequences without changing the
alignment.


=head2 map_chars

 Title     : map_chars
 Usage     : $ali->map_chars('\.','-')
 Function  : 

             Does a s/$arg1/$arg2/ on the sequences. Useful for gap
             characters

             Notice that the from (arg1) is interpreted as a regex,
             so be careful about quoting meta characters (eg
             $ali->map_chars('.','-') wont do what you want)

 Returns   : None
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
 Returns   : 
 Argument  : 

=cut

sub uppercase {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 match_line

 Title    : match_line()
 Usage    : $align->match_line()
 Function : Generates a match line - much like consensus string
            except that a line indicating the '*' for a match.
 Argument : (optional) Match line characters ('*' by default)
            (optional) Strong match char (':' by default)
            (optional) Weak match char ('.' by default)

=cut

sub match_line {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 match

 Title     : match()
 Usage     : $ali->match()
 Function  : 

             Goes through all columns and changes residues that are
             identical to residue in first sequence to match '.'
             character. Sets match_char.

             USE WITH CARE: Most MSE formats do not support match
             characters in sequences, so this is mostly for output
             only. NEXUS format (Bio::AlignIO::nexus) can handle
             it.

 Returns   : 1
 Argument  : a match character, optional, defaults to '.'

=cut

sub match {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 unmatch

 Title     : unmatch()
 Usage     : $ali->unmatch()
 Function  : 

             Undoes the effect of method match. Unsets match_char.

 Returns   : 1
 Argument  : a match character, optional, defaults to '.'

=cut

sub unmatch {
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head1 MSE attibutes

Methods for setting and reading the MSE attributes. 

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

=head2 missing_char

 Title     : missing_char
 Usage     : $myalign->missing_char("?")
 Function  : Gets/sets the missing_char attribute of the alignment
             It is generally recommended to set it to 'n' or 'N' 
             for nucleotides and to 'X' for protein. 
 Returns   : An missing_char string,
 Argument  : An missing_char string (optional)

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
 Argument  : An match_char string (optional)

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
 Argument  : An gap_char string (optional)

=cut

sub gap_char {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 symbol_chars

 Title   : symbol_chars
 Usage   : my @symbolchars = $aln->symbol_chars;
 Function: Returns all the seen symbols (other than gaps)
 Returns : array of characters that are the seen symbols
 Argument: boolean to include the gap/missing/match characters

=cut

sub symbol_chars{
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head1 Alignment descriptors

These read only methods describe the MSE in various ways. 


=head2 consensus_string

 Title     : consensus_string
 Usage     : $str = $ali->consensus_string($threshold_percent)
 Function  : Makes a strict consensus 
 Returns   : consensus string
 Argument  : Optional threshold ranging from 0 to 100.
             The consensus residue has to appear at least threshold %
             of the sequences at a given location, otherwise a '?'
             character will be placed at that location.
             (Default value = 0%)

=cut

sub consensus_string {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 consensus_iupac

 Title     : consensus_iupac
 Usage     : $str = $ali->consensus_iupac()
 Function  : 

             Makes a consensus using IUPAC ambiguity codes from DNA
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

=head2 is_flush

 Title     : is_flush
 Usage     : if( $ali->is_flush() )
           : 
           :
 Function  : Tells you whether the alignment 
           : is flush, ie all of the same length
           : 
           :
 Returns   : 1 or 0
 Argument  : 

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
 Returns   : integer
 Argument  : 

=cut

sub length {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 maxname_length

 Title     : maxname_length
 Usage     : $ali->maxname_length()
 Function  : 

             Gets the maximum length of the displayname in the
             alignment. Used in writing out various MSE formats.

 Returns   : integer
 Argument  : 

=cut

sub maxname_length {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 num_residues

 Title     : num_residues
 Usage     : $no = $ali->num_residues
 Function  : number of residues in total in the alignment
 Returns   : integer
 Argument  :
 Note      : replaces no_residues

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
 Argument  : None
 Note      : replaces no_sequences

=cut

sub num_sequences {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 percentage_identity

 Title   : percentage_identity
 Usage   : $id = $align->percentage_identity
 Function: The function calculates the percentage identity of the alignment
 Returns : The percentage identity of the alignment (as defined by the 
	   implementation)
 Argument: None

=cut

sub percentage_identity{
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 overall_percentage_identity

 Title   : overall_percentage_identity
 Usage   : $id = $align->overall_percentage_identity
 Function: The function calculates the percentage identity of 
           the conserved columns
 Returns : The percentage identity of the conserved columns
 Args    : None

=cut

sub overall_percentage_identity{
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

=cut

sub average_percentage_identity{
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head1 Alignment positions

Methods to map a sequence position into an alignment column and back.
column_from_residue_number() does the former. The latter is really a
property of the sequence object and can done using
L<Bio::LocatableSeq::location_from_column>:

    # select somehow a sequence from the alignment, e.g.
    my $seq = $aln->get_seq_by_pos(1);
    #$loc is undef or Bio::LocationI object
    my $loc = $seq->location_from_column(5);


=head2 column_from_residue_number

 Title   : column_from_residue_number
 Usage   : $col = $ali->column_from_residue_number( $seqname, $resnumber)
 Function:

           This function gives the position in the alignment
           (i.e. column number) of the given residue number in the
           sequence with the given name. For example, for the
           alignment

  	     Seq1/91-97 AC..DEF.GH
  	     Seq2/24-30 ACGG.RTY..
  	     Seq3/43-51 AC.DDEFGHI

           column_from_residue_number( "Seq1", 94 ) returns 6.
           column_from_residue_number( "Seq2", 25 ) returns 2.
           column_from_residue_number( "Seq3", 50 ) returns 9.

           An exception is thrown if the residue number would lie
           outside the length of the alignment
           (e.g. column_from_residue_number( "Seq2", 22 )

	  Note: If the parent sequence is represented by more than one
	  alignment sequence and the residue number is present in
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
           :
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
 Function  : 

             Sets the names to be name_# where # is the number of
             times this name has been used.

 Returns   : None 
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
 Returns   : None
 Argument  : None

=cut

sub set_displayname_normal {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head1 Deprecated methods

=head2 no_residues

 Title     : no_residues
 Usage     : $no = $ali->no_residues
 Function  : number of residues in total in the alignment
 Returns   : integer
 Argument  : 
 Note      : deprecated in favor of num_residues()

=cut

sub no_residues {
    # immediate deprecation
    shift->deprecated();
}

=head2 no_sequences

 Title     : no_sequences
 Usage     : $depth = $ali->no_sequences
 Function  : number of sequence in the sequence alignment
 Returns   : integer
 Argument  : None
 Note      : deprecated in favor of num_sequences() 

=cut

sub no_sequences {
    # immediate deprecation
    shift->deprecated();
}

1;
