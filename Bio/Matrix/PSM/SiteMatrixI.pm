
=head1 NAME

Bio::Matrix::PSM::SiteMatrixI - SiteMatrixI implementation, holds a
position scoring matrix (or position weight matrix) and log-odds

=head1 SYNOPSIS

  # You cannot use this module directly; see Bio::Matrix::PSM::SiteMatrix
  # for an example implementation

=head1 DESCRIPTION

SiteMatrix is designed to provide some basic methods when working with position
scoring (weight) matrices, such as transcription factor binding sites for
example. A DNA PSM consists of four vectors with frequencies {A,C,G,T}. This is
the minimum information you should provide to construct a PSM object. The
vectors can be provided as strings with frequenciesx10 rounded to an int, going
from {0..a} and 'a' represents the maximum (10). This is like MEME's compressed
representation of a matrix and it is quite useful when working with relational
DB. If arrays are provided as an input (references to arrays actually) they can
be any number, real or integer (frequency or count).

When creating the object you can ask the constructor to make a simple pseudo
count correction by adding a number (typically 1) to all positions (with the
-correction option). After adding the number the frequencies will be
calculated. Only use correction when you supply counts, not frequencies.

Throws an exception if: You mix as an input array and string (for example A
matrix is given as array, C - as string). The position vector is (0,0,0,0). One
of the probability vectors is shorter than the rest.

Summary of the methods I use most frequently (details bellow):

  iupac - return IUPAC compliant consensus as a string
  score - Returns the score as a real number
  IC - information content. Returns a real number
  id - identifier. Returns a string
  accession - accession number. Returns a string
  next_pos - return the sequence probably for each letter, IUPAC
      symbol, IUPAC probability and simple sequence
  consenus letter for this position. Rewind at the end. Returns a hash.
  pos - current position get/set. Returns an integer.
  regexp - construct a regular expression based on IUPAC consensus.
      For example AGWV will be [Aa][Gg][AaTt][AaCcGg]
  width - site width
  get_string - gets the probability vector for a single base as a string.
  get_array - gets the probability vector for a single base as an array.
  get_logs_array - gets the log-odds vector for a single base as an array.

New methods, which might be of interest to anyone who wants to store PSM in a relational
database without creating an entry for each position is the ability to compress the
PSM vector into a string with losing usually less than 1% of the data.
this can be done with:

  my $str=$matrix->get_compressed_freq('A');

or

  my $str=$matrix->get_compressed_logs('A');

Loading from a database should be done with new, but is not yest implemented.
However you can still uncompress such string with:

  my @arr=Bio::Matrix::PSM::_uncompress_string ($str,1,1); for PSM

or

  my @arr=Bio::Matrix::PSM::_uncompress_string ($str,1000,2); for log odds

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Stefan Kirov

Email skirov@utk.edu

=head1 APPENDIX

=cut


# Let the code begin...

package Bio::Matrix::PSM::SiteMatrixI;

# use strict;
use base qw(Bio::Root::RootI);

=head2 calc_weight

 Title   : calc_weight
 Usage   : $self->calc_weight({A=>0.2562,C=>0.2438,G=>0.2432,T=>0.2568});
 Function: Recalculates the PSM (or weights) based on the PFM (the frequency matrix)
           and user supplied background model.
 Throws  : if no model is supplied
 Example :
 Returns :
 Args    : reference to a hash with background frequencies for A,C,G and T

=cut

sub calc_weight {
  my $self = shift;
  $self->throw_not_implemented();
}


=head2 next_pos

 Title   : next_pos
 Usage   : my %base=$site->next_pos;
 Function: 

           Retrieves the next position features: frequencies and weights for
           A,C,G,T, the main letter (as in consensus) and the
           probabilty for this letter to occur at this position and
           the current position

 Throws  :
 Example :
 Returns : hash (pA,pC,pG,pT,lA,lC,lG,lT,base,prob,rel)
 Args    : none


=cut

sub next_pos {
  my $self = shift;
  $self->throw_not_implemented();
}

=head2 curpos

 Title   : curpos
 Usage   : my $pos=$site->curpos;
 Function: Gets/sets the current position. Converts to 0 if argument is minus and
            to width if greater than width
 Throws  :
 Example :
 Returns : integer
 Args    : integer

=cut

sub curpos {
    my $self = shift;
   $self->throw_not_implemented();
}

=head2 e_val

 Title   : e_val
 Usage   : my $score=$site->e_val;
 Function: Gets/sets the e-value
 Throws  :
 Example :
 Returns : real number
 Args    : real number

=cut

sub e_val {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 consensus

 Title   : consensus
 Usage   :
 Function: Returns the consensus
 Returns : string
 Args    : (optional) threshold value 1 to 10, default 5
           '5' means the returned characters had a 50% or higher presence at
           their position

=cut

sub consensus {
  my $self = shift;
  $self->throw_not_implemented();
}

=head2 accession_number

 Title   : accession_number
 Usage   :
 Function: accession number, this will be unique id for the SiteMatrix object as
 			well for any other object, inheriting from SiteMatrix
 Throws  :
 Example :
 Returns : string
 Args    : string

=cut

sub accession_number {
  my $self = shift;
  $self->throw_not_implemented();
}


=head2 width

 Title   : width
 Usage   : my $width=$site->width;
 Function: Returns the length of the site
 Throws  :
 Example :
 Returns : number
 Args    :

=cut

sub width {
  my $self = shift;
  $self->throw_not_implemented();
}

=head2 IUPAC

 Title   : IUPAC
 Usage   : my $iupac_consensus=$site->IUPAC;
 Function: Returns IUPAC compliant consensus
 Throws  :
 Example :
 Returns : string
 Args    :

=cut

sub IUPAC {
  my $self = shift;
  $self->throw_not_implemented();
}

=head2 IC

 Title   : IC
 Usage   : my $ic=$site->IC;
 Function: Information content
 Throws  :
 Example :
 Returns : real number
 Args    : none

=cut

sub IC {
my $self=shift;
$self->throw_not_implemented();
}

=head2 get_string

 Title   : get_string
 Usage   : my $freq_A=$site->get_string('A');
 Function: Returns given probability vector as a string. Useful if you want to
           store things in a rel database, where arrays are not first choice
 Throws  : If the argument is outside {A,C,G,T}
 Example :
 Returns : string
 Args    : character {A,C,G,T}

=cut

sub get_string {
 my $self=shift;
 $self->throw_not_implemented();
}

=head2 id

 Title   : id
 Usage   : my $id=$site->id;
 Function: Gets/sets the site id
 Throws  :
 Example :
 Returns : string
 Args    : string

=cut

sub id {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 regexp

 Title   : regexp
 Usage   : my $regexp=$site->regexp;
 Function: Returns a regular expression which matches the IUPAC convention.
           N will match X, N, - and .
 Throws  :
 Example :
 Returns : string
 Args    :

=cut

sub regexp {
 my $self=shift;
 $self->throw_not_implemented();
}

=head2 regexp_array

 Title   : regexp_array
 Usage   : my @regexp=$site->regexp;
 Function: Returns a regular expression which matches the IUPAC convention.
           N will match X, N, - and .
 Throws  :
 Example :
 Returns : array
 Args    :
 To do   : I have separated regexp and regexp_array, but
           maybe they can be rewritten as one - just check what
           should be returned

=cut

sub regexp_array {
 my $self=shift;
 $self->throw_not_implemented();
}

=head2 get_array

 Title   : get_array
 Usage   : my @freq_A=$site->get_array('A');
 Function: Returns an array with frequencies for a specified base
 Throws  :
 Example :
 Returns : array
 Args    : char

=cut

sub get_array {
  my $self=shift;
  $self->throw_not_implemented();
}


=head2 _to_IUPAC

 Title   : _to_IUPAC
 Usage   :
 Function: Converts a single position to IUPAC compliant symbol and
            returns its probability. For rules see the implementation.
 Throws  :
 Example :
 Returns : char, real number
 Args    : real numbers for A,C,G,T (positional)

=cut

sub _to_IUPAC {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 _to_cons

 Title   : _to_cons
 Usage   :
 Function: Converts a single position to simple consensus character and
            returns its probability. For rules see the implementation,
 Throws  :
 Example :
 Returns : char, real number
 Args    : real numbers for A,C,G,T (positional)

=cut

sub _to_cons {
    my $self = shift;
    $self->throw_not_implemented();
}


=head2 _calculate_consensus

 Title   : _calculate_consensus
 Usage   :
 Function: Internal stuff
 Throws  :
 Example :
 Returns :
 Args    :

=cut

sub _calculate_consensus {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 _compress_array

 Title   : _compress_array
 Usage   :
 Function:  Will compress an array of real signed numbers to a string (ie vector of bytes)
 			-127 to +127 for bi-directional(signed) and 0..255 for unsigned ;
 Throws  :
 Example :  Internal stuff
 Returns :  String
 Args    :  array reference, followed by an max value and
 			direction (optional, default 1-unsigned),1 unsigned, any other is signed.

=cut

sub _compress_array {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 _uncompress_string

 Title   : _uncompress_string
 Usage   :
 Function:  Will uncompress a string (vector of bytes) to create an array of real
            signed numbers (opposite to_compress_array)
 Throws  :
 Example :  Internal stuff
 Returns :  string, followed by an max value and
 			direction (optional, default 1-unsigned), 1 unsigned, any other is signed.
 Args    :  array

=cut

sub _uncompress_string {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 get_compressed_freq

 Title   : get_compressed_freq
 Usage   :
 Function:  A method to provide a compressed frequency vector. It uses one byte to
 			code the frequence for one of the probability vectors for one position.
			Useful for relational database. Improvment of the previous 0..a coding.
 Throws  :
 Example :  my $strA=$self->get_compressed_freq('A');
 Returns :  String
 Args    :  char

=cut

sub get_compressed_freq {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 get_compressed_logs

 Title   : get_compressed_logs
 Usage   :
 Function:  A method to provide a compressed log-odd vector. It uses one byte to
 			code the log value for one of the log-odds vectors for one position.
 Throws  :
 Example :  my $strA=$self->get_compressed_logs('A');
 Returns :  String
 Args    :  char

=cut

sub get_compressed_logs {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 sequence_match_weight

 Title   : sequence_match_weight
 Usage   :
 Function:  This method will calculate the score of a match, based on the PWM
            if such is associated with the matrix object. Returns undef if no
             PWM data is available.
 Throws  :   if the length of the sequence is different from the matrix width
 Example :  my $score=$matrix->sequence_match_weight('ACGGATAG');
 Returns :  Floating point
 Args    :  string

=cut

sub sequence_match_weight {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 get_all_vectors

 Title   : get_all_vectors
 Usage   :
 Function:  returns all possible sequence vectors to satisfy the PFM under
            a given threshold
 Throws  :  If threshold outside of 0..1 (no sense to do that)
 Example :  my @vectors=$self->get_all_vectors(4);
 Returns :  Array of strings
 Args    :  (optional) floating

=cut

sub get_all_vectors {
 my $self = shift;
    $self->throw_not_implemented();
}
1;
