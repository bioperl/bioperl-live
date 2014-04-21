#---------------------------------------------------------

=head1 NAME

Bio::Matrix::PSM::ProtMatrix - SiteMatrixI implementation, holds a
position scoring matrix (or position weight matrix) with log-odds scoring
information.

=head1 SYNOPSIS

   use Bio::Matrix::PSM::ProtMatrix;
   # Create from memory by supplying probability matrix hash both as strings or
   # arrays where the frequencies   Hash entries of the form lN refer to an array
   # of position-specific log-odds scores for amino acid N. Hash entries of the
   # form pN represent the position-specific probability of finding amino acid N.

   my %param = (
             'id' => 'A. thaliana protein atp1',
             '-e_val' => $score,
             'lS' => [ '-2', '3', '-3', '2', '-3', '1', '1', '3' ],
             'lF' => [ '-1', '-4', '0', '-5', '0', '-5', '-4', '-4' ],
             'lT' => [ '-1', '1', '0', '1', '-2', '-1', '0', '1' ],
             'lN' => [ '-3', '-1', '-2', '3', '-5', '5', '-2', '0' ],
             'lK' => [ '-2', '0', '-3', '2', '-3', '2', '-3', '-1' ],
             'lY' => [ '-2', '-3', '-3', '-4', '-3', '-4', '-4', '-4' ],
             'lE' => [ '-3', '4', '-3', '2', '-4', '-2', '-3', '2' ],
             'lV' => [ '0', '-2', '1', '-4', '1', '-4', '-1', '-3' ],
             'lQ' => [ '-1', '0', '-2', '3', '-4', '1', '-3', '0' ],
             'lM' => [ '8', '-3', '8', '-3', '1', '-3', '-3', '-3' ],
             'lC' => [ '-2', '-3', '-3', '-4', '-3', '-4', '-3', '-3' ],
             'lL' => [ '1', '-3', '1', '-4', '3', '-4', '-2', '-4' ],
             'lA' => [ '-2', '1', '-2', '0', '-2', '-2', '2', '2' ],
             'lW' => [ '-2', '-4', '-3', '-5', '-4', '-5', '-5', '-5' ],
             'lP' => [ '-3', '-2', '-4', '-3', '-1', '-3', '6', '-3' ],
             'lH' => [ '-2', '-2', '-3', '-2', '-5', '-2', '-2', '-3' ],
             'lD' => [ '-4', '-1', '-3', '1', '-3', '-1', '-3', '4' ],
             'lR' => [ '-2', '-1', '-3', '0', '-4', '4', '-4', '-3' ],
             'lI' => [ '0', '-3', '0', '-4', '6', '-4', '-2', '-2' ],
             'lG' => [ '-4', '-2', '-4', '-2', '-5', '-3', '-1', '-2' ],
             'pS' => [ '0', '33', '0', '16', '1', '12', '11', '25' ],
             'pF' => [ '0', '0', '2', '0', '3', '0', '0', '0' ],
             'pT' => [ '0', '8', '7', '10', '1', '2', '7', '8' ],
             'pN' => [ '0', '0', '2', '13', '0', '36', '1', '4' ],
             'pK' => [ '0', '5', '0', '13', '1', '15', '0', '2' ],
             'pY' => [ '0', '0', '0', '0', '0', '0', '0', '0' ],
             'pE' => [ '0', '41', '1', '12', '0', '0', '0', '15' ],
             'pV' => [ '0', '3', '9', '0', '2', '0', '3', '1' ],
             'pQ' => [ '0', '0', '0', '15', '0', '4', '0', '3' ],
             'pM' => [ '100', '0', '66', '0', '2', '0', '0', '0' ],
             'pC' => [ '0', '0', '0', '0', '0', '0', '0', '0' ],
             'pL' => [ '0', '0', '8', '0', '25', '0', '4', '0' ],
             'pA' => [ '0', '10', '1', '9', '2', '0', '22', '16' ],
             'pW' => [ '0', '0', '0', '0', '0', '0', '0', '0' ],
             'pP' => [ '0', '0', '0', '0', '3', '1', '45', '0' ],
             'pH' => [ '0', '0', '0', '0', '0', '0', '1', '0' ],
             'pD' => [ '0', '0', '1', '7', '2', '2', '0', '22' ],
             'pR' => [ '0', '0', '0', '3', '0', '27', '0', '0' ],
             'pI' => [ '0', '0', '3', '0', '59', '1', '2', '3' ],
             'pG' => [ '0', '0', '0', '1', '0', '0', '4', '1' ],
   );

   my $matrix = Bio::Matrix::PSM::ProtMatrix( %param );


   my $site = Bio::Matrix::PSM::ProtMatrix->new(%param);
   # Or get it from a file:
   use Bio::Matrix::PSM::IO;
   my $psmIO = Bio::Matrix::PSM::IO->new(-file => $file, -format => 'psi-blast');
   while (my $psm = $psmIO->next_psm) {
      #Now we have a Bio::Matrix::PSM::Psm object, 
      # see Bio::Matrix::PSM::PsmI for details
      #This is a Bio::Matrix::PSM::ProtMatrix object now
      my $matrix = $psm->matrix;   
   }

   # Get a simple consensus, where alphabet is:
   # {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V,}
   # choosing the highest probability or N if prob is too low
   my $consensus = $site->consensus;

   # Retrieving and using regular expressions:
   my $regexp = $site->regexp;
   my $count = grep($regexp,$seq);
   my $count = ($seq=~ s/$regexp/$1/eg);
   print "Motif $mid is present $count times in this sequence\n";

=head1 DESCRIPTION

ProtMatrix is designed to provide some basic methods when working with
position scoring (weight) matrices related to protein sequences.  A
protein PSM consists of 20 vectors with 20 frequencies (one per amino
acid per position).  This is the minimum information you should
provide to construct a PSM object.  The vectors can be provided as
strings with frequencies where the frequency is {0..a} and a=1. This
is the way MEME compressed representation of a matrix and it is quite
useful when working with relational DB.  If arrays are provided as an
input (references to arrays actually) they can be any number, real or
integer (frequency or count).

When creating the object the constructor will check for positions that
equal 0.  If such is found it will increase the count for all
positions by one and recalculate the frequency.  Potential bug - if
you are using frequencies and one of the positions is 0 it will change
significantly.  However, you should never have frequency that equals
0.

Throws an exception if: You mix as an input array and string (for
example A matrix is given as array, C - as string).  The position
vector is (0,0,0,0).  One of the probability vectors is shorter than
the rest.

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

New methods, which might be of interest to anyone who wants to store
PSM in a relational database without creating an entry for each
position is the ability to compress the PSM vector into a string with
losing usually less than 1% of the data.  this can be done with:

   my $str=$matrix->get_compressed_freq('A');
or

   my $str=$matrix->get_compressed_logs('A');

Loading from a database should be done with new, but is not yet implemented.
However you can still uncompress such string with:

   my @arr=Bio::Matrix::PSM::_uncompress_string ($str,1,1); for PSM

or

   my @arr=Bio::Matrix::PSM::_uncompress_string ($str,1000,2); for log odds

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.   Your participation is much appreciated.

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

=head1 AUTHOR - James Thompson

Email tex@biosysadmin.com

=head1 APPENDIX

=cut


# Let the code begin...
package Bio::Matrix::PSM::ProtMatrix;
use strict;

use base qw(Bio::Root::Root Bio::Matrix::PSM::SiteMatrixI);

=head2 new

 Title    : new
 Usage    : my $site = Bio::Matrix::PSM::ProtMatrix->new( 
               %probs,
               %logs,
               -IC    => $ic,
               -e_val => $score, 
               -id    => $mid
               -model => \%model
            );
 Function : Creates a new Bio::Matrix::PSM::ProtMatrix object from memory
 Throws   : If inconsistent data for all vectors (all 20 amino acids) is
               provided, if you mix input types (string vs array) or if a
               position freq is 0.
 Example  :
 Returns  : Bio::Matrix::PSM::ProtMatrix object
 Args     : Hash references to log-odds scores and probabilities for
            position-specific scoring info, e-value (optional), information
            content (optional), id (optional), model for background distribution
            of proteins (optional).

=cut

sub new {
   my ($class, @args) = @_;
   my $self = $class->SUPER::new(@args);
   my $consensus;
   #Too many things to rearrange, and I am creating simultanuously >500 
   # such objects routinely, so this becomes performance issue
   my %input;
   while( @args ) {
     (my $key = shift @args) =~ s/-//gi; #deletes all dashes (only dashes)!
     $input{$key} = shift @args;
   }

   # get a protein alphabet for processing log-odds scores and probabilities
   # maybe change this later on to allow for non-standard aa lists? 
   my @alphabet = qw/A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V/;

   foreach my $aa (@alphabet) {
      $self->{"log$aa"}  = defined($input{"l$aa"}) ? $input{"l$aa"} 
         : $self->throw("Error: No log-odds information for $aa!");
      $self->{"prob$aa"}  = defined($input{"p$aa"}) ? $input{"p$aa"} 
         : $self->throw("Error: No probability information for $aa!");
   }
   
   $self->{_position} = 0;
   $self->{IC}        = $input{IC};
   $self->{e_val}     = $input{e_val};
   $self->{sites}     = $input{sites};
   $self->{width}     = $input{width};
   $self->{accession_number} = $input{accession_number};
   $self->{_correction} = defined($input{correction}) ?
      $input{correction} : 1 ; # Correction might be unwanted- supply your own
   # No id provided, null for the sake of rel db
   $self->{id} = defined($input{id}) ? $input{id} : 'null'; 
   $self->{_alphabet} = \@alphabet;

   #Make consensus, throw if any one of the vectors is shorter
   $self = _calculate_consensus($self,$input{model});
   return $self;
}

=head2 alphabet

 Title    : Returns an array (or array reference if desired) to the alphabet 
 Usage    :
 Function : Returns an array (or array reference) containing all of the
            allowable characters for this matrix.
 Throws   :
 Example  :
 Returns  : Array or arrary reference.
 Args     :

=cut

sub alphabet {
   my $self = shift;
   if ( wantarray ) {
      return $self->{_alphabet};
   } else {
      return @{$self->{_alphabet}};
   }
}

=head2 _calculate_consensus

 Title    : _calculate_consensus
 Usage    :
 Function : Calculates the consensus sequence for this matrix. 
 Throws   :
 Example  :
 Returns  :
 Args     :

=cut

sub _calculate_consensus {
   my $self   = shift;
   my $thresh = shift;
   
   # verify that all of the array lengths in %probs are the same
   my @lengths = map { scalar(@$_) } map {$self->{"prob$_"}} @{ $self->{_alphabet} };
   my $len = shift @lengths;
   for ( @lengths ) {
      if ( $_ ne $len ) { $self->throw( "Probability matrix is damaged!\n" ) };
   }

   # iterate over probs, generate the most likely sequence and put it into 
   # $self->{seq}. Put the probability of this sequence into $self->{seqp}.
   for ( my $i = 0; $i < $len; $i++ ) {
      # get a list of all the probabilities at position $i, ordered by $self->{_alphabet}
      my @probs = map { ${$self->{"prob$_"}}[$i] } @{ $self->{_alphabet} };
      # calculate the consensus of @probs, put sequence into seqp and probabilities into seqp
      (${$self->{seq}}[$i],${$self->{seqp}}[$i]) = $self->_to_cons( @probs, $thresh );
   }

   return $self;
}

=head2 next_pos

 Title    : next_pos
 Usage    :
 Function : Retrives the next position features: frequencies for all 20 amino
            acids, log-odds scores for all 20 amino acids at this position,
            the main (consensus) letter at this position, the probability
            for the consensus letter to occur at this position and the relative
            current position as an integer.
 Throws   :
 Example  :
 Returns  : hash (or hash reference) (pA,pR,pN,pD,...,logA,logR,logN,logD,aa,prob,rel)
            - pN entries represent the probability for amino acid N
            to be at this position
            - logN entries represent the log-odds score for having amino acid
            N at this position
            - aa is the consensus amino acid
            - prob is the probability for the consensus amino acid to be at this 
            position
            - rel is the relative index of the current position (integer)
 Args      : none


=cut

sub next_pos {
      my $self = shift;
         $self->throw("instance method called on class") unless ref $self;

      my $len = @{$self->{seq}};
      my $pos = $self->{_position};

      # return a PSM if we're still within range
      if ($pos<$len) {

         my %probs = map { ("p$_", ${$self->{"prob$_"}}[$pos]) } @{$self->{_alphabet}};
         my %logs  = map { ("l$_", ${$self->{"log$_"}}[$pos]) } @{$self->{_alphabet}};
         my $base  = ${$self->{seq}}[$pos];
         my $prob  = ${$self->{seqp}}[$pos];

         $self->{_position}++;
         my %hash = ( %probs, %logs, base => $base, rel => $pos, prob => $prob );
         
         # decide whether to return the hash or a reference to it
         if ( wantarray ) {
            return %hash;
         } else {
            return \%hash;
         }
   } else { # otherwise, reset $self->{_position} and return nothing
      $self->{_position} = 0; 
      return;
   }
}


=head2 curpos

 Title    : curpos
 Usage    :
 Function : Gets/sets the current position. 
 Throws   :
 Example  :
 Returns  : Current position (integer).
 Args     : New position (integer).

=cut

sub curpos {
      my $self = shift;
      if (@_) { $self->{_position} = shift; }
      return $self->{_position};
}


=head2 e_val

 Title    : e_val
 Usage    :
 Function : Gets/sets the e-value
 Throws   :
 Example  :
 Returns  : 
 Args     : real number

=cut

sub e_val {
      my $self = shift;
      if (@_) { $self->{e_val} = shift; }
      return $self->{e_val};
}


=head2 IC

 Title    : IC
 Usage    :
 Function : Position-specific information content.
 Throws   :
 Example  :
 Returns  : Information content for current position.
 Args     : Information content for current position.

=cut

sub IC {
      my $self = shift;
      if (@_) { $self->{IC} = shift; }
      return $self->{IC};
}

=head2 accession_number

 Title    : accession_number
 Usage    :
 Function: accession number, this will be unique id for the ProtMatrix object as
            well for any other object, inheriting from ProtMatrix.
 Throws   :
 Example  :
 Returns  : New accession number (string)
 Args     : Accession number (string)

=cut

sub accession_number {
      my $self = shift;
      if (@_) { $self->{accession_number} = shift; }
      return $self->{accession_number};
}

=head2 consensus

 Title    : consensus
 Usage    :
 Function : Returns the consensus sequence for this PSM.
 Throws   : if supplied with thresold outisde 5..10 range
 Example  :
 Returns  : string
 Args     : (optional) threshold value 5 to 10 (corresponds to 50-100% at each position

=cut

sub consensus {
   my $self = shift;
   my $thresh=shift;
   $self->_calculate_consensus($thresh) if ($thresh); #Change of threshold
   my $consensus='';

   foreach my $letter (@{$self->{seq}}) {
       $consensus .= $letter;
   }

   return $consensus;
}

sub IUPAC {
   my $self = shift;
   return $self->consensus;
}


=head2 get_string

 Title   : get_string
 Usage   :
 Function: Returns given probability vector as a string. Useful if you want to
            store things in a rel database, where arrays are not first choice
 Throws  : If the argument is outside {A,C,G,T}
 Example :
 Returns : string
 Args    : character {A,C,G,T}

=cut

sub get_string {
   my $self = shift;
   my $base = shift;
   my $string = '';

   my @prob = @{$self->{"prob$base"}};
   if ( ! @prob ) {
      $self->throw( "No such base: $base\n");
   }

   foreach my $prob (@prob) {
      my $corrected = $prob*10;
      my $next = sprintf("%.0f",$corrected);
      $next = 'a' if ($next eq '10');
      $string .= $next;
   }
   return $string;
}



=head2 width

 Title    : width
 Usage    :
 Function : Returns the length of the site
 Throws   :
 Example  :
 Returns  : number
 Args     :

=cut

sub width {
   my $self = shift;
   my $width = @{$self->{probA}};
   return $width;
}

=head2 get_array

 Title    : get_array
 Usage    :
 Function : Returns an array with frequencies for a specified amino acid.
 Throws   :
 Example  :
 Returns  : Array representing frequencies for specified amino acid.
 Args     : Single amino acid (character).

=cut

sub get_array {
   my $self = shift;
   my $letter = uc(shift);

   $self->throw ("No such base: $letter!\n") unless grep { /$letter/ } @{$self->{_alphabet}};

   return @{$self->{"prob$letter"}}; 
}


=head2 get_logs_array

 Title    : get_logs_array
 Usage    :
 Function : Returns an array with log_odds for a specified base
 Throws   :
 Example  :
 Returns  : Array representing log-odds scores for specified amino acid.
 Args     : Single amino acid (character).

=cut

sub get_logs_array {
   my $self = shift;
   my $letter = uc(shift);

   $self->throw ("No such base: $letter!\n") unless grep { /$letter/ } @{$self->{_alphabet}};

   return @{$self->{"log$letter"}}; 
}

=head2 id

 Title    : id
 Usage    :
 Function : Gets/sets the site id
 Throws   :
 Example  :
 Returns  : string
 Args     : string

=cut

sub id {
      my $self = shift;
      if (@_) { $self->{id} = shift; }
      return $self->{id};
}

=head2 regexp

 Title    : regexp
 Usage    :
 Function : Returns a case-insensitive regular expression which matches the
            IUPAC convention.  X's in consensus sequence will match anything.     
 Throws   :
 Example  :
 Returns  : string
 Args     : Threshold for calculating consensus sequence (number in range 0-100
            representing a percentage). Threshold defaults to 20.

=cut

sub regexp {
   my $self = shift;
   my $threshold = 20;
   if ( @_ ) { my $threshold = shift };

   my @alphabet = @{$self->{_alphabet}};
   my $width = $self->width;
   my (@regexp, $i);
   for ( $i = 0; $i < $width; $i++ ) {
      # get an array of the residues at this position with p > $threshold
      my @letters = map { uc($_).lc($_) } grep { $self->{"prob$_"}->[$i] >= $threshold } @alphabet;

      my $reg;
      if ( scalar(@letters) == 0 ) {
         $reg = '\.';
      } else {
         $reg = '['.join('',@letters).']';
      }
      push @regexp, $reg;
   }

   if ( wantarray ) { 
      return @regexp;
   } else {
      return join '', @regexp;
   }
}


=head2 regexp_array

 Title    : regexp_array
 Usage    :
 Function : Returns an array of position-specific regular expressions.
             X's in consensus sequence will match anything.      
 Throws   :
 Example  :
 Returns  : Array of position-specific regular expressions.
 Args     : Threshold for calculating consensus sequence (number in range 0-100
            representing a percentage). Threshold defaults to 20.
 Notes    : Simply calls regexp method in list context.

=cut

sub regexp_array {
   my $self = shift;
   
   return @{ $self->regexp };
}


=head2 _compress_array

 Title    : _compress_array
 Usage    :
 Function :  Will compress an array of real signed numbers to a string (ie vector of bytes)
             -127 to +127 for bi-directional(signed) and 0..255 for unsigned ;
 Throws   :
 Example  :  Internal stuff
 Returns  :  String
 Args     :  array reference, followed by max value and direction (optional, defaults to 1),
             direction of 1 is unsigned, anything else is signed. 

=cut

sub _compress_array {
   my ($array,$lm,$direct)=@_;
   my $str;
   return unless(($array) && ($lm));
   $direct=1 unless ($direct);
   my $k1= ($direct==1) ? (255/$lm) : (127/$lm);
   foreach my $c (@{$array}) {
      $c=$lm if ($c>$lm);
      $c=-$lm if (($c<-$lm) && ($direct !=1));
      $c=0 if (($c<0) && ($direct ==1));
      my $byte=int($k1*$c);
      $byte=127+$byte if ($direct !=1);#Clumsy, should be really shift the bits
      my $char=chr($byte);
      $str.=$char;
   }
   return $str;
}

=head2 _uncompress_string

 Title    : _uncompress_string
 Usage    :
 Function :   Will uncompress a string (vector of bytes) to create an array of real
                  signed numbers (opposite to_compress_array)
 Throws   :
 Example  :   Internal stuff
 Returns  :   string, followed by max value and direction (optional, defaults to 1),
              direction of 1 is unsigned, anything else is signed.
 Args     :   array

=cut

sub _uncompress_string {
   my ($str,$lm,$direct)=@_;
   my @array;
   return unless(($str) && ($lm));
   $direct=1 unless ($direct);
   my $k1= ($direct==1) ? (255/$lm) : (127/$lm);
   while (my $c=chop($str)) {
      my $byte=ord($c);
      $byte=$byte-127 if ($direct !=1);#Clumsy, should be really shift the bits
      my $num=$byte/$k1;
      unshift @array,$num;
   }

   return @array;
}

=head2 get_compressed_freq

 Title    : get_compressed_freq
 Usage    :
 Function:   A method to provide a compressed frequency vector. It uses one byte to
             code the frequence for one of the probability vectors for one position.
             Useful for relational database. Improvment of the previous 0..a coding.
 Throws   :
 Example  :   my $strA=$self->get_compressed_freq('A');
 Returns  :   String
 Args     :   char 

=cut

sub get_compressed_freq {
   my $self=shift;
   my $base=shift;
   my $string='';
   my @prob;
   BASE: {
      if ($base eq 'A') {
         @prob = @{$self->{probA}} unless (!defined($self->{probA}));
         last BASE;
      }
         if ($base eq 'G') {
         @prob = @{$self->{probG}} unless (!defined($self->{probG}));
         last BASE;
      }
         if ($base eq 'C') {
         @prob = @{$self->{probC}} unless (!defined($self->{probC}));
         last BASE;
      }
         if ($base eq 'T') {
         @prob = @{$self->{probT}} unless (!defined($self->{probT}));
         last BASE;
      }
      $self->throw ("No such base: $base!\n");
   }
   my $str= _compress_array(\@prob,1,1);
   return $str;
}

=head2 sequence_match_weight

 Title    : sequence_match_weight
 Usage    :
 Function :   This method will calculate the score of a match, based on the PSM
              if such is associated with the matrix object. Returns undef if no
              PSM data is available.
 Throws   :   if the length of the sequence is different from the matrix width
 Example  :   my $score=$matrix->sequence_match_weight('ACGGATAG');
 Returns  :   Floating point
 Args     :   string

=cut

sub sequence_match_weight {
   my ($self,$seq)=@_;
   return unless ($self->{logA});

   my $seqlen = length($seq);
   my $width  = $self->width;
   $self->throw("Error: Input sequence size ($seqlen) not equal to PSM size ($width)!\n")
      unless (length($seq) == $self->width);

   my ($score,$i) = (0,0);
   foreach my $letter ( split //, $seq ) {
      # add up the score for this position
      $score += $self->{"log$letter"}->[$i];
      $i++;
   }
   return $score;
}


=head2 _to_IUPAC

 Title   : _to_IUPAC
 Usage   :
 Function: Converts a single position to IUPAC compliant symbol and returns its probability.
            Currently returns the most likely amino acid/probability combination.
 Throws  :
 Example :
 Returns : char, real number representing an amino acid and a probability.
 Args    : real numbers for all 20 amino acids (ordered by alphabet contained
            in $self->{_alphabet}, minimum probability threshold.

=cut

sub _to_IUPAC {
   my ($self,@probs,$thresh) = @_;

   # provide a default threshold of 5, corresponds to 5% threshold for 
   # inferring that the aa at any position is the true aa
   $thresh = 5 unless ( defined $thresh );

   my ($IUPAC_aa,$max_prob) = ('X',$thresh);
   for my $aa ( @{$self->{_alphabet}} ) {
      my $prob = shift @probs;
      if ( $prob > $max_prob ) {
         $IUPAC_aa = $aa;
         $max_prob = $prob;
      }
   }
   
   return $IUPAC_aa, $max_prob;
}

=head2 _to_cons

 Title   : _to_cons
 Usage   :
 Function: Converts a single position to simple consensus character and returns
            its probability. Currently just calls the _to_IUPAC subroutine. 
 Throws  :
 Example :
 Returns : char, real number
 Args    : real numbers for A,C,G,T (positional)

=cut

sub _to_cons {
   return _to_IUPAC( @_ );
}

=head2 get_all_vectors

 Title    : get_all_vectors
 Usage    :
 Function :  returns all possible sequence vectors to satisfy the PFM under
             a given threshold
 Throws   :  If threshold outside of 0..1 (no sense to do that)
 Example  :  my @vectors = $self->get_all_vectors(4);
 Returns  :  Array of strings
 Args     :  (optional) floating

=cut

#sub get_all_vectors {
#   my $self = shift;
#   my $thresh = shift;
#
#   $self->throw("Out of range. Threshold should be >0 and 1<.\n") if (($thresh<0) || ($thresh>1));
#
#   my @seq = split(//,$self->consensus($thresh*10));
#   my @perm;
#   for my $i (0..@{$self->{probA}}) {
#      push @{$perm[$i]},'A' if ($self->{probA}->[$i]>$thresh);
#      push @{$perm[$i]},'C' if ($self->{probC}->[$i]>$thresh);
#      push @{$perm[$i]},'G' if ($self->{probG}->[$i]>$thresh);
#      push @{$perm[$i]},'T' if ($self->{probT}->[$i]>$thresh);
#      push @{$perm[$i]},'N' if  ($seq[$i] eq 'N');
#   }
#   my $fpos=shift @perm;
#   my @strings=@$fpos;
#   foreach my $pos (@perm) {
#      my @newstr;
#      foreach my $let (@$pos) {
#         foreach my $string (@strings) {
#           my $newstring = $string . $let;
#           push @newstr,$newstring;
#         }
#      }
#      @strings=@newstr;
#   }
#   return @strings;
#}

1;
