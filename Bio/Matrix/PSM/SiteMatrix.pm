# $Id$
#---------------------------------------------------------

=head1 NAME

Bio::Matrix::PSM::SiteMatrix - SiteMatrixI implementation, holds a
position scoring matrix (or position weight matrix) and log-odds

=head1 SYNOPSIS

  use Bio::Matrix::PSM::SiteMatrix;
  # Create from memory by supplying probability matrix hash
  # both as strings or arrays
  # where the frequencies  $a,$c,$g and $t are supplied either as
  # arrayref or string. Accordingly, lA, lC, lG and lT are the log
  # odds (only as arrays, no checks done right now)
  my ($a,$c,$g,$t,$score,$ic, $mid)=@_; 
  #or
  my ($a,$c,$g,$t,$score,$ic,$mid)=('05a011','110550','400001',
                                    '100104',0.001,19.2,'CRE1');
  #Where a stands for all (this frequency=1), see explanation bellow
  my %param=(-pA=>$a,-pC=>$c,-pG=>$g,-pT=>$t,
             -lA=>$la, -lC=>$lc,-lG=>$lg,-lT=>$l,
             -IC=>$ic,-e_val=>$score, -id=>$mid);
  my $site=new Bio::Matrix::PSM::SiteMatrix(%param);
  #Or get it from a file:
  use Bio::Matrix::PSM::IO;
  my $psmIO= new Bio::Matrix::PSM::IO(-file=>$file, -format=>'transfac');
  while (my $psm=$psmIO->next_psm) {
    #Now we have a Bio::Matrix::PSM::Psm object, 
    # see Bio::Matrix::PSM::PsmI for details
    #This is a Bio::Matrix::PSM::SiteMatrix object now
    my $matrix=$psm->matrix;  
  }

  # Get a simple consensus, where alphabet is {A,C,G,T,N}, 
  # choosing the highest probability or N if prob is too low
  my $consensus=$site->consensus;

  #Getting/using regular expression
  my $regexp=$site->regexp;
  my $count=grep($regexp,$seq);
  my $count=($seq=~ s/$regexp/$1/eg);
  print "Motif $mid is present $count times in this sequence\n";

=head1 DESCRIPTION

SiteMatrix is designed to provide some basic methods when working with
position scoring (weight) matrices, such as transcription factor
binding sites for example.  A DNA PSM consists of four vectors with
frequencies {A,C,G,T). This is the minimum information you should
provide to construct a PSM object. The vectors can be provided as
strings with frequencies where the frequency is {0..a} and a=1. This
is the way MEME compressed representation of a matrix and it is quite
useful when working with relational DB.  If arrays are provided as an
input (references to arrays actually) they can be any number, real or
integer (frequency or count).

When creating the object the constructor will check for positions that
equal 0.  If such is found it will increase the count for all
positions by one and recalculate the frequency.  Potential bug- if you
are using frequencies and one of the positions is 0 it will change
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

  bioperl-l@bioperl.org                 - General discussion
  http://bio.perl.org/MailList.html     - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Stefan Kirov

Email skirov@utk.edu


=head1 APPENDIX

=cut


# Let the code begin...
package Bio::Matrix::PSM::SiteMatrix;
use Bio::Matrix::PSM::SiteMatrixI;
use Bio::Root::Root;
use vars qw(@ISA);
use strict;

@ISA=qw(Bio::Root::Root Bio::Matrix::PSM::SiteMatrixI);

=head2 new

 Title   : new
 Usage   : my $site=new Bio::Matrix::PSM::SiteMatrix(-pA=>$a,-pC=>$c,
						     -pG=>$g,-pT=>$t,
						     -IC=>$ic,
						     -e_val=>$score, 
						     -id=>$mid);
 Function:  Creates a new Bio::Matrix::PSM::SiteMatrix object from memory
 Throws : If inconsistent data for all vectors (A,C,G and T) is
          provided, if you mix input types (string vs array) or if a
          position freq is 0.
 Example :
 Returns :  Bio::Matrix::PSM::SiteMatrix object
 Args    :  hash


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
    $self->{_position}   = 0;
    $self->{IC}     = $input{IC};
    $self->{e_val}  = $input{e_val};
    $self->{sites}  = $input{sites};
    $self->{width}  = $input{width};
    $self->{accession_number}=$input{accession_number};
    $self->{_correction}   =  defined($input{correction}) ? 
	$input{correction} : 1 ; # Correction might be unwanted- supply your own
    # No id provided, null for the sake of rel db
	$self->{logA}   = $input{lA};
	$self->{logC}   = $input{lC};
	$self->{logG}   = $input{lG};
	$self->{logT}   = $input{lT};
    $self->{id}= defined($input{id}) ? $input{id} : 'null'; 
	return $self unless (defined($input{pA}) && defined($input{pC}) && defined($input{pG}) && defined($input{pT}));
#This should go to _initialize?
#Check for input type- no mixing alllowed, throw ex
    if (ref($input{pA}) =~ /ARRAY/i ) {
	$self->throw("Mixing matrix data types not allowed: C is not reference") unless(ref($input{pC}));
	$self->throw("Mixing matrix data types not allowed: G is not reference") unless (ref($input{pG}));
	$self->throw("Mixing matrix data types not allowed: T is not reference") unless (ref($input{pT}));
	$self->{probA}   = $input{pA};
	$self->{probC}   = $input{pC};
	$self->{probG}   = $input{pG};
	$self->{probT}   = $input{pT};
    }
    else {
	$self->throw("Mixing matrix data types not allowed: C is reference") if (ref($input{pC}));
	$self->throw("Mixing matrix data types not allowed: G is reference") if (ref($input{pG}));
	$self->throw("Mixing matrix data types not allowed: T is reference") if (ref($input{pT}));
	$self->{probA}   = [split(//,$input{pA})];
	$self->{probC}   = [split(//,$input{pC})];
	$self->{probG}   = [split(//,$input{pG})];
	$self->{probT}   = [split(//,$input{pT})];
	for (my $i=0; $i<@{$self->{probA}}+1; $i++) {
	    ${$self->{probA}}[$i]='10' if ( ${$self->{probA}}[$i] and ${$self->{probA}}[$i] eq 'a');
	    ${$self->{probC}}[$i]='10' if ( ${$self->{probC}}[$i] and ${$self->{probC}}[$i] eq 'a');
	    ${$self->{probG}}[$i]='10' if ( ${$self->{probG}}[$i] and ${$self->{probG}}[$i] eq 'a');
	    ${$self->{probT}}[$i]='10' if ( ${$self->{probT}}[$i] and ${$self->{probT}}[$i] eq 'a');
	}
#If this is MEME like output(probabilities, rather than count) here is the place for a check
    }
#Check for position with 0 for all bases, throw exception if so
#Correct 0 positions- inc by 1
    for (my $i=0;$i<$#{$self->{probA}}+1;$i++) {
	$self->throw("Position meaningless-all frequencies are 0") if ((${$self->{probA}}[$i]+${$self->{probC}}[$i]+${$self->{probG}}[$i]+${$self->{probT}}[$i])==0);
	$self->{_corrected}= ((${$self->{probA}}[$i]==0) || 
			      (${$self->{probG}}[$i]==0) || 
			      (${$self->{probC}}[$i]==0) || 
			      (${$self->{probT}}[$i]==0));
	if ($self->{_corrected}) {
	    ${$self->{probA}}[$i] += $self->{_correction};
	    ${$self->{probC}}[$i] += $self->{_correction};
	    ${$self->{probG}}[$i] += $self->{_correction};
	    ${$self->{probT}}[$i] += $self->{_correction};
	}
	my $div= ${$self->{probA}}[$i]+ ${$self->{probC}}[$i]+ ${$self->{probG}}[$i]+ ${$self->{probT}}[$i];
	${$self->{probA}}[$i]=${$self->{probA}}[$i]/$div;
	${$self->{probC}}[$i]=${$self->{probC}}[$i]/$div;
	${$self->{probG}}[$i]=${$self->{probG}}[$i]/$div;
	${$self->{probT}}[$i]=${$self->{probT}}[$i]/$div;
    }
#Make consensus, throw if any one of the vectors is shorter
    $self=_calculate_consensus($self);
    return $self;
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
    my $self=shift;
    my $thresh=shift;
    my ($lc,$lt,$lg)=($#{$self->{probC}},$#{$self->{probT}},$#{$self->{probG}});
    my $len=$#{$self->{probA}};
    $self->throw("Probability matrix is damaged for C: $len vs $lc") if ($len != $lc);
    $self->throw("Probability matrix is damaged for T: $len vs $lt") if ($len != $lt);
    $self->throw("Probability matrix is damaged for G: $len vs $lg") if ($len != $lg);
    for (my $i=0; $i<$len+1; $i++) {
	(${$self->{IUPAC}}[$i],${$self->{IUPACp}}[$i])=_to_IUPAC(${$self->{probA}}[$i],${$self->{probC}}[$i],${$self->{probG}}[$i],${$self->{probT}}[$i],$thresh);
	(${$self->{seq}}[$i],${$self->{seqp}}[$i])=_to_cons(${$self->{probA}}[$i],${$self->{probC}}[$i],${$self->{probG}}[$i],${$self->{probT}}[$i],$thresh);
    }
    return $self;
}

=head2 next_pos

 Title   : next_pos
 Usage   :
 Function: Retrives the next position features: frequencies for A,C,G,T, the main
            letter (as in consensus) and the probabilty for this letter to occur at this position
            and the current position
 Throws  :
 Example :
 Returns : hash (pA,pC,pG,pT,base,prob,rel)
 Args    : none


=cut

sub next_pos {
    my $self = shift;
    die "instance method called on class" unless ref $self;
    my $len=@{$self->{seq}};
    my $pos=$self->{_position};
    # End reached?
    if ($self->{_position}<$len) {
	my $pA=${$self->{probA}}[$pos];
	my $pC=${$self->{probC}}[$pos];
	my $pG=${$self->{probG}}[$pos];
	my $pT=${$self->{probT}}[$pos];
	my $base=${$self->{seq}}[$pos];
	my $prob=${$self->{seqp}}[$pos];
	$self->{_position}++;
	my %seq=(pA=>$pA,pT=>$pT,pC=>$pC,pG=>$pG, base=>$base,rel=>$pos, prob=>$prob);
	return %seq;
    }
    else {$self->{_position}=0; return undef;}
}


=head2 curpos

 Title   : curpos
 Usage   :
 Function: Gets/sets the current position. Converts to 0 if argument is minus and
            to width if greater than width
 Throws  :
 Example :
 Returns : integer
 Args    : integer

=cut

sub curpos {
    my $self = shift;
    my $prev = $self->{_position};
    if (@_) { $self->{_position} = shift; }
    return $prev;
}


=head2 e_val

 Title   : e_val
 Usage   :
 Function: Gets/sets the e-value
 Throws  :
 Example :
 Returns : real number
 Args    : real number

=cut

sub e_val {
    my $self = shift;
    my $prev = $self->{e_val};
    if (@_) { $self->{e_val} = shift; }
    return $prev;
}


=head2 IC

 Title   : IC
 Usage   :
 Function: Information content
 Throws  :
 Example :
 Returns : real number
 Args    : none

=cut

sub IC {
    my $self = shift;
    my $prev = $self->{IC};
    if (@_) { $self->{IC} = shift; }
    return $prev;
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
    my $prev = $self->{accession_number};
    if (@_) { $self->{accession_number} = shift; }
    return $prev;
}

=head2 consensus

 Title   : consensus
 Usage   :
 Function: Returns the consensus
 Throws  : if supplied with thresold outisde 5..10 range
 Example :
 Returns : string
 Args    : (optional) threshold value 5 to 10

=cut

sub consensus {
  my $self = shift;
  my $thresh=shift;
  _calculate_consensus($self,$thresh)  if ($thresh); #Change of threshold
  my $consensus='';
  foreach my $letter (@{$self->{seq}}) {
     $consensus .= $letter;
  }
  return $consensus;
}


=head2 width

 Title   : width
 Usage   :
 Function: Returns the length of the site
 Throws  :
 Example :
 Returns : number
 Args    :

=cut

sub width {
  my $self = shift;
  my $width=@{$self->{probA}};
  return $width;
}

=head2 IUPAC

 Title   : IUPAC
 Usage   :
 Function: Returns IUPAC compliant consensus
 Throws  :
 Example :
 Returns : string
 Args    :

=cut

sub IUPAC {
	my $self = shift;
	my $iu=$self->{IUPAC};
	my $iupac='';
	foreach my $let (@{$iu}) {
		$iupac .= $let;
	}
return $iupac;
}

=head2 _to_IUPAC

 Title   : _to_IUPAC
 Usage   :
 Function: Converts a single position to IUPAC compliant symbol and returns its probability.
            For rules see the implementation
 Throws  :
 Example :
 Returns : char, real number
 Args    : real numbers for A,C,G,T (positional)

=cut

sub _to_IUPAC {
	my $A=shift;
	my $C=shift;
	my $G=shift;
	my $T=shift;
	my $all=$A+$G+$C+$T;
	my $a=$A/$all;
	my $g=$G/$all;
	my $c=$C/$all;
	my $t=$T/$all;
	my $single=0.7*$all;
	my $double=0.8*$all;
	my $triple=0.9*$all;
	return 'A',$a if ($a>$single);
	return 'G',$g if ($g>$single);
	return 'C',$c if ($c>$single);
	return 'T',$t if ($t>$single);
	my $r=$g+$a;
	return 'R',$r if ($r>$double);
	my $y=$t+$c;
	return 'Y',$y if ($y>$double);
	my $m=$a+$c;
	return 'M',$m if ($m>$double);
	my $k=$g+$t;
	return 'K',$k if ($k>$double);
	my $s=$g+$c;
	return 'S',$s if ($s>$double);
	my $w=$a+$t;
	return 'W',$w if ($w>$double);
	my $d=$r+$t;
	return 'D',$d if ($d>$triple);
	my $v=$r+$c;
	return 'V',$v if ($v>$triple);
	my $b=$y+$g;
	return 'B',$b if ($b>$triple);
	my $h=$y+$a;
	return 'H',$h if ($h>$triple);
	return 'N',0;
}

=head2 _to_cons

 Title   : _to_cons
 Usage   :
 Function: Converts a single position to simple consensus character and returns its probability.
            For rules see the implementation
 Throws  :
 Example :
 Returns : char, real number
 Args    : real numbers for A,C,G,T (positional)

=cut

sub _to_cons {
	my $A=shift;
	my $C=shift;
	my $G=shift;
	my $T=shift;
  my $thresh=shift;
	my $all=$A+$G+$C+$T;
	my $a=$A*10/$all;
	my $g=$G*10/$all;
	my $c=$C*10/$all;
	my $t=$T*10/$all;
 #Check for user supplied threshold
  $thresh=5 unless ($thresh);  #default threshold
  return 'N',10 if (($a<$thresh) &&  ($c<$thresh) &&  ($g<$thresh) && ($t<$thresh));
	return 'A',$a if ($a>$thresh);
	return 'G',$g if ($g>$thresh);
	return 'C',$c if ($c>$thresh);
	return 'T',$t if ($t>$thresh);
	return 'N',10 if (($a==$t) && ($a==$c) && ($a==$g));
	return 'A',$a if (($a>$t) &&($a>$c) && ($a>$g)); #Is this a good idea?
	return 'C',$c if (($c>$t) &&($c>$a) && ($c>$g));
	return 'G',$g if (($g>$t) &&($g>$c) && ($g>$a));
	return 'N',10;
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
	my $self=shift;
	my $base=shift;
	my $string='';
	my @prob;
	BASE: {
		if ($base eq 'A') {@prob= @{$self->{probA}}; last BASE; }
		if ($base eq 'C') {@prob= @{$self->{probC}}; last BASE; }
		if ($base eq 'G') {@prob= @{$self->{probG}}; last BASE; }
		if ($base eq 'T') {@prob= @{$self->{probT}}; last BASE; }
		$self->throw ("No such base: $base!\n");
	}
foreach  my $prob (@prob) {
	my $corrected=$prob*10;
	my $next=sprintf("%.0f",$corrected);
	$next='a' if ($next eq '10');
	$string .= $next;
}
return $string;
}


=head2 get_array

 Title   : get_array
 Usage   :
 Function: Returns an array with frequencies for a specified base
 Throws  :
 Example :
 Returns : array
 Args    : char

=cut

sub get_array {
	my $self=shift;
	my $base=uc(shift);
	return  @{$self->{probA}} if ($base eq 'A');
	return  @{$self->{probC}} if ($base eq 'C');
	return  @{$self->{probG}} if ($base eq 'G');
	return  @{$self->{probT}} if ($base eq 'T');
	$self->throw ("No such base: $base!\n");
}


=head2 get_logs_array

 Title   : get_logs_array
 Usage   :
 Function: Returns an array with log_odds for a specified base
 Throws  :
 Example :
 Returns : array
 Args    : char

=cut

sub get_logs_array {
	my $self=shift;
	my $base=uc(shift);
	return  @{$self->{logA}} if (($base eq 'A')  && ($self->{logA}));
	return  @{$self->{logC}} if (($base eq 'C')  && ($self->{logC}));
	return  @{$self->{logG}} if (($base eq 'G')  && ($self->{logG}));
	return  @{$self->{logT}} if (($base eq 'T')  && ($self->{logT}));
	$self->throw ("No such base: $base!\n") if (!grep(/$base/,qw(A C G T)));
  return ();
}

=head2 id

 Title   : id
 Usage   :
 Function: Gets/sets the site id
 Throws  :
 Example :
 Returns : string
 Args    : string

=cut

sub id {
    my $self = shift;
    my $prev = $self->{id};
    if (@_) { $self->{id} = shift; }
    return $prev;
}


=head2 regexp

 Title   : regexp
 Usage   :
 Function: Returns a regular expression which matches the IUPAC convention.
            N will match X, N, - and .
 Throws  :
 Example :
 Returns : string
 Args    :

=cut

sub regexp {
	my $self=shift;
	my $regexp;
	foreach my $letter (@{$self->{IUPAC}}) {
		my $reg;
		LETTER: {
			if ($letter eq 'A') { $reg='Aa'; last LETTER; }
			if ($letter eq 'C') { $reg='Cc'; last LETTER; }
			if ($letter eq 'G') { $reg='Gg'; last LETTER; }
			if ($letter eq 'T') { $reg='Tt'; last LETTER; }
			if ($letter eq 'M') { $reg='AaCc'; last LETTER; }
			if ($letter eq 'R') { $reg='AaGg'; last LETTER; }
			if ($letter eq 'W') { $reg='AaTt'; last LETTER; }
			if ($letter eq 'S') { $reg='CcGg'; last LETTER; }
			if ($letter eq 'Y') { $reg='CcTt'; last LETTER; }
			if ($letter eq 'K') { $reg='GgTt'; last LETTER; }
			if ($letter eq 'V') { $reg='AaCcGg'; last LETTER; }
			if ($letter eq 'H') { $reg='AaCcTt'; last LETTER; }
			if ($letter eq 'D') { $reg='AaGgTt'; last LETTER; }
			if ($letter eq 'B') { $reg='CcGgTt'; last LETTER; }
			 $reg="\.";
		}
		$regexp .= "[$reg]";
	}
return $regexp;
}


=head2 regexp_array

 Title   : regexp_array
 Usage   :
 Function: Returns a regular expression which matches the IUPAC convention.
            N will match X, N, - and .
 Throws  :
 Example :
 Returns : array
 Args    :
 To do   : I have separated regexp and regexp_array, but
           maybe they can be rewritten as one - just check what should be returned

=cut

sub regexp_array {
	my $self=shift;
	my @regexp;
	foreach my $letter (@{$self->{IUPAC}}) {
		my $reg;
		LETTER: {
			if ($letter eq 'A') { $reg='Aa'; last LETTER; }
			if ($letter eq 'C') { $reg='Cc'; last LETTER; }
			if ($letter eq 'G') { $reg='Gg'; last LETTER; }
			if ($letter eq 'T') { $reg='Tt'; last LETTER; }
			if ($letter eq 'M') { $reg='AaCc'; last LETTER; }
			if ($letter eq 'R') { $reg='AaGg'; last LETTER; }
			if ($letter eq 'W') { $reg='AaTt'; last LETTER; }
			if ($letter eq 'S') { $reg='CcGg'; last LETTER; }
			if ($letter eq 'Y') { $reg='CcTt'; last LETTER; }
			if ($letter eq 'K') { $reg='GgTt'; last LETTER; }
			if ($letter eq 'V') { $reg='AaCcGg'; last LETTER; }
			if ($letter eq 'H') { $reg='AaCcTt'; last LETTER; }
			if ($letter eq 'D') { $reg='AaGgTt'; last LETTER; }
			if ($letter eq 'B') { $reg='CcGgTt'; last LETTER; }
			$reg="\."; 
		}
		push @regexp,$reg;
	}
return @regexp;
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
	my ($array,$lm,$direct)=@_;
	my $str;
	return undef unless(($array) && ($lm));
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
	my ($str,$lm,$direct)=@_;
	my @array;
	return undef unless(($str) && ($lm));
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
	my $self=shift;
	my $base=shift;
	my $string='';
	my @prob;
	BASE: {
		if ($base eq 'A') {
      @prob= @{$self->{probA}} unless (!defined($self->{probA}));
      last BASE;
    }
  		if ($base eq 'G') {
      @prob= @{$self->{probG}} unless (!defined($self->{probG}));
      last BASE;
    }
  		if ($base eq 'C') {
      @prob= @{$self->{probC}} unless (!defined($self->{probC}));
      last BASE;
    }
  		if ($base eq 'T') {
      @prob= @{$self->{probT}} unless (!defined($self->{probT}));
      last BASE;
    }
		$self->throw ("No such base: $base!\n");
	}
	my $str= _compress_array(\@prob,1,1);
  return $str;
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
	my $self=shift;
	my $base=shift;
	my $string='';
	my @prob;
	BASE: {
		if ($base eq 'A') {@prob= @{$self->{logA}} unless (!defined($self->{logA})); last BASE; }
		if ($base eq 'C') {@prob= @{$self->{logC}} unless (!defined($self->{logC})); last BASE; }
		if ($base eq 'G') {@prob= @{$self->{logG}} unless (!defined($self->{logG})); last BASE; }
		if ($base eq 'T') {@prob= @{$self->{logT}} unless (!defined($self->{logT})); last BASE; }
		$self->throw ("No such base: $base!\n");
	}
	return _compress_array(\@prob,1000,2);
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
my ($self,$seq)=@_;
return undef unless ($self->{logA});
$self->throw ("I can calculate the score only for sequence which are exactly my size\n") unless (length($seq)==@{$self->{logA}});
my @seq=split(//,$seq);
my ($score,$i);
foreach my $pos (@seq) {
  my $tv='log' . $pos;
  $score+=$self->{$tv}->[$i];
  $i++;
}
return $score;
}


=head2 get_all_vectors

 Title   : get_all_vectors
 Usage   :
 Function:  returns all possible sequence vectors to satisfy the PFM under
            a given threshold
 Throws  :  If threshold outside of 3..7 (no sense to do that)
 Example :  my @vectors=$self->get_all_vectors(4);
 Returns :  Array of strings
 Args    :  (optional) floating

=cut

sub get_all_vectors {
	my $self=shift;
	my $thresh=shift;
  $self->throw("Out of range. Threshold should be >3 and 7<.\n") if (($thresh<3) || ($thresh>7));
  my @seq=split(//,$self->consensus($thresh));
  my @perm;
  $thresh=$thresh/10;
  for my $i (0..@{$self->{probA}}) {
    push @{$perm[$i]},'A' if ($self->{probA}->[$i]>$thresh);
    push @{$perm[$i]},'C' if ($self->{probC}->[$i]>$thresh);
    push @{$perm[$i]},'G' if ($self->{probG}->[$i]>$thresh);
    push @{$perm[$i]},'T' if ($self->{probT}->[$i]>$thresh);
    push @{$perm[$i]},'N' if  ($seq[$i] eq 'N');
  }
  my $fpos=shift @perm;
  my @strings=@$fpos;
  foreach my $pos (@perm) {
    my @newstr;
    foreach my $let (@$pos) {
      foreach my $string (@strings) {
        my $newstring = $string . $let;
        push @newstr,$newstring;
      }
    }
    @strings=@newstr;
  }
	return @strings;
}


1;
