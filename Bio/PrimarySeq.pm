#
# bioperl module for Bio::PrimarySeq
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PrimarySeq - Bioperl lightweight Sequence Object

=head1 SYNOPSIS

  # Bio::SeqIO for file reading, Bio::DB::GenBank for
  # database reading

  use Bio::Seq;
  use Bio::SeqIO;
  use Bio::DB::GenBank;

  # make from memory

  $seqobj = Bio::PrimarySeq->new ( -seq => 'ATGGGGTGGGCGGTGGGTGGTTTG',
                                   -id  => 'GeneFragment-12',
                                   -accession_number => 'X78121',
                                   -alphabet => 'dna',
                                   -is_circular => 1 );
  print "Sequence ", $seqobj->id(), " with accession ",
    $seqobj->accession_number, "\n";

  # read from file

  $inputstream = Bio::SeqIO->new(-file => "myseq.fa",
                                 -format => 'Fasta');
  $seqobj = $inputstream->next_seq();
  print "Sequence ", $seqobj->id(), " and desc ", $seqobj->desc, "\n";

  # to get out parts of the sequence.

  print "Sequence ", $seqobj->id(), " with accession ",
    $seqobj->accession_number, " and desc ", $seqobj->desc, "\n";

  $string  = $seqobj->seq();
  $string2 = $seqobj->subseq(1,40);

=head1 DESCRIPTION

PrimarySeq is a lightweight Sequence object, storing the sequence, its
name, a computer-useful unique name, and other fundamental attributes.
It does not contain sequence features or other information.  To have a
sequence with sequence features you should use the Seq object which uses
this object.

Although new users will use Bio::PrimarySeq a lot, in general you will
be using it from the Bio::Seq object. For more information on Bio::Seq
see L<Bio::Seq>. For interest you might like to know that
Bio::Seq has-a Bio::PrimarySeq and forwards most of the function calls
to do with sequence to it (the has-a relationship lets us get out of a
otherwise nasty cyclical reference in Perl which would leak memory).

Sequence objects are defined by the Bio::PrimarySeqI interface, and this
object is a pure Perl implementation of the interface. If that's
gibberish to you, don't worry. The take home message is that this
object is the bioperl default sequence object, but other people can
use their own objects as sequences if they so wish. If you are
interested in wrapping your own objects as compliant Bioperl sequence
objects, then you should read the Bio::PrimarySeqI documentation

The documentation of this object is a merge of the Bio::PrimarySeq and
Bio::PrimarySeqI documentation.  This allows all the methods which you can
call on sequence objects here.

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

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PrimarySeq;
use vars qw($MATCHPATTERN $GAP_SYMBOLS);
use strict;

$MATCHPATTERN = 'A-Za-z\-\.\*\?=~';
$GAP_SYMBOLS = '-~';

use base qw(Bio::Root::Root Bio::PrimarySeqI
            Bio::IdentifiableI Bio::DescribableI);

#
# setup the allowed values for alphabet()
#

my %valid_type = map {$_, 1} qw( dna rna protein );

=head2 new

 Title   : new
 Usage   : $seq    = Bio::PrimarySeq->new( -seq => 'ATGGGGGTGGTGGTACCCT',
                                           -id  => 'human_id',
                                           -accession_number => 'AL000012',
                                           );

 Function: Returns a new primary seq object from
           basic constructors, being a string for the sequence
           and strings for id and accession_number.

           Note that you can provide an empty sequence string. However, in
           this case you MUST specify the type of sequence you wish to
           initialize by the parameter -alphabet. See alphabet() for possible
           values.
 Returns : a new Bio::PrimarySeq object
 Args    : -seq         => sequence string
           -display_id  => display id of the sequence (locus name)
           -accession_number => accession number
           -primary_id  => primary id (Genbank id)
           -version     => version number
           -namespace   => the namespace for the accession
           -authority   => the authority for the namespace
           -description => description text
           -desc        => alias for description
           -alphabet    => sequence type (alphabet) (dna|rna|protein)
           -id          => alias for display id
           -is_circular => boolean field for whether or not sequence is circular
           -direct      => boolean field for directly setting sequence (requires alphabet also set)
           -ref_to_seq  => boolean field indicating the sequence is a reference (?!?)
           -nowarnonempty => boolean field for whether or not to warn when sequence is empty

=cut


sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my($seq,$id,$acc,$pid,$ns,$auth,$v,$oid,
       $desc,$description,
       $alphabet,$given_id,$is_circular,$direct,$ref_to_seq,$len,$nowarnonempty) =
      $self->_rearrange([qw(SEQ
                            DISPLAY_ID
                            ACCESSION_NUMBER
                            PRIMARY_ID
                            NAMESPACE
                            AUTHORITY
                            VERSION
                            OBJECT_ID
                            DESC
                            DESCRIPTION
                            ALPHABET
                            ID
                            IS_CIRCULAR
                            DIRECT
                            REF_TO_SEQ
                            LENGTH
                            NOWARNONEMPTY
                            )],
                            @args);
  
    # private var _nowarnonempty, need to be set before calling _guess_alphabet
    $self->{'_nowarnonempty'} = $nowarnonempty; 

    if( defined $id && defined $given_id ) {
      if( $id ne $given_id ) {
        $self->throw("Provided both id and display_id constructor ".
            "functions. [$id] [$given_id]");
      }
    }
    if( defined $given_id ) { $id = $given_id; }

    # let's set the length before the seq -- if there is one, this length is
    # going to be invalidated
    defined $len && $self->length($len);

    # if alphabet is provided we set it first, so that it won't be guessed
    # when the sequence is set
    $alphabet && $self->alphabet($alphabet);

    # bernd's idea: define ids so that invalid sequence messages
    # can be more informative...
    defined $id  && $self->display_id($id);
    $acc         && $self->accession_number($acc);
    defined $pid && $self->primary_id($pid);

    # if there is an alphabet, and direct is passed in, assume the alphabet
    # and sequence is ok
    
    if( $direct && $ref_to_seq) {
      $self->{'seq'} = $$ref_to_seq;
        if( ! $alphabet ) {
          $self->_guess_alphabet();
        } # else it has been set already above
    } else {
        # print STDERR "DEBUG: setting sequence to [$seq]\n";
        # note: the sequence string may be empty
        $self->seq($seq) if defined($seq);
    }

    $desc        && $self->desc($desc);
    $description && $self->description($description);
    $is_circular && $self->is_circular($is_circular);
    $ns          && $self->namespace($ns);
    $auth        && $self->authority($auth);
    defined($v)  && $self->version($v);
    defined($oid) && $self->object_id($oid);


    return $self;
}

sub direct_seq_set {
    my $obj = shift;
    return $obj->{'seq'} = shift if @_;
    return;
}


=head2 seq

 Title   : seq
 Usage   : $string    = $obj->seq()
 Function: Returns the sequence as a string of letters. The
           case of the letters is left up to the implementer.
           Suggested cases are upper case for proteins and lower case for
           DNA sequence (IUPAC standard), but you should not rely on this.
 Returns : A scalar
 Args    : Optionally on set the new value (a string). An optional second
           argument presets the alphabet (otherwise it will be guessed).

=cut

sub seq {
   my ($obj,@args) = @_;

   if( scalar(@args) == 0 ) {
       return $obj->{'seq'};
   }

   my ($value,$alphabet) = @args;

   if(@args) {
       if(defined($value) && (! $obj->validate_seq($value))) {
           $obj->throw("Attempting to set the sequence '".(defined($obj->id) ||
               "[unidentified sequence]")."' to [$value] which does not look healthy");
       }
       # if a sequence was already set we make sure that we re-adjust the
       # alphabet, otherwise we skip guessing if alphabet is already set
       # note: if the new seq is empty or undef, we don't consider that a
       # change (we wouldn't have anything to guess on anyway)
       my $is_changed_seq =
         exists($obj->{'seq'}) && (CORE::length($value || '') > 0);
       $obj->{'seq'} = $value;
       # new alphabet overridden by arguments?
       if($alphabet) {
           # yes, set it no matter what
           $obj->alphabet($alphabet);
       } elsif ($is_changed_seq || (! defined($obj->alphabet()))) {
           # if we changed a previous sequence to a new one or if there is no
           # alphabet yet at all, we need to guess the (possibly new) alphabet
           $obj->_guess_alphabet();
       } # else (seq not changed and alphabet was defined) do nothing
       # if the seq is changed, make sure we unset a possibly set length
       $obj->length(undef) if $is_changed_seq || $obj->{'seq'};
   }
   return $obj->{'seq'};
}

=head2 validate_seq

 Title   : validate_seq
 Usage   : if(! $seq->validate_seq($seq_str) ) {
                print "sequence $seq_str is not valid for an object of
                alphabet ",$seq->alphabet, "\n";
           }
 Function: Validates a given sequence string. A validating sequence string
           must be accepted by seq(). A string that does not validate will
           lead to an exception if passed to seq().

           The implementation provided here does not take alphabet() into
           account. Allowed are all letters (A-Z) and '-','.','*','?','=',
           and '~'.

 Example :
 Returns : 1 if the supplied sequence string is valid for the object, and
           0 otherwise.
 Args    : The sequence string to be validated.


=cut

sub validate_seq {
    my ($self,$seqstr) = @_;
    if( ! defined $seqstr ){ $seqstr = $self->seq(); }
    return 0 unless( defined $seqstr);
    if((CORE::length($seqstr) > 0) &&
        ($seqstr !~ /^([$MATCHPATTERN]+)$/)) {
            $self->warn("sequence '".(defined($self->id) || "[unidentified sequence]").
            "' doesn't validate, mismatch is " .
            join(",",($seqstr =~ /([^$MATCHPATTERN]+)/g)));
        return 0;
    }
    return 1;
}

=head2 subseq

 Title   : subseq
 Usage   : $substring = $obj->subseq(10,40);
           $substring = $obj->subseq(10,40,NOGAP)
           $substring = $obj->subseq(-START=>10,-END=>40,-REPLACE_WITH=>'tga')
 Function: returns the subseq from start to end, where the first sequence
           character has coordinate 1 number is inclusive, ie 1-2 are the 
           first two characters of the sequence
 Returns : a string
 Args    : integer for start position
           integer for end position
                 OR
           Bio::LocationI location for subseq (strand honored)
           Specify -NOGAP=>1 to return subseq with gap characters removed
           Specify -REPLACE_WITH=>$new_subseq to replace the subseq returned
           with $new_subseq in the sequence object

=cut

sub subseq {
   my $self = shift;
   my @args = @_;
   my ($start,$end,$nogap,$replace) = $self->_rearrange([qw(START 
                                                            END
                                                            NOGAP
                                                            REPLACE_WITH)],@args);
   
   # if $replace is specified, have the constructor validate it as seq
   my $dummy = new Bio::PrimarySeq(-seq=>$replace, -alphabet=>$self->alphabet) if defined($replace);

   if( ref($start) && $start->isa('Bio::LocationI') ) {
       my $loc = $start;
       my $seq = "";
       foreach my $subloc ($loc->each_Location()) {
           my $piece = $self->subseq(-START=>$subloc->start(),
                                     '-END'=>$subloc->end(), 
                                     -REPLACE_WITH=>$replace,
                                     -NOGAP=>$nogap);
           $piece =~ s/[$GAP_SYMBOLS]//g if $nogap;
           if($subloc->strand() < 0) {
               $piece = Bio::PrimarySeq->new('-seq' => $piece)->revcom()->seq();
           }
           $seq .= $piece;
       }
       return $seq;
   } elsif(  defined  $start && defined $end ) {
       if( $start > $end ){
           $self->throw("Bad start,end parameters. Start [$start] has to be ".
             "less than end [$end]");
       }
       if( $start <= 0 ) {
           $self->throw("Bad start parameter ($start). Start must be positive.");
       }

       # remove one from start, and then length is end-start
       $start--;
       my @ss_args = map { eval "defined $_"  ? $_ : () } qw( $self->{seq} $start $end-$start $replace);
       my $seqstr = eval join( '', "substr(", join(',',@ss_args), ")");

       if( $end > $self->length) {
           if ($self->is_circular) {
               my $start = 0;
               my $end = $end - $self->length;
               my @ss_args = map { eval "defined $_"  ? $_ : () } qw( $self->{seq} $start $end-$start $replace);
               my $appendstr = eval join( '', "substr(", join(',',@ss_args), ")");
               $seqstr .= $appendstr;
       } else {
           $self->throw("Bad end parameter ($end). End must be less than the total length of sequence (total=".$self->length.")")
           }
       } 

       $seqstr =~ s/[$GAP_SYMBOLS]//g if ($nogap);
       return $seqstr;

   } else {
       $self->warn("Incorrect parameters to subseq - must be two integers or a Bio::LocationI object. Got:", $self,$start,$end,$replace,$nogap);
       return;
   }
}

=head2 length

 Title   : length
 Usage   : $len = $seq->length();
 Function: Get the length of the sequence in number of symbols (bases
           or amino acids).

           You can also set this attribute, even to a number that does
           not match the length of the sequence string. This is useful
           if you don''t want to set the sequence too, or if you want
           to free up memory by unsetting the sequence. In the latter
           case you could do e.g.

               $seq->length($seq->length);
               $seq->seq(undef);

           Note that if you set the sequence to a value other than
           undef at any time, the length attribute will be
           invalidated, and the length of the sequence string will be
           reported again. Also, we won''t let you lie about the length.

 Example :
 Returns : integer representing the length of the sequence.
 Args    : Optionally, the value on set

=cut

sub length {
    my $self = shift;
    my $len = CORE::length($self->seq() || '');

    if(@_) {
        my $val = shift;
        if(defined($val) && $len && ($len != $val)) {
            $self->throw("You're trying to lie about the length: ".
                "is $len but you say ".$val);
        }
        $self->{'_seq_length'} = $val;
    } elsif(defined($self->{'_seq_length'})) {
        return $self->{'_seq_length'};
    }
    return $len;
}

=head2 display_id

 Title   : display_id or display_name
 Usage   : $id_string = $obj->display_id();
 Function: returns the display id, aka the common name of the Sequence object.

           The semantics of this is that it is the most likely string to
           be used as an identifier of the sequence, and likely to have
           "human" readability.  The id is equivalent to the ID field of
           the GenBank/EMBL databanks and the id field of the
           Swissprot/sptrembl database. In fasta format, the >(\S+) is
           presumed to be the id, though some people overload the id to
           embed other information. Bioperl does not use any embedded
           information in the ID field, and people are encouraged to use
           other mechanisms (accession field for example, or extending
           the sequence object) to solve this.

           With the new Bio::DescribeableI interface, display_name aliases
           to this method.

 Returns : A string
 Args    : None


=cut

sub display_id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'display_id'} = $value;
    }
    return $obj->{'display_id'};
}

=head2 accession_number

 Title   : accession_number or object_id
 Usage   : $unique_key = $obj->accession_number;
 Function: Returns the unique biological id for a sequence, commonly
           called the accession_number. For sequences from established
           databases, the implementors should try to use the correct
           accession number. Notice that primary_id() provides the
           unique id for the implemetation, allowing multiple objects
           to have the same accession number in a particular implementation.

           For sequences with no accession number, this method should
           return "unknown".

           [Note this method name is likely to change in 1.3]

           With the new Bio::IdentifiableI interface, this is aliased
           to object_id

 Returns : A string
 Args    : A string (optional) for setting

=cut

sub accession_number {
    my( $obj, $acc ) = @_;

    if (defined $acc) {
        $obj->{'accession_number'} = $acc;
    } else {
        $acc = $obj->{'accession_number'};
        $acc = 'unknown' unless defined $acc;
    }
    return $acc;
}

=head2 primary_id

 Title   : primary_id
 Usage   : $unique_key = $obj->primary_id;
 Function: Returns the unique id for this object in this
           implementation. This allows implementations to manage their
           own object ids in a way the implementaiton can control
           clients can expect one id to map to one object.

           For sequences with no natural primary id, this method
           should return a stringified memory location.

 Returns : A string
 Args    : A string (optional, for setting)

=cut

sub primary_id {
    my $obj = shift;

    if(@_) {
        $obj->{'primary_id'} = shift;
    }
    if( ! defined($obj->{'primary_id'}) ) {
        return "$obj";
    }
    return $obj->{'primary_id'};
}


=head2 alphabet

 Title   : alphabet
 Usage   : if( $obj->alphabet eq 'dna' ) { /Do Something/ }
 Function: Get/Set the alphabet of sequence, one of
           'dna', 'rna' or 'protein'. This is case sensitive.

           This is not called <type> because this would cause
           upgrade problems from the 0.5 and earlier Seq objects.

 Returns : a string either 'dna','rna','protein'. NB - the object must
           make a call of the type - if there is no alphabet specified it
           has to guess.
 Args    : optional string to set : 'dna' | 'rna' | 'protein'


=cut

sub alphabet {
    my ($obj,$value) = @_;
    if (defined $value) {
        $value = lc $value;
        unless ( $valid_type{$value} ) {
        $obj->throw("Alphabet '$value' is not a valid alphabet (".
            join(',', map "'$_'", sort keys %valid_type) .") lowercase");
        }
        $obj->{'alphabet'} = $value;
    }
    return $obj->{'alphabet'};
}

=head2 desc

 Title   : desc or description
 Usage   : $obj->desc($newval)
 Function: Get/set description of the sequence.

           'description' is an alias for this for compliance with the
           Bio::DescribeableI interface.

 Example :
 Returns : value of desc (a string)
 Args    : newvalue (a string or undef, optional)


=cut

sub desc{
    my $self = shift;

    return $self->{'desc'} = shift if @_;
    return $self->{'desc'};
}

=head2 can_call_new

 Title   : can_call_new
 Usage   :
 Function:
 Example :
 Returns : true
 Args    :


=cut

sub can_call_new {
   my ($self) = @_;

   return 1;
}

=head2 id

 Title   : id
 Usage   : $id = $seq->id()
 Function: This is mapped on display_id
 Example :
 Returns :
 Args    :


=cut

sub  id {
   return shift->display_id(@_);
}

=head2 is_circular

 Title   : is_circular
 Usage   : if( $obj->is_circular) { /Do Something/ }
 Function: Returns true if the molecule is circular
 Returns : Boolean value
 Args    : none

=cut

sub is_circular{
    my $self = shift;

    return $self->{'is_circular'} = shift if @_;
    return $self->{'is_circular'};
}


=head1 Methods for Bio::IdentifiableI compliance

=cut

=head2 object_id

 Title   : object_id
 Usage   : $string    = $obj->object_id()
 Function: A string which represents the stable primary identifier
           in this namespace of this object. For DNA sequences this
           is its accession_number, similarly for protein sequences.

           This is aliased to accession_number().
 Returns : A scalar


=cut

sub object_id {
    return shift->accession_number(@_);
}

=head2 version

 Title   : version
 Usage   : $version    = $obj->version()
 Function: A number which differentiates between versions of
           the same object. Higher numbers are considered to be
           later and more relevant, but a single object described
           the same identifier should represent the same concept.

 Returns : A number

=cut

sub version{
    my ($self,$value) = @_;
    if( defined $value) {
        $self->{'_version'} = $value;
    }
    return $self->{'_version'};
}


=head2 authority

 Title   : authority
 Usage   : $authority    = $obj->authority()
 Function: A string which represents the organisation which
           granted the namespace, written as the DNS name for
           organisation (eg, wormbase.org).

 Returns : A scalar

=cut

sub authority {
    my ($obj,$value) = @_;
    if( defined $value) {
        $obj->{'authority'} = $value;
    }
    return $obj->{'authority'};
}

=head2 namespace

 Title   : namespace
 Usage   : $string    = $obj->namespace()
 Function: A string representing the name space this identifier
           is valid in, often the database name or the name
           describing the collection.

 Returns : A scalar


=cut

sub namespace{
    my ($self,$value) = @_;
    if( defined $value) {
        $self->{'namespace'} = $value;
    }
    return $self->{'namespace'} || "";
}

=head1 Methods for Bio::DescribableI compliance

This comprises of display_name and description.

=cut

=head2 display_name

 Title   : display_name
 Usage   : $string    = $obj->display_name()
 Function: A string which is what should be displayed to the user.
           The string should have no spaces (ideally, though a cautious
           user of this interface would not assumme this) and should be
           less than thirty characters (though again, double checking
           this is a good idea).

           This is aliased to display_id().
 Returns : A scalar

=cut

sub display_name {
    return shift->display_id(@_);
}

=head2 description

 Title   : description
 Usage   : $string    = $obj->description()
 Function: A text string suitable for displaying to the user a
           description. This string is likely to have spaces, but
           should not have any newlines or formatting - just plain
           text. The string should not be greater than 255 characters
           and clients can feel justified at truncating strings at 255
           characters for the purposes of display.

           This is aliased to desc().
 Returns : A scalar

=cut

sub description {
    return shift->desc(@_);
}

=head1 Methods Inherited from Bio::PrimarySeqI

These methods are available on Bio::PrimarySeq, although they are
actually implemented on Bio::PrimarySeqI

=head2 revcom

 Title   : revcom
 Usage   : $rev = $seq->revcom()
 Function: Produces a new Bio::SeqI implementing object which
           is the reversed complement of the sequence. For protein
           sequences this throws an exception of
           "Sequence is a protein. Cannot revcom".

           The id is the same id as the orginal sequence, and the
           accession number is also indentical. If someone wants to
           track that this sequence has be reversed, it needs to
           define its own extensions.

           To do an inplace edit of an object you can go:

           $seqobj = $seqobj->revcom();

           This of course, causes Perl to handle the garbage
           collection of the old object, but it is roughly speaking as
           efficient as an inplace edit.

 Returns : A new (fresh) Bio::SeqI object
 Args    : none

=cut

=head2 trunc

 Title   : trunc
 Usage   : $subseq = $myseq->trunc(10,100);
 Function: Provides a truncation of a sequence,

 Example :
 Returns : A fresh Bio::SeqI implementing object.
 Args    :


=cut

=head1 Internal methods

These are internal methods to PrimarySeq

=cut

=head2 _guess_alphabet

 Title   : _guess_alphabet
 Usage   :
 Function: Automatically guess and set the type of sequence: dna, rna or protein
 Example :
 Returns : one of strings 'dna', 'rna' or 'protein'.
 Args    : none


=cut

sub _guess_alphabet {
    my ($self) = @_;
    my $type;

    # Remove char's that clearly don't denote nucleic or amino acids
    my $str = $self->seq();
    $str =~ s/[-.?]//gi;

    # Check for sequences without valid letters
    my $total = CORE::length($str);
    if( $total == 0 ) {
        if (!$self->{'_nowarnonempty'}) {
            $self->warn("Got a sequence with no letters in it ".
               "cannot guess alphabet");
        }
        return '';
    }

    if ($str =~ m/[EFIJLOPQXZ]/i) {
        # Start with a safe method to find proteins.
        # Unambiguous IUPAC letters for proteins are: E,F,I,J,L,O,P,Q,X,Z
        $type = 'protein';
    } else {
        # Alphabet is unsure, could still be DNA, RNA or protein.
        # DNA and RNA contain mostly A, T, U, G, C and N, but the other letters
        # they use are also among the 15 valid letters that a protein sequence
        # can contain at this stage. Make our best guess based on sequence
        # composition. If it contains over 70% of ACGTUN, it is likely nucleic.
        if( ($str =~ tr/ATUGCNatugcn//) / $total > 0.7 ) {
            if ( $str =~ m/U/i ) {
                $type = 'rna';
            } else {
                $type = 'dna';
            }
        } else {
            $type = 'protein';
        }
    }

    $self->alphabet($type);
    return $type;
}

############################################################################
# aliases due to name changes or to compensate for our lack of consistency #
############################################################################

sub accession {
    my $self = shift;

    $self->warn(ref($self)."::accession is deprecated, ".
        "use accession_number() instead");
    return $self->accession_number(@_);
}

1;

