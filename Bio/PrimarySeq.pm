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

Bio::PrimarySeq - Bioperl lightweight sequence object

=head1 SYNOPSIS

  # Bio::SeqIO for file reading, Bio::DB::GenBank for
  # database reading

  use Bio::Seq;
  use Bio::SeqIO;
  use Bio::DB::GenBank;

  # make from memory

  $seqobj = Bio::PrimarySeq->new (
      -seq              => 'ATGGGGTGGGCGGTGGGTGGTTTG',
      -id               => 'GeneFragment-12',
      -accession_number => 'X78121',
      -alphabet         => 'dna',
      -is_circular      => 1,
  );
  print "Sequence ", $seqobj->id(), " with accession ",
    $seqobj->accession_number, "\n";

  # read from file

  $inputstream = Bio::SeqIO->new(
      -file   => "myseq.fa",
      -format => 'Fasta',
  );
  $seqobj = $inputstream->next_seq();
  print "Sequence ", $seqobj->id(), " and desc ", $seqobj->desc, "\n";

  # to get out parts of the sequence.

  print "Sequence ", $seqobj->id(), " with accession ",
    $seqobj->accession_number, " and desc ", $seqobj->desc, "\n";

  $string  = $seqobj->seq();
  $string2 = $seqobj->subseq(1,40);

=head1 DESCRIPTION

PrimarySeq is a lightweight sequence object, storing the sequence, its
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

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut



package Bio::PrimarySeq;

use strict;

our $MATCHPATTERN = 'A-Za-z\-\.\*\?=~';
our $GAP_SYMBOLS  = '-~';

use base qw(Bio::Root::Root Bio::PrimarySeqI
            Bio::IdentifiableI Bio::DescribableI);


# Setup the allowed values for alphabet()
my %valid_type = map {$_, 1} qw( dna rna protein );


=head2 new

 Title   : new
 Usage   : $seqobj = Bio::PrimarySeq->new( -seq => 'ATGGGGGTGGTGGTACCCT',
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
 Args    : -seq              => sequence string
           -ref_to_seq       => ... or reference to a sequence string
           -display_id       => display id of the sequence (locus name)
           -accession_number => accession number
           -primary_id       => primary id (Genbank id)
           -version          => version number
           -namespace        => the namespace for the accession
           -authority        => the authority for the namespace
           -description      => description text
           -desc             => alias for description
           -alphabet         => skip alphabet guess and set it to dna, rna or protein
           -id               => alias for display id
           -is_circular      => boolean to indicate that sequence is circular
           -direct           => boolean to directly set sequences. The next time -seq,
                                seq() or -ref_to_seq is use, the sequence will not be
                                validated. Be careful with this...
           -nowarnonempty    => boolean to avoid warning when sequence is empty

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($seq, $id, $acc, $pid, $ns, $auth, $v, $oid, $desc, $description,
        $alphabet, $given_id, $is_circular, $direct, $ref_to_seq, $len,
        $nowarnonempty) =
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

    # Private var _nowarnonempty, needs to be set before calling _guess_alphabet
    $self->{'_nowarnonempty'} = $nowarnonempty;
    $self->{'_direct'} = $direct;

    if( defined $id && defined $given_id ) {
        if( $id ne $given_id ) {
          $self->throw("Provided both id and display_id constructors: [$id] [$given_id]");
        }
    }
    if( defined $given_id ) { $id = $given_id; }

    # Bernd's idea: set ids now for more informative invalid sequence messages
    defined $id  && $self->display_id($id);
    $acc         && $self->accession_number($acc);
    defined $pid && $self->primary_id($pid);

    # Set alphabet now to avoid guessing it later, when sequence is set
    $alphabet && $self->alphabet($alphabet);

    # Set the length before the seq. If there is a seq, length will be updated later
    $self->{'length'} = $len || 0;

    # Set the sequence (but also alphabet and length)
    if ($ref_to_seq) {
        $self->_set_seq_by_ref($ref_to_seq, $alphabet);
    } else {
        if (defined $seq) {
            # Note: the sequence string may be empty
            $self->seq($seq);
        }
    }

    $desc                  && $self->desc($desc);
    $description           && $self->description($description);
    $ns                    && $self->namespace($ns);
    $auth                  && $self->authority($auth);
    # Any variable that can have a value "0" must be tested with defined
    # or it will fail to be added to the new object
    defined($v)            && $self->version($v);
    defined($oid)          && $self->object_id($oid);
    defined($is_circular)  && $self->is_circular($is_circular);

    return $self;
}


=head2 seq

 Title   : seq
 Usage   : $string = $seqobj->seq();
 Function: Get or set  the sequence as a string of letters. The case of
           the letters is left up to the implementer. Suggested cases are
           upper case for proteins and lower case for DNA sequence (IUPAC
           standard), but you should not rely on this. An error is thrown if
           the sequence contains invalid characters: see validate_seq().
 Returns : A scalar
 Args    : - Optional new sequence value (a string) to set
           - Optional alphabet (it is guessed by default)

=cut

sub seq {
    my ($self, @args) = @_;

    if( scalar @args == 0 ) {
        return $self->{'seq'};
    }

    my ($seq_str, $alphabet) = @args;
    if (@args) {
        $self->_set_seq_by_ref(\$seq_str, $alphabet);
    }

    return $self->{'seq'};
}


sub _set_seq_by_ref {
    # Set a sequence by reference. A reference is used to avoid the cost of
    # copying the sequence (which can be very large) between functions.
    my ($self, $seq_str_ref, $alphabet) = @_;

    # Validate sequence if sequence is not empty and we are not in direct mode
    if ( (! $self->{'_direct'}) && (defined $$seq_str_ref) ) {
        $self->validate_seq($$seq_str_ref, 1);
    }
    delete $self->{'_direct'}; # next sequence will have to be validated

    # Record sequence length
    my $len = CORE::length($$seq_str_ref || '');
    my $is_changed_seq = (exists $self->{'seq'}) && ($len > 0);
    # Note: if the new seq is empty or undef, this is not considered a change
    delete $self->{'_freeze_length'} if $is_changed_seq;
    $self->{'length'} = $len if not exists $self->{'_freeze_length'};

    # Set sequence
    $self->{'seq'} = $$seq_str_ref;

    # Set or guess alphabet
    if ($alphabet) {
        # Alphabet specified, set it no matter what
        $self->alphabet($alphabet);
    } elsif ($is_changed_seq || (! defined($self->alphabet()))) {
        # If we changed a previous sequence to a new one or if there is no
        # alphabet yet at all, we need to guess the (possibly new) alphabet
        $self->_guess_alphabet();
    } # else (seq not changed and alphabet was defined) do nothing

    return 1;
}


=head2 validate_seq

 Title   : validate_seq
 Usage   : if(! $seqobj->validate_seq($seq_str) ) {
                print "sequence $seq_str is not valid for an object of
                alphabet ",$seqobj->alphabet, "\n";
           }
 Function: Test that the given sequence is valid, i.e. contains only valid
           characters. The allowed characters are all letters (A-Z) and '-','.',
           '*','?','=' and '~'. Spaces are not valid. Note that this
           implementation does not take alphabet() into account and that empty
           sequences are considered valid.
 Returns : 1 if the supplied sequence string is valid, 0 otherwise.
 Args    : - Sequence string to be validated
           - Boolean to optionally throw an error if the sequence is invalid

=cut

sub validate_seq {
    my ($self, $seqstr, $throw) = @_;
    if ( (defined $seqstr                ) &&
         ($seqstr !~ /^[$MATCHPATTERN]*$/) ) {
        if ($throw) {
            $self->throw("Failed validation of sequence '".(defined($self->id) ||
                '[unidentified sequence]')."'. Invalid characters were: " .
                join('',($seqstr =~ /[^$MATCHPATTERN]/g)));
        }
        return 0;
    }
    return 1;
}


=head2 subseq

 Title   : subseq
 Usage   : $substring = $seqobj->subseq(10,40);
           $substring = $seqobj->subseq(10,40,'nogap');
           $substring = $seqobj->subseq(-start=>10, -end=>40, -replace_with=>'tga');
           $substring = $seqobj->subseq($location_obj);
           $substring = $seqobj->subseq($location_obj, -nogap => 1);
 Function: Return the subseq from start to end, where the first sequence
           character has coordinate 1 number is inclusive, ie 1-2 are the
           first two characters of the sequence. The given start coordinate
           has to be larger than the end, even if the sequence is circular.
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
    my ($start, $end, $nogap, $replace) = $self->_rearrange([qw(START
                                                             END
                                                             NOGAP
                                                             REPLACE_WITH)], @args);

    # If -replace_with is specified, validate the replacement sequence
    if (defined $replace) {
        $self->validate_seq( $replace ) ||
            $self->throw("Replacement sequence does not look valid");
    }

    if( ref($start) && $start->isa('Bio::LocationI') ) {
        my $loc = $start;
        my $seq = '';

        # For Split objects if Guide Strand is negative,
        # pass the sublocations in reverse
        my $order = 0;
        if ($loc->isa('Bio::Location::SplitLocationI')) {
            # guide_strand can return undef, so don't compare directly
            # to avoid 'uninitialized value' warning
            my $guide_strand = defined ($loc->guide_strand) ? ($loc->guide_strand) : 0;
            $order = ($guide_strand == -1) ? -1 : 0;
        }
        # Reversing order using ->each_Location(-1) does not work well for
        # cut by origin-splits (like "complement(join(1900..END,START..50))"),
        # so use "reverse" instead
        my @sublocs = ($order == -1) ? reverse $loc->each_Location(): $loc->each_Location;
        foreach my $subloc (@sublocs) {
            my $piece = $self->subseq(-start        => $subloc->start(),
                                      -end          => $subloc->end(),
                                      -replace_with => $replace,
                                      -nogap        => $nogap);
            $piece =~ s/[$GAP_SYMBOLS]//g if $nogap;

            # strand can return undef, so don't compare directly
            # to avoid 'uninitialized value' warning
            my $strand = defined ($subloc->strand) ? ($subloc->strand) : 0;
            if ($strand < 0) {
                $piece = $self->_revcom_from_string($piece, $self->alphabet);
            }
            $seq .= $piece;
        }
        return $seq;
    } elsif( defined $start && defined $end ) {
        if( $start > $end ){
            $self->throw("Bad start,end parameters. Start [$start] has to be ".
              "less than end [$end]");
        }
        if( $start <= 0 ) {
            $self->throw("Bad start parameter ($start). Start must be positive.");
        }

        # Remove one from start, and then length is end-start
        $start--;

        my $seqstr;
        if (defined $replace) {
            $seqstr = substr $self->{seq}, $start, $end-$start, $replace;
        } else {
            $seqstr = substr $self->{seq}, $start, $end-$start;
        }


        if ($end > $self->length) {
            if ($self->is_circular) {
                my $start = 0;
                my $end = $end - $self->length;

                my $appendstr;
                if (defined $replace) {
                    $appendstr = substr $self->{seq}, $start, $end-$start, $replace;
                } else {
                    $appendstr = substr $self->{seq}, $start, $end-$start;
                }

                $seqstr .= $appendstr;
            } else {
                $self->throw("Bad end parameter ($end). End must be less than ".
                    "the total length of sequence (total=".$self->length.")")
            }
        }

        $seqstr =~ s/[$GAP_SYMBOLS]//g if ($nogap);
        return $seqstr;

    } else {
        $self->warn("Incorrect parameters to subseq - must be two integers or ".
            "a Bio::LocationI object. Got:", $self,$start,$end,$replace,$nogap);
        return;
    }
}


=head2 length

 Title   : length
 Usage   : $len = $seqobj->length();
 Function: Get the stored length of the sequence in number of symbols (bases
           or amino acids). In some circumstances, you can also set this attribute:

           1. For empty sequences, you can set the length to anything you want:
              my $seqobj = Bio::PrimarySeq->new( -length => 123 );
              my $len = $seqobj->len; # 123
           2. To save memory when using very long sequences, you can set the
              length of the sequence to the length of the sequence (and nothing
              else):
              my $seqobj = Bio::PrimarySeq->new( -seq => 'ACGT...' ); # 1 Mbp sequence
              # process $seqobj... then after you're done with it
              $seqobj->length($seqobj->length);
              $seqobj->seq(undef); # free memory!
              my $len = $seqobj->len; # 1 Mbp

           Note that if you set seq() to a value other than undef at any time,
           the length attribute will be reset.
 Returns : integer representing the length of the sequence.
 Args    : Optionally, the value on set

=cut

sub length {
    my ($self, $val) = @_;
    if (defined $val) {
        my $len = $self->{'length'};
        if ($len && ($len != $val)) {
            $self->throw("Can not set the length to $val, current length value is $len");
        }
        $self->{'length'} = $val;
        $self->{'_freeze_length'} = undef;
    }
    return $self->{'length'};
}


=head2 display_id

 Title   : display_id or display_name
 Usage   : $id_string = $seqobj->display_id();
 Function: Get or set the display id, aka the common name of the sequence object.

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
 Returns : A string for the display ID
 Args    : Optional string for the display ID to set

=cut

sub display_id {
    my ($self, $value) = @_;
    if( defined $value) {
        $self->{'display_id'} = $value;
    }
    return $self->{'display_id'};
}


=head2 accession_number

 Title   : accession_number or object_id
 Usage   : $unique_key = $seqobj->accession_number;
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
    my( $self, $acc ) = @_;
    if (defined $acc) {
        $self->{'accession_number'} = $acc;
    } else {
        $acc = $self->{'accession_number'};
        $acc = 'unknown' unless defined $acc;
    }
    return $acc;
}


=head2 primary_id

 Title   : primary_id
 Usage   : $unique_key = $seqobj->primary_id;
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
    my $self = shift;

    if(@_) {
        $self->{'primary_id'} = shift;
    }
    if( ! defined($self->{'primary_id'}) ) {
        return "$self";
    }
    return $self->{'primary_id'};
}


=head2 alphabet

 Title   : alphabet
 Usage   : if( $seqobj->alphabet eq 'dna' ) { # Do something }
 Function: Get/set the alphabet of sequence, one of
           'dna', 'rna' or 'protein'. This is case sensitive.

           This is not called <type> because this would cause
           upgrade problems from the 0.5 and earlier Seq objects.
 Returns : a string either 'dna','rna','protein'. NB - the object must
           make a call of the type - if there is no alphabet specified it
           has to guess.
 Args    : optional string to set : 'dna' | 'rna' | 'protein'


=cut

sub alphabet {
    my ($self,$value) = @_;
    if (defined $value) {
        $value = lc $value;
        unless ( $valid_type{$value} ) {
        $self->throw("Alphabet '$value' is not a valid alphabet (".
            join(',', map "'$_'", sort keys %valid_type) .") lowercase");
        }
        $self->{'alphabet'} = $value;
    }
    return $self->{'alphabet'};
}


=head2 desc

 Title   : desc or description
 Usage   : $seqobj->desc($newval);
 Function: Get/set description of the sequence.

           'description' is an alias for this for compliance with the
           Bio::DescribeableI interface.
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
 Usage   : $id = $seqobj->id();
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
 Usage   : if( $seqobj->is_circular) { # Do something }
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

=head2 object_id

 Title   : object_id
 Usage   : $string = $seqobj->object_id();
 Function: Get or set a string which represents the stable primary identifier
           in this namespace of this object. For DNA sequences this
           is its accession_number, similarly for protein sequences.

           This is aliased to accession_number().
 Returns : A scalar
 Args    : Optional object ID to set.

=cut

sub object_id {
    return shift->accession_number(@_);
}


=head2 version

 Title   : version
 Usage   : $version = $seqobj->version();
 Function: Get or set a number which differentiates between versions of
           the same object. Higher numbers are considered to be
           later and more relevant, but a single object described
           the same identifier should represent the same concept.
 Returns : A number
 Args    : Optional version to set.

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
 Usage   : $authority = $seqobj->authority();
 Function: Get or set a string which represents the organisation which
           granted the namespace, written as the DNS name of the
           organisation (eg, wormbase.org).
 Returns : A scalar
 Args    : Optional authority to set.

=cut

sub authority {
    my ($self, $value) = @_;
    if( defined $value) {
        $self->{'authority'} = $value;
    }
    return $self->{'authority'};
}


=head2 namespace

 Title   : namespace
 Usage   : $string = $seqobj->namespace();
 Function: Get or set a string representing the name space this identifier
           is valid in, often the database name or the name describing the
           collection.
 Returns : A scalar
 Args    : Optional namespace to set.

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

=head2 display_name

 Title   : display_name
 Usage   : $string = $seqobj->display_name();
 Function: Get or set a string which is what should be displayed to the user.
           The string should have no spaces (ideally, though a cautious
           user of this interface would not assumme this) and should be
           less than thirty characters (though again, double checking
           this is a good idea).

           This is aliased to display_id().
 Returns : A string for the display name
 Args    : Optional string for the display name to set.

=cut

sub display_name {
    return shift->display_id(@_);
}


=head2 description

 Title   : description
 Usage   : $string = $seqobj->description();
 Function: Get or set a text string suitable for displaying to the user a
           description. This string is likely to have spaces, but
           should not have any newlines or formatting - just plain
           text. The string should not be greater than 255 characters
           and clients can feel justified at truncating strings at 255
           characters for the purposes of display.

           This is aliased to desc().
 Returns : A string for the description
 Args    : Optional string for the description to set.

=cut

sub description {
    return shift->desc(@_);
}


=head1 Methods Inherited from Bio::PrimarySeqI

These methods are available on Bio::PrimarySeq, although they are
actually implemented on Bio::PrimarySeqI

=head2 revcom

 Title   : revcom
 Usage   : $rev = $seqobj->revcom();
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

=head2 trunc

 Title   : trunc
 Usage   : $subseq = $myseq->trunc(10,100);
 Function: Provides a truncation of a sequence,
 Returns : A fresh Bio::SeqI implementing object.
 Args    : Numbers for the start and end positions

=head1 Internal methods

These are internal methods to PrimarySeq

=head2 _guess_alphabet

 Title   : _guess_alphabet
 Usage   :
 Function: Automatically guess and set the type of sequence: dna, rna, protein
           or '' if the sequence was empty. This method first removes dots (.),
           dashes (-) and question marks (?) before guessing the alphabet
           using the IUPAC conventions for ambiguous residues. Since the DNA and
           RNA characters are also valid characters for proteins, there is
           no foolproof way of determining the right alphabet. This is our best
           guess only!
 Returns : string 'dna', 'rna', 'protein' or ''.
 Args    : none

=cut

sub _guess_alphabet {
    my ($self) = @_;
    # Guess alphabet
    my $alphabet = $self->_guess_alphabet_from_string($self->seq, $self->{'_nowarnonempty'});
    # Set alphabet unless it is unknown
    $self->alphabet($alphabet) if $alphabet;
    return $alphabet;
}


sub _guess_alphabet_from_string {
    # Get the alphabet from a sequence string
    my ($self, $str, $nowarnonempty) = @_;

    $nowarnonempty = 0 if not defined $nowarnonempty;

    # Remove chars that clearly don't denote nucleic or amino acids
    $str =~ s/[-.?]//gi;

    # Check for sequences without valid letters
    my $alphabet;
    my $total = CORE::length($str);
    if( $total == 0 ) {
        if (not $nowarnonempty) {
            $self->warn("Got a sequence without letters. Could not guess alphabet");
        }
        $alphabet = '';
    }

    # Determine alphabet now
    if (not defined $alphabet) {
        if ($str =~ m/[EFIJLOPQXZ]/i) {
            # Start with a safe method to find proteins.
            # Unambiguous IUPAC letters for proteins are: E,F,I,J,L,O,P,Q,X,Z
            $alphabet = 'protein';
        } else {
            # Alphabet is unsure, could still be DNA, RNA or protein
            # DNA and RNA contain mostly A, T, U, G, C and N, but the other
            # letters they use are also among the 15 valid letters that a
            # protein sequence can contain at this stage. Make our best guess
            # based on sequence composition. If it contains over 70% of ACGTUN,
            # it is likely nucleic.
            if( ($str =~ tr/ATUGCNWSKMatugcnwskm//) / $total > 0.7 ) {
                if ( $str =~ m/U/i ) {
                    $alphabet = 'rna';
                } else {
                    $alphabet = 'dna';
                }
            } else {
                $alphabet = 'protein';
            }
        }
    }

    return $alphabet;
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
