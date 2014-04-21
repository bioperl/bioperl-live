#
# BioPerl module for Bio::Seq
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

Bio::Seq - Sequence object, with features

=head1 SYNOPSIS

    # This is the main sequence object in Bioperl

    # gets a sequence from a file
    $seqio  = Bio::SeqIO->new( '-format' => 'embl' , -file => 'myfile.dat');
    $seqobj = $seqio->next_seq();

    # SeqIO can both read and write sequences; see Bio::SeqIO
    # for more information and examples

    # get from database
    $db = Bio::DB::GenBank->new();
    $seqobj = $db->get_Seq_by_acc('X78121');

    # make from strings in script
    $seqobj = Bio::Seq->new( -display_id => 'my_id',
                             -seq => $sequence_as_string);

    # gets sequence as a string from sequence object
    $seqstr   = $seqobj->seq(); # actual sequence as a string
    $seqstr   = $seqobj->subseq(10,50); # slice in biological coordinates

    # retrieves information from the sequence
    # features must implement Bio::SeqFeatureI interface

    @features = $seqobj->get_SeqFeatures(); # just top level
    foreach my $feat ( @features ) {
        print "Feature ",$feat->primary_tag," starts ",$feat->start," ends ",
            $feat->end," strand ",$feat->strand,"\n";

        # features retain link to underlying sequence object
        print "Feature sequence is ",$feat->seq->seq(),"\n"
    }

    # sequences may have a species

    if( defined $seq->species ) {
        print "Sequence is from ",$species->binomial," [",$species->common_name,"]\n";
    }

    # annotation objects are Bio::AnnotationCollectionI's
    $ann      = $seqobj->annotation(); # annotation object

    # references is one type of annotations to get. Also get
    # comment and dblink. Look at Bio::AnnotationCollection for
    # more information

    foreach my $ref ( $ann->get_Annotations('reference') ) {
        print "Reference ",$ref->title,"\n";
    }

    # you can get truncations, translations and reverse complements, these
    # all give back Bio::Seq objects themselves, though currently with no
    # features transfered

    my $trunc = $seqobj->trunc(100,200);
    my $rev   = $seqobj->revcom();

    # there are many options to translate - check out the docs
    my $trans = $seqobj->translate();

    # these functions can be chained together

    my $trans_trunc_rev = $seqobj->trunc(100,200)->revcom->translate();



=head1 DESCRIPTION

A Seq object is a sequence with sequence features placed on it. The
Seq object contains a PrimarySeq object for the actual sequence and
also implements its interface.

In Bioperl we have 3 main players that people are going to use frequently

  Bio::PrimarySeq  - just the sequence and its names, nothing else.
  Bio::SeqFeatureI - a feature on a sequence, potentially with a sequence
                     and a location and annotation.
  Bio::Seq         - A sequence and a collection of sequence features
                     (an aggregate) with its own annotation.

Although Bioperl is not tied heavily to file formats these distinctions do
map to file formats sensibly and for some bioinformaticians this might help

  Bio::PrimarySeq  - Fasta file of a sequence
  Bio::SeqFeatureI - A single entry in an EMBL/GenBank/DDBJ feature table
  Bio::Seq         - A single EMBL/GenBank/DDBJ entry

By having this split we avoid a lot of nasty circular references
(sequence features can hold a reference to a sequence without the sequence
holding a reference to the sequence feature). See L<Bio::PrimarySeq> and
L<Bio::SeqFeatureI> for more information.

Ian Korf really helped in the design of the Seq and SeqFeature system.

=head2 Examples

A simple and fundamental block of code:

  use Bio::SeqIO;

  my $seqIOobj = Bio::SeqIO->new(-file=>"1.fa"); # create a SeqIO object
  my $seqobj = $seqIOobj->next_seq;              # get a Seq object

With the Seq object in hand one has access to a powerful set of Bioperl
methods and related Bioperl objects. This next script will take a file of sequences
in EMBL format and create a file of the reverse-complemented sequences
in Fasta format using Seq objects. It also prints out details about the
exons it finds as sequence features in Genbank Flat File format.

  use Bio::Seq;
  use Bio::SeqIO;

  $seqin = Bio::SeqIO->new( -format => 'EMBL' , -file => 'myfile.dat');
  $seqout= Bio::SeqIO->new( -format => 'Fasta', -file => '>output.fa');

  while((my $seqobj = $seqin->next_seq())) {
      print "Seen sequence ",$seqobj->display_id,", start of seq ",
            substr($seqobj->seq,1,10),"\n";
      if( $seqobj->alphabet eq 'dna') {
            $rev = $seqobj->revcom;
            $id  = $seqobj->display_id();
            $id  = "$id.rev";
            $rev->display_id($id);
            $seqout->write_seq($rev);
      }

      foreach $feat ( $seqobj->get_SeqFeatures() ) {
           if( $feat->primary_tag eq 'exon' ) {
              print STDOUT "Location ",$feat->start,":",
                    $feat->end," GFF[",$feat->gff_string,"]\n";
           }
      }
  }

Let's examine the script. The lines below import the Bioperl modules.
Seq is the main Bioperl sequence object and SeqIO is the Bioperl support
for reading sequences from files and to files

  use Bio::Seq;
  use Bio::SeqIO;

These two lines create two SeqIO streams: one for reading in sequences
and one for outputting sequences:

  $seqin = Bio::SeqIO->new( -format => 'EMBL' , -file => 'myfile.dat');
  $seqout= Bio::SeqIO->new( -format => 'Fasta', -file => '>output.fa');

Notice that in the "$seqout" case there is a greater-than sign,
indicating the file is being opened for writing.

Using the

  '-argument' => value

syntax is common in Bioperl. The file argument is like an argument
to open() . You can also pass in filehandles or FileHandle objects by
using the -fh argument (see L<Bio::SeqIO> documentation for details).
Many formats in Bioperl are handled, including Fasta, EMBL, GenBank,
Swissprot (swiss), PIR, and GCG.

  $seqin = Bio::SeqIO->new( -format => 'EMBL' , -file => 'myfile.dat');
  $seqout= Bio::SeqIO->new( -format => 'Fasta', -file => '>output.fa');

This is the main loop which will loop progressively through sequences
in a file, and each call to $seqio-E<gt>next_seq() provides a new Seq
object from the file:

  while((my $seqobj = $seqio->next_seq())) {

This print line below accesses fields in the Seq object directly. The
$seqobj-E<gt>display_id is the way to access the display_id attribute
of the Seq object. The $seqobj-E<gt>seq method gets the actual
sequence out as string. Then you can do manipulation of this if
you want to (there are however easy ways of doing truncation,
reverse-complement and translation).

  print "Seen sequence ",$seqobj->display_id,", start of seq ",
               substr($seqobj->seq,1,10),"\n";

Bioperl has to guess the alphabet of the sequence, being either 'dna',
'rna', or 'protein'. The alphabet attribute is one of these three
possibilities.

  if( $seqobj->alphabet eq 'dna') {

The $seqobj-E<gt>revcom method provides the reverse complement of the Seq
object as another Seq object. Thus, the $rev variable is a reference to
another Seq object. For example, one could repeat the above print line
for this Seq object (putting $rev in place of $seqobj). In this
case we are going to output the object into the file stream we built
earlier on.

  $rev = $seqobj->revcom;

When we output it, we want the id of the outputted object
to be changed to "$id.rev", ie, with .rev on the end of the name. The
following lines retrieve the id of the sequence object, add .rev
to this and then set the display_id of the rev sequence object to
this. Notice that to set the display_id attribute you just need
call the same method, display_id(), with the new value as an argument.
Getting and setting values with the same method is common in Bioperl.

  $id  = $seqobj->display_id();
  $id  = "$id.rev";
  $rev->display_id($id);

The write_seq method on the SeqIO output object, $seqout, writes the
$rev object to the filestream we built at the top of the script.
The filestream knows that it is outputting in fasta format, and
so it provides fasta output.

  $seqout->write_seq($rev);

This block of code loops over sequence features in the sequence
object, trying to find ones who have been tagged as 'exon'.
Features have start and end attributes and can be outputted
in Genbank Flat File format, GFF, a standarized format for sequence
features.

  foreach $feat ( $seqobj->get_SeqFeatures() ) {
      if( $feat->primary_tag eq 'exon' ) {
          print STDOUT "Location ",$feat->start,":",
             $feat->end," GFF[",$feat->gff_string,"]\n";
      }
  }

The code above shows how a few Bio::Seq methods suffice to read, parse,
reformat and analyze sequences from a file. A full list of methods
available to Bio::Seq objects is shown below. Bear in mind that some of
these methods come from PrimarySeq objects, which are simpler
than Seq objects, stripped of features (see L<Bio::PrimarySeq> for
more information).

  # these methods return strings, and accept strings in some cases:

  $seqobj->seq();              # string of sequence
  $seqobj->subseq(5,10);       # part of the sequence as a string
  $seqobj->accession_number(); # when there, the accession number
  $seqobj->alphabet();         # one of 'dna','rna',or 'protein'
  $seqobj->version()           # when there, the version
  $seqobj->keywords();         # when there, the Keywords line
  $seqobj->length()            # length
  $seqobj->desc();             # description
  $seqobj->primary_id();       # a unique id for this sequence regardless
                               # of its display_id or accession number
  $seqobj->display_id();       # the human readable id of the sequence

Some of these values map to fields in common formats. For example, The
display_id() method returns the LOCUS name of a Genbank entry,
the (\S+) following the E<gt> character in a Fasta file, the ID from
a SwissProt file, and so on. The desc() method will return the DEFINITION
line of a Genbank file, the description following the display_id in a
Fasta file, and the DE field in a SwissProt file.

  # the following methods return new Seq objects, but
  # do not transfer features across to the new object:

  $seqobj->trunc(5,10)  # truncation from 5 to 10 as new object
  $seqobj->revcom       # reverse complements sequence
  $seqobj->translate    # translation of the sequence

  # if new() can be called this method returns 1, else 0

  $seqobj->can_call_new

  # the following method determines if the given string will be accepted
  # by the seq() method - if the string is acceptable then validate()
  # returns 1, or 0 if not

  $seqobj->validate_seq($string)

  # the following method returns or accepts a Species object:

  $seqobj->species();

Please see L<Bio::Species> for more information on this object.

  # the following method returns or accepts an Annotation object
  # which in turn allows access to Annotation::Reference
  # and Annotation::Comment objects:

  $seqobj->annotation();

These annotations typically refer to entire sequences, unlike
features.  See L<Bio::AnnotationCollectionI>,
L<Bio::Annotation::Collection>, L<Bio::Annotation::Reference>, and
L<Bio::Annotation::Comment> for details.

It is also important to be able to describe defined portions of a
sequence. The combination of some description and the corresponding
sub-sequence is called a feature - an exon and its coordinates within
a gene is an example of a feature, or a domain within a protein.

  # the following methods return an array of SeqFeatureI objects:

  $seqobj->get_SeqFeatures # The 'top level' sequence features
  $seqobj->get_all_SeqFeatures # All sequence features, including sub-seq
                               # features, such as features in an exon

  # to find out the number of features use:

  $seqobj->feature_count

Here are just some of the methods available to SeqFeatureI objects:

  # these methods return numbers:

  $feat->start          # start position (1 is the first base)
  $feat->end            # end position (2 is the second base)
  $feat->strand         # 1 means forward, -1 reverse, 0 not relevant

  # these methods return or accept strings:

  $feat->primary_tag    # the name of the sequence feature, eg
                        # 'exon', 'glycoslyation site', 'TM domain'
  $feat->source_tag     # where the feature comes from, eg, 'EMBL_GenBank',
                        # or 'BLAST'

  # this method returns the more austere PrimarySeq object, not a
  # Seq object - the main difference is that PrimarySeq objects do not
  # themselves contain sequence features

  $feat->seq            # the sequence between start,end on the
                        # correct strand of the sequence

See L<Bio::PrimarySeq> for more details on PrimarySeq objects.

  # useful methods for feature comparisons, for start/end points

  $feat->overlaps($other)  # do $feat and $other overlap?
  $feat->contains($other)  # is $other completely within $feat?
  $feat->equals($other)    # do $feat and $other completely agree?

  # one can also add features

  $seqobj->add_SeqFeature($feat)     # returns 1 if successful

  # sub features. For complex join() statements, the feature
  # is one sequence feature with many sub SeqFeatures

  $feat->sub_SeqFeature  # returns array of sub seq features

Please see L<Bio::SeqFeatureI> and L<Bio::SeqFeature::Generic>,
for more information on sequence features.

It is worth mentioning that one can also retrieve the start and end
positions of a feature using a Bio::LocationI object:

  $location = $feat->location # $location is a Bio::LocationI object
  $location->start;           # start position
  $location->end;             # end position

This is useful because one needs a Bio::Location::SplitLocationI object
in order to retrieve the coordinates inside the Genbank or EMBL join()
statements (e.g. "CDS    join(51..142,273..495,1346..1474)"):

  if ( $feat->location->isa('Bio::Location::SplitLocationI') &&
               $feat->primary_tag eq 'CDS' )  {
    foreach $loc ( $feat->location->sub_Location ) {
      print $loc->start . ".." . $loc->end . "\n";
    }
  }

See L<Bio::LocationI> and L<Bio::Location::SplitLocationI> for more
information.

=head1 Implemented Interfaces

This class implements the following interfaces.

=over 4

=item Bio::SeqI

Note that this includes implementing Bio::PrimarySeqI.

=item Bio::IdentifiableI

=item Bio::DescribableI

=item Bio::AnnotatableI

=item Bio::FeatureHolderI

=back

=head1 FEEDBACK


=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

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

=head1 AUTHOR - Ewan Birney, inspired by Ian Korf objects

Email birney@ebi.ac.uk

=head1 CONTRIBUTORS

 Jason Stajich E<lt>jason@bioperl.orgE<gt>
 Mark A. Jensen maj -at- fortinbras -dot- us

=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut

#'
# Let the code begin...


package Bio::Seq;
use strict;

use Bio::Annotation::Collection;
use Bio::PrimarySeq;

use base qw(Bio::Root::Root Bio::SeqI Bio::IdentifiableI Bio::DescribableI Bio::AnnotatableI Bio::FeatureHolderI Bio::AnnotationCollectionI);

=head2 new

 Title   : new
 Usage   : $seq = Bio::Seq->new( -seq => 'ATGGGGGTGGTGGTACCCT',
                                 -id  => 'human_id',
                                 -accession_number => 'AL000012',
                               );

 Function: Returns a new Seq object from
           basic constructors, being a string for the sequence
           and strings for id and accession_number
 Returns : a new Bio::Seq object

=cut

sub new {
    my($caller,@args) = @_;

    if( $caller ne 'Bio::Seq') {
        $caller = ref($caller) if ref($caller);
    }

    # we know our inherietance hierarchy
    my $self = Bio::Root::Root->new(@args);
    bless $self,$caller;

    # this is way too sneaky probably. We delegate the construction of
    # the Seq object onto PrimarySeq and then pop primary_seq into
    # our primary_seq slot

    my $pseq = Bio::PrimarySeq->new(@args);

    # as we have just made this, we know it is ok to set hash directly
    # rather than going through the method

    $self->{'primary_seq'} = $pseq;

    # setting this array is now delayed until the final
    # moment, again speed ups for non feature containing things
    # $self->{'_as_feat'} = [];


    my ($ann, $pid,$feat,$species) = &Bio::Root::RootI::_rearrange($self,[qw(ANNOTATION PRIMARY_ID FEATURES SPECIES)], @args);

    # for a number of cases - reading fasta files - these are never set. This
    # gives a quick optimisation around testing things later on

    if( defined $ann || defined $pid || defined $feat || defined $species ) {
        $pid && $self->primary_id($pid);
        $species && $self->species($species);
        $ann && $self->annotation($ann);
        
        if( defined $feat ) {
            if( ref($feat) !~ /ARRAY/i ) {
                if( ref($feat) && $feat->isa('Bio::SeqFeatureI') ) {
                    $self->add_SeqFeature($feat);
                } else {
                    $self->warn("Must specify a valid Bio::SeqFeatureI or ArrayRef of Bio::SeqFeatureI's with the -features init parameter for ".ref($self));
                }
            } else {
                foreach my $feature ( @$feat ) {
                    $self->add_SeqFeature($feature);
                }
            }
        }
    }

    return $self;
}


=head1 PrimarySeq interface


The PrimarySeq interface provides the basic sequence getting
and setting methods for on all sequences.

These methods implement the Bio::PrimarySeq interface by delegating
to the primary_seq inside the object. This means that you
can use a Seq object wherever there is a PrimarySeq, and
of course, you are free to use these functions anyway.

=cut

=head2 seq

 Title   : seq
 Usage   : $string = $obj->seq()
 Function: Get/Set the sequence as a string of letters. The
           case of the letters is left up to the implementer.
           Suggested cases are upper case for proteins and lower case for
           DNA sequence (IUPAC standard),
           but implementations are suggested to keep an open mind about
           case (some users... want mixed case!)
 Returns : A scalar
 Args    : Optionally on set the new value (a string). An optional second
           argument presets the alphabet (otherwise it will be guessed).
           Both parameters may also be given in named parameter style
           with -seq and -alphabet being the names.

=cut

sub seq {
    return shift->primary_seq()->seq(@_);
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
           implementation does not take alphabet() into account.
 Returns : 1 if the supplied sequence string is valid, 0 otherwise.
 Args    : - Sequence string to be validated
           - Boolean to throw an error if the sequence is invalid

=cut

sub validate_seq {
    return shift->primary_seq()->validate_seq(@_);
}


=head2 length

 Title   : length
 Usage   : $len = $seq->length()
 Function:
 Example :
 Returns : Integer representing the length of the sequence.
 Args    : None

=cut

sub length {
    return shift->primary_seq()->length(@_);
}


=head1 Methods from the Bio::PrimarySeqI interface

=head2 subseq

 Title   : subseq
 Usage   : $substring = $obj->subseq(10,40);
 Function: Returns the subseq from start to end, where the first base
           is 1 and the number is inclusive, ie 1-2 are the first two
           bases of the sequence

           Start cannot be larger than end but can be equal

 Returns : A string
 Args    : 2 integers


=cut

sub subseq {
    return shift->primary_seq()->subseq(@_);
}


=head2 display_id

 Title   : display_id
 Usage   : $id = $obj->display_id or $obj->display_id($newid);
 Function: Gets or sets the display id, also known as the common name of
           the Seq object.

           The semantics of this is that it is the most likely string
           to be used as an identifier of the sequence, and likely to
           have "human" readability.  The id is equivalent to the LOCUS
           field of the GenBank/EMBL databanks and the ID field of the
           Swissprot/sptrembl database. In fasta format, the >(\S+) is
           presumed to be the id, though some people overload the id
           to embed other information. Bioperl does not use any
           embedded information in the ID field, and people are
           encouraged to use other mechanisms (accession field for
           example, or extending the sequence object) to solve this.

           Notice that $seq->id() maps to this function, mainly for
           legacy/convenience issues.
 Returns : A string
 Args    : None or a new id

=cut

sub display_id {
   return shift->primary_seq->display_id(@_);
}


=head2 accession_number

 Title   : accession_number
 Usage   : $unique_biological_key = $obj->accession_number;
 Function: Returns the unique biological id for a sequence, commonly
           called the accession_number. For sequences from established
           databases, the implementors should try to use the correct
           accession number. Notice that primary_id() provides the
           unique id for the implemetation, allowing multiple objects
           to have the same accession number in a particular implementation.

           For sequences with no accession number, this method should return
           "unknown".

           Can also be used to set the accession number.
 Example : $key = $seq->accession_number or $seq->accession_number($key)
 Returns : A string
 Args    : None or an accession number

=cut

sub accession_number {
   return shift->primary_seq->accession_number(@_);
}


=head2 desc

 Title   : desc
 Usage   : $seqobj->desc($string) or $seqobj->desc()
 Function: Sets or gets the description of the sequence
 Example :
 Returns : The description
 Args    : The description or none

=cut

sub desc {
   return shift->primary_seq->desc(@_);
}


=head2 primary_id

 Title   : primary_id
 Usage   : $unique_implementation_key = $obj->primary_id;
 Function: Returns the unique id for this object in this
           implementation. This allows implementations to manage
           their own object ids in a way the implementation can control
           clients can expect one id to map to one object.

           For sequences with no natural id, this method should return
           a stringified memory location.

           Can also be used to set the primary_id (or unset to undef).

           [Note this method name is likely to change in 1.3]

 Example : $id = $seq->primary_id or $seq->primary_id($id)
 Returns : A string
 Args    : None or an id, or undef to unset the primary id.

=cut

sub primary_id {
    # Note: this used to not delegate to the primary seq. This is
    # really bad in very subtle ways. E.g., if you created the object
    # with a primary id given to the constructor and then later you
    # change the primary id, if this method wouldn't delegate you'd
    # have different values for primary id in the PrimarySeq object
    # compared to this instance. Not good.

    # I can't remember why not delegating was ever deemed
    # advantageous, but I hereby claim that its problems far outweigh
    # its advantages, if there are any. Convince me otherwise if you
    # disagree. HL 2004/08/05

    return shift->primary_seq->primary_id(@_);
}


=head2 can_call_new

 Title   : can_call_new
 Usage   : if ( $obj->can_call_new ) {
               $newobj = $obj->new( %param );
           }
 Function: can_call_new returns 1 or 0 depending
           on whether an implementation allows new
           constructor to be called. If a new constructor
           is allowed, then it should take the followed hashed
           constructor list.

           $myobject->new( -seq => $sequence_as_string,
                           -display_id  => $id
                           -accession_number => $accession
                           -alphabet => 'dna',
                         );
 Example :
 Returns : 1 or 0
 Args    : None

=cut

sub can_call_new {
    return 1;
}


=head2 alphabet

 Title   : alphabet
 Usage   : if ( $obj->alphabet eq 'dna' ) { /Do Something/ }
 Function: Get/Set the type of sequence being one of
           'dna', 'rna' or 'protein'. This is case sensitive.

           This is not called <type> because this would cause
           upgrade problems from the 0.5 and earlier Seq objects.

 Returns : A string either 'dna','rna','protein'. NB - the object must
           make a call of the type - if there is no type specified it
           has to guess.
 Args    : optional string to set : 'dna' | 'rna' | 'protein'

=cut

sub alphabet {
   my $self = shift;
   return $self->primary_seq->alphabet(@_) if @_ && defined $_[0];
   return $self->primary_seq->alphabet();
}


=head2 is_circular

 Title   : is_circular
 Usage   : if( $obj->is_circular) { /Do Something/ }
 Function: Returns true if the molecule is circular
 Returns : Boolean value
 Args    : none

=cut

sub is_circular {
    return shift->primary_seq()->is_circular(@_);
}


=head1 Methods for Bio::IdentifiableI compliance

=head2 object_id

 Title   : object_id
 Usage   : $string    = $obj->object_id()
 Function: a string which represents the stable primary identifier
           in this namespace of this object. For DNA sequences this
           is its accession_number, similarly for protein sequences

           This is aliased to accession_number().
 Returns : A scalar

=cut

sub object_id {
    return shift->accession_number(@_);
}


=head2 version

 Title   : version
 Usage   : $version    = $obj->version()
 Function: a number which differentiates between versions of
           the same object. Higher numbers are considered to be
           later and more relevant, but a single object described
           the same identifier should represent the same concept

 Returns : A number

=cut

sub version{
    return shift->primary_seq->version(@_);
}


=head2 authority

 Title   : authority
 Usage   : $authority    = $obj->authority()
 Function: a string which represents the organisation which
           granted the namespace, written as the DNS name for
           organisation (eg, wormbase.org)

 Returns : A scalar

=cut

sub authority {
    return shift->primary_seq()->authority(@_);
}


=head2 namespace

 Title   : namespace
 Usage   : $string    = $obj->namespace()
 Function: A string representing the name space this identifier
           is valid in, often the database name or the name
           describing the collection

 Returns : A scalar

=cut

sub namespace{
    return shift->primary_seq()->namespace(@_);
}


=head1 Methods for Bio::DescribableI compliance

=head2 display_name

 Title   : display_name
 Usage   : $string    = $obj->display_name()
 Function: A string which is what should be displayed to the user
           the string should have no spaces (ideally, though a cautious
           user of this interface would not assumme this) and should be
           less than thirty characters (though again, double checking
           this is a good idea)

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
           characters for the purposes of display

           This is aliased to desc().
 Returns : A scalar

=cut

sub description {
    return shift->desc(@_);
}


=head1 Methods for implementing Bio::AnnotatableI

=head2 annotation

 Title   : annotation
 Usage   : $ann = $seq->annotation or 
           $seq->annotation($ann)
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


=head1 Methods for delegating Bio::AnnotationCollectionI

=head2 get_Annotations()

 Usage   : my @annotations = $seq->get_Annotations('key')
 Function: Retrieves all the Bio::AnnotationI objects for a specific key
           for this object
 Returns : list of Bio::AnnotationI - empty if no objects stored for a key
 Args    : string which is key for annotations

=cut

sub get_Annotations { shift->annotation->get_Annotations(@_); }


=head2 add_Annotation()

 Usage   : $seq->add_Annotation('reference',$object);
           $seq->add_Annotation($object,'Bio::MyInterface::DiseaseI');
           $seq->add_Annotation($object);
           $seq->add_Annotation('disease',$object,'Bio::MyInterface::DiseaseI');
 Function: Adds an annotation for a specific key for this sequence object.

           If the key is omitted, the object to be added must provide a value
           via its tagname().

           If the archetype is provided, this and future objects added under
           that tag have to comply with the archetype and will be rejected
           otherwise.

 Returns : none
 Args    : annotation key ('disease', 'dblink', ...)
           object to store (must be Bio::AnnotationI compliant)
           [optional] object archetype to map future storage of object
           of these types to

=cut

sub add_Annotation { shift->annotation->add_Annotation(@_) }


=head2 remove_Annotations()

 Usage   : $seq->remove_Annotations()
 Function: Remove the annotations for the specified key from this sequence 
           object
 Returns : an list of Bio::AnnotationI compliant objects which were stored
           under the given key(s) for this sequence object
 Args    : the key(s) (tag name(s), one or more strings) for which to
           remove annotations (optional; if none given, flushes all
           annotations)

=cut

sub remove_Annotations { shift->annotation->remove_Annotations(@_) }


=head2 get_num_of_annotations()

 Usage   : my $count = $seq->get_num_of_annotations()
 Alias   : num_Annotations
 Function: Returns the count of all annotations stored for this sequence
           object      
 Returns : integer
 Args    : none

=cut

sub get_num_of_annotations { shift->annotation->get_num_of_annotations(@_) }
sub num_Annotations { shift->get_num_of_annotations }; #DWYM


=head1 Methods to implement Bio::FeatureHolderI

This includes methods for retrieving, adding, and removing features.

=cut

=head2 get_SeqFeatures

 Title   : get_SeqFeatures
 Usage   :
 Function: Get the feature objects held by this feature holder.

           Features which are not top-level are subfeatures of one or
           more of the returned feature objects, which means that you
           must traverse the subfeature arrays of each top-level
           feature object in order to traverse all features associated
           with this sequence.

           Top-level features can be obtained by tag, specified in 
           the argument.

           Use get_all_SeqFeatures() if you want the feature tree
           flattened into one single array.

 Example :
 Returns : an array of Bio::SeqFeatureI implementing objects
 Args    : [optional] scalar string (feature tag)

=cut

sub get_SeqFeatures{
   my $self = shift;
   my $tag = shift;

   if( !defined $self->{'_as_feat'} ) {
       $self->{'_as_feat'} = [];
   }
   if ($tag) {
       return map { $_->primary_tag eq $tag ? $_ : () } @{$self->{'_as_feat'}};
   }
   else {
       return @{$self->{'_as_feat'}};
   }
}


=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : @feat_ary = $seq->get_all_SeqFeatures();
 Function: Returns the tree of feature objects attached to this
           sequence object flattened into one single array. Top-level
           features will still contain their subfeature-arrays, which
           means that you will encounter subfeatures twice if you
           traverse the subfeature tree of the returned objects.

           Use get_SeqFeatures() if you want the array to contain only
           the top-level features.

 Returns : An array of Bio::SeqFeatureI implementing objects.
 Args    : None

=cut

# this implementation is inherited from FeatureHolderI

=head2 feature_count

 Title   : feature_count
 Usage   : $seq->feature_count()
 Function: Return the number of SeqFeatures attached to a sequence
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


=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   : $seq->add_SeqFeature($feat);
 Function: Adds the given feature object to the feature array of this
           sequence. The object passed is required to implement the
           Bio::SeqFeatureI interface.
           The 'EXPAND' qualifier (see L<Bio::FeatureHolderI>) is supported, but
           has no effect,
 Returns : 1 on success
 Args    : A Bio::SeqFeatureI implementing object.

=cut

sub add_SeqFeature {
   my ($self, @feat) = @_;

   $self->{'_as_feat'} = [] unless $self->{'_as_feat'};

   if (scalar @feat > 1) {
      $self->deprecated(
         -message => 'Providing an array of features to Bio::Seq add_SeqFeature()'.
                     ' is deprecated and will be removed in a future version. '.
                     'Add a single feature at a time instead.',
         -warn_version    => 1.007,
         -throw_version   => 1.009,
      );
   }

   for my $feat ( @feat ) {

       next if $feat eq 'EXPAND'; # Need to support it for FeatureHolderI compliance

       if( !$feat->isa("Bio::SeqFeatureI") ) {
           $self->throw("Expected a Bio::SeqFeatureI object, but got a $feat.");
       }

       # make sure we attach ourselves to the feature if the feature wants it
       my $aseq = $self->primary_seq;
       $feat->attach_seq($aseq) if $aseq;

       push(@{$self->{'_as_feat'}},$feat);
   }
   return 1;
}


=head2 remove_SeqFeatures

 Title   : remove_SeqFeatures
 Usage   : $seq->remove_SeqFeatures();
 Function: Flushes all attached SeqFeatureI objects.

           To remove individual feature objects, delete those from the returned
           array and re-add the rest.
 Example :
 Returns : The array of Bio::SeqFeatureI objects removed from this seq.
 Args    : None

=cut

sub remove_SeqFeatures {
    my $self = shift;

    return () unless $self->{'_as_feat'};
    my @feats = @{$self->{'_as_feat'}};
    $self->{'_as_feat'} = [];
    return @feats;
}


=head1 Methods provided in the Bio::PrimarySeqI interface

These methods are inherited from the PrimarySeq interface
and work as one expects, building new Bio::Seq objects
or other information as expected. See L<Bio::PrimarySeq>
for more information.

Sequence Features are B<not> transferred to the new objects.
This is possibly a mistake. Anyone who feels the urge in
dealing with this is welcome to give it a go.

=head2 revcom

 Title   : revcom
 Usage   : $rev = $seq->revcom()
 Function: Produces a new Bio::Seq object which
           is the reversed complement of the sequence. For protein
           sequences this throws an exception of "Sequence is a protein.
           Cannot revcom"

           The id is the same id as the original sequence, and the
           accession number is also identical. If someone wants to track
           that this sequence has be reversed, it needs to define its own
           extensions

           To do an in-place edit of an object you can go:

           $seq = $seq->revcom();

           This of course, causes Perl to handle the garbage collection of
           the old object, but it is roughly speaking as efficient as an
           in-place edit.

 Returns : A new (fresh) Bio::Seq object
 Args    : None

=head2 trunc

 Title   : trunc
 Usage   : $subseq = $myseq->trunc(10,100);
 Function: Provides a truncation of a sequence

 Example :
 Returns : A fresh Seq object
 Args    : A Seq object

=head2 id

 Title   : id
 Usage   : $id = $seq->id()
 Function: This is mapped on display_id
 Returns : value of display_id()
 Args    : [optional] value to update display_id

=cut

sub  id {
    return shift->display_id(@_);
}


=head1 Seq only methods

These methods are specific to the Bio::Seq object, and not
found on the Bio::PrimarySeq object

=head2 primary_seq

 Title   : primary_seq
 Usage   : $seq->primary_seq or $seq->primary_seq($newval)
 Function: Get or set a PrimarySeq object
 Example :
 Returns : PrimarySeq object
 Args    : None or PrimarySeq object

=cut

sub primary_seq {
   my ($obj,$value) = @_;

   if( defined $value) {
       if( ! ref $value || ! $value->isa('Bio::PrimarySeqI') ) {
           $obj->throw("$value is not a Bio::PrimarySeq compliant object");
       }

       $obj->{'primary_seq'} = $value;
       # descend down over all seqfeature objects, seeing whether they
       # want an attached seq.

       foreach my $sf ( $obj->get_SeqFeatures() ) {
           $sf->attach_seq($value);
       }

   }
   return $obj->{'primary_seq'};

}


=head2 species

 Title   : species
 Usage   : $species = $seq->species() or $seq->species($species)
 Function: Gets or sets the species
 Returns : L<Bio::Species> object
 Args    : None or L<Bio::Species> object

See L<Bio::Species> for more information

=cut

sub species {
    my ($self, $species) = @_;
    if ($species) {
        $self->{'species'} = $species;
    } else {
        return $self->{'species'};
    }
}


# Internal methods follow...

# keep AUTOLOAD happy
sub DESTROY { }

############################################################################
# aliases due to name changes or to compensate for our lack of consistency #
############################################################################

# in all other modules we use the object in the singular --
# lack of consistency sucks
*flush_SeqFeature = \&remove_SeqFeatures;
*flush_SeqFeatures = \&remove_SeqFeatures;

# this is now get_SeqFeatures() (from FeatureHolderI)
*top_SeqFeatures = \&get_SeqFeatures;

# this is now get_all_SeqFeatures() in FeatureHolderI
sub all_SeqFeatures{
    return shift->get_all_SeqFeatures(@_);
}

sub accession {
    my $self = shift;
    $self->warn(ref($self)."::accession is deprecated, ".
                "use accession_number() instead");
    return $self->accession_number(@_);
}

1;
