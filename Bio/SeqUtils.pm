# BioPerl module for Bio::SeqUtils
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqUtils - Additional methods for PrimarySeq objects

=head1 SYNOPSIS

    use Bio::SeqUtils;
    # get a Bio::PrimarySeqI compliant object, $seq, somehow
    $util = Bio::SeqUtils->new();
    $polypeptide_3char = $util->seq3($seq);
    # or
    $polypeptide_3char = Bio::SeqUtils->seq3($seq);

    # set the sequence string (stored in one char code in the object)
    Bio::SeqUtils->seq3($seq, $polypeptide_3char);

    # translate a sequence in all six frames
    @seqs = Bio::SeqUtils->translate_6frames($seq);

    # inplace editing of the sequence
    Bio::SeqUtils->mutate($seq,
                          Bio::LiveSeq::Mutation->new(-seq => 'c',
                                                      -pos => 3
                                                     ));
    # mutate a sequence to desired similarity%
    $newseq = Bio::SeqUtils-> evolve
        ($seq, $similarity, $transition_transversion_rate);

    # concatenate two or more sequences with annotations and features,
    # the first sequence will be modified
    Bio::SeqUtils->cat(@seqs);
    my $catseq=$seqs[0];

    # truncate a sequence, retaining features and adjusting their
    # coordinates if necessary
    my $truncseq = Bio::SeqUtils->trunc_with_features($seq, 100, 200);

    # reverse complement a sequence and its features
    my $revcomseq = Bio::SeqUtils->revcom_with_features($seq);

    # simulate cloning of a fragment into a vector. Cut the vector at
    # positions 1000 and 1100 (deleting postions 1001 to 1099) and
    # "ligate" a fragment into the sites. The fragment is
    # reverse-complemented in this example (option "flip"). 
    # All features of the vector and fragment are preserved and 
    # features that are affected by the deletion/insertion are 
    # modified accordingly.
    # $vector and $fragment must be Bio::SeqI compliant objects 
    my $new_molecule = Bio::Sequtils->ligate(
      -vector => $vector, 
      -fragment => $fragment,
      -left => 1000,
      -right => 1100,
      -flip => 1 
    );

    # delete a segment of a sequence (from pos 1000 to 1100, inclusive),
    # again preserving features and annotations
    my $new_molecule = Bio::SeqUtils->cut( $seq, 1000, 1100 );

    # insert a fragment into a recipient between positions 1000 and
    # 1001. $recipient is a Bio::SeqI compliant object
    my $new_molecule =  Bio::SeqUtils::PbrTools->insert( 
      $recipient_seq, 
      $fragment_seq,
      1000
    );

=head1 DESCRIPTION

This class is a holder of methods that work on Bio::PrimarySeqI-
compliant sequence objects, e.g. Bio::PrimarySeq and
Bio::Seq. These methods are not part of the Bio::PrimarySeqI
interface and should in general not be essential to the primary function
of sequence objects. If you are thinking of adding essential
functions, it might be better to create your own sequence class.
See L<Bio::PrimarySeqI>, L<Bio::PrimarySeq>, and L<Bio::Seq> for more.

The methods take as their first argument a sequence object. It is
possible to use methods without first creating a SeqUtils object,
i.e. use it as an anonymous hash.

The first two methods, seq3() and seq3in(), give out or read in protein
sequences coded in three letter IUPAC amino acid codes.

The next two methods, translate_3frames() and translate_6frames(), wrap
around the standard translate method to give back an array of three
forward or all six frame translations.

The mutate() method mutates the sequence string with a mutation
description object.

The cat() method concatenates two or more sequences. The first sequence
is modified by addition of the remaining sequences. All annotations and
sequence features will be transferred.

The revcom_with_features() and trunc_with_features() methods are similar
to the revcom() and trunc() methods from Bio::Seq, but also adjust any
features associated with the sequence as appropriate.

There are also methods that simulate molecular cloning with rich
sequence objects. 
The delete() method cuts a segment out of a sequence and re-joins the
left and right fragments (like splicing or digesting and re-ligating a
molecule).  Positions (and types) of sequence features are adjusted
accordingly: 
Features that span the deleted segment are converted to split featuress
to indicate the disruption. (Sub)Features that extend into the deleted
segment are truncated.
A new molecule is created and returned.

The insert() method inserts a fragment (which can be a rich Bio::Seq
object) into another sequence object adding all annotations and
features to the final product.  
Features that span the insertion site are converted to split features
to indicate the disruption. 
A new feature is added to indicate the inserted fragment itself.  
A new molecule is created and returned.

The ligate() method simulates digesting a recipient (vector) and
ligating a fragment into it, which can also be flipped if needed. It
is simply a combination of a deletion and an insertion step and
returns a new molecule. The rules for modifying feature locations
outlined above are also used here, e.g. features that span the cut
sites are converted to split features with truncated sub-locations.


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

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki-at-bioperl-dot-org

=head1 CONTRIBUTORS

Roy R. Chaudhuri - roy.chaudhuri at gmail.com
Frank Schwach - frank.schwach@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqUtils;
use strict;
use warnings;
use Scalar::Util qw(blessed);
use parent qw(Bio::Root::Root);

# new inherited from RootI

our %ONECODE = (
    'Ala' => 'A',
    'Asx' => 'B',
    'Cys' => 'C',
    'Asp' => 'D',
    'Glu' => 'E',
    'Phe' => 'F',
    'Gly' => 'G',
    'His' => 'H',
    'Ile' => 'I',
    'Lys' => 'K',
    'Leu' => 'L',
    'Met' => 'M',
    'Asn' => 'N',
    'Pro' => 'P',
    'Gln' => 'Q',
    'Arg' => 'R',
    'Ser' => 'S',
    'Thr' => 'T',
    'Val' => 'V',
    'Trp' => 'W',
    'Xaa' => 'X',
    'Tyr' => 'Y',
    'Glx' => 'Z',
    'Ter' => '*',
    'Sec' => 'U',
    'Pyl' => 'O',
    'Xle' => 'J'
);

our %THREECODE = (
    'A' => 'Ala',
    'B' => 'Asx',
    'C' => 'Cys',
    'D' => 'Asp',
    'E' => 'Glu',
    'F' => 'Phe',
    'G' => 'Gly',
    'H' => 'His',
    'I' => 'Ile',
    'K' => 'Lys',
    'L' => 'Leu',
    'M' => 'Met',
    'N' => 'Asn',
    'P' => 'Pro',
    'Q' => 'Gln',
    'R' => 'Arg',
    'S' => 'Ser',
    'T' => 'Thr',
    'V' => 'Val',
    'W' => 'Trp',
    'Y' => 'Tyr',
    'Z' => 'Glx',
    'X' => 'Xaa',
    '*' => 'Ter',
    'U' => 'Sec',
    'O' => 'Pyl',
    'J' => 'Xle'
);

=head2 seq3

 Title   : seq3
 Usage   : $string = Bio::SeqUtils->seq3($seq)
 Function: Read only method that returns the amino acid sequence as a
           string of three letter codes. alphabet has to be
           'protein'. Output follows the IUPAC standard plus 'Ter' for
           terminator. Any unknown character, including the default
           unknown character 'X', is changed into 'Xaa'. A noncoded
           aminoacid selenocystein is recognized (Sec, U).

 Returns : A scalar
 Args    : character used for stop in the protein sequence optional,
           defaults to '*' string used to separate the output amino
           acid codes, optional, defaults to ''

=cut

sub seq3 {
    my ( $self, $seq, $stop, $sep ) = @_;

    $seq->isa('Bio::PrimarySeqI')
      || $self->throw('Not a Bio::PrimarySeqI object but [$self]');
    $seq->alphabet eq 'protein'
      || $self->throw('Not a protein sequence');

    if ( defined $stop ) {
        length $stop != 1
          and $self->throw('One character stop needed, not [$stop]');
        $THREECODE{$stop} = "Ter";
    }
    $sep ||= '';

    my $aa3s;
    foreach my $aa ( split //, uc $seq->seq ) {
        $THREECODE{$aa} and $aa3s .= $THREECODE{$aa} . $sep, next;
        $aa3s .= 'Xaa' . $sep;
    }
    $sep and substr( $aa3s, -( length $sep ), length $sep ) = '';
    return $aa3s;
}

=head2 seq3in

 Title   : seq3in
 Usage   : $seq = Bio::SeqUtils->seq3in($seq, 'MetGlyTer')
 Function: Method for changing of the sequence of a
           Bio::PrimarySeqI sequence object. The three letter amino
           acid input string is converted into one letter code.  Any
           unknown character triplet, including the default 'Xaa', is
           converted into 'X'.

 Returns : Bio::PrimarySeq object
 Args    : sequence string
           optional character to be used for stop in the protein sequence,
              defaults to '*'
           optional character to be used for unknown in the protein sequence,
              defaults to 'X'

=cut

sub seq3in {
    my ( $self, $seq, $string, $stop, $unknown ) = @_;

    $seq->isa('Bio::PrimarySeqI')
      || $self->throw("Not a Bio::PrimarySeqI object but [$self]");
    $seq->alphabet eq 'protein'
      || $self->throw('Not a protein sequence');

    if ( defined $stop ) {
        length $stop != 1
          and $self->throw("One character stop needed, not [$stop]");
        $ONECODE{'Ter'} = $stop;
    }
    if ( defined $unknown ) {
        length $unknown != 1
          and $self->throw("One character stop needed, not [$unknown]");
        $ONECODE{'Xaa'} = $unknown;
    }

    my ( $aas, $aa3 );
    my $length = ( length $string ) - 2;
    for ( my $i = 0 ; $i < $length ; $i += 3 ) {
        $aa3 = substr( $string, $i, 3 );
        $aa3 = ucfirst( lc($aa3) );
        $ONECODE{$aa3} and $aas .= $ONECODE{$aa3}, next;
        $aas .= $ONECODE{'Xaa'};
    }
    $seq->seq($aas);
    return $seq;
}

=head2 translate_3frames

 Title   : translate_3frames
 Usage   : @prots = Bio::SeqUtils->translate_3frames($seq)
 Function: Translate a nucleotide sequence in three forward frames.
           The IDs of the sequences are appended with '-0F', '-1F', '-2F'.
 Returns : An array of seq objects
 Args    : sequence object
           same arguments as to Bio::PrimarySeqI::translate

=cut

sub translate_3frames {
    my ( $self, $seq, @args ) = @_;

    $self->throw( 'Object [$seq] '
          . 'of class ['
          . ref($seq)
          . ']  can not be translated.' )
      unless $seq->can('translate');

    my ( $stop, $unknown, $frame, $tableid, $fullCDS, $throw ) = @args;
    my @seqs;
    my $f = 0;
    while ( $f != 3 ) {
        my $translation =
          $seq->translate( $stop, $unknown, $f, $tableid, $fullCDS, $throw );
        $translation->id( $seq->id . "-" . $f . "F" );
        push @seqs, $translation;
        $f++;
    }

    return @seqs;
}

=head2 translate_6frames

 Title   : translate_6frames
 Usage   : @prots = Bio::SeqUtils->translate_6frames($seq)
 Function: translate a nucleotide sequence in all six frames
           The IDs of the sequences are appended with '-0F', '-1F', '-2F',
           '-0R', '-1R', '-2R'.
 Returns : An array of seq objects
 Args    : sequence object
           same arguments as to Bio::PrimarySeqI::translate

=cut

sub translate_6frames {
    my ( $self, $seq, @args ) = @_;

    my @seqs  = $self->translate_3frames( $seq,         @args );
    my @seqs2 = $self->translate_3frames( $seq->revcom, @args );
    foreach my $seq2 (@seqs2) {
        my ($tmp) = $seq2->id;
        $tmp =~ s/F$/R/g;
        $seq2->id($tmp);
    }
    return @seqs, @seqs2;
}

=head2 valid_aa

 Title   : valid_aa
 Usage   : my @aa = $table->valid_aa
 Function: Retrieves a list of the valid amino acid codes.
           The list is ordered so that first 21 codes are for unique
           amino acids. The rest are ['B', 'Z', 'X', '*'].
 Returns : array of all the valid amino acid codes
 Args    : [optional] $code => [0 -> return list of 1 letter aa codes,
                                1 -> return list of 3 letter aa codes,
                                2 -> return associative array of both ]

=cut

sub valid_aa {
    my ( $self, $code ) = @_;

    if ( !$code ) {
        my @codes;
        foreach my $c ( sort values %ONECODE ) {
            push @codes, $c unless ( $c =~ /[BZX\*]/ );
        }
        push @codes, qw(B Z X *);    # so they are in correct order ?
        return @codes;
    }
    elsif ( $code == 1 ) {
        my @codes;
        foreach my $c ( sort keys %ONECODE ) {
            push @codes, $c unless ( $c =~ /(Asx|Glx|Xaa|Ter)/ );
        }
        push @codes, ( 'Asx', 'Glx', 'Xaa', 'Ter' );
        return @codes;
    }
    elsif ( $code == 2 ) {
        my %codes = %ONECODE;
        foreach my $c ( keys %ONECODE ) {
            my $aa = $ONECODE{$c};
            $codes{$aa} = $c;
        }
        return %codes;
    }
    else {
        $self->warn(
            "unrecognized code in " . ref($self) . " method valid_aa()" );
        return ();
    }
}

=head2 mutate

 Title   : mutate
 Usage   : Bio::SeqUtils->mutate($seq,$mutation1, $mutation2);
 Function: Inplace editing of the sequence.

           The second argument can be a Bio::LiveSeq::Mutation object
           or an array of them. The mutations are applied sequentially
           checking only that their position is within the current
           sequence.  Insertions are inserted before the given
           position.

 Returns : boolean
 Args    : sequence object
           mutation, a Bio::LiveSeq::Mutation object, or an array of them

See L<Bio::LiveSeq::Mutation>.

=cut

sub mutate {
    my ( $self, $seq, @mutations ) = @_;

    $self->throw( 'Object [$seq] '
          . 'of class ['
          . ref($seq)
          . '] should be a Bio::PrimarySeqI ' )
      unless $seq->isa('Bio::PrimarySeqI');
    $self->throw( 'Object [$mutations[0]] '
          . 'of class ['
          . ref( $mutations[0] )
          . '] should be a Bio::LiveSeq::Mutation' )
      unless $mutations[0]->isa('Bio::LiveSeq::Mutation');

    foreach my $mutation (@mutations) {
        $self->throw('Attempting to mutate sequence beyond its length')
          unless $mutation->pos - 1 <= $seq->length;

        my $string = $seq->seq;
        substr $string, $mutation->pos - 1, $mutation->len, $mutation->seq;
        $seq->seq($string);
    }
    1;
}

=head2 cat

  Title   : cat
  Usage   : Bio::SeqUtils->cat(@seqs);
            my $catseq=$seqs[0];
  Function: Concatenates a list of Bio::Seq objects, adding them all on to the
            end of the first sequence. Annotations and sequence features are
            copied over from any additional objects, and the coordinates of any
            copied features are adjusted appropriately.
  Returns : a boolean
  Args    : array of sequence objects

Note that annotations have no sequence locations. If you concatenate
sequences with the same annotations they will all be added.

=cut

sub cat {
    my ( $self, $seq, @seqs ) = @_;
    $self->throw( 'Object [$seq] '
          . 'of class ['
          . ref($seq)
          . '] should be a Bio::PrimarySeqI ' )
      unless $seq->isa('Bio::PrimarySeqI');

    for my $catseq (@seqs) {
        $self->throw( 'Object [$catseq] '
              . 'of class ['
              . ref($catseq)
              . '] should be a Bio::PrimarySeqI ' )
          unless $catseq->isa('Bio::PrimarySeqI');

        $self->throw(
                'Trying to concatenate sequences with different alphabets: '
              . $seq->display_id . '('
              . $seq->alphabet
              . ') and '
              . $catseq->display_id . '('
              . $catseq->alphabet
              . ')' )
          unless $catseq->alphabet eq $seq->alphabet;

        my $length = $seq->length;
        $seq->seq( $seq->seq . $catseq->seq );

        # move annotations
        if (    $seq->isa("Bio::AnnotatableI")
            and $catseq->isa("Bio::AnnotatableI") )
        {
            foreach my $key ( $catseq->annotation->get_all_annotation_keys() ) {

                foreach my $value ( $catseq->annotation->get_Annotations($key) )
                {
                    $seq->annotation->add_Annotation( $key, $value );
                }
            }
        }

        # move SeqFeatures
        if ( $seq->isa('Bio::SeqI') and $catseq->isa('Bio::SeqI') ) {
            for my $feat ( $catseq->get_SeqFeatures ) {
                $seq->add_SeqFeature( $self->_coord_adjust( $feat, $length ) );
            }
        }

    }
    1;
}

=head2 trunc_with_features

 Title   : trunc_with_features
 Usage   : $trunc=Bio::SeqUtils->trunc_with_features($seq, $start, $end);
 Function: Like Bio::Seq::trunc, but keeps features (adjusting coordinates
           where necessary. Features that partially overlap the region have
           their location changed to a Bio::Location::Fuzzy.
 Returns : A new sequence object
 Args    : A sequence object, start coordinate, end coordinate (inclusive)


=cut

sub trunc_with_features {
    use Bio::Range;
    my ( $self, $seq, $start, $end ) = @_;
    $self->throw( 'Object [$seq] '
          . 'of class ['
          . ref($seq)
          . '] should be a Bio::SeqI ' )
      unless $seq->isa('Bio::SeqI');
    my $trunc = $seq->trunc( $start, $end );
    my $truncrange =
      Bio::Range->new( -start => $start, -end => $end, -strand => 0 );

    # make sure that there is no annotation or features in $trunc
    # (->trunc() now clone objects except for Bio::Seq::LargePrimarySeq)
    $trunc->annotation->remove_Annotations;
    $trunc->remove_SeqFeatures;

    # move annotations
    foreach my $key ( $seq->annotation->get_all_annotation_keys() ) {
        foreach my $value ( $seq->annotation->get_Annotations($key) ) {
            $trunc->annotation->add_Annotation( $key, $value );
        }
    }

    # move features
    foreach (
        grep {
            $_ = $self->_coord_adjust( $_, 1 - $start, $end + 1 - $start )
              if $_->overlaps($truncrange)
        } $seq->get_SeqFeatures
      )
    {
        $trunc->add_SeqFeature($_);
    }
    return $trunc;
}

=head2 delete
 
 Title   : delete
 Function: cuts a segment out of a sequence and re-joins the left and right fragments
           (like splicing or digesting and re-ligating a molecule).
           Positions (and types) of sequence features are adjusted accordingly:
           Features that span the cut site are converted to split featuress to
           indicate the disruption. 
           Features that extend into the cut-out fragment are truncated.
           A new molecule is created and returned.
 Usage   : my $cutseq =  Bio::SeqUtils::PbrTools->cut( $seq, 1000, 1100 );
 Args    : a Bio::PrimarySeqI compliant object to cut,
           first nt of the segment to be deleted 
           last nt of the segment to be deleted  
           optional:
           hash-ref of options:
            clone_obj: if true, clone the input sequence object rather
                       than calling "new" on the object's class 

 Returns : a new Bio::Seq object 

=cut

sub delete {
    my $self = shift;
    my ( $seq, $left, $right, $opts_ref ) = @_;
    $self->throw( 'was expecting 3-4 paramters but got ' . @_ )
      unless @_ == 3 || @_ == 4;

    $self->throw(
        'Object of class [' . ref($seq) . '] should be a Bio::PrimarySeqI ' )
      unless blessed($seq) && $seq->isa('Bio::PrimarySeqI');

    $self->throw("Left coordinate ($left) must be >= 1") if $left < 1;
    if ( $right > $seq->length ) {
        $self->throw( "Right coordinate ($right) must be less than "
              . 'sequence length ('
              . $seq->length
              . ')' );
    }

    # piece together the sequence string of the remaining fragments
    my $left_seq  = $seq->subseq( 1,          $left - 1 );
    my $right_seq = $seq->subseq( $right + 1, $seq->length );
    if ( !$left_seq || !$right_seq ) {
        $self->throw(
'could not assemble sequences. At least one of the fragments is empty'
        );
    }
    my $seq_str = $left_seq . $right_seq;

    # create the new seq object with the same class as the recipient
    # or (if requested), make a clone of the existing object. In the
    # latter case we need to remove sequence features from the cloned
    # object instead of copying them
    my $product;
    if ( $opts_ref->{clone_obj} ) {
        $product = $self->_new_seq_via_clone( $seq, $seq_str );
    }
    else {
        $product = $self->_new_seq_from_old( $seq, { seq => $seq_str } );
    }

    # move sequence features
    if ( $product->isa('Bio::SeqI') && $seq->isa('Bio::SeqI') ) {
        for my $feat ( $seq->get_SeqFeatures ) {
            my $adjfeat = $self->_coord_adjust_deletion( $feat, $left, $right );
            $product->add_SeqFeature($adjfeat) if $adjfeat;
        }
    }

    # add a feature to annotatde the deletion
    my $deletion_feature = Bio::SeqFeature::Generic->new(
        -primary_tag => 'misc_feature',
        -tag      => { note => 'deletion of ' . ( $right - $left + 1 ) . 'bp' },
        -location => Bio::Location::Simple->new(
            -start         => $left - 1,
            -end           => $left,
            -location_type => 'IN-BETWEEN'
        )
    );
    $product->add_SeqFeature($deletion_feature);
    return $product;
}

=head2 insert
 
 Title   : insert
 Function: inserts a fragment (a Bio::Seq object) into a nother sequence object
           adding all annotations and features to the final product.
           Features that span the insertion site are converted to split
           features to indicate the disruption.
           A new feature is added to indicate the inserted fragment itself.
           A new molecule is created and returned.
 Usage   : # insert a fragment after pos 1000
           my $insert_seq =  Bio::SeqUtils::PbrTools->insert( 
             $recipient_seq, 
             $fragment_seq,
             1000
           );
 Args    : recipient sequence (a Bio::PrimarySeqI compliant object),
           a fragmetn to insert (Bio::PrimarySeqI compliant object), 
           insertion position (fragment is inserted to the right of this pos)
            pos=0 will prepend the fragment to the recipient
           optional:
           hash-ref of options:
            clone_obj: if true, clone the input sequence object rather
                       than calling "new" on the object's class 
 Returns : a new Bio::Seq object 

=cut

sub insert {
    my $self = shift;
    my ( $recipient, $fragment, $insert_pos, $opts_ref ) = @_;
    $self->throw( 'was expecting 3-4 paramters but got ' . @_ )
      unless @_ == 3 || @_ == 4;

    $self->throw( 'Recipient object of class ['
          . ref($recipient)
          . '] should be a Bio::PrimarySeqI ' )
      unless blessed($recipient) && $recipient->isa('Bio::PrimarySeqI');

    $self->throw( 'Fragment object of class ['
          . ref($fragment)
          . '] should be a Bio::PrimarySeqI ' )
      unless blessed($fragment) && $fragment->isa('Bio::PrimarySeqI');

    $self->throw( 'Can\'t concatenate sequences with different alphabets: '
          . 'recipient is '
          . $recipient->alphabet
          . ' and fragment is '
          . $fragment->alphabet )
      unless $recipient->alphabet eq $fragment->alphabet;

    if ( $insert_pos < 0 or $insert_pos > $recipient->length ) {
        $self->throw( "insertion position ($insert_pos) must be between 0 and "
              . 'recipient sequence length ('
              . $recipient->length
              . ')' );
    }

    if ( $fragment->can('is_circular') && $fragment->is_circular ) {
        $self->throw('Can\'t insert circular fragments');
    }

    if ( !$recipient->seq ) {
        $self->throw(
            'Recipient has no sequence, can not insert into this object');
    }

    # construct raw sequence of the new molecule
    my $left_seq =
        $insert_pos > 0
      ? $recipient->subseq( 1, $insert_pos )
      : '';
    my $mid_seq = $fragment->seq;
    my $right_seq =
        $insert_pos < $recipient->length
      ? $recipient->subseq( $insert_pos + 1, $recipient->length )
      : '';
    my $seq_str = $left_seq . $mid_seq . $right_seq;

    # create the new seq object with the same class as the recipient
    # or (if requested), make a clone of the existing object. In the
    # latter case we need to remove sequence features from the cloned
    # object instead of copying them
    my $product;
    if ( $opts_ref->{clone_obj} ) {
        $product = $self->_new_seq_via_clone( $recipient, $seq_str );
    }
    else {
        my @desc;
        push @desc, 'Inserted fragment: ' . $fragment->desc
          if defined $fragment->desc;
        push @desc, 'Recipient: ' . $recipient->desc
          if defined $recipient->desc;
        $product = $self->_new_seq_from_old(
            $recipient,
            {
                seq              => $seq_str,
                display_id       => $recipient->display_id,
                accession_number => $recipient->accession_number || '',
                alphabet         => $recipient->alphabet,
                desc             => join( '; ', @desc ),
                verbose          => $recipient->verbose || $fragment->verbose,
                is_circular      => $recipient->is_circular || 0,
            }
        );

    }    # if clone_obj

    # move annotations from fragment to product
    if (   $product->isa("Bio::AnnotatableI")
        && $fragment->isa("Bio::AnnotatableI") )
    {
        foreach my $key ( $fragment->annotation->get_all_annotation_keys ) {
            foreach my $value ( $fragment->annotation->get_Annotations($key) ) {
                $product->annotation->add_Annotation( $key, $value );
            }
        }
    }

    # move sequence features to product with adjusted coordinates
    if ( $product->isa('Bio::SeqI') ) {

        # for the fragment, just shift the features to new position
        if ( $fragment->isa('Bio::SeqI') ) {
            for my $feat ( $fragment->get_SeqFeatures ) {
                my $adjfeat = $self->_coord_adjust( $feat, $insert_pos );
                $product->add_SeqFeature($adjfeat) if $adjfeat;
            }
        }

        # for recipient, shift and modify features according to insertion.
        if ( $recipient->isa('Bio::SeqI') ) {
            for my $feat ( $recipient->get_SeqFeatures ) {
                my $adjfeat =
                  $self->_coord_adjust_insertion( $feat, $insert_pos,
                    $fragment->length );
                $product->add_SeqFeature($adjfeat) if $adjfeat;
            }
        }
    }

    # add a feature to annotate the insertion
    my $insertion_feature = Bio::SeqFeature::Generic->new(
        -start       => $insert_pos + 1,
        -end         => $insert_pos + $fragment->length,
        -primary_tag => 'misc_feature',
        -tag         => { note => 'inserted fragment' },
    );
    $product->add_SeqFeature($insertion_feature);

    return $product;
}

=head2 ligate
 
 title   : ligate
 function: pastes a fragment (which can also have features) into a recipient 
           sequence between two "cut" sites, preserving features and adjusting 
           their locations.
           This is a shortcut for deleting a segment from a sequence object followed
           by an insertion of a fragmnet and is supposed to be used to simulate
           in-vitro cloning where a recipient (a vector) is digested and a fragment 
           is then ligated into the recipient molecule. The fragment can be flipped
           (reverse-complemented with all its features).
           A new sequence object is returned to represent the product of the reaction.
           Features and annotations are transferred from the insert to the product
           and features on the recipient are adjusted according to the methods 
           L</"delete"> amd L</"insert">:
           Features spanning the insertion site will be split up into two sub-locations.
           (Sub-)features in the deleted region are themselves deleted.
           (Sub-)features that extend into the deleted region are truncated. 
           The class of the product object depends on the class of the recipient (vector)
           sequence object. if it is not possible to instantiate a new
           object of that class, a Bio::Primaryseq object is created instead.
 usage   : # insert the flipped fragment between positions 1000 and 1100 of the 
           # vector, i.e. everything between these two positions is deleted and
           # replaced by the fragment
           my $new_molecule = Bio::Sequtils::Pbrtools->ligate(
             -recipient => $vector, 
             -fragment => $fragment,
             -left => 1000,
             -right => 1100,
             -flip      => 1, 
             -clone_obj => 1 
           );
 args    : recipient: the recipient/vector molecule
           fragment: molecule that is to be ligated into the vector
           left: left cut site (fragment will be inserted to the right of 
                 this position)
           optional:
            right: right cut site (fragment will be inseterted to the 
                   left of this position). defaults to left+1
            flip: boolean, if true, the fragment is reverse-complemented 
                  (including features) before inserting
            clone_obj: if true, clone the recipient object to create the product
                       instead of calling "new" on its class
 returns : a new Bio::Seq object of the ligated fragments

=cut

sub ligate {
    my $self = shift;
    my ( $recipient, $fragment, $left, $right, $flip, $clone_obj ) =
      $self->_rearrange( [qw(RECIPIENT FRAGMENT LEFT RIGHT FLIP CLONE_OBJ )],
        @_ );
    $self->throw("missing required parameter 'recipient'") unless $recipient;
    $self->throw("missing required parameter 'fragment'")  unless $fragment;
    $self->throw("missing required parameter 'left'")      unless defined $left;

    $right ||= $left + 1;

    $self->throw(
        "Fragment must be a Bio::PrimarySeqI compliant object but it is a "
          . ref($fragment) )
      unless blessed($fragment) && $fragment->isa('Bio::PrimarySeqI');

    $fragment = $self->revcom_with_features($fragment) if $flip;

    my $opts_ref = {};
    $opts_ref->{clone_obj} = 1 if $clone_obj;

    # clone in two steps: first delete between the insertion sites,
    # then insert the fragment. Step 1 is skipped if insert positions
    # are adjacent (no deletion)
    my ( $product1, $product2 );
    eval {
        if ( $right == $left + 1 ) {
            $product1 = $recipient;
        }
        else {
            $product1 =
              $self->delete( $recipient, $left + 1, $right - 1, $opts_ref );
        }
    };
    $self->throw( "Failed in step 1 (cut recipient): " . $@ ) if $@;
    eval { $product2 = $self->insert( $product1, $fragment, $left, $opts_ref ) };
    $self->throw( "Failed in step 2 (insert fragment): " . $@ ) if $@;

    return $product2;

}

=head2 _coord_adjust_deletion
 
 title   : _coord_adjust_deletion
 function: recursively adjusts coordinates of seqfeatures on a molecule
           where a segment has been deleted.
           (sub)features that span the deletion site become split features.
           (sub)features that extend into the deletion site are truncated.
           A note is added to the feature to inform about the size and
           position of the deletion.
 usage   : my $adjusted_feature = Bio::Sequtils::_coord_adjust_deletion( 
             $feature,
             $start,
             $end
           );
 args    : a Bio::SeqFeatureI compliant object,
           start (inclusive) position of the deletion site,
           end (inclusive) position of the deletion site
 returns : a Bio::SeqFeatureI compliant object 

=cut

sub _coord_adjust_deletion {
    my ( $self, $feat, $left, $right ) = @_;

    $self->throw( 'object [$feat] '
          . 'of class ['
          . ref($feat)
          . '] should be a Bio::SeqFeatureI ' )
      unless $feat->isa('Bio::SeqFeatureI');
    $self->throw('missing coordinates: need a left and a right position')
      unless defined $left && defined $right;

    if ( $left > $right ) {
        if ( $feat->can('is_circular') && $feat->is_circular ) {

            # todo handle circular molecules
            $self->throw(
'can not yet handle deletions in circular molecules if deletion spans origin'
            );
        }
        else {
            $self->throw(
                    "left coordinate ($left) must be less than right ($right)"
                  . " but it was greater" );
        }
    }
    my $deletion = Bio::Location::Simple->new(
        -start => $left,
        -end   => $right,
    );
    my $del_length = $right - $left + 1;

    my @adjsubfeat;
    for my $subfeat ( $feat->get_SeqFeatures ) {
        my $adjsubfeat =
          $self->_coord_adjust_deletion( $subfeat, $left, $right );
        push @adjsubfeat, $adjsubfeat if $adjsubfeat;
    }

    my @loc;
    my $note;
    for ( $feat->location->each_Location ) {
        next if $deletion->contains($_);    # this location will be deleted;
        my $strand     = $_->strand;
        my $type       = $_->location_type;
        my $start      = $_->start;
        my $start_type = $_->can('start_pos_type') ? $_->start_pos_type : undef;
        my $end        = $_->end;
        my $end_type   = $_->can('end_pos_type') ? $_->end_pos_type : undef;
        my @newcoords  = ();
        if ( $start < $deletion->start && $end > $deletion->end )
        {                                   # split the feature
            @newcoords = (
                [ $start, ( $deletion->start - 1 ), $start_type, $end_type ],
                [
                    ( $deletion->start ), $end - $del_length,
                    $start_type,          $end_type
                ]
            );
            $note =
                $del_length
              . 'bp internal deletion between pos '
              . ( $deletion->start - 1 ) . ' and '
              . $deletion->start;
        }
        elsif ( $_->start < $deletion->start && $_->end >= $deletion->start )
        {    # truncate feature end
            @newcoords =
              ( [ $start, ( $deletion->start - 1 ), $start_type, $end_type ] );
            $note =
              ( $end - $deletion->start + 1 ) . 'bp deleted from feature ';
            if ( $feat->strand ) {
                $note .= $feat->strand == 1 ? "3' " : "5' ";
            }
            $note .= 'end';
        }
        elsif ( $_->start <= $deletion->end && $_->end > $deletion->end )
        {    # truncate feature start and shift end
            @newcoords = (
                [
                    ( $deletion->start ), $end - $del_length,
                    $start_type,          $end_type
                ]
            );
            $note =
              ( $deletion->end - $start + 1 ) . 'bp deleted from feature ';
            if ( $feat->strand ) {
                $note .= $feat->strand == 1 ? "5' end" : "3' end";
            }
            else {
                $note .= 'start';
            }
        }
        elsif ( $start >= $deletion->end ) {    # just shift entire location
            @newcoords = (
                [
                    $start - $del_length, $end - $del_length,
                    $start_type,          $end_type
                ]
            );
        }
        else {                                  # not affected by deletion
            @newcoords = ( [ $start, $end, $start_type, $end_type ] );
        }

        # if we have no coordinates, we return nothing
        # the feature is deleted
        return unless @newcoords;

        my @subloc =
          $self->_location_objects_from_coordinate_list( \@newcoords, $strand,
            $type );
        push @loc, $self->_single_loc_object_from_collection(@subloc);
    }    # each location

    # create new feature based on original one and move annotation across
    my $newfeat =
      Bio::SeqFeature::Generic->new( -primary => $feat->primary_tag );
    foreach my $key ( $feat->annotation->get_all_annotation_keys() ) {
        foreach my $value ( $feat->annotation->get_Annotations($key) ) {
            $newfeat->annotation->add_Annotation( $key, $value );
        }
    }
    foreach my $key ( $feat->get_all_tags() ) {
        $newfeat->add_tag_value( $key, $feat->get_tag_values($key) );
    }

    # If we have a note about the deleted bases, add it
    if ($note) {
        $newfeat->add_tag_value( 'note', $note );
    }

    # set modified location(s) for the new feature and
    # add its subfeatures if any
    my $loc = $self->_single_loc_object_from_collection(@loc);
    $loc ? $newfeat->location($loc) : return;
    $newfeat->add_SeqFeature($_) for @adjsubfeat;

    return $newfeat;

}

=head2 _coord_adjust_insertion
 
 title   : _coord_adjust_insertion
 function: recursively adjusts coordinates of seqfeatures on a molecule
           where another sequence has been inserted.
           (sub)features that span the insertion site become split features
           and a note is added about the size and positin of the insertion.
           Features with an IN-BETWEEN location at the insertion site
           are lost (such features can only exist between adjacent bases)
 usage   : my $adjusted_feature = Bio::Sequtils::_coord_adjust_insertion( 
             $feature,
             $insert_pos,
             $insert_length
           );
 args    : a Bio::SeqFeatureI compliant object,
           insertion position (insert to the right of this position)
           length of inserted fragment
 returns : a Bio::SeqFeatureI compliant object 

=cut

sub _coord_adjust_insertion {
    my ( $self, $feat, $insert_pos, $insert_len ) = @_;

    $self->throw( 'object [$feat] '
          . 'of class ['
          . ref($feat)
          . '] should be a Bio::SeqFeatureI ' )
      unless $feat->isa('Bio::SeqFeatureI');
    $self->throw('missing insert position') unless defined $insert_pos;
    $self->throw('missing insert length')   unless defined $insert_len;

    my @adjsubfeat;
    for my $subfeat ( $feat->get_SeqFeatures ) {
        push @adjsubfeat,
          $self->_coord_adjust_insertion( $subfeat, $insert_pos, $insert_len );
    }

    my @loc;
    my $note;
    for ( $feat->location->each_Location ) {

        # loose IN-BETWEEN features at the insertion site
        next
          if ( $_->location_type eq 'IN-BETWEEN' && $_->start == $insert_pos );
        my $strand     = $_->strand;
        my $type       = $_->location_type;
        my $start      = $_->start;
        my $start_type = $_->can('start_pos_type') ? $_->start_pos_type : undef;
        my $end        = $_->end;
        my $end_type   = $_->can('end_pos_type') ? $_->end_pos_type : undef;
        my @newcoords  = ();
        if ( $start <= $insert_pos && $end > $insert_pos ) { # split the feature
            @newcoords = (
                [ $start, $insert_pos, $start_type, $end_type ],
                [
                    ( $insert_pos + 1 + $insert_len ), $end + $insert_len,
                    $start_type,                       $end_type
                ]
            );
            $note =
                $insert_len
              . 'bp internal insertion between pos '
              . $insert_pos . ' and '
              . ( $insert_pos + $insert_len + 1 );

        }
        elsif ( $start > $insert_pos ) {    # just shift entire location
            @newcoords = (
                [
                    $start + $insert_len, $end + $insert_len,
                    $start_type,          $end_type
                ]
            );
        }
        else {                              # not affected
            @newcoords = ( [ $start, $end, $start_type, $end_type ] );
        }

        # if we have deleted all coordinates, return nothing
        # (possible if all locations are IN-BETWEEN)
        return unless @newcoords;

        my @subloc =
          $self->_location_objects_from_coordinate_list( \@newcoords, $strand,
            $type );

        # put together final location which could be a split now
        push @loc, $self->_single_loc_object_from_collection(@subloc);
    }    # each location

    # create new feature based on original one and move annotation across
    my $newfeat =
      Bio::SeqFeature::Generic->new( -primary => $feat->primary_tag );
    foreach my $key ( $feat->annotation->get_all_annotation_keys() ) {
        foreach my $value ( $feat->annotation->get_Annotations($key) ) {
            $newfeat->annotation->add_Annotation( $key, $value );
        }
    }
    foreach my $key ( $feat->get_all_tags() ) {
        $newfeat->add_tag_value( $key, $feat->get_tag_values($key) );
    }

    # If we have a note about the inserted bases, add it
    if ($note) {
        $newfeat->add_tag_value( 'note', $note );
    }

    # set modified location(s) for the new feature and
    # add its subfeatures if any
    my $loc = $self->_single_loc_object_from_collection(@loc);
    $loc ? $newfeat->location($loc) : return;
    $newfeat->add_SeqFeature($_) for @adjsubfeat;

    return $newfeat;

}

=head2 _single_loc_object_from_collection
 
 Title   : _single_loc_object_from_collection
 Function: takes an array of location objects. Returns either a split 
           location object if there are more than one locations in the 
           array or returns the single location if there is only one
 Usage   : my $loc = _single_loc_object_from_collection( @sublocs );
 Args    : array of Bio::Location objects
 Returns : a single Bio:;Location object containing all locations

=cut

sub _single_loc_object_from_collection {
    my ( $self, @locs ) = @_;
    my $loc;
    if ( @locs > 1 ) {
        $loc = Bio::Location::Split->new;
        $loc->add_sub_Location(@locs);
    }
    elsif ( @locs == 1 ) {
        $loc = shift @locs;
    }
    return $loc;
}    # _single_loc_object_from_collection

=head2 _location_objects_from_coordinate_list
 
 Title   : _location_objects_from_coordinate_list
 Function: takes an array-ref of start/end coordinates, a strand and a 
           type and returns a list of Bio::Location objects (Fuzzy by 
           default, Simple in case of in-between coordinates).
           If location type is not "IN-BETWEEN", individual types may be
           passed in for start and end location as per Bio::Location::Fuzzy
           documentation.
 Usage   : my @loc_objs = $self->_location_objects_from_coordinate_list( 
             \@coords, 
             $strand, 
             $type
           );
 Args    : array-ref of array-refs each containing:
           start, end [, start-type, end-type]   
             where types are optional. If given, must be
             a one of ('BEFORE', 'AFTER', 'EXACT','WITHIN', 'BETWEEN')
           strand (all locations must be on same strand)
           location-type (EXACT, IN-BETWEEN etc)
 Returns : list of Bio::Location objects

=cut

sub _location_objects_from_coordinate_list {
    my $self = shift;
    my ( $coords_ref, $strand, $type ) = @_;
    $self->throw( 'expected 3 parameters but got ' . @_ ) unless @_ == 3;
    $self->throw('first argument must be an ARRAY reference#')
      unless ref($coords_ref) eq 'ARRAY';

    my @loc;
    foreach my $coords_set (@$coords_ref) {
        my ( $start, $end, $start_type, $end_type ) = @$coords_set;

        # taken from Bio::SeqUtils::_coord_adjust
        if ( $type ne 'IN-BETWEEN' ) {
            my $loc = Bio::Location::Fuzzy->new(
                -start         => $start,
                -end           => $end,
                -strand        => $strand,
                -location_type => $type
            );
            $loc->start_pos_type($start_type) if $start_type;
            $loc->end_pos_type($end_type)     if $end_type;
            push @loc, $loc;
        }
        else {
            push @loc,
              Bio::Location::Simple->new(
                -start         => $start,
                -end           => $end,
                -strand        => $strand,
                -location_type => $type
              );
        }
    }    # each coords_set
    return @loc;
}    # _location_objects_from_coordinate_list

=head2 _new_seq_via_clone
 
 Title   : _new_seq_via_clone
 Function: clone a sequence object using Bio::Root::Root::clone and set the new sequence string
           sequence features are removed.
 Usage   : my $new_seq = $self->_new_seq_via_clone( $seq_obj, $seq_str );
 Args    : original seq object [, new sequence string]
 Returns : a clone of the original sequence object, optionally with new sequence string

=cut

sub _new_seq_via_clone {
    my ( $self, $in_seq_obj, $seq_str ) = @_;
    my $out_seq_obj = $in_seq_obj->clone;
    $out_seq_obj->remove_SeqFeatures if $out_seq_obj->can('remove_SeqFeatures');
    if ( blessed $out_seq_obj->seq
        && $out_seq_obj->seq->isa('Bio::PrimarySeq') )
    {
        $out_seq_obj->seq->seq($seq_str);
    }
    else {
        $out_seq_obj->seq($seq_str);
    }
    return $out_seq_obj;

}    # _new_seq_via_clone

=head2 _new_seq_from_old
 
 Title   : _new_seq_from_old
 Function: creates a new sequence obejct, if possible of the same class as the old and adds 
           attributes to it. Also copies annotation across to the new object.
 Usage   : my $new_seq = $self->_new_seq_from_old( $seq_obj, { seq => $seq_str, display_id => 'some_ID'});
 Args    : old sequence object
           hashref of attributes for the new sequence (sequence string etc.)
 Returns : a new Bio::Seq object

=cut

sub _new_seq_from_old {
    my ( $self, $in_seq_obj, $attr ) = @_;
    $self->throw('attributes must be a hashref')
      if $attr && ref($attr) ne 'HASH';

    my $seqclass;
    if ( $in_seq_obj->can_call_new ) {
        $seqclass = ref($in_seq_obj);
    }
    else {
        $seqclass = 'Bio::Primaryseq';
        $self->_attempt_to_load_seq;
    }

    my $out_seq_obj = $seqclass->new(
        -seq        => $attr->{seq}        || $in_seq_obj->seq,
        -display_id => $attr->{display_id} || $in_seq_obj->display_id,
        -accession_number => $attr->{accession_number}
          || $in_seq_obj->accession_number
          || '',
        -alphabet    => $in_seq_obj->alphabet,
        -desc        => $attr->{desc} || $in_seq_obj->desc,
        -verbose     => $attr->{verbose} || $in_seq_obj->verbose,
        -is_circular => $attr->{is_circular} || $in_seq_obj->is_circular || 0,
    );

    # move the annotation across to the product
    if (   $out_seq_obj->isa("Bio::AnnotatableI")
        && $in_seq_obj->isa("Bio::AnnotatableI") )
    {
        foreach my $key ( $in_seq_obj->annotation->get_all_annotation_keys ) {
            foreach my $value ( $in_seq_obj->annotation->get_Annotations($key) )
            {
                $out_seq_obj->annotation->add_Annotation( $key, $value );
            }
        }
    }
    return $out_seq_obj;
}    # _new_seq_from_old

=head2 _coord_adjust

  Title   : _coord_adjust
  Usage   : my $newfeat=Bio::SeqUtils->_coord_adjust($feature, 100, $seq->length);
  Function: Recursive subroutine to adjust the coordinates of a feature
            and all its subfeatures. If a sequence length is specified, then
            any adjusted features that have locations beyond the boundaries
            of the sequence are converted to Bio::Location::Fuzzy objects.

  Returns : A Bio::SeqFeatureI compliant object.
  Args    : A Bio::SeqFeatureI compliant object,
            the number of bases to add to the coordinates
            (optional) the length of the parent sequence


=cut

sub _coord_adjust {
    my ( $self, $feat, $add, $length ) = @_;
    $self->throw( 'Object [$feat] '
          . 'of class ['
          . ref($feat)
          . '] should be a Bio::SeqFeatureI ' )
      unless $feat->isa('Bio::SeqFeatureI');
    my @adjsubfeat;
    for my $subfeat ( $feat->get_SeqFeatures ) {
        push @adjsubfeat, $self->_coord_adjust( $subfeat, $add, $length );
    }
    my @loc;
    for ( $feat->location->each_Location ) {
        my @coords = ( $_->start, $_->end );
        my $strand = $_->strand;
        my $type   = $_->location_type;
        foreach (@coords) {
            $self->throw("can not handle negative feature positions (got: $_)")
              if $_ < 0;
            if ( $add + $_ < 1 ) {
                $_ = '<1';
            }
            elsif ( defined $length and $add + $_ > $length ) {
                $_ = ">$length";
            }
            else {
                $_ = $add + $_;
            }
        }
        push @loc,
          $self->_location_objects_from_coordinate_list( [ \@coords ],
            $strand, $type );
    }
    my $newfeat =
      Bio::SeqFeature::Generic->new( -primary => $feat->primary_tag );
    foreach my $key ( $feat->annotation->get_all_annotation_keys() ) {
        foreach my $value ( $feat->annotation->get_Annotations($key) ) {
            $newfeat->annotation->add_Annotation( $key, $value );
        }
    }
    foreach my $key ( $feat->get_all_tags() ) {
        $newfeat->add_tag_value( $key, $feat->get_tag_values($key) );
    }
    my $loc = $self->_single_loc_object_from_collection(@loc);
    $loc ? $newfeat->location($loc) : return;
    $newfeat->add_SeqFeature($_) for @adjsubfeat;
    return $newfeat;
}

=head2 revcom_with_features

 Title   : revcom_with_features
 Usage   : $revcom=Bio::SeqUtils->revcom_with_features($seq);
 Function: Like Bio::Seq::revcom, but keeps features (adjusting coordinates
           as appropriate.
 Returns : A new sequence object
 Args    : A sequence object


=cut

sub revcom_with_features {
    my ( $self, $seq ) = @_;
    $self->throw( 'Object [$seq] '
          . 'of class ['
          . ref($seq)
          . '] should be a Bio::SeqI ' )
      unless $seq->isa('Bio::SeqI');
    my $revcom = $seq->revcom;

    # make sure that there is no annotation or features in $trunc
    # (->revcom() now clone objects except for Bio::Seq::LargePrimarySeq)
    $revcom->annotation->remove_Annotations;
    $revcom->remove_SeqFeatures;

    #move annotations
    foreach my $key ( $seq->annotation->get_all_annotation_keys() ) {
        foreach my $value ( $seq->annotation->get_Annotations($key) ) {
            $revcom->annotation->add_Annotation( $key, $value );
        }
    }

    #move features
    for ( map { $self->_feature_revcom( $_, $seq->length ) }
        reverse $seq->get_SeqFeatures )
    {
        $revcom->add_SeqFeature($_);
    }
    return $revcom;
}

=head2 _feature_revcom

  Title   : _feature_revcom
  Usage   : my $newfeat=Bio::SeqUtils->_feature_revcom($feature, $seq->length);
  Function: Recursive subroutine to reverse complement a feature and
            all its subfeatures. The length of the parent sequence must be
            specified.

  Returns : A Bio::SeqFeatureI compliant object.
  Args    : A Bio::SeqFeatureI compliant object,
            the length of the parent sequence


=cut

sub _feature_revcom {
    my ( $self, $feat, $length ) = @_;
    $self->throw( 'Object [$feat] '
          . 'of class ['
          . ref($feat)
          . '] should be a Bio::SeqFeatureI ' )
      unless $feat->isa('Bio::SeqFeatureI');
    my @adjsubfeat;
    for my $subfeat ( $feat->get_SeqFeatures ) {
        push @adjsubfeat, $self->_feature_revcom( $subfeat, $length );
    }
    my @loc;
    for ( $feat->location->each_Location ) {
        my $type = $_->location_type;
        my $strand;
        if    ( $_->strand == -1 ) { $strand = 1 }
        elsif ( $_->strand == 1 )  { $strand = -1 }
        else                       { $strand = $_->strand }
        my $newend =
          $self->_coord_revcom( $_->start, $_->start_pos_type, $length );
        my $newstart =
          $self->_coord_revcom( $_->end, $_->end_pos_type, $length );
        my $newstart_type = $_->end_pos_type;
        $newstart_type = 'BEFORE' if $_->end_pos_type eq 'AFTER';
        $newstart_type = 'AFTER'  if $_->end_pos_type eq 'BEFORE';
        my $newend_type = $_->start_pos_type;
        $newend_type = 'BEFORE' if $_->start_pos_type eq 'AFTER';
        $newend_type = 'AFTER'  if $_->start_pos_type eq 'BEFORE';
        push @loc,
          $self->_location_objects_from_coordinate_list(
            [ [ $newstart, $newend, $newstart_type, $newend_type ] ],
            $strand, $type );
    }
    my $newfeat =
      Bio::SeqFeature::Generic->new( -primary => $feat->primary_tag );
    foreach my $key ( $feat->annotation->get_all_annotation_keys() ) {
        foreach my $value ( $feat->annotation->get_Annotations($key) ) {
            $newfeat->annotation->add_Annotation( $key, $value );
        }
    }
    foreach my $key ( $feat->get_all_tags() ) {
        $newfeat->add_tag_value( $key, $feat->get_tag_values($key) );
    }

    my $loc = $self->_single_loc_object_from_collection(@loc);
    $loc ? $newfeat->location($loc) : return;

    $newfeat->add_SeqFeature($_) for @adjsubfeat;
    return $newfeat;
}

sub _coord_revcom {
    my ( $self, $coord, $type, $length ) = @_;
    if ( $type eq 'BETWEEN' or $type eq 'WITHIN' ) {
        $coord =~ s/(\d+)(\D*)(\d+)/$length+1-$3.$2.$length+1-$1/ge;
    }
    else {
        $coord =~ s/(\d+)/$length+1-$1/ge;
        $coord =~ tr/<>/></;
        $coord = '>' . $coord
          if $type eq 'BEFORE' and substr( $coord, 0, 1 ) ne '>';
        $coord = '<' . $coord
          if $type eq 'AFTER' and substr( $coord, 0, 1 ) ne '<';
    }
    return $coord;
}

=head2 evolve

  Title   : evolve
  Usage   : my $newseq = Bio::SeqUtils->
                evolve($seq, $similarity, $transition_transversion_rate);
  Function: Mutates the sequence by point mutations until the similarity of
            the new sequence has decreased to the required level.
            Transition/transversion rate is adjustable.
  Returns : A new Bio::PrimarySeq object
  Args    : sequence object
            percentage similarity (e.g. 80)
            tr/tv rate, optional, defaults to 1 (= 1:1)

Set the verbosity of the Bio::SeqUtils object to positive integer to
see the mutations as they happen.

This method works only on nucleotide sequences. It prints a warning if
you set the target similarity to be less than 25%.

Transition/transversion ratio is an observed attribute of an sequence
comparison. We are dealing here with the transition/transversion rate
that we set for our model of sequence evolution.

=cut

sub evolve {
    my ( $self, $seq, $sim, $rate ) = @_;
    $rate ||= 1;

    $self->throw( 'Object [$seq] '
          . 'of class ['
          . ref($seq)
          . '] should be a Bio::PrimarySeqI ' )
      unless $seq->isa('Bio::PrimarySeqI');

    $self->throw(
        "[$sim] " . ' should be a positive integer or float under 100' )
      unless $sim =~ /^[+\d.]+$/ and $sim <= 100;

    $self->warn(
        "Nucleotide sequences are 25% similar by chance.
        Do you really want to set similarity to [$sim]%?\n"
    ) unless $sim > 25;

    $self->throw('Only nucleotide sequences are supported')
      if $seq->alphabet eq 'protein';

    # arrays of possible changes have transitions as first items
    my %changes;
    $changes{'a'} = [ 't', 'c', 'g' ];
    $changes{'t'} = [ 'a', 'c', 'g' ];
    $changes{'c'} = [ 'g', 'a', 't' ];
    $changes{'g'} = [ 'c', 'a', 't' ];

    # given the desired rate, find out where cut off points need to be
    # when random numbers are generated from 0 to 100
    # we are ignoring identical mutations (e.g. A->A) to speed things up
    my $bin_size   = 100 / ( $rate + 2 );
    my $transition = 100 - ( 2 * $bin_size );
    my $first_transversion = $transition + $bin_size;

    # unify the look of sequence strings
    my $string = lc $seq->seq;    # lower case
    $string =~
      s/u/t/;    # simplyfy our life; modules should deal with the change anyway
                 # store the original sequence string
    my $oristring = $string;
    my $length    = $seq->length;

    # stop evolving if the limit has been reached
    until ( $self->_get_similarity( $oristring, $string ) <= $sim ) {

        # find the location in the string to change
        my $loc = int( rand $length ) + 1;

        # nucleotide to change
        my $oldnuc = substr $string, $loc - 1, 1;
        my $newnuc;

        # nucleotide it is changed to
        my $choose = rand(100);
        if ( $choose < $transition ) {
            $newnuc = $changes{$oldnuc}[0];
        }
        elsif ( $choose < $first_transversion ) {
            $newnuc = $changes{$oldnuc}[1];
        }
        else {
            $newnuc = $changes{$oldnuc}[2];
        }

        # do the change
        substr $string, $loc - 1, 1, $newnuc;

        $self->debug("$loc$oldnuc>$newnuc\n");
    }

    return new Bio::PrimarySeq(
        -id          => $seq->id . "-$sim",
        -description => $seq->description,
        -seq         => $string
    );
}

sub _get_similarity {
    my ( $self, $oriseq, $seq ) = @_;

    my $len = length($oriseq);
    my $c;

    for ( my $i = 0 ; $i < $len ; $i++ ) {
        $c++ if substr( $oriseq, $i, 1 ) eq substr( $seq, $i, 1 );
    }
    return 100 * $c / $len;
}

1;
