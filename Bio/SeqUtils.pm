#
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

    # truncate a sequence, retaining features and adjusting their
    # coordinates if necessary
    my $truncseq = Bio::SeqUtils->trunc_with_features($seq, 100, 200);

    # reverse complement a sequence and its features
    my $revcomseq = Bio::SeqUtils->revcom_with_features($seq);

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

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki-at-bioperl-dot-org

=head1 CONTRIBUTORS

Roy R. Chaudhuri - roy.chaudhuri at gmail.com

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqUtils;
use vars qw(%ONECODE %THREECODE);
use strict;
use Carp;

use base qw(Bio::Root::Root);
# new inherited from RootI

BEGIN {
    # Note : Ambiguity code 'J' = I/L (used for ambiguities in mass-spec data)
    %ONECODE =
    ('Ala' => 'A', 'Asx' => 'B', 'Cys' => 'C', 'Asp' => 'D',
     'Glu' => 'E', 'Phe' => 'F', 'Gly' => 'G', 'His' => 'H',
     'Ile' => 'I', 'Lys' => 'K', 'Leu' => 'L', 'Met' => 'M',
     'Asn' => 'N', 'Pro' => 'P', 'Gln' => 'Q', 'Arg' => 'R',
     'Ser' => 'S', 'Thr' => 'T', 'Val' => 'V', 'Trp' => 'W',
     'Xaa' => 'X', 'Tyr' => 'Y', 'Glx' => 'Z', 'Ter' => '*',
     'Sec' => 'U', 'Pyl' => 'O', 'Xle' => 'J'
     );

    %THREECODE =
    ('A' => 'Ala', 'B' => 'Asx', 'C' => 'Cys', 'D' => 'Asp',
     'E' => 'Glu', 'F' => 'Phe', 'G' => 'Gly', 'H' => 'His',
     'I' => 'Ile', 'K' => 'Lys', 'L' => 'Leu', 'M' => 'Met',
     'N' => 'Asn', 'P' => 'Pro', 'Q' => 'Gln', 'R' => 'Arg',
     'S' => 'Ser', 'T' => 'Thr', 'V' => 'Val', 'W' => 'Trp',
     'Y' => 'Tyr', 'Z' => 'Glx', 'X' => 'Xaa', '*' => 'Ter',
     'U' => 'Sec', 'O' => 'Pyl', 'J' => 'Xle'
     );
}

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
   my ($self, $seq, $stop, $sep ) = @_;

   $seq->isa('Bio::PrimarySeqI') ||
       $self->throw('Not a Bio::PrimarySeqI object but [$self]');
   $seq->alphabet eq 'protein' ||
       $self->throw('Not a protein sequence');

   if (defined $stop) {
       length $stop != 1 and $self->throw('One character stop needed, not [$stop]');
       $THREECODE{$stop} = "Ter";
   }
   $sep ||= '';

   my $aa3s;
   foreach my $aa  (split //, uc $seq->seq) {
       $THREECODE{$aa} and $aa3s .= $THREECODE{$aa}. $sep, next;
       $aa3s .= 'Xaa'. $sep;
   }
   $sep and substr($aa3s, -(length $sep), length $sep) = '' ;
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
   my ($self, $seq, $string, $stop, $unknown) = @_;

   $seq->isa('Bio::PrimarySeqI') ||
	  $self->throw("Not a Bio::PrimarySeqI object but [$self]");
   $seq->alphabet eq 'protein' ||
	  $self->throw('Not a protein sequence');

   if (defined $stop) {
		length $stop != 1 and $self->throw("One character stop needed, not [$stop]");
		$ONECODE{'Ter'} = $stop;
   }
   if (defined $unknown) {
		length $unknown != 1 and $self->throw("One character stop needed, not [$unknown]");
		$ONECODE{'Xaa'} = $unknown;
   }

   my ($aas, $aa3);
   my $length = (length $string) - 2;
   for (my $i = 0 ; $i < $length ; $i += 3)  {
		$aa3 = substr($string, $i, 3);
		$aa3 = ucfirst(lc($aa3)); 
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
    my ($self, $seq, @args ) = @_;

    $self->throw('Object [$seq] '. 'of class ['. ref($seq).  ']  can not be translated.')
	unless $seq->can('translate');

    my ($stop, $unknown, $frame, $tableid, $fullCDS, $throw) = @args;
    my @seqs;
    my $f = 0;
    while ($f != 3) {
        my $translation = $seq->translate($stop, $unknown,$f,$tableid, $fullCDS, $throw );
	$translation->id($seq->id. "-". $f. "F");
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
    my ($self, $seq, @args ) = @_;

    my @seqs = $self->translate_3frames($seq, @args);
    my @seqs2 = $self->translate_3frames($seq->revcom, @args);
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

sub valid_aa{
   my ($self,$code) = @_;

   if( ! $code ) { 
       my @codes;
       foreach my $c ( sort values %ONECODE ) {
	   push @codes, $c unless ( $c =~ /[BZX\*]/ );
       }
       push @codes, qw(B Z X *); # so they are in correct order ?
       return @codes;
  }
   elsif( $code == 1 ) { 
       my @codes;
       foreach my $c ( sort keys %ONECODE ) {
	   push @codes, $c unless ( $c =~ /(Asx|Glx|Xaa|Ter)/ );
       }
       push @codes, ('Asx', 'Glx', 'Xaa', 'Ter' );
       return @codes;
   }
   elsif( $code == 2 ) { 
       my %codes = %ONECODE;
       foreach my $c ( keys %ONECODE ) {
	   my $aa = $ONECODE{$c};
	   $codes{$aa} = $c;
       }
       return %codes;
   } else {
       $self->warn("unrecognized code in ".ref($self)." method valid_aa()");
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
    my ($self, $seq, @mutations ) = @_;

    $self->throw('Object [$seq] '. 'of class ['. ref($seq).
                 '] should be a Bio::PrimarySeqI ')
	unless $seq->isa('Bio::PrimarySeqI');
    $self->throw('Object [$mutations[0]] '. 'of class ['. ref($mutations[0]).
                 '] should be a Bio::LiveSeq::Mutation')
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
  Usage   : my $catseq = Bio::SeqUtils->cat(@seqs)
  Function: Concatenates an array of Bio::Seq objects, using the first sequence
            as a target. Annotations and sequence features are copied over 
            from any additional objects. Adjusts the coordinates of copied 
            features.
  Returns : a boolean
  Args    : array of sequence objects

Note that annotations have no sequence locations. If you concatenate
sequences with the same annotations they will all be added.

=cut

sub cat {
    my ($self, $seq, @seqs) = @_;
    $self->throw('Object [$seq] '. 'of class ['. ref($seq).
                 '] should be a Bio::PrimarySeqI ')
        unless $seq->isa('Bio::PrimarySeqI');
    

    for my $catseq (@seqs) {
        $self->throw('Object [$catseq] '. 'of class ['. ref($catseq).
                     '] should be a Bio::PrimarySeqI ')
            unless $catseq->isa('Bio::PrimarySeqI');

        $self->throw('Trying to concatenate sequences with different alphabets: '.
                     $seq->display_id. '('. $seq->alphabet. ') and '. $catseq->display_id.
                     '('. $catseq->alphabet. ')')
            unless $catseq->alphabet eq $seq->alphabet;


        my $length=$seq->length;
        $seq->seq($seq->seq.$catseq->seq);

        # move annotations
        if ($seq->isa("Bio::AnnotatableI") and $catseq->isa("Bio::AnnotatableI")) {
            foreach my $key ( $catseq->annotation->get_all_annotation_keys() ) {

                foreach my $value ( $catseq->annotation->get_Annotations($key) ) {
                    $seq->annotation->add_Annotation($key, $value);
                }
            } 
        }
        
        # move SeqFeatures
        if ( $seq->isa('Bio::SeqI') and $catseq->isa('Bio::SeqI')) {
            for my $feat ($catseq->get_SeqFeatures) {
                $seq->add_SeqFeature($self->_coord_adjust($feat, $length));
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

sub trunc_with_features{
    use Bio::Range;
    my ($self,$seq,$start,$end) = @_;
    $self->throw('Object [$seq] '. 'of class ['. ref($seq).
                 '] should be a Bio::SeqI ')
    unless $seq->isa('Bio::SeqI');
    my $trunc=$seq->trunc($start, $end);
    my $truncrange=Bio::Range->new(-start=>$start, -end=>$end, -strand=>0);
    #move annotations
    foreach my $key ( $seq->annotation->get_all_annotation_keys() ) {
	foreach my $value ( $seq->annotation->get_Annotations($key) ) {
	    $trunc->annotation->add_Annotation($key, $value);
	}
    } 
    
    #move features
    $trunc->add_SeqFeature(grep {$_=$self->_coord_adjust($_, 1-$start, $end+1-$start) if $_->overlaps($truncrange)} $seq->get_SeqFeatures);
    return $trunc;
}



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
    my ($self, $feat, $add, $length)=@_;
    $self->throw('Object [$feat] '. 'of class ['. ref($feat).
                 '] should be a Bio::SeqFeatureI ')
        unless $feat->isa('Bio::SeqFeatureI');
    my @adjsubfeat;
    for my $subfeat ($feat->remove_SeqFeatures) {
        push @adjsubfeat, $self->_coord_adjust($subfeat, $add, $length);
    }
    my @loc;
    for ($feat->location->each_Location) {
        my @coords=($_->start, $_->end);
        my $strand=$_->strand;
	my $type=$_->location_type;
        map s/(\d+)/if ($add+$1<1) {'<1'} elsif (defined $length and $add+$1>$length) {">$length"} else {$add+$1}/ge, @coords;
	my($newstart,$newend)=@coords;
	unless ($type eq 'IN-BETWEEN') {
	    push @loc, Bio::Location::Fuzzy->new(-start=>$newstart,
						 -end=>$newend,
						 -strand=>$strand,
						 -location_type=>$type
						);
	} else {
	    push @loc, Bio::Location::Simple->new(-start=>$newstart,
					  -end=>$newend,
					  -strand=>$strand,
					  -location_type=>$type
					 );
	}
    }
    my $newfeat=Bio::SeqFeature::Generic->new(-primary=>$feat->primary_tag);
    foreach my $key ( $feat->annotation->get_all_annotation_keys() ) {
	foreach my $value ( $feat->annotation->get_Annotations($key) ) {
	    $newfeat->annotation->add_Annotation($key, $value);
	}
    } 
    foreach my $key ( $feat->get_all_tags() ) {
	    $newfeat->add_tag_value($key, $feat->get_tag_values($key));
    } 
    if (@loc==1) {
        $newfeat->location($loc[0])
    } else {
        my $loc=Bio::Location::Split->new;
        $loc->add_sub_Location(@loc);
        $newfeat->location($loc);
    }
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

sub revcom_with_features{
    my ($self,$seq) = @_;
    $self->throw('Object [$seq] '. 'of class ['. ref($seq).
                 '] should be a Bio::SeqI ')
    unless $seq->isa('Bio::SeqI');
    my $revcom=$seq->revcom;
    
    #move annotations
    foreach my $key ( $seq->annotation->get_all_annotation_keys() ) {
	foreach my $value ( $seq->annotation->get_Annotations($key) ) {
	    $revcom->annotation->add_Annotation($key, $value);
	}
    } 

    #move features
    $revcom->add_SeqFeature(map {$self->_feature_revcom($_, $seq->length)} reverse $seq->get_SeqFeatures);
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
    my ($self, $feat, $length)=@_;
    $self->throw('Object [$feat] '. 'of class ['. ref($feat).
                 '] should be a Bio::SeqFeatureI ')
        unless $feat->isa('Bio::SeqFeatureI');
    my @adjsubfeat;
    for my $subfeat ($feat->remove_SeqFeatures) {
        push @adjsubfeat, $self->_feature_revcom($subfeat, $length);
    }
    my @loc;
    for ($feat->location->each_Location) {
	my $type=$_->location_type;
        my $strand;
	if ($_->strand==-1) {$strand=1}
	elsif ($_->strand==1) {$strand=-1}
	else {$strand=$_->strand}
	my $newend=$self->_coord_revcom($_->start,
					$_->start_pos_type,
					$length);
	my $newstart=$self->_coord_revcom($_->end,
					  $_->end_pos_type,
					  $length);
	unless ($type eq 'IN-BETWEEN') {
	    push @loc, Bio::Location::Fuzzy->new(-start=>$newstart,
						 -end=>$newend,
						 -strand=>$strand,
						 -location_type=>$type
						);
	} else {
	    push @loc, Bio::Location::Simple->new(-start=>$newstart,
						  -end=>$newend,
						  -strand=>$strand,
						  -location_type=>$type
						 );
	}
    }
    my $newfeat=Bio::SeqFeature::Generic->new(-primary=>$feat->primary_tag);
    foreach my $key ( $feat->annotation->get_all_annotation_keys() ) {
	foreach my $value ( $feat->annotation->get_Annotations($key) ) {
	    $newfeat->annotation->add_Annotation($key, $value);
	}
    } 
    foreach my $key ( $feat->get_all_tags() ) {
	    $newfeat->add_tag_value($key, $feat->get_tag_values($key));
    } 
    if (@loc==1) {
        $newfeat->location($loc[0])
    } else {
        my $loc=Bio::Location::Split->new;
        $loc->add_sub_Location(@loc);
        $newfeat->location($loc);
    }
    $newfeat->add_SeqFeature($_) for @adjsubfeat;
    return $newfeat;
}

sub _coord_revcom {
    my ($self, $coord, $type, $length)=@_;
    if ($type eq 'BETWEEN' or $type eq 'WITHIN') {
	$coord=~s/(\d+)(.*)(\d+)/$length+1-$3.$2.$length+1-$1/ge;
    } else {
	$coord=~s/(\d+)/$length+1-$1/ge;
	$coord='>'.$coord if $type eq 'BEFORE';
	$coord='<'.$coord if $type eq 'AFTER';
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
    my ($self, $seq, $sim, $rate) = @_;
    $rate ||= 1;

    $self->throw('Object [$seq] '. 'of class ['. ref($seq).
                     '] should be a Bio::PrimarySeqI ')
            unless $seq->isa('Bio::PrimarySeqI');
    
    $self->throw("[$sim] ". ' should be a positive integer or float under 100')
            unless $sim =~ /^[+\d.]+$/ and $sim <= 100;

    $self->warn("Nucleotide sequences are 25% similar by chance.
        Do you really want to set similarity to [$sim]%?\n")
            unless $sim >25 ;

    $self->throw('Only nucleotide sequences are supported')
            if $seq->alphabet eq 'protein';


    # arrays of possible changes have transitions as first items
    my %changes;
    $changes{'a'} = ['t', 'c', 'g'];
    $changes{'t'} = ['a', 'c', 'g'];
    $changes{'c'} = ['g', 'a', 't'];
    $changes{'g'} = ['c', 'a', 't'];


    # given the desired rate, find out where cut off points need to be
    # when random numbers are generated from 0 to 100
    # we are ignoring identical mutations (e.g. A->A) to speed things up
    my $bin_size = 100/($rate + 2);  
    my $transition = 100 - (2*$bin_size);
    my $first_transversion = $transition + $bin_size;

    # unify the look of sequence strings
    my $string = lc $seq->seq; # lower case
    $string =~ s/u/t/; # simplyfy our life; modules should deal with the change anyway
    # store the original sequence string
    my $oristring = $string;
    my $length = $seq->length;

    while (1) {
        # find the location in the string to change
        my $loc = int (rand $length) + 1;


        # nucleotide to change
        my $oldnuc = substr $string, $loc-1, 1;
        my $newnuc;

        # nucleotide it is changed to
        my $choose = rand(100);
        if ($choose < $transition ) {
            $newnuc =  $changes{$oldnuc}[0];
        }
        elsif ($choose < $first_transversion ) {
            $newnuc =  $changes{$oldnuc}[1];
        } else {
            $newnuc =  $changes{$oldnuc}[2];
        }

        # do the change
        substr $string, $loc-1, 1 , $newnuc;

        $self->debug("$loc$oldnuc>$newnuc\n");

        # stop evolving if the limit has been reached
        last if $self->_get_similarity($oristring, $string) <= $sim;

    }

    return new Bio::PrimarySeq(-id => $seq->id. "-$sim",
                               -description => $seq->description,
                               -seq => $string
                              )
}


sub _get_similarity  {
    my ($self, $oriseq, $seq) = @_;

    my $len = length($oriseq);
    my $c;

    for (my $i = 0; $i< $len; $i++ ) {
        $c++ if substr($oriseq, $i, 1) eq substr($seq, $i, 1);
    }
    return 100 * $c/$len;
}

1;

