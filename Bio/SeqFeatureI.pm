# $Id$
#
# BioPerl module for Bio::SeqFeatureI
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeatureI - Abstract interface of a Sequence Feature

=head1 SYNOPSIS

    # get a seqfeature somehow, eg,

    foreach $feat ( $seq->top_SeqFeatures() ) {
            print "Feature from ", $feat->start, "to ", 
	          $feat->end, " Primary tag  ", $feat->primary_tag, 
	          ", produced by ", $feat->source_tag(), "\n";

            if( $feat->strand == 0 ) {
		print "Feature applicable to either strand\n";
            } else {
                print "Feature on strand ", $feat->strand,"\n"; # -1,1
            }

            foreach $tag ( $feat->all_tags() ) {
		print "Feature has tag ", $tag, "with values, ",
		      join(' ',$feat->each_tag_value($tag)), "\n";
            }
	    print "new feature\n" if $feat->has_tag('new');
	    # features can have sub features
	    my @subfeat = $feat->get_SeqFeatures();
	}

=head1 DESCRIPTION

This interface is the functions one can expect for any Sequence
Feature, whatever its implementation or whether it is a more complex
type (eg, a Gene). This object doesn\'t actually provide any
implemention, it just provides the definitions of what methods one can
call. See Bio::SeqFeature::Generic for a good standard implementation
of this object

=head1 FEEDBACK

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeatureI;
use vars qw(@ISA);
use strict;

use Bio::RangeI;
use Bio::Seq;

use Carp;

@ISA = qw(Bio::RangeI);

=head1 SeqFeatureI specific methods

New method interfaces.

=cut

=head2 get_SeqFeatures

 Title   : get_SeqFeatures
 Usage   : @feats = $feat->get_SeqFeatures();
 Function: Returns an array of sub Sequence Features
 Returns : An array
 Args    : none


=cut

sub get_SeqFeatures{
   my ($self,@args) = @_;

   $self->throw_not_implemented();
}

=head2 display_name

 Title   : display_name
 Usage   : $name = $feat->display_name()
 Function: Returns the human-readable name of the feature for displays.
 Returns : a string
 Args    : none

=cut

sub display_name { 
    shift->throw_not_implemented();
}

=head2 primary_tag

 Title   : primary_tag
 Usage   : $tag = $feat->primary_tag()
 Function: Returns the primary tag for a feature,
           eg 'exon'
 Returns : a string 
 Args    : none


=cut

sub primary_tag{
   my ($self,@args) = @_;

   $self->throw_not_implemented();

}

=head2 source_tag

 Title   : source_tag
 Usage   : $tag = $feat->source_tag()
 Function: Returns the source tag for a feature,
           eg, 'genscan' 
 Returns : a string 
 Args    : none


=cut

sub source_tag{
   my ($self,@args) = @_;

   $self->throw_not_implemented();
}

=head2 has_tag

 Title   : has_tag
 Usage   : $tag_exists = $self->has_tag('some_tag')
 Function: 
 Returns : TRUE if the specified tag exists, and FALSE otherwise
 Args    :


=cut

sub has_tag{
   my ($self,@args) = @_;

   $self->throw_not_implemented();

}

=head2 get_tag_values

 Title   : get_tag_values
 Usage   : @values = $self->get_tag_values('some_tag')
 Function: 
 Returns : An array comprising the values of the specified tag.
 Args    :


=cut

sub get_tag_values {
    shift->throw_not_implemented();
}

=head2 get_all_tags

 Title   : get_all_tags
 Usage   : @tags = $feat->get_all_tags()
 Function: gives all tags for this feature
 Returns : an array of strings
 Args    : none


=cut

sub get_all_tags{
    shift->throw_not_implemented();
}

=head2 attach_seq

 Title   : attach_seq
 Usage   : $sf->attach_seq($seq)
 Function: Attaches a Bio::Seq object to this feature. This
           Bio::Seq object is for the *entire* sequence: ie
           from 1 to 10000

           Note that it is not guaranteed that if you obtain a feature from
           an object in bioperl, it will have a sequence attached. Also,
           implementors of this interface can choose to provide an empty
           implementation of this method. I.e., there is also no guarantee 
           that if you do attach a sequence, seq() or entire_seq() will not
           return undef.

           The reason that this method is here on the interface is to enable
           you to call it on every SeqFeatureI compliant object, and
           that it will be implemented in a useful way and set to a useful 
           value for the great majority of use cases. Implementors who choose
           to ignore the call are encouraged to specifically state this in
           their documentation.

 Example :
 Returns : TRUE on success
 Args    : a Bio::PrimarySeqI compliant object


=cut

sub attach_seq {
    shift->throw_not_implemented();
}

=head2 seq

 Title   : seq
 Usage   : $tseq = $sf->seq()
 Function: returns the truncated sequence (if there is a sequence attached) 
           for this feature
 Example :
 Returns : sub seq (a Bio::PrimarySeqI compliant object) on attached sequence
           bounded by start & end, or undef if there is no sequence attached
 Args    : none


=cut

sub seq {
    shift->throw_not_implemented();
}

=head2 entire_seq

 Title   : entire_seq
 Usage   : $whole_seq = $sf->entire_seq()
 Function: gives the entire sequence that this seqfeature is attached to
 Example :
 Returns : a Bio::PrimarySeqI compliant object, or undef if there is no
           sequence attached
 Args    : none


=cut

sub entire_seq {
    shift->throw_not_implemented();
}


=head2 seq_id

 Title   : seq_id
 Usage   : $obj->seq_id($newval)
 Function: There are many cases when you make a feature that you
           do know the sequence name, but do not know its actual
           sequence. This is an attribute such that you can store
           the ID (e.g., display_id) of the sequence.

           This attribute should *not* be used in GFF dumping, as
           that should come from the collection in which the seq
           feature was found.
 Returns : value of seq_id
 Args    : newvalue (optional)


=cut

sub seq_id {
    shift->throw_not_implemented();
}

=head2 unique_id

 Title   : unique_id
 Usage   : $obj->unique_id($newval)
 Function: Get/set the unique-identify string for this feature
 Returns : value of seq_id
 Args    : newvalue (optional)


=cut

sub unique_id {
    shift->throw_not_implemented();
}

=head2 gff_string

 Title   : gff_string
 Usage   : $str = $feat->gff_string;
           $str = $feat->gff_string($gff_formatter);
 Function: Provides the feature information in GFF format.

           The implementation provided here returns GFF2 by default. If you
           want a different version, supply an object implementing a method
           gff_string() accepting a SeqFeatureI object as argument. E.g., to
           obtain GFF1 format, do the following:

                my $gffio = Bio::Tools::GFF->new(-gff_version => 1);
                $gff1str = $feat->gff_string($gff1io);

 Returns : A string
 Args    : Optionally, an object implementing gff_string().


=cut

sub gff_string{
   my ($self,$formatter) = @_;

   $formatter = $self->_static_gff_formatter unless $formatter;
   return $formatter->gff_string($self);
}

my $static_gff_formatter = undef;

=head2 _static_gff_formatter

 Title   : _static_gff_formatter
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _static_gff_formatter{
   my ($self,@args) = @_;

   if( !defined $static_gff_formatter ) {
       $static_gff_formatter = Bio::Tools::GFF->new('-gff_version' => 2);
   }
   return $static_gff_formatter;
}

=head1 Bio::RangeI methods

List of interfaces inherited from Bio::RangeI (see L<Bio::RangeI>
for details).

=cut

=head2 start

 Title   : start
 Usage   : $start = $feat->start
 Function: Returns the start coordinate of the feature
 Returns : integer
 Args    : none


=head2 end

 Title   : end
 Usage   : $end = $feat->end
 Function: Returns the end coordinate of the feature
 Returns : integer
 Args    : none

=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
 Function: Returns strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none


=cut

=head1 Decorating methods

These methods have an implementation provided by Bio::SeqFeatureI,
but can be validly overwritten by subclasses

=head2 spliced_seq

  Title   : spliced_seq

  Usage   : $seq = $feature->spliced_seq()
            $seq = $feature_with_remote_locations->spliced_seq($db_for_seqs)

  Function: Provides a sequence of the feature which is the most
            semantically "relevant" feature for this sequence. A default
            implementation is provided which for simple cases returns just
            the sequence, but for split cases, loops over the split location
            to return the sequence. In the case of split locations with
            remote locations, eg

            join(AB000123:5567-5589,80..1144)

            in the case when a database object is passed in, it will attempt
            to retrieve the sequence from the database object, and "Do the right thing",
            however if no database object is provided, it will generate the correct
            number of N's (DNA) or X's (protein, though this is unlikely).

            This function is deliberately "magical" attempting to second guess
            what a user wants as "the" sequence for this feature

            Implementing classes are free to override this method with their
            own magic if they have a better idea what the user wants

  Args    : [optional] A Bio::DB::RandomAccessI compliant object
  Returns : A Bio::Seq

=cut

sub spliced_seq {
    my ($self,$db) = shift;

    if( ! $self->location->isa("Bio::Location::SplitLocationI") ) {
	return $self->seq(); # nice and easy!
    }

    # redundant test, but the above ISA is probably not ideal.
    if( ! $self->location->isa("Bio::Location::SplitLocationI") ) {
	$self->throw("not atomic, not split, yikes, in trouble!");
    }

    my $seqstr;
    my $seqid = $self->entire_seq->display_id;
    # This is to deal with reverse strand features
    # so we are really sorting features 5' -> 3' on their strand
    # i.e. rev strand features will be sorted largest to smallest
    # as this how revcom CDSes seem to be annotated in genbank.
    # Might need to eventually allow this to be programable?    
    # (can I mention how much fun this is NOT! --jason)
    
    my ($mixed,$fstrand) = (0);
    if( $self->isa('Bio::Das::SegmentI') &&
	! $self->absolute ) { 
	$self->warn("Calling spliced_seq with a Bio::Das::SegmentI which does have absolute set to 1 -- be warned you may not be getting things on the correct strand");
    }
    
    my @locs = map { $_->[0] }
    # sort so that most negative is first basically to order
    # the features on the opposite strand 5'->3' on their strand
    # rather than they way most are input which is on the fwd strand

    sort { $a->[1] <=> $b->[1] } # Yes Tim, Schwartzian transformation
    map { 
	$fstrand = $_->strand unless defined $fstrand;
	$mixed = 1 if defined $_->strand && $fstrand != $_->strand;
	[ $_, $_->start* ($_->strand || 1)];	    
    } $self->location->each_Location; 
    
    if ( $mixed ) { 
	$self->warn("Mixed strand locations, spliced seq using the input order rather than trying to sort");    
	@locs = $self->location->each_Location; 
    }

    foreach my $loc ( @locs  ) {
	if( ! $loc->isa("Bio::Location::Atomic") ) {
	    $self->throw("Can only deal with one level deep locations");
	}
	my $called_seq;
	if( $fstrand != $loc->strand ) {
	    $self->warn("feature strand is different from location strand!");
	}
	# deal with remote sequences

	if( $loc->seq_id ne $seqid ) {
	    if( defined $db ) {
		my $sid = $loc->seq_id;
		$sid =~ s/\.\d+//g;
		eval {
		    $called_seq = $db->get_Seq_by_acc($sid);
		};
		if( $@ ) {
		    $self->warn("In attempting to join a remote location, sequence $sid was not in database. Will provide padding N's. Full exception \n\n$@");
		    $called_seq = undef;
		}
	    } else {
		$called_seq = undef;
	    }
	    if( !defined $called_seq ) {
		$seqstr .= 'N' x $self->length;
		next;
	    }
	} else {
	    $called_seq = $self->entire_seq;
	}
	
	if( $self->isa('Bio::Das::SegmentI') ) {
	    my ($s,$e) = ($loc->start,$loc->end);	    
	    $seqstr .= $called_seq->subseq($s,$e)->seq();
	} else { 
	    # This is dumb subseq should work on locations...
	    if( $loc->strand == 1 ) {
		$seqstr .= $called_seq->subseq($loc->start,$loc->end);
	    } else {
		$seqstr .= $called_seq->trunc($loc->start,$loc->end)->revcom->seq();
	    }
	}
    }
    my $out = Bio::Seq->new( -id => $self->entire_seq->display_id . "_spliced_feat",
				      -seq => $seqstr);
    
    return $out;
}

=head1 RangeI methods

These methods are inherited from RangeI and can be used
directly from a SeqFeatureI interface. Remember that a 
SeqFeature is-a RangeI, and so wherever you see RangeI you
can use a feature ($r in the below documentation).

=head2 overlaps

  Title   : overlaps
  Usage   : if($feat->overlaps($r)) { do stuff }
            if($feat->overlaps(200)) { do stuff }
  Function: tests if $feat overlaps $r
  Args    : a RangeI to test for overlap with, or a point
  Returns : true if the Range overlaps with the feature, false otherwise


=head2 contains

  Title   : contains
  Usage   : if($feat->contains($r) { do stuff }
  Function: tests whether $feat totally contains $r
  Args    : a RangeI to test for being contained
  Returns : true if the argument is totaly contained within this range


=head2 equals

  Title   : equals
  Usage   : if($feat->equals($r))
  Function: test whether $feat has the same start, end, strand as $r
  Args    : a RangeI to test for equality
  Returns : true if they are describing the same range


=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
triplets (start, stop, strand) from which new ranges could be built.

=head2 intersection

  Title   : intersection
  Usage   : ($start, $stop, $strand) = $feat->intersection($r)
  Function: gives the range that is contained by both ranges
  Args    : a RangeI to compare this one to
  Returns : nothing if they do not overlap, or the range that they do overlap

=head2 union

  Title   : union
  Usage   : ($start, $stop, $strand) = $feat->union($r);
          : ($start, $stop, $strand) = Bio::RangeI->union(@ranges);
  Function: finds the minimal range that contains all of the ranges
  Args    : a range or list of ranges to find the union of
  Returns : the range containing all of the ranges

=cut

=head2 location

 Title   : location
 Usage   : my $location = $seqfeature->location()
 Function: returns a location object suitable for identifying location 
	   of feature on sequence or parent feature  
 Returns : Bio::LocationI object
 Args    : none


=cut

sub location {
   my ($self) = @_;

   $self->throw_not_implemented();
}


1;
