# $Id$
#
# bioperl module for Bio::PrimaryQual
#
# Cared for by Chad Matsalla <bioinformatics@dieselwurks.com>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::PrimaryQual - Bioperl lightweight Quality Object

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bio.perl.org/MailList.html - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Chad Matsalla

Email bioinformatics@dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Seq::PrimaryQual;
use vars qw(@ISA %valid_type);
use strict;

use Bio::Root::RootI;
use Bio::Seq::QualI;
@ISA = qw(Bio::Root::RootI Bio::Seq::QualI);


=head2 new()

 Title   : new()
 Usage   : $qual = Bio::Seq::PrimaryQual->new( -qual => '10 20 30 40 50 50 20 10',
                                           -id  => 'human_id',
                                           -accession_number => 'AL000012',
                                           );

 Function: Returns a new Bio::Seq::PrimaryQual object from basic 
	constructors, being a string _or_ a reference to an array for the
	sequence and strings for id and accession_number. Note that you
	can provide an empty quality string.
 Returns : a new Bio::Seq::PrimaryQual object

=cut




sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my($qual,$id,$acc,$pid,$desc,$given_id) =
        $self->_rearrange([qw(QUAL
                              DISPLAY_ID
                              ACCESSION_NUMBER
                              PRIMARY_ID
                              DESC
                              ID
                              )],
                          @args);
 
    if( defined $id && defined $given_id ) {
        if( $id ne $given_id ) {
            $self->throw("Provided both id and display_id constructor functions. [$id] [$given_id]");   
        }
    }
    if( defined $given_id ) { $id = $given_id; }

    # if moltype is provided we set it first, so that it won't be guessed
    # when the sequence is set

	# what is the best way to deal with moltype?
	# IGNORE it. Bwahahahaha.
	# $moltype && $self->moltype($moltype);

	    # note: the sequence string may be empty

    $self->qual($qual);
	# if defined($qual);
    $id      && $self->display_id($id);
    $acc     && $self->accession_number($acc);
    $pid     && $self->primary_id($pid);
    $desc    && $self->desc($desc);

    return $self;
}



sub moltype {
	return undef;
}

=head2 qual()

 Title   : qual()
 Usage   : @quality_values  = @{$obj->qual()};
 Function: Returns the quality as a reference to an array containing the
        quality values. The individual elements of the quality array are
        not validated and can be any numeric value.
 Returns : A reference to an array.
 Status  :

=cut

sub qual {
	my ($obj,$value) = @_;
	if( defined $value) {
		if (!$value) {
			$obj->warn("Setting quality to an empty value.");
			my @foo = [];
			$obj->{'qual'} = \@foo;
		}
		elsif(! $obj->validate_qual($value)) {
			$obj->throw("Attempting to set the quality to [$value] which does not look healthy");
		}
		elsif (ref($value) eq "ARRAY") {
			# if the user passed in a reference to an array
			$obj->{'qual'} = $value;
		}
		else {
			my @quality_array = split/\s+/,$value;
			$obj->{'qual'} = \@quality_array;
		}
			# if(($is_changed_qual && (CORE::length($value) > 0)) ||
			# (! defined($obj->moltype()))) {
			#	$obj->_guess_type();
			# }
	}
   return $obj->{'qual'};
}

=head2 validate_qual($qualstring)

 Title	 : validate_qual($qualstring)
 Usage	 : print("Valid.") if { &validate_qual($self,$qualities); }
 Function: Make sure that the quality, if it has length > 0, contains at
	least one digit. Note that quality strings are parsed into arrays
	using split/\d+/,$quality_string, so make sure that your quality
	scalar looks like this if you want it to be parsed properly.
 Returns : 1 for a valid sequence (WHY? Shouldn't it return 0? <boggle>)
 Args    : a scalar (any scalar, why PrimarySeq author?) and a scalar
	containing the string to validate.

=cut

sub validate_qual {
	# how do I validate quality values?
	# \d+\s+\d+..., I suppose
   my ($self,$qualstr) = @_;
	# why the CORE??
	if (CORE::length($qualstr) == 0) {
		return 0;
	}
   if((CORE::length($qualstr) > 0) && $qualstr =~ /\d/) {
	# print("The length of this thing is greater then zero and contains at least one digit <boggle>!\n");
       return 1;
   }
   return 0;
}

=head2 subqual($start,$end)

 Title   : subqual($start,$end)
 Usage   : @subset_of_quality_values = @{$obj->subqual(10,40)};
 Function: returns the quality values from $start to $end, where the
        first value is 1 and the number is inclusive, ie 1-2 are the first
        two bases of the sequence. Start cannot be larger than end but can
        be equal.
 Returns : A reference to an array.
 Args    : a start position and an end position

=cut


sub subqual {
   my ($self,$start,$end) = @_;

   if( $start > $end ){
       $self->throw("in subqual, start [$start] has to be greater than end [$end]");
   }

   if( $start <= 0 || $end > $self->length ) {
       $self->throw("You have to have start positive and length less than the total length of sequence [$start:$end] Total ".$self->length."");
   }

   # remove one from start, and then length is end-start

   $start--;
	$end--;
	my @sub_qual_array = @{$self->{qual}}[$start..$end];

 	#   return substr $self->seq(), $start, ($end-$start);
	return \@sub_qual_array;

}

=head2 display_id()

 Title   : display_id()
 Usage   : $id_string = $obj->display_id();
 Function: returns the display id, aka the common name of the Quality
        object.
        The semantics of this is that it is the most likely string to be
        used as an identifier of the quality sequence, and likely to have
        "human" readability.  The id is equivalent to the ID field of the
        GenBank/EMBL databanks and the id field of the Swissprot/sptrembl
        database. In fasta format, the >(\S+) is presumed to be the id,
        though some people overload the id to embed other information.
        Bioperl does not use any embedded information in the ID field,
        and people are encouraged to use other mechanisms (accession field
        for example, or extending the sequence object) to solve this.
        Notice that $seq->id() maps to this function, mainly for
        legacy/convience issues
 Returns : A string
 Args    : None
 Status  : Virtual

=cut

sub display_id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'display_id'} = $value;
    }
    return $obj->{'display_id'};

}

=head2 accession_number()

 Title   : accession_number()
 Usage   : $unique_biological_key = $obj->accession_number();
 Function: Returns the unique biological id for a sequence, commonly
        called the accession_number. For sequences from established
        databases, the implementors should try to use the correct
        accession number. Notice that primary_id() provides the unique id
        for the implemetation, allowing multiple objects to have the same
        accession number in a particular implementation. For sequences
        with no accession number, this method should return "unknown".
 Returns : A string
 Args    : None
 Status  : Virtual


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

=head2 primary_id()

 Title   : primary_id()
 Usage   : $unique_implementation_key = $obj->primary_id();
 Function: Returns the unique id for this object in this implementation.
        This allows implementations to manage their own object ids in a
        way the implementaiton can control clients can expect one id to
        map to one object. For sequences with no accession number, this
        method should return a stringified memory location.
 Returns : A string
 Args    : None

=cut

sub primary_id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'primary_id'} = $value;
    }
   if( ! exists $obj->{'primary_id'} ) {
       return "$obj";
   }
   return $obj->{'primary_id'};

}

=head2 desc()

 Title   : desc()
 Usage   : $qual->desc($newval);
           $description = $qual->desc();
 Function: Get/set description text for a qual object
 Example :
 Returns : Value of desc
 Args    : newvalue (optional)

=cut

sub desc {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'desc'} = $value;
    }
    return $obj->{'desc'};

}

=head2 id()

 Title   : id()
 Usage   : $id = $qual->id();
 Function: Return the ID of the quality. This should normally be (and
	actually is in the implementation provided here) just a synonym
	for display_id().
 Returns : A string.
 Args    : None.

=cut

sub id {
   my ($self,$value) = @_;
   if( defined $value ) {
        return $self->display_id($value);
   }
   return $self->display_id();
}

=head2 length()

 Title	 : length()
 Usage	 : $length = $qual->length();
 Function: Return the length of the array holding the quality values.
	Under most circumstances, this should match the number of quality
	values but no validation is done when the PrimaryQual object is
	constructed and non-digits could be put into this array. Is this a
	bug? Just enough rope...
 Returns : A scalar (the number of elements in the quality array).
 Args    : None.

=cut

sub length {
	my $self = shift;
	if (ref($self->{qual}) ne "ARRAY") {
		$self->warn("{qual} is not an array here. Why? It appears to be ".ref($self->{qual})."(".$self->{qual}."). Good thing this can _never_ happen.");
		::dumpValue($self);
	}
	return scalar(@{$self->{qual}});
}

=head2 qualat($position)

 Title   : qualat($position)
 Usage   : $quality = $obj->qualat(10);
 Function: Return the quality value at the given location, where the
        first value is 1 and the number is inclusive, ie 1-2 are the first
        two bases of the sequence. Start cannot be larger than end but can
        be equal.
 Returns : A scalar.
 Args    : A position.

=cut

sub qualat {
        my ($self,$val) = @_;
        my @qualat = @{$self->subqual($val,$val)};
        if (scalar(@qualat) == 1) {
                return $qualat[0];
        }
        else {
                $self->throw("AAAH! qualat provided more then one quality.");
        }
} 

1;
