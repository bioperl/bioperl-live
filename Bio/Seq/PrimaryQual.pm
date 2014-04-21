#
# bioperl module for Bio::PrimaryQual
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

 use Bio::Seq::PrimaryQual;

 # you can use either a space-delimited string for quality

 my $string_quals = "10 20 30 40 50 40 30 20 10";
 my $qualobj = Bio::Seq::PrimaryQual->new(
     -qual             => $string_quals,
     -id               => 'QualityFragment-12',
     -accession_number => 'X78121',
 );

 # _or_ you can use an array of quality values

 my @q2 = split/ /,$string_quals;
 $qualobj = Bio::Seq::PrimaryQual->new(
     -qual              => \@q2,
     -primary_id        => 'chads primary_id',
     -desc              => 'chads desc',
     -accession_number  => 'chads accession_number',
     -id                => 'chads id'
 );

 # to get the quality values out:

 my @quals = @{$qualobj->qual()};

 # to give _new_ quality values

 my $newqualstring = "50 90 1000 20 12 0 0";
 $qualobj->qual($newqualstring);


=head1 DESCRIPTION

This module provides a mechanism for storing quality
values. Much more useful as part of
Bio::Seq::SeqWithQuality where these quality values
are associated with the sequence information.

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

=head1 AUTHOR - Chad Matsalla

Email bioinformatics@dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


package Bio::Seq::PrimaryQual;

use strict;

use base qw(Bio::Root::Root Bio::Seq::QualI);

our $MATCHPATTERN = '0-9eE\.\s+-';


=head2 new()

 Title   : new()
 Usage   : $qual = Bio::Seq::PrimaryQual->new(
               -qual             => '10 20 30 40 50 50 20 10',
               -id               => 'human_id',
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

    # default: turn ON the warnings (duh)
    my($qual,$id,$acc,$pid,$desc,$given_id,$header) =
        $self->_rearrange([qw(QUAL
                              DISPLAY_ID
                              ACCESSION_NUMBER
                              PRIMARY_ID
                              DESC
                              ID
                              HEADER
                              )],
                          @args);
    if( defined $id && defined $given_id ) {
        if( $id ne $given_id ) {
            $self->throw("Provided both id and display_id constructor functions. [$id] [$given_id]");   
        }
    }
    if( defined $given_id ) { $id = $given_id; }
    
    # note: the sequence string may be empty
    $self->qual(defined($qual) ? $qual : []);
     $header && $self->header($header);
    $id      && $self->display_id($id);
    $acc     && $self->accession_number($acc);
    $pid     && $self->primary_id($pid);
    $desc    && $self->desc($desc);

    return $self;
}


=head2 qual()

 Title   : qual()
 Usage   : @quality_values  = @{$obj->qual()};
 Function: Get or set the quality as a reference to an array containing the
           quality values. An error is generated if the quality scores are
           invalid, see validate_qual().
 Returns : A reference to an array.

=cut

sub qual {
    my ($self,$value) = @_;

    if( ! defined $value || length($value) == 0 ) { 
        $self->{'qual'} ||= [];
    } elsif( ref($value) =~ /ARRAY/i ) {
        # if the user passed in a reference to an array
        $self->{'qual'} = $value;
    } else {
        $self->validate_qual($value, 1);
        $value =~ s/^\s+//;
        $self->{'qual'} = [split(/\s+/,$value)];
    }
    
    return $self->{'qual'};
}


=head2 seq()

 Title    : seq()
 Usager   : $sequence = $obj->seq();
 Function : Returns the quality numbers as a space-separated string.
 Returns  : Single string.
 Args     : None.

=cut

sub seq {
    return join ' ', @{ shift->qual };
}


=head2 validate_qual($qualstring)

 Title   : validate_qual($qualstring)
 Usage   : print("Valid.") if { &validate_qual($self, $quality_string); }
 Function: Test that the given quality string is valid. It is expected to
           contain space-delimited numbers that can be parsed using split /\d+/.
           However, this validation takes shortcuts and only tests that the
           string contains characters valid in numbers: 0-9 . eE +-
           Note that empty quality strings are valid too.
 Returns : 1 for a valid sequence, 0 otherwise
 Args    : - Scalar containing the quality string to validate.
           - Boolean to optionally throw an error if validation failed

=cut

sub validate_qual {
    my ($self, $qualstr, $throw) = @_;
    if ( (defined $qualstr                ) &&
         ($qualstr !~ /^[$MATCHPATTERN]*$/) ) {
        if ($throw) {
            $self->throw("Failed validation of quality score from  '".
               (defined($self->id)||'[unidentified sequence]')."'. No numeric ".
               "value found.\n");
        }
        return 0;
    }
    return 1;
}


=head2 subqual($start,$end)

 Title   : subqual($start,$end)
 Usage   : @subset_of_quality_values = @{$obj->subqual(10,40)};
 Function: returns the quality values from $start to $end, where the
           first value is 1 and the number is inclusive, ie 1-2 are the
           first two bases of the sequence. Start cannot be larger than
           end but can be equal.
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
           and people are encouraged to use other mechanisms (accession
           field for example, or extending the sequence object) to solve
           this. Notice that $seq->id() maps to this function, mainly for
           legacy/convience issues
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


=head2 header()

 Title   : header()
 Usage   : $header = $obj->header();
 Function: Get/set the header that the user wants printed for this
     quality object.
 Returns : A string
 Args    : None

=cut

sub header {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'header'} = $value;
    }
    return $obj->{'header'};

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

 Title   : length()
 Usage   : $length = $qual->length();
 Function: Return the length of the array holding the quality values.
           Under most circumstances, this should match the number of quality
           values but no validation is done when the PrimaryQual object is
           constructed and non-digits could be put into this array. Is this
           a bug? Just enough rope...
 Returns : A scalar (the number of elements in the quality array).
 Args    : None.

=cut

sub length {
    my $self = shift;
    if (ref($self->{qual}) ne "ARRAY") {
        $self->warn("{qual} is not an array here. Why? It appears to be ".ref($self->{qual})."(".$self->{qual}."). Good thing this can _never_ happen.");
    }
    return scalar(@{$self->{qual}});
}


=head2 qualat()

 Title   : qualat
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
    } else {
        $self->throw("qualat() provided more than one quality.");
    }
} 

=head2 to_string()

 Title   : to_string()
 Usage   : $quality = $obj->to_string();
 Function: Return a textual representation of what the object contains.
           For this module, this function will return:
                qual
                display_id
                accession_number
                primary_id
                desc
                id
                length
 Returns : A scalar.
 Args    : None.

=cut

sub to_string {
    my ($self,$out,$result) = shift;
    $out = "qual: ".join(',',@{$self->qual()});
    foreach (qw(display_id accession_number primary_id desc id length)) {
        $result = $self->$_();
        if (!$result) { $result = "<unset>"; }
        $out .= "$_: $result\n";
    }
    return $out;
}


sub to_string_automatic {
    my ($self,$sub_result,$out) = shift;
    foreach (sort keys %$self) {
        print("Working on $_\n");
        eval { $self->$_(); };
        if ($@) { $sub_result = ref($_); }
        elsif (!($sub_result = $self->$_())) {
            $sub_result = "<unset>";
        }
        if (ref($sub_result) eq "ARRAY") {
            print("This thing ($_) is an array!\n");
            $sub_result = join(',',@$sub_result);
        }
        $out .= "$_: ".$sub_result."\n";
    }
    return $out;
} 

1;
