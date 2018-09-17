#
# BioPerl module for Bio::Seq::Meta
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::Meta - Generic superclass for sequence objects with
residue-based meta information

=head1 SYNOPSIS

  use Bio::LocatableSeq;
  use Bio::Seq::Meta;
  use Bio::Tools::OddCodes;
  use Bio::SeqIO;

  my $seq = Bio::Seq::Meta->new(-id=>'test',
                                   -seq=>'ACTGCTAGCT',
                                   -start=>2434,
                                   -end=>2443,
                                   -strand=>1,
                                   -verbose=>1, # to see warnings
                                  );

  # the existing sequence object can be a Bio::PrimarySeq, too

  # to test this is a meta seq object
  $seq->isa("Bio::Seq::Meta")
      || $seq->throw("$seq is not a Bio::Seq::Meta");


  $seq->meta('1234567890');
  $seq = Bio::Seq::Meta->new(-id=>'test',
                             -seq=>'HACILMIFGT',
                             -start=>2434,
                             -end=>2443,
                             -strand=>1,
                             -meta=>'1234567890',
                             -verbose=>1, # to see warnings
                            );

  # accessors
  $string     = $seq->meta_text();
  $substring  = $seq->submeta_text(2,5);
  $unique_key = $seq->accession_number();

  # storing output from Bio::Tools::OddCodes as meta data
  my $protcodes = Bio::Tools::OddCodes->new(-seq => $seq);
  my @codes = qw(structural chemical functional charge hydrophobic);
  map { $seq->named_meta($_, ${$protcodes->$_($seq) } )} @codes;

  my $out = Bio::SeqIO->new(-format=>'metafasta');
  $out->write_seq($seq);

=head1 DESCRIPTION

This class implements generic methods for sequences with residue-based
meta information. Meta sequences with meta data are Bio::LocatableSeq
objects with additional methods to store that meta information. See
L<Bio::LocatableSeq> and L<Bio::Seq::MetaI>.

The meta information in this class is always one character per residue
long and blank values are space characters (ASCII 32).

After the latest rewrite, the meta information no longer covers all
the residues automatically. Methods to check the length of meta
information (L<meta_length>)and to see if the ends are flushed to the
sequence have been added (L<is_flush>). To force the old
functionality, set L<force_flush> to true.

It is assumed that meta data values do not depend on the nucleotide
sequence strand value.

Application specific implementations should inherit from this class to
override and add to these methods.

L<Bio::Seq::Meta::Array> allows for more complex meta values (scalars
or objects) to be used.

=head2 Method naming

Character based meta data is read and set by method meta() and its
variants. These are the suffixes and prefixes used in the variants:

    [named_] [sub] meta [_text]

=over 3

=item _text

Suffix B<_text> guaranties that output is a string. Note that it does
not limit the input.

In this implementation, the output is always text, so these methods
are redundant.

=item sub

Prefix B<sub>, like in subseq(), means that the method applies to sub
region of the sequence range and takes start and end as arguments.
Unlike subseq(), these methods are able to set values.  If the range
is not defined, it defaults to the complete sequence.

=item named

Prefix B<named_> in method names allows the used to attach multiple
meta strings to one sequence by explicitly naming them. The name is
always the first argument to the method. The "unnamed" methods use the
class wide default name for the meta data and are thus special cases
"named" methods.

Note that internally names are keys in a hash and any misspelling of a
name will silently store the data under a wrong name. The used names
(keys) can be retrieved using method meta_names(). See L<meta_names>.

=back

=head1 NOTE

This Bio::Seq::MetaI implementation inherits from Bio::LocatableSeq, which
itself inherits from Bio::PrimarySeq. It is not a Bio::SeqI, so bless-ing
objects of this class into a Bio::SeqI or vice versa and will not work as
expected (see bug 2262). This may be addressed in a future refactor of
Bio::LocatableSeq.


=head1 SEE ALSO

L<Bio::LocatableSeq>, 
L<Bio::Seq::MetaI>, 
L<Bio::Seq::Meta::Array>

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

Email heikki-at-bioperl-dot-org

=head1 CONTRIBUTORS

Chad Matsalla, bioinformatics@dieselwurks.com

Aaron Mackey, amackey@virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::Meta;
use vars qw($DEFAULT_NAME $GAP $META_GAP);
use strict;

#use overload '""' => \&to_string;

use base qw(Bio::LocatableSeq Bio::Seq::MetaI);


BEGIN {

    $DEFAULT_NAME = 'DEFAULT';
    $GAP = '-';
    $META_GAP = ' ';
}

=head2 new

 Title   : new
 Usage   : $metaseq = Bio::Seq::Meta->new
	        ( -meta => 'aaaaaaaabbbbbbbb',
                  -seq =>  'TKLMILVSHIVILSRM'
	          -id  => 'human_id',
	          -accession_number => 'S000012',
	        );
 Function: Constructor for Bio::Seq::Meta class, meta data being in a
           string. Note that you can provide an empty quality string.
 Returns : a new Bio::Seq::Meta object

=cut


sub new {
    my ($class, @args) = @_;

    my $self = $class->SUPER::new(@args);

    my($meta, $forceflush, $nm) =
        $self->_rearrange([qw(META
                              FORCE_FLUSH
                              NAMED_META)],
                          @args);

    #$self->{'_meta'} = {};
    $self->{'_meta'}->{$DEFAULT_NAME} = "";

    $meta && $self->meta($meta);
    if ($nm && ref($nm) eq 'HASH') {
        while (my ($name, $meta) = each %$nm) {
            $self->named_meta($name, $meta);
        }
    }
    $forceflush && $self->force_flush($forceflush);

    return $self;
}


=head2 meta

 Title   : meta
 Usage   : $meta_values  = $obj->meta($values_string);
 Function:

           Get and set method for the meta data starting from residue
           position one. Since it is dependent on the length of the
           sequence, it needs to be manipulated after the sequence.

           The length of the returned value always matches the length
           of the sequence, if force_flush() is set. See L<force_flush>.

 Returns : meta data in a string
 Args    : new value, string, optional

=cut

sub meta {
   shift->named_meta($DEFAULT_NAME, shift);
}

=head2 meta_text

 Title   : meta_text
 Usage   : $meta_values  = $obj->meta_text($values_arrayref);
 Function: Variant of meta() guarantied to return a textual
           representation  of meta data. For details, see L<meta>.
 Returns : a string
 Args    : new value, optional

=cut

sub meta_text {
    shift->meta(shift);
}

=head2 named_meta

 Title   : named_meta()
 Usage   : $meta_values  = $obj->named_meta($name, $values_arrayref);
 Function: A more general version of meta(). Each meta data set needs
           to be named. See also L<meta_names>.
 Returns : a string
 Args    : scalar, name of the meta data set
           new value, optional

=cut

sub named_meta {
   my ($self, $name, $value) = @_;

   $name ||= $DEFAULT_NAME;
   if( defined $value) {

       $self->throw("I need a scalar value, not [". ref($value). "]")
	   if ref($value);

       # test for length
       my $diff = $self->length - CORE::length($value);
       if ($diff > 0) {
           $value .= (" " x $diff);
       }

       $self->{'_meta'}->{$name} = $value;

       #$self->_test_gap_positions($name) if $self->verbose > 0;
   }

   return " " x $self->length 
    if $self->force_flush && not defined $self->{'_meta'}->{$name};


   $self->_do_flush if $self->force_flush;

   return $self->{'_meta'}->{$name};
}

=head2 _test_gap_positions

 Title   : _test_gap_positions
 Usage   : $meta_values  = $obj->_test_gap_positions($name);
 Function: Internal test for correct position of gap characters.
           Gap being only '-' this time.

           This method is called from named_meta() when setting meta
           data but only if verbose is positive as this can be an
           expensive process on very long sequences. Set verbose(1) to
           see warnings when gaps do not align in sequence and meta
           data and turn them into errors by setting verbose(2).

 Returns : true on success, prints warnings
 Args    : none

=cut

sub _test_gap_positions {
    my $self = shift;
    my $name = shift;
    my $success = 1;

    $self->seq || return $success;
    my $len = CORE::length($self->seq);
    for (my $i=0; $i < $len; $i++) {
        my $s = substr $self->{seq}, $i, 1;
        my $m = substr $self->{_meta}->{$name}, $i, 1;
        $self->warn("Gap mismatch [$m/$s] in column [". ($i+1). "] of [$name] meta data in seq [". $self->id. "]")
            and $success = 0
                if ($s eq $META_GAP) && $s ne $m;
    }
    return $success;
}

=head2 named_meta_text

 Title   : named_meta_text()
 Usage   : $meta_values  = $obj->named_meta_text($name, $values_arrayref);
 Function: Variant of named_meta() guarantied to return a textual
           representation  of the named meta data.
           For details, see L<meta>.
 Returns : a string
 Args    : scalar, name of the meta data set
           new value, optional

=cut

sub named_meta_text {
    shift->named_meta(@_);
}

=head2 submeta

 Title   : submeta
 Usage   : $subset_of_meta_values = $obj->submeta(10, 20, $value_string);
           $subset_of_meta_values = $obj->submeta(10, undef, $value_string);
 Function:

           Get and set method for meta data for subsequences.

           Numbering starts from 1 and the number is inclusive, ie 1-2
           are the first two residue of the sequence. Start cannot be
           larger than end but can be equal.

           If the second argument is missing the returned values
           should extend to the end of the sequence.

           The return value may be a string or an array reference,
           depending on the implementation. If in doubt, use
           submeta_text() which is a variant guarantied to return a
           string.  See L<submeta_text>.

 Returns : A reference to an array or a string
 Args    : integer, start position
           integer, end position, optional when a third argument present
           new value, optional

=cut

sub submeta {
   shift->named_submeta($DEFAULT_NAME, @_);
}

=head2 submeta_text

 Title   : submeta_text
 Usage   : $meta_values  = $obj->submeta_text(20, $value_string);
 Function: Variant of submeta() guarantied to return a textual 
           representation  of meta data. For details, see L<meta>.
 Returns : a string
 Args    : new value, optional


=cut

sub submeta_text {
    shift->submeta(@_);
}

=head2 named_submeta

 Title   : named_submeta
 Usage   : $subset_of_meta_values = $obj->named_submeta($name, 10, 20, $value_string);
           $subset_of_meta_values = $obj->named_submeta($name, 10);
 Function: Variant of submeta() guarantied to return a textual
           representation  of meta data. For details, see L<meta>.
 Returns : A reference to an array or a string
 Args    : scalar, name of the meta data set
           integer, start position
           integer, end position, optional when a third argument present
           new value, optional

=cut

sub named_submeta {
    my ($self, $name, $start, $end, $value) = @_;

    $name ||= $DEFAULT_NAME;
    $start ||=1;


    $start =~ /^[+]?\d+$/ and $start > 0 or
        $self->throw("Need at least a positive integer start value");

    if ($value) {
        $end ||= $start+length($value)-1;
        $self->warn("You are setting meta values beyond the length of the sequence\n".
                    "[$start > ". length($self->seq)."] in sequence ". $self->id)
            if $start > length $self->seq;

        # pad meta data if needed
        $self->{_meta}->{$name} = () unless defined $self->{_meta}->{$name};
        if (length($self->{_meta}->{$name}) < $start) {
            $self->{'_meta'}->{$name} .=  " " x ( $start - length($self->{'_meta'}->{$name}) -1);
        }

        my $tail = '';
        $tail = substr ($self->{_meta}->{$name}, $start-1+length($value))
            if length($self->{_meta}->{$name}) >= $start-1+length($value);
        
        substr ($self->{_meta}->{$name}, --$start) = $value;
        $self->{_meta}->{$name} .= $tail;

        return substr ($self->{_meta}->{$name}, $start, $end - $start + 1);

    } else {

        $end or $end = length $self->seq;

        # pad meta data if needed
        if (length($self->{_meta}->{$name}) < $end) {
            $self->{'_meta'}->{$name} .=  " " x ( $start - length($self->{'_meta'}->{$name}));
        }

        return substr ($self->{_meta}->{$name}, $start-1, $end - $start + 1)
    }
}


=head2 named_submeta_text

 Title   : named_submeta_text
 Usage   : $meta_values  = $obj->named_submeta_text($name, 20, $value_string);
 Function: Variant of submeta() guarantied to return a textual
           representation  of meta data. For details, see L<meta>.
 Returns : a string
 Args    : scalar, name of the meta data
 Args    : integer, start position, optional
           integer, end position, optional
           new value, optional

=cut

sub named_submeta_text {
    shift->named_submeta(@_);
}

=head2 meta_names

 Title   : meta_names
 Usage   : @meta_names  = $obj->meta_names()
 Function: Retrieves an array of meta data set names. The default
           (unnamed) set name is guarantied to be the first name.
 Returns : an array of names
 Args    : none

=cut

sub meta_names {
    my ($self) = @_;

    my @r;
    foreach  ( sort keys %{$self->{'_meta'}} ) {
        push (@r, $_) unless $_ eq $DEFAULT_NAME;
    }
    unshift @r, $DEFAULT_NAME if $self->{'_meta'}->{$DEFAULT_NAME};
    return @r;
}


=head2 meta_length

 Title   : meta_length()
 Usage   : $meeta_len  = $obj->meta_length();
 Function: return the number of elements in the meta set
 Returns : integer
 Args    : -

=cut

sub meta_length {
   my ($self) = @_;
   return $self->named_meta_length($DEFAULT_NAME);
}


=head2 named_meta_length

 Title   : named_meta_length()
 Usage   : $meta_len  = $obj->named_meta_length($name);
 Function: return the number of elements in the named meta set
 Returns : integer
 Args    : -

=cut

sub named_meta_length {
   my ($self, $name) = @_;
   $name ||= $DEFAULT_NAME;
   return length ($self->{'_meta'}->{$name});
}


=head2 force_flush

 Title   : force_flush()
 Usage   : $force_flush = $obj->force_flush(1);
 Function: Automatically pad with empty values or truncate meta values
           to sequence length. Not done by default.
 Returns : boolean 1 or 0
 Args    : optional boolean value

Note that if you turn this forced padding off, the previously padded
values are not changed.

=cut

sub force_flush {
    my ($self, $value) = @_;

    if (defined $value) {
        if ($value) {
            $self->{force_flush} = 1;
            $self->_do_flush;
        } else {
            $self->{force_flush} = 0;
        }
    }

    return $self->{force_flush};
}


=head2 _do_flush

 Title   : _do_flush
 Usage   : 
 Function: internal method to do the force that meta values are same 
           length as the sequence . Called from L<force_flush>
 Returns : 
 Args    : 

=cut


sub _do_flush {
    my ($self) = @_;

    foreach my $name ( ('DEFAULT', $self->meta_names) ) {

        # elongnation
        if ($self->length > $self->named_meta_length($name)) {
            $self->{'_meta'}->{$name} .= $META_GAP x ($self->length - $self->named_meta_length($name)) ;
        }
        # truncation
        elsif ( $self->length < $self->named_meta_length($name) ) {
            $self->{_meta}->{$name} = substr($self->{_meta}->{$name}, 0, $self->length-1);
        }
    }

}


=head2 is_flush

 Title   : is_flush
 Usage   : $is_flush  = $obj->is_flush()
           or  $is_flush = $obj->is_flush($my_meta_name)
 Function: Boolean to tell if all meta values are in
           flush with the sequence length.
           Returns true if force_flush() is set
           Set verbosity to a positive value to see failed meta sets
 Returns : boolean 1 or 0
 Args    : optional name of the meta set

=cut

sub is_flush {

    my ($self, $name) = shift;

    return 1 if $self->force_flush;

    my $sticky = '';


    if ($name) {
        $sticky .= "$name " if $self->length != $self->named_meta_length($name);
    } else {
        foreach my $m ($self->meta_names) {
            $sticky .= "$m " if ($self->named_meta_length($m) > 0) && ($self->length != $self->named_meta_length($m));
        }
    }

    if ($sticky) {
        print "These meta set are not flush: $sticky\n" if $self->verbose; 
        return 0;
    }

    return 1;
}


=head1 Bio::PrimarySeqI methods

=head2 revcom

 Title   : revcom
 Usage   : $newseq = $seq->revcom();
 Function: Produces a new Bio::Seq::MetaI implementing object where
           the order of residues and their meta information is reversed.
 Returns : A new (fresh) Bio::Seq::Meta object
 Args    : none
 Throws  : if the object returns false on is_flush()

Note: The method does nothing to meta values, it reorders them, only.

=cut

sub revcom {
    my $self = shift;

    $self->throw("Can not get a reverse complement. The object is not flush.")
        unless $self->is_flush;

    my $new = $self->SUPER::revcom;
    foreach (keys %{$self->{_meta}}) {
        $new->named_meta($_, scalar reverse $self->{_meta}->{$_} );
    };
    return $new;
}

=head2 trunc

 Title   : trunc
 Usage   : $subseq = $seq->trunc(10,100);
 Function: Provides a truncation of a sequence together with meta data
 Returns : a fresh Bio::Seq::Meta implementing object
 Args    : Two integers denoting first and last residue of the sub-sequence.

=cut

sub trunc {
    my ($self, $start, $end) = @_;

    # test arguments
    $start =~ /^[+]?\d+$/ and $start > 0 or
        $self->throw("Need at least a positive integer start value as start");
    $end =~ /^[+]?\d+$/ and $end > 0 or
        $self->throw("Need at least a positive integer start value as end");
    $end >= $start or
        $self->throw("End position has to be larger or equal to start");
    $end <= $self->length or
        $self->throw("End position can not be larger than sequence length");

    my $new = $self->SUPER::trunc($start, $end);
    $start--;
    foreach (keys %{$self->{_meta}}) {
        $new->named_meta($_,
                         substr($self->{_meta}->{$_}, $start, $end - $start)
                        );
    };
    return $new;
}


sub to_string {
    my ($self) = @_;
    my $out = Bio::SeqIO->new(-format=>'metafasta');
    $out->write_seq($self);
    return 1;
}

1;
