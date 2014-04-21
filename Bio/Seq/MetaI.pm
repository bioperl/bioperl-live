#
# BioPerl module for Bio::Seq::MetaI
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

Bio::Seq::MetaI - Interface for sequence objects with residue-based
meta information

=head1 SYNOPSIS

  # get a Bio::Seq::MetaI compliant object somehow

  # to test this is a meta seq object
  $obj->isa("Bio::Seq::MetaI")
     || $obj->throw("$obj not a Bio::Seq::MetaI");

  # accessors
  $string     = $obj->meta;
  $string     = $obj->meta_text;
  $substring  = $obj->submeta(12,50);
  $unique_key = $obj->accession_number();


=head1 DESCRIPTION

This class defines an abstract interface for basic residue-based meta
information. Examples of this kind of meta data are secondary
structures (RNA and protein), protein hydrophobicity assignments, or
other alternative alphabets for polypeptides, sequence quality data
and nucleotide alignments with translations.

The length of the meta data sequence is not dependent on the amount of
the meta information. The meta information always covers all the
residues, but a blank value is used to denote unavailable
information. If necessary the implementation quietly truncates or
extends meta information with blank values. Definition of blank is
implementation dependent. Gaps in MSAs should not have meta
information.

At this point a residue in a sequence object can have only one meta
value. If you need more, use multiple copies of the sequence object.

Meta data storage can be implemented in various ways, e.g: string,
array of scalars, array of hashes, array of objects.

If the implementation so chooses, there can be more then one meta
values associated to each residue. See L<named_meta> and
L<names_submeta>. Note that use of arbitrary names is very prone to
typos leading to creation of additional copies of meta data sets.

Bio::Seq::Meta provides basic, pure perl implementation of sequences
with meta information. See L<Bio::Seq::Meta>. Application specific
implementations will override and add to these methods.

=head2 Method naming

Character based meta data is read and set by method meta() and its
variants. These are the suffixes and prefixes used in the variants:

    [named_] [sub] meta [_text]

=over 3

=item _text

Suffix B<_text> guaranties that output is a string. Note that it does
not limit the input.

=item sub

Prefix B<sub>, like in subseq(), means that the method applies to sub
region of the sequence range and takes start and end as arguments.
Unlike subseq(), these methods are able to set values.  If the range
is not defined, it defaults to the complete sequence.

=item named_

Prefix B<named_> in method names allows the used to attach multiple
meta strings to one sequence by explicitly naming them. The name is
always the first argument to the method. The "unnamed" methods use the
class wide default name for the meta data and are thus special cases
"named" methods.

Note that internally names are keys in a hash and any misspelling of a
name will silently store the data under a wrong name. The used names
(keys) can be retrieved using method meta_names(). See L<meta_names>.

=back


=head1 SEE ALSO

L<Bio::Seq::Meta>, 
L<Bio::Seq::Meta::Array>, 
L<Bio::Seq::EncodedSeq>, 
L<Bio::Tools::OddCodes>, 
L<Bio::Seq::Quality>


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

Chad Matsalla, bioinformatics@dieselwurks.com;
Aaron Mackey, amackey@virginia.edu;
Peter Schattner schattner@alum.mit.edu;
Richard Adams, Richard.Adams@ed.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::MetaI;
use strict;

use base qw(Bio::Root::RootI);

=head2 meta

 Title   : meta
 Usage   : $meta_values  = $obj->meta($values_string);
 Function:

           Get and set method for the unnamed meta data starting from
           residue position one. Since it is dependent on the length
           of the sequence, it needs to be manipulated after the
           sequence.

           The implementation may choose to accept argument values in
           a string or in an array (reference) or in a hash
           (reference).

           The return value may be a string or an array reference,
           depending on the implentation. If in doubt, use meta_text()
           which is a variant guarantied to return a string.  See
           L<meta_text>.

           The length of the returned value always matches the length
           of the sequence.

 Returns : A reference to an array or a string
 Args    : new value, optional

=cut

sub meta { shift->throw_not_implemented }

=head2 meta_text

 Title   : meta_text()
 Usage   : $meta_values  = $obj->meta_text($values_arrayref);
 Function: Variant of meta() guarantied to return a textual
           representation of the meta data. For details, see L<meta>.
 Returns : a string
 Args    : new value, optional

=cut

sub meta_text { shift->throw_not_implemented }

=head2 named_meta

 Title   : named_meta()
 Usage   : $meta_values  = $obj->named_meta($name, $values_arrayref);
 Function: A more general version of meta(). Each meta data set needs
           to be named. See also L<meta_names>.
 Returns : a string
 Args    : scalar, name of the meta data set
           new value, optional

=cut

sub named_meta { shift->throw_not_implemented }

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

sub named_meta_text { shift->throw_not_implemented }

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

           If implementation tries to set values beyond the current
           sequence, they should be ignored.

           The return value may be a string or an array reference,
           depending on the implentation. If in doubt, use
           submeta_text() which is a variant guarantied to return a
           string.  See L<submeta_text>.

 Returns : A reference to an array or a string
 Args    : integer, start position, optional
           integer, end position, optional
           new value, optional

=cut

sub submeta { shift->throw_not_implemented }

=head2 submeta_text

 Title   : submeta_text
 Usage   : $meta_values  = $obj->submeta_text(20, $value_string);
 Function: Variant of submeta() guarantied to return a textual
           representation  of meta data. For details, see L<meta>.
 Returns : a string
 Args    : integer, start position, optional
           integer, end position, optional
           new value, optional

=cut

sub submeta_text { shift->throw_not_implemented }

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

sub named_submeta { shift->throw_not_implemented }

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

sub named_submeta_text { shift->throw_not_implemented }

=head2 meta_names

 Title   : meta_names
 Usage   : @meta_names  = $obj->meta_names()
 Function: Retrives an array of meta data set names. The default (unnamed)
           set name is guarantied to be the first name.
 Returns : an array of names
 Args    : none

=cut

sub meta_names { shift->throw_not_implemented }


=head2 force_flush

 Title   : force_flush()
 Usage   : $force_flush = $obj->force_flush(1);
 Function: Automatically pad with empty values or truncate meta values to
            sequence length
 Returns : boolean 1 or 0
 Args    : optional boolean value

=cut

sub force_flush { shift->throw_not_implemented }


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

sub is_flush { shift->throw_not_implemented }


=head2 meta_length

 Title   : meta_length()
 Usage   : $meeta_len  = $obj->meta_length();
 Function: return the number of elements in the meta set
 Returns : integer
 Args    : -

=cut

sub meta_length { shift->throw_not_implemented }

=head2 named_meta_length

 Title   : named_meta_length()
 Usage   : $meeta_len  = $obj->named_meta_length($name);
 Function: return the number of elements in the named meta set
 Returns : integer
 Args    : -

=cut

sub named_meta_length { shift->throw_not_implemented }


=head1 Bio::PrimarySeqI methods

Implemeting classes will need to rewrite these Bio::PrimaryI methods.

=cut

=head2 revcom

 Title   : revcom
 Usage   : $newseq = $seq->revcom();
 Function: Produces a new Bio::Seq::MetaI implementing object where
           the order of residues and their meta information is reversed.
 Returns : A new (fresh) Bio::Seq::MetaI object
 Args    : none

=cut

sub revcom { shift->throw_not_implemented }

=head2 trunc

 Title   : trunc
 Usage   : $subseq = $myseq->trunc(10,100);
 Function: Provides a truncation of a sequence
 Returns : a fresh Bio::Seq::MetaI implementing object
 Args    : Two integers denoting first and last residue of the sub-sequence.

=cut

sub trunc { shift->throw_not_implemented }


1;
