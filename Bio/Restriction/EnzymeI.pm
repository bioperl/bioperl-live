#------------------------------------------------------------------
#
# BioPerl module Bio::Restriction::EnzymeI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho, heikki-at-bioperl-dot-org
#
# You may distribute this module under the same terms as perl itself
#------------------------------------------------------------------

## POD Documentation:

=head1 NAME

Bio::Restriction::EnzymeI - Interface class for restriction endonuclease

=head1 SYNOPSIS

  # do not run this class directly

=head1 DESCRIPTION

This module defines methods for a single restriction endonuclease.  For an
implementation, see L<Bio::Restriction::Enzyme>.

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
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 CONTRIBUTORS

Rob Edwards, redwards@utmem.edu

=head1 SEE ALSO

L<Bio::Restriction::Enzyme>

=head1 APPENDIX

Methods beginning with a leading underscore are considered private and
are intended for internal use by this module. They are not considered
part of the public interface and are described here for documentation
purposes only.

=cut

package Bio::Restriction::EnzymeI;
use strict;

use base qw(Bio::Root::RootI);

=head1 Essential methods

=cut

=head2 name

 Title    : name
 Usage    : $re->name($newval)
 Function : Gets/Sets the restriction enzyme name
 Example  : $re->name('EcoRI')
 Returns  : value of name
 Args     : newvalue (optional)

This will also clean up the name. I have added this because some
people get confused about restriction enzyme names.  The name should
be One upper case letter, and two lower case letters (because it is
derived from the organism name, eg.  EcoRI is from E. coli). After
that it is all confused, but the numbers should be roman numbers not
numbers, therefore we'll correct those. At least this will provide
some standard, I hope.

=cut

sub name {  shift->throw_not_implemented; }

=head2 site

 Title     : site
 Usage     : $re->site();
 Function  : Gets/sets the recognition sequence for the enzyme.
 Example   : $seq_string = $re->site();
 Returns   : String containing recognition sequence indicating
           : cleavage site as in  'G^AATTC'.
 Argument  : n/a
 Throws    : n/a

Side effect: the sequence is always converted to upper case.

The cut site can also be set by using methods L<cut|cut> and
L<complementary_cut|complementary_cut>.

This will pad out missing sequence with N's. For example the enzyme
Acc36I cuts at ACCTGC(4/8). This will be returned as ACCTGCNNNN^

Note that the common notation ACCTGC(4/8) means that the forward
strand cut is four nucleotides after the END of the recognition
site. The forward cut() in the coordinates used here in Acc36I
ACCTGC(4/8) is at 6+4 i.e. 10.

** This is the main setable method for the recognition site.

=cut

sub site {  shift->throw_not_implemented; }

=head2 revcom_site

 Title     : revcom_site
 Usage     : $re->revcom_site();
 Function  : Gets/sets the complementary recognition sequence for the enzyme.
 Example   : $seq_string = $re->revcom_site();
 Returns   : String containing recognition sequence indicating
           : cleavage site as in  'G^AATTC'.
 Argument  : Sequence of the site
 Throws    : n/a

This is the same as site, except it returns the revcom site. For
palindromic enzymes these two are identical. For non-palindromic
enzymes they are not!

See also L<site|site> above.

=cut

sub cuts_after {  shift->throw_not_implemented; }

=head2 cut

 Title     : cut
 Usage     : $num = $re->cut(1);
 Function  : Sets/gets an integer indicating the position of cleavage
             relative to the 5' end of the recognition sequence in the
             forward strand.

             For type II enzymes, sets the symmetrically positioned
             reverse strand cut site by calling complementary_cut().

 Returns   : Integer, 0 if not set
 Argument  : an integer for the forward strand cut site (optional)


Note that the common notation ACCTGC(4/8) means that the forward
strand cut is four nucleotides after the END of the recognition
site. The forwad cut in the coordinates used here in Acc36I
ACCTGC(4/8) is at 6+4 i.e. 10.

Note that REBASE uses notation where cuts within symmetic sites are
marked by '^' within the forward sequence but if the site is
asymmetric the parenthesis syntax is used where numbering ALWAYS
starts from last nucleotide in the forward strand. That's why AciI has
a site usually written as CCGC(-3/-1) actualy cuts in

  C^C G C
  G G C^G

In our notation, these locations are 1 and 3.

The cuts locations in the notation used are relative to the first
(non-N) nucleotide of the reported forward strand of the recognition
sequence. The following diagram numbers the phosphodiester bonds
(marked by + ) which can be cut by the restriction enzymes:

                           1   2   3   4   5   6   7   8  ...
     N + N + N + N + N + G + A + C + T + G + G + N + N + N
  ... -5  -4  -3  -2  -1

=cut

sub cut {  shift->throw_not_implemented; }

=head2 complementary_cut

 Title     : complementary_cut
 Usage     : $num = $re->complementary_cut('1');
 Function  : Sets/Gets an integer indicating the position of cleavage
           : on the reverse strand of the restriction site.
 Returns   : Integer
 Argument  : An integer (optional)
 Throws    : Exception if argument is non-numeric.

This method determines the cut on the reverse strand of the sequence.
For most enzymes this will be within the sequence, and will be set
automatically based on the forward strand cut, but it need not be.

B<Note> that the returned location indicates the location AFTER the
first non-N site nucleotide in the FORWARD strand.

=cut

sub complementary_cut {  shift->throw_not_implemented; }

=head1 Read only (usually) recognition site descriptive methods

=cut

=head2 type

 Title     : type
 Usage     : $re->type();
 Function  : Get/set the restriction system type
 Returns   : 
 Argument  : optional type: ('I'|II|III)

Restriction enzymes have been catezorized into three types. Some
REBASE formats give the type, but the following rules can be used to
classify the known enzymes:

=over 4

=item 1

Bipartite site (with 6-8 Ns in the middle and the cut site
is E<gt> 50 nt away) =E<gt> type I

=item 2

Site length E<lt> 3  =E<gt> type I

=item 3

5-6 asymmetric site and cuts E<gt>20 nt away =E<gt> type III

=item 4

All other  =E<gt> type II

=back

There are some enzymes in REBASE which have bipartite recognition site
and cat far from the site but are still classified as type I. I've no
idea if this is really so.

=cut

sub type {  shift->throw_not_implemented; }

=head2 seq

 Title     : seq
 Usage     : $re->seq();
 Function  : Get the Bio::PrimarySeq.pm object representing
           : the recognition sequence
 Returns   : A Bio::PrimarySeq object representing the
             enzyme recognition site
 Argument  : n/a
 Throws    : n/a


=cut

sub seq {  shift->throw_not_implemented; }

=head2 string

 Title     : string
 Usage     : $re->string();
 Function  : Get a string representing the recognition sequence.
 Returns   : String. Does NOT contain a  '^' representing the cut location
             as returned by the site() method.
 Argument  : n/a
 Throws    : n/a

=cut

sub string {  shift->throw_not_implemented; }

=head2 revcom

 Title     : revcom
 Usage     : $re->revcom();
 Function  : Get a string representing the reverse complement of
           : the recognition sequence.
 Returns   : String
 Argument  : n/a
 Throws    : n/a

=cut

sub revcom {  shift->throw_not_implemented; }

=head2 recognition_length

 Title     : recognition_length
 Usage     : $re->recognition_length();
 Function  : Get the length of the RECOGNITION sequence.
             This is the total recognition sequence,
             inluding the ambiguous codes.
 Returns   : An integer
 Argument  : Nothing

See also: L<non_ambiguous_length>

=cut

sub recognition_length {  shift->throw_not_implemented; }

=head2 non_ambiguous_length

 Title     : non_ambiguous_length
 Usage     : $re->non_ambiguous_length();
 Function  : Get the nonambiguous length of the RECOGNITION sequence.
             This is the total recognition sequence,
             excluding the ambiguous codes.
 Returns   : An integer
 Argument  : Nothing

See also: L<non_ambiguous_length>

=cut

sub non_ambiguous_length {  shift->throw_not_implemented; }

=head2 cutter

 Title    : cutter
 Usage    : $re->cutter
 Function : Returns the "cutter" value of the recognition site.

            This is a value relative to site length and lack of
            ambiguity codes. Hence: 'RCATGY' is a five (5) cutter site
            and 'CCTNAGG' a six cutter

            This measure correlates to the frequency of the enzyme
            cuts much better than plain recognition site length.

 Example  : $re->cutter
 Returns  : integer or float number
 Args     : none

Why is this better than just stripping the ambiguous codes? Think about
it like this: You have a random sequence; all nucleotides are equally
probable. You have a four nucleotide re site. The probability of that
site finding a match is one out of 4^4 or 256, meaning that on average
a four cutter finds a match every 256 nucleotides. For a six cutter,
the average fragment length is 4^6 or 4096. In the case of ambiguity
codes the chances are finding the match are better: an R (A|T) has 1/2
chance of finding a match in a random sequence. Therefore, for RGCGCY
the probability is one out of (2*4*4*4*4*2) which exactly the same as
for a five cutter! Cutter, although it can have non-integer values
turns out to be a useful and simple measure.

From bug 2178: VHDB are ambiguity symbols that match three different
nucleotides, so they contribute less to the effective recognition sequence
length than e.g. Y which matches only two nucleotides. A symbol which matches n
of the 4 nucleotides has an effective length of 1 - log(n) / log(4).

=cut

sub cutter {  shift->throw_not_implemented; }

=head2 is_palindromic

 Title     : is_palindromic
 Usage     : $re->is_palindromic();
 Function  : Determines if the recognition sequence is palindromic
           : for the current restriction enzyme.
 Returns   : Boolean
 Argument  : n/a
 Throws    : n/a

A palindromic site (EcoRI):

  5-GAATTC-3
  3-CTTAAG-5

=cut

sub is_palindromic {  shift->throw_not_implemented; }

=head2 overhang

 Title     : overhang
 Usage     : $re->overhang();
 Function  : Determines the overhang of the restriction enzyme
 Returns   : "5'", "3'", "blunt" of undef
 Argument  : n/a
 Throws    : n/a

A blunt site in SmaI returns C<blunt>

  5' C C C^G G G 3'
  3' G G G^C C C 5'

A 5' overhang in EcoRI returns C<5'>

  5' G^A A T T C 3'
  3' C T T A A^G 5'

A 3' overhang in KpnI returns C<3'>

  5' G G T A C^C 3'
  3' C^C A T G G 5'

=cut

sub overhang {  shift->throw_not_implemented; }

=head2 overhang_seq

 Title     : overhang_seq
 Usage     : $re->overhang_seq();
 Function  : Determines the overhang sequence of the restriction enzyme
 Returns   : a Bio::LocatableSeq
 Argument  : n/a
 Throws    : n/a

I do not think it is necessary to create a seq object of these. (Heikki)

Note: returns empty string for blunt sequences and undef for ones that
we don't know.  Compare these:

A blunt site in SmaI returns empty string

  5' C C C^G G G 3'
  3' G G G^C C C 5'

A 5' overhang in EcoRI returns C<AATT>

  5' G^A A T T C 3'
  3' C T T A A^G 5'

A 3' overhang in KpnI returns C<GTAC>

  5' G G T A C^C 3'
  3' C^C A T G G 5'

Note that you need to use method L<overhang|overhang> to decide
whether it is a 5' or 3' overhang!!!

Note: The overhang stuff does not work if the site is asymmetric! Rethink! 

=cut

sub overhang_seq {  shift->throw_not_implemented; }

=head2 compatible_ends

 Title     : compatible_ends
 Usage     : $re->compatible_ends($re2);
 Function  : Determines if the two restriction enzyme cut sites
              have compatible ends.
 Returns   : 0 if not, 1 if only one pair ends match, 2 if both ends.
 Argument  : a Bio::Restriction::Enzyme
 Throws    : unless the argument is a Bio::Resriction::Enzyme and
             if there are Ns in the ovarhangs

In case of type II enzymes which which cut symmetrically, this
function can be considered to return a boolean value.

=cut

sub compatible_ends {shift->throw_not_implemented;}

=head2 is_ambiguous

 Title     : is_ambiguous
 Usage     : $re->is_ambiguous();
 Function  : Determines if the restriction enzyme contains ambiguous sequences
 Returns   : Boolean
 Argument  : n/a
 Throws    : n/a

=cut

sub is_ambiguous {  shift->throw_not_implemented; }

=head2 Additional methods from Rebase

=cut


=head2 is_prototype

 Title    : is_prototype
 Usage    : $re->is_prototype
 Function : Get/Set method for finding out if this enzyme is a prototype
 Example  : $re->is_prototype(1)
 Returns  : Boolean
 Args     : none

Prototype enzymes are the most commonly available and usually first
enzymes discoverd that have the same recognition site. Using only
prototype enzymes in restriciton analysis avoids redundacy and
speeds things up.

=cut

sub is_prototype {  shift->throw_not_implemented; }

=head2 prototype_name

 Title    : prototype_name
 Usage    : $re->prototype_name
 Function : Get/Set method for the name of prototype for
            this enzyme's recognition site
 Example  : $re->prototype_name(1)
 Returns  : prototype enzyme name string or an empty string
 Args     : optional prototype enzyme name string

If the enzyme itself is the protype, its own name is returned.  Not to
confuse the negative result with an unset value, use method
L<is_prototype|is_prototype>.

This method is called I<prototype_name> rather than I<prototype>,
because it returns a string rather than on object.

=cut

sub prototype_name {  shift->throw_not_implemented; }

=head2 isoschizomers

 Title     : isoschizomers
 Usage     : $re->isoschizomers(@list);
 Function  : Gets/Sets a list of known isoschizomers (enzymes that
             recognize the same site, but don't necessarily cut at
             the same position).
 Arguments : A reference to an array that contains the isoschizomers
 Returns   : A reference to an array of the known isoschizomers or 0
             if not defined.

Added for compatibility to REBASE

=cut

sub isoschizomers {  shift->throw_not_implemented; }

=head2 purge_isoschizomers

 Title     : purge_isoschizomers
 Usage     : $re->purge_isoschizomers();
 Function  : Purges the set of isoschizomers for this enzyme
 Arguments : 
 Returns   : 1

=cut

sub purge_isoschizomers {  shift->throw_not_implemented; }

=head2 methylation_sites

 Title     : methylation_sites
 Usage     : $re->methylation_sites(\%sites);
 Function  : Gets/Sets known methylation sites (positions on the sequence
             that get modified to promote or prevent cleavage).
 Arguments : A reference to a hash that contains the methylation sites
 Returns   : A reference to a hash of the methylation sites or
             an empty string if not defined.

There are three types of methylation sites:

=over 3

=item *  (6) = N6-methyladenosine

=item *  (5) = 5-methylcytosine

=item *  (4) = N4-methylcytosine

=back

These are stored as 6, 5, and 4 respectively.  The hash has the
sequence position as the key and the type of methylation as the value.
A negative number in the sequence position indicates that the DNA is
methylated on the complementary strand.

Note that in REBASE, the methylation positions are given 
Added for compatibility to REBASE.

=cut

sub methylation_sites {  shift->throw_not_implemented; }

=head2 purge_methylation_sites

 Title     : purge_methylation_sites
 Usage     : $re->purge_methylation_sites();
 Function  : Purges the set of methylation_sites for this enzyme
 Arguments : 
 Returns   : 

=cut

sub purge_methylation_sites {  shift->throw_not_implemented; }

=head2 microbe

 Title     : microbe
 Usage     : $re->microbe($microbe);
 Function  : Gets/Sets microorganism where the restriction enzyme was found
 Arguments : A scalar containing the microbes name
 Returns   : A scalar containing the microbes name or 0 if not defined

Added for compatibility to REBASE

=cut

sub microbe {  shift->throw_not_implemented; }

=head2 source

 Title     : source
 Usage     : $re->source('Rob Edwards');
 Function  : Gets/Sets the person who provided the enzyme
 Arguments : A scalar containing the persons name
 Returns   : A scalar containing the persons name or 0 if not defined

Added for compatibility to REBASE

=cut

sub source {  shift->throw_not_implemented; }

=head2 vendors

 Title     : vendors
 Usage     : $re->vendor(@list_of_companies);
 Function  : Gets/Sets the a list of companies that you can get the enzyme from.
             Also sets the commercially_available boolean
 Arguments : A reference to an array containing the names of companies
             that you can get the enzyme from
 Returns   : A reference to an array containing the names of companies
             that you can get the enzyme from

Added for compatibility to REBASE

=cut

sub vendors {  shift->throw_not_implemented; }

=head2 purge_vendors

 Title     : purge_vendors
 Usage     : $re->purge_references();
 Function  : Purges the set of references for this enzyme
 Arguments : 
 Returns   : 

=cut

sub purge_vendors {  shift->throw_not_implemented; }

=head2 vendor

 Title     : vendor
 Usage     : $re->vendor(@list_of_companies);
 Function  : Gets/Sets the a list of companies that you can get the enzyme from.
             Also sets the commercially_available boolean
 Arguments : A reference to an array containing the names of companies
             that you can get the enzyme from
 Returns   : A reference to an array containing the names of companies
             that you can get the enzyme from

Added for compatibility to REBASE

=cut

sub vendor {  shift->throw_not_implemented; }

=head2 references

 Title     : references
 Usage     : $re->references(string);
 Function  : Gets/Sets the references for this enzyme
 Arguments : an array of string reference(s) (optional)
 Returns   : an array of references

Use L<purge_references|purge_references> to reset the list of references

This should be a L<Bio::Biblio> or L<Bio::Annotation::Reference> object, but its not (yet)

=cut

sub references {  shift->throw_not_implemented; }

=head2 purge_references

 Title     : purge_references
 Usage     : $re->purge_references();
 Function  : Purges the set of references for this enzyme
 Arguments : 
 Returns   : 1

=cut

sub purge_references {  shift->throw_not_implemented; }

=head2 clone

 Title     : clone
 Usage     : $re->clone
 Function  : Deep copy of the object
 Arguments : -
 Returns   : new Bio::Restriction::EnzymeI object

This works as long as the object is a clean in-memory object using
scalars, arrays and hashes. You have been warned.

If you have module Storable, it is used, otherwise local code is used.
Todo: local code cuts circular references.

=cut

sub clone {  shift->throw_not_implemented; }

1;

