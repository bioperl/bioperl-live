#
# BioPerl module for Bio::Search::Tiling::TilingI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Mark A. Jensen <maj@fortinbras.us>
#
# Copyright Mark A. Jensen
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Tiling::TilingI - Abstract interface for an HSP tiling module

=head1 SYNOPSIS

Not used directly. Useful POD here for developers, however.

The interface is desgined to make the following code conversion as
simple as possible:

From:

 # Bio::Search::SearchUtils-based
 while ( local $_ = $result->next_hit ) {
    printf( "E-value: %g; Fraction aligned: %f; Number identical: %d\n",
      $hit->significance, $hit->frac_aligned_query, $hit->num_identical);
 }

To:

 # TilingI-based
 while ( local $_ = $result->next_hit ) {
    my $tiling = Bio::Search::Tiling::MyTiling($_);
    printf( "E-value: %g; Fraction aligned: %f; Number identical: %d\n",
      $hit->significance, $tiling->frac_aligned_query, $tiling->num_identical);
 }



=head1 DESCRIPTION

This module provides strong suggestions for any intended HSP tiling
object implementation. An object subclassing TilingI should override
the methods defined here according to their descriptions below.

See the section STATISTICS METHODS for hints on implementing methods
that are valid across different algorithms and report types.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Mark A. Jensen

Email maj@fortinbras.us

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::Tiling::TilingI;
use strict;
use warnings;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;

use base qw(Bio::Root::Root);

=head2 STATISTICS METHODS

The tiling statistics can be thought of as global counterparts to
similar statistics defined for the individual HSPs. We therefore
prescribe definitions for many of the synonymous methods defined in
L<Bio::Search::HSP::HSPI>.

The tiling statistics must be able to keep track of the coordinate
systems in which both the query and subject sequences exist; i.e.,
either nucleotide or amino acid. This information is typically
inferred from the name of the algorithm used to perform the original
search (contained in C<$hit_object-E<gt>algorithm>). Here is a table
of algorithm information that may be useful (if you trust us).

 algorithm   query on hit   coordinates(q/h)
 ---------   ------------   ---------------
  blastn      dna on dna         dna/dna
  blastp      aa  on aa           aa/aa
  blastx      xna on aa          dna/aa
 tblastn      aa  on xna          aa/dna
 tblastx      xna on xna         dna/dna
   fasta      dna on dna         dna/dna
   fasta      aa  on aa           aa/aa
   fastx      xna on aa          dna/aa
   fasty      xna on aa          dna/aa
  tfasta      aa  on xna          aa/dna
  tfasty      aa  on xna          aa/dna
 megablast    dna on dna         dna/dna

  xna: translated nucleotide data

Statistics methods must also be aware of differences in reporting
among the algorithms. Hit attributes are not necessarily normalized
over all algorithms. Devs, please feel free to add examples to the
list below.

=over

=item NCBI BLAST vs WU-BLAST (AB-BLAST) lengths

The total length of the alignment is reported differently between these two flavors. C<$hit_object-E<gt>length()> will contain the number in the denominator of the stats line; i.e., 120 in 

 Identical = 34/120 Positives = 67/120

NCBI BLAST uses the total length of the query sequence as input by the user (a.k.a. "with gaps"). WU-BLAST uses the length of the query sequence actually aligned by the algorithm (a.k.a. "without gaps").

=back

Finally, developers should remember that sequence data may or may not
be associated with the HSPs contained in the hit object. This will
typically depend on whether a full report (e.g, C<blastall -m0>) or a
summary (e.g., C<blastall -m8>) was parsed. Statistics methods that
depend directly on the sequence data will need to check that
that data is present.

=head2 identities

 Title   : identities
 Alias   : num_identical
 Usage   : $num_identities = $tiling->identities()
 Function: Return the estimated or exact number of identities in the
           tiling, accounting for overlapping HSPs
 Example : 
 Returns : number of identical residue pairs
 Args    :

=cut

sub identities{
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

#HSPI synonym
sub num_identical { shift->identities( @_ ) }

=head2 conserved

 Title   : conserved
 Alias   : num_conserved
 Usage   : $num_conserved = $tiling->conserved()
 Function: Return the estimated or exact number of conserved sites in the 
           tiling, accounting for overlapping HSPs
 Example : 
 Returns : number of conserved residue pairs
 Args    :

=cut

sub conserved{
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

#HSPI synonym
sub num_conserved { shift->conserved( @_ ) }

=head2 length

 Title   : length
 Usage   : $max_length = $tiling->length($type)
 Function: Return the total number of residues of the subject or query
           sequence covered by the tiling
 Returns : number of "logical" residues covered
 Args    : scalar $type, one of 'hit', 'subject', 'query'

=cut

sub length{
    my ($self, $type, @args) = @_;
    $self->throw_not_implemented;
}

=head2 frac_identical

 Title   : frac_identical
 Usage   : $tiling->frac_identical($type)
 Function: Return the fraction of sequence length consisting
           of identical pairs
 Returns : scalar float
 Args    : scalar $type, one of 'hit', 'subject', 'query'
 Note    : This method must take account of the $type coordinate
           system and the length reporting method (see STATISTICS
           METHODS above)

=cut

sub frac_identical {
    my ($self, $type, @args) = @_;
    $self->throw_not_implemented;
}

=head2 percent_identity

 Title   : percent_identity
 Usage   : $tiling->percent_identity($type)
 Function: Return the fraction of sequence length consisting
           of identical pairs as a percentage
 Returns : scalar float
 Args    : scalar $type, one of 'hit', 'subject', 'query'

=cut

sub percent_identity {
    my ($self, $type, @args) = @_;
    return $self->frac_identical($type, @args) * 100;
}

=head2 frac_conserved

 Title   : frac_conserved
 Usage   : $tiling->frac_conserved($type)
 Function: Return the fraction of sequence length consisting
           of conserved pairs
 Returns : scalar float
 Args    : scalar $type, one of 'hit', 'subject', 'query'
 Note    : This method must take account of the $type coordinate
           system and the length reporting method (see STATISTICS
           METHODS above)

=cut

sub frac_conserved{
    my ($self, $type, @args) = @_;
    $self->throw_not_implemented;
}

=head2 percent_conserved

 Title   : percent_conserved
 Usage   : $tiling->percent_conserved($type)
 Function: Return the fraction of sequence length consisting
           of conserved pairs as a percentage
 Returns : scalar float
 Args    : scalar $type, one of 'hit', 'subject', 'query'

=cut

sub percent_conserved {
    my ($self, $type, @args) = @_;
    return $self->frac_conserved($type, @args) * 100;
}

=head2 frac_aligned

 Title   : frac_aligned
 Usage   : $tiling->frac_aligned($type)
 Function: Return the fraction of B<input> sequence length consisting
           that was aligned by the algorithm
 Returns : scalar float
 Args    : scalar $type, one of 'hit', 'subject', 'query'
 Note    : This method must take account of the $type coordinate
           system and the length reporting method (see STATISTICS
           METHODS above)

=cut

sub frac_aligned{
    my ($self, $type, @args) = @_;
    $self->throw_not_implemented;
}

# aliases for back compat
sub frac_aligned_query { shift->frac_aligned('query', @_) }
sub frac_aligned_hit { shift->frac_aligned('hit', @_) }

=head2 range

 Title   : range
 Usage   : $tiling->range($type)
 Function: Returns the extent of the longest tiling
           as ($min_coord, $max_coord)
 Returns : array of two scalar integers
 Args    : scalar $type, one of 'hit', 'subject', 'query'

=cut

sub range {
    my ($self, $type, @args) = @_;
    $self->throw_not_implemented;
}

=head1 TILING ITERATORS

=head2 next_tiling

 Title   : next_tiling
 Usage   : @hsps = $self->next_tiling($type);
 Function: Obtain a tiling of HSPs over the $type ('hit', 'subject',
           'query') sequence
 Example :
 Returns : an array of HSPI objects
 Args    : scalar $type: one of 'hit', 'subject', 'query', with
           'subject' an alias for 'hit'

=cut

sub next_tiling{
    my ($self,$type,@args) = @_;
    $self->throw_not_implemented;
}

=head2 rewind_tilings

 Title   : rewind_tilings
 Usage   : $self->rewind_tilings($type)
 Function: Reset the next_tilings($type) iterator
 Example :
 Returns : True on success
 Args    : scalar $type: one of 'hit', 'subject', 'query', with
           'subject' an alias for 'hit'

=cut

sub rewind_tilings{
    my ($self, $type, @args) = @_;
    $self->throw_not_implemented;
}

#alias
sub rewind { shift->rewind_tilings(@_) }

=head1 INFORMATIONAL ACCESSORS

=head2 algorithm

 Title   : algorithm
 Usage   : $tiling->algorithm
 Function: Retrieve the algorithm name associated with the 
           invocant's hit object
 Returns : scalar string 
 Args    : 

=cut

sub algorithm{
    my ($self, @args) = @_;
    $self->throw_not_implemented;
}

1;
