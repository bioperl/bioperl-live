#
# BioPerl module for Bio::Search::Hit::ModelHit
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Chris Fields <cjfields at bioperl dot org>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Hit::ModelHit - A model-based implementation of the Bio::Search::Hit::HitI interface

=head1 SYNOPSIS

    use Bio::Search::Hit::ModelHit;
    my $hit = Bio::Search::Hit::ModelHit->new(-algorithm => 'rnamotif');

    # typically one gets HitI objects from a SearchIO stream via a ResultI
    use Bio::SearchIO;
    my $parser = Bio::SearchIO->new(-format => 'infernal', -file => 'trap.inf');

    my $result = $parser->next_result;
    my $hit    = $result->next_hit;

=head1 DESCRIPTION

This object handles the hit data from a database search using models or
descriptors instead of sequences, such as Infernal, HMMER, RNAMotif, etc.

Unless you're writing a parser, you won't ever need to create a ModelHit or
any other HitI-implementing object. If you use the SearchIO system, HitI objects
are created automatically from a SearchIO stream which returns
Bio::Search::Hit::HitI objects.

Note that several HitI-based methods have been overridden from ModelHit due to
their unreliability when dealing with queries that aren't sequence-based. It may
be possible to reimplement these at a later point, but for the time being they
will throw warnings and return w/o results.

For documentation on what you can do with ModelHit (and other HitI objects),
please see the API documentation in
L<Bio::Search::Hit::HitI|Bio::Search::Hit::HitI>.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Chris Fields

Email cjfields at bioperl dot org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::Hit::ModelHit;

use strict;

use base qw(Bio::Search::Hit::GenericHit);

=head1 HitI methods implemented in parent class Bio::Search::Hit::ModelHit

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::Hit::ModelHit->new();
 Function: Builds a new Bio::Search::Hit::ModelHit object
 Returns : Bio::Search::Hit::ModelHit
 Args    : -name         => Name of Hit (required)
           -description  => Description (optional)
           -accession    => Accession number (optional)
           -ncbi_gi      => NCBI GI UID (optional)
           -length       => Length of the Hit (optional)
           -score        => Raw Score for the Hit (optional)
           -bits         => Bit Score for the Hit (optional)
           -significance => Significance value for the Hit (optional)
           -algorithm    => Algorithm used (BLASTP, FASTX, etc...)
           -hsps         => Array ref of HSPs for this Hit.
           -found_again  => boolean, true if hit appears in a
                            "previously found" section of a PSI-Blast report.
           -hsp_factory  => Bio::Factory::ObjectFactoryI able to create HSPI
                            objects.

=cut

=head2 add_hsp

 Title   : add_hsp
 Usage   : $hit->add_hsp($hsp)
 Function: Add a HSP to the collection of HSPs for a Hit
 Returns : number of HSPs in the Hit
 Args    : Bio::Search::HSP::HSPI object, OR hash ref containing data suitable
           for creating a HSPI object (&hsp_factory must be set to get it back)

=cut

=head2 hsp_factory

 Title   : hsp_factory
 Usage   : $hit->hsp_factory($hsp_factory)
 Function: Get/set the factory used to build HSPI objects if necessary.
 Returns : Bio::Factory::ObjectFactoryI
 Args    : Bio::Factory::ObjectFactoryI

=cut

=head2 Bio::Search::Hit::HitI methods

Implementation of Bio::Search::Hit::HitI methods

=head2 name

 Title   : name
 Usage   : $hit_name = $hit->name();
 Function: returns the name of the Hit sequence
 Returns : a scalar string
 Args    : [optional] scalar string to set the name

=cut

=head2 accession

 Title   : accession
 Usage   : $acc = $hit->accession();
 Function: Retrieve the accession (if available) for the hit
 Returns : a scalar string (empty string if not set)
 Args    : none

=cut

=head2 description

 Title   : description
 Usage   : $desc = $hit->description();
 Function: Retrieve the description for the hit
 Returns : a scalar string
 Args    : [optional] scalar string to set the descrition

=cut

=head2 length

 Title   : length
 Usage   : my $len = $hit->length
 Function: Returns the length of the hit
 Returns : integer
 Args    : [optional] integer to set the length

=cut

=head2 algorithm

 Title   : algorithm
 Usage   : $alg = $hit->algorithm();
 Function: Gets the algorithm specification that was used to obtain the hit
           For BLAST, the algorithm denotes what type of sequence was aligned
           against what (BLASTN: dna-dna, BLASTP prt-prt, BLASTX translated
           dna-prt, TBLASTN prt-translated dna, TBLASTX translated
           dna-translated dna).
 Returns : a scalar string
 Args    : [optional] scalar string to set the algorithm

=cut

=head2 raw_score

 Title   : raw_score
 Usage   : $score = $hit->raw_score();
 Function: Gets the "raw score" generated by the algorithm.  What
           this score is exactly will vary from algorithm to algorithm,
           returning undef if unavailable.
 Returns : a scalar value
 Args    : [optional] scalar value to set the raw score

=cut

=head2 score

Equivalent to L<raw_score()|raw_score>

=cut

=head2 significance

 Title   : significance
 Usage   : $significance = $hit->significance();
 Function: Used to obtain the E or P value of a hit, i.e. the probability that
           this particular hit was obtained purely by random chance.  If
           information is not available (nor calculatable from other
           information sources), return undef.
 Returns : a scalar value or undef if unavailable
 Args    : [optional] scalar value to set the significance

=cut

=head2 bits

 Usage     : $hit_object->bits();
 Purpose   : Gets the bit score of the best HSP for the current hit.
 Example   : $bits = $hit_object->bits();
 Returns   : Integer or undef if bit score is not set
 Argument  : n/a
 Comments  : For BLAST1, the non-bit score is listed in the summary line.

See Also   : L<score()|score>

=cut

=head2 next_hsp

 Title    : next_hsp
 Usage    : while( $hsp = $obj->next_hsp()) { ... }
 Function : Returns the next available High Scoring Pair
 Example  :
 Returns  : Bio::Search::HSP::HSPI object or null if finished
 Args     : none

=cut

=head2 hsps

 Usage     : $hit_object->hsps();
 Purpose   : Get a list containing all HSP objects.
           : Get the numbers of HSPs for the current hit.
 Example   : @hsps = $hit_object->hsps();
           : $num  = $hit_object->hsps();  # alternatively, use num_hsps()
 Returns   : Array context : list of Bio::Search::HSP::BlastHSP.pm objects.
           : Scalar context: integer (number of HSPs).
           :                 (Equivalent to num_hsps()).
 Argument  : n/a. Relies on wantarray
 Throws    : Exception if the HSPs have not been collected.

See Also   : L<hsp()|hsp>, L<num_hsps()|num_hsps>

=cut

=head2 num_hsps

 Usage     : $hit_object->num_hsps();
 Purpose   : Get the number of HSPs for the present hit.
 Example   : $nhsps = $hit_object->num_hsps();
 Returns   : Integer or '-' if HSPs have not been callected
 Argument  : n/a

See Also   : L<hsps()|hsps>

=cut

=head2 rewind

 Title   : rewind
 Usage   : $hit->rewind;
 Function: Allow one to reset the HSP iterator to the beginning
           Since this is an in-memory implementation
 Returns : none
 Args    : none

=cut

=head2 ambiguous_aln

 Usage     : $ambig_code = $hit_object->ambiguous_aln();
 Purpose   : Sets/Gets ambiguity code data member.
 Example   : (see usage)
 Returns   : String = 'q', 's', 'qs', '-'
           :   'q'  = query sequence contains overlapping sub-sequences
           :          while sbjct does not.
           :   's'  = sbjct sequence contains overlapping sub-sequences
           :          while query does not.
           :   'qs' = query and sbjct sequence contains overlapping sub-sequences
           :          relative to each other.
           :   '-'  = query and sbjct sequence do not contains multiple domains
           :          relative to each other OR both contain the same distribution
           :          of similar domains.
 Argument  : n/a
 Throws    : n/a
 Comment   : Note: "sbjct" is synonymous with "hit"

=cut

=head2 overlap

See documentation in L<Bio::Search::Hit::HitI::overlap()|Bio::Search::Hit::HitI>

=cut

sub overlap {
    my $self = shift;
    $self->{'_overlap'} = shift if @_;
    return exists $self->{'_overlap'} ? $self->{'_overlap'} : 0;
}


=head2 n

 Usage     : $hit_object->n();
 Purpose   : Gets the N number for the current hit.
           : This is the number of HSPs in the set which was ascribed
           : the lowest P-value (listed on the description line).
           : This number is not the same as the total number of HSPs.
           : To get the total number of HSPs, use num_hsps().
 Example   : $n = $hit_object->n();
 Returns   : Integer
 Argument  : n/a
 Throws    : Exception if HSPs have not been set.
 Comments  : Calling n() on such reports will result in a call to num_hsps().
           : The num_hsps() method will count the actual number of
           : HSPs in the alignment listing, which may exceed N in
           : some cases.

See Also   : L<num_hsps()|num_hsps>

=cut

sub n {
    my $self = shift;
    $self->{'_n'} = shift if @_;
    return exists $self->{'_n'} ? $self->{'_n'} : $self->num_hsps;
}


=head2 p

 Usage     : $hit_object->p( [format] );
 Purpose   : Get the P-value for the best HSP
 Example   : $p =  $sbjct->p;
           : $p =  $sbjct->p('exp');  # get exponent only.
           : ($num, $exp) =  $sbjct->p('parts');  # split sci notation into parts
 Returns   : Float or scientific notation number (the raw P-value, DEFAULT).
           : Integer if format == 'exp' (the magnitude of the base 10 exponent).
           : 2-element list (float, int) if format == 'parts' and P-value
           :                is in scientific notation (See Comments).
 Argument  : format: string of 'raw' | 'exp' | 'parts'
           :    'raw' returns value given in report. Default. (1.2e-34)
           :    'exp' returns exponent value only (34)
           :    'parts' returns the decimal and exponent as a
           :            2-element list (1.2, -34) (See Comments).
 Throws    : Warns if no P-value is defined. Uses expect instead.
 Comments  : Using the 'parts' argument is not recommended since it will not
           : work as expected if the P-value is not in scientific notation.
           : That is, floats are not converted into sci notation before
           : splitting into parts.

See Also   : L<expect()|expect>, L<signif()|signif>, L<Bio::Search::SearchUtils::get_exponent()|Bio::Search::SearchUtils>

=cut

=head2 hsp

 Usage     : $hit_object->hsp( [string] );
 Purpose   : Get a single HSPI object for the present HitI object.
 Example   : $hspObj  = $hit_object->hsp;  # same as 'best'
           : $hspObj  = $hit_object->hsp('best');
           : $hspObj  = $hit_object->hsp('worst');
 Returns   : Object reference for a Bio::Search::HSP::BlastHSP.pm object.
 Argument  : String (or no argument).
           :   No argument (default) = highest scoring HSP (same as 'best').
           :   'best' or 'first' = highest scoring HSP.
           :   'worst' or 'last' = lowest scoring HSP.
 Throws    : Exception if the HSPs have not been collected.
           : Exception if an unrecognized argument is used.

See Also   : L<hsps()|hsps>, L<num_hsps>()

=cut

sub hsp {
    my( $self, $option ) = @_;
    $option ||= 'best';

    if (not ref $self->{'_hsps'}) {
        $self->throw("Can't get HSPs: data not collected.");
    }

    my @hsps = $self->hsps;

    return $hsps[0]      if $option =~ /best|first|1/i;
    return $hsps[$#hsps] if $option =~ /worst|last/i;

    $self->throw("Can't get HSP for: $option\n" .
                 "Valid arguments: 'best', 'worst'");
}

=head2 rank

 Title   : rank
 Usage   : $obj->rank($newval)
 Function: Get/Set the rank of this Hit in the Query search list
           i.e. this is the Nth hit for a specific query
 Returns : value of rank
 Args    : newvalue (optional)


=cut

sub rank {
    my $self = shift;
    return $self->{'_rank'} = shift if @_;
    return $self->{'_rank'} || 1;
}

=head2 locus

 Title   : locus
 Usage   : $locus = $hit->locus();
 Function: Retrieve the locus (if available) for the hit
 Returns : a scalar string (empty string if not set)
 Args    : none

=cut

sub locus {
    my ($self,$value) = @_;
    my $previous = $self->{'_locus'};
    if( defined $value || ! defined $previous ) {
      unless (defined $value) {
        if ($self->{'_name'} =~/(gb|emb|dbj|ref)\|(.*)\|(.*)/) {
                  $value = $previous = $3;
                } else {
          $value = $previous = '';
        }
      }
          $self->{'_locus'} = $value;
    }
        return $previous;
}

=head2 each_accession_number

 Title   : each_accession_number
 Usage   : @each_accession_number = $hit->each_accession_number();
 Function: Get each accession number listed in the description of the hit.
           If there are no alternatives, then only the primary accession will
           be given
 Returns : list of all accession numbers in the description
 Args    : none

=cut

sub each_accession_number {
    my ($self,$value) = @_;
    my $desc = $self->{'_description'};
    #put primary accnum on the list
    my @accnums;
    push (@accnums,$self->{'_accession'});
    if( defined $desc )  {
      while ($desc =~ /(\b\S+\|\S*\|\S*\s?)/g) {
        my $id = $1;
        my ($acc, $version);
	if ($id =~ /(gb|emb|dbj|sp|pdb|bbs|ref|tp[gde])\|(.*)\|(.*)/) {
	    ($acc, $version) = split /\./, $2;
	} elsif ($id =~ /(pir|prf|pat|gnl)\|(.*)\|(.*)/) {
	    ($acc, $version) = split /\./, $3;
	} elsif( $id =~ /(gim|gi|bbm|bbs|lcl)\|(\d*)/) {
	    $acc = $id;
	} elsif( $id =~ /(oth)\|(.*)\|(.*)\|(.*)/ ) { # discontinued...
	    ($acc,$version) = ($2);
	} else {
                     #punt, not matching the db's at ftp://ftp.ncbi.nih.gov/blast/db/README
                     #Database Name                     Identifier Syntax
          #============================      ========================
          #GenBank                           gb|accession|locus
          #EMBL Data Library                 emb|accession|locus
          #DDBJ, DNA Database of Japan       dbj|accession|locus
          #NBRF PIR                          pir||entry
          #Protein Research Foundation       prf||name
          #SWISS-PROT                        sp|accession|entry name
          #Brookhaven Protein Data Bank      pdb|entry|chain
          #Patents                           pat|country|number
          #GenInfo Backbone Id               bbs|number
          #General database identifier           gnl|database|identifier
          #NCBI Reference Sequence           ref|accession|locus
          #Local Sequence identifier         lcl|identifier
              $acc=$id;
            }
            push(@accnums, $acc);
          }
    }
    return @accnums;
}

=head2 tiled_hsps

See documentation in L<Bio::Search::SearchUtils::tile_hsps()|Bio::Search::SearchUtils>

=cut

=head2 query_length

 Title   : query_length
 Usage   : $obj->query_length($newval)
 Function: Get/Set the query_length
 Returns : value of query_length (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub query_length {
    my $self = shift;

    return $self->{'_query_length'} = shift if @_;
    return $self->{'_query_length'};
}

=head2 ncbi_gi

 Title   : ncbi_gi
 Usage   : $acc = $hit->ncbi_gi();
 Function: Retrieve the NCBI Unique ID (aka the GI #),
           if available, for the hit
 Returns : a scalar string (empty string if not set)
 Args    : none

=cut

sub ncbi_gi {
    my ($self,$value) = @_;
    my $previous = $self->{'_ncbi_gi'};
    if( defined $value || ! defined $previous ) {
        $value = $previous = '' unless defined $value;
        $self->{'_ncbi_gi'} = $value;
    }
        return $previous;
}

=head1 ModelHit methods overridden in ModelHit

The following methods have been overridden due to their current reliance on
sequence-based queries. They may be implemented in future versions of this class.

=head2 length_aln

=cut

sub length_aln {
    my $self = shift;
    $self->warn('$hit->length_aln not implemented for Model-based searches');
    return;
}

=head2 gaps

=cut

sub gaps {
    my $self = shift;
    $self->warn('$hit->gaps not implemented for Model-based searches');
    return;
}

=head2 matches

=cut

sub matches {
    my $self = shift;
    $self->warn('$hit->matches not implemented for Model-based searches');
    return;
}

=head2 start

=cut

sub start {
    my $self = shift;
    $self->warn('$hit->start not implemented for Model-based searches');
    return;
}


=head2 end

=cut

sub end {
    my $self = shift;
    $self->warn('$hit->end not implemented for Model-based searches');
    return;
}

=head2 range

=cut

sub range {
    my $self = shift;
    $self->warn('$hit->range not implemented for Model-based searches');
    return;
}

=head2 frac_identical

=cut

sub frac_identical {
    my $self = shift;
    $self->warn('$hit->frac_identical not implemented for Model-based searches');
    return;
}

=head2 frac_conserved

=cut

sub frac_conserved {
    my $self = shift;
    $self->warn('$hit->frac_conserved not implemented for Model-based searches');
    return;
}

=head2 frac_aligned_query

=cut

sub frac_aligned_query {
    my $self = shift;
    $self->warn('$hit->frac_aligned_query not implemented for Model-based searches');
    return;
}

=head2 frac_aligned_hit

=cut

sub frac_aligned_hit {
    my $self = shift;
    $self->warn('$hit->frac_aligned_hit not implemented for Model-based searches');
    return;
}

=head2 num_unaligned_hit

=cut

*num_unaligned_sbjct = \&num_unaligned_hit;

sub num_unaligned_hit {
    my $self = shift;
    $self->warn('$hit->num_unaligned_hit/num_unaligned_sbjct not implemented for Model-based searches');
    return;
}

=head2 num_unaligned_query

=cut

sub num_unaligned_query {
    my $self = shift;
    $self->warn('$hit->num_unaligned_query not implemented for Model-based searches');
    return;
}

=head2 seq_inds

=cut

sub seq_inds {
    my $self = shift;
    $self->warn('$hit->seq_inds not implemented for Model-based searches');
    return;
}

=head2 strand

=cut

sub strand {
    my $self = shift;
    $self->warn('$hit->strand not implemented for Model-based searches');
    return;
}

=head2 frame

=cut

sub frame {
    my $self = shift;
    $self->warn('$hit->frame not implemented for Model-based searches');
    return;
}

=head2 logical_length

=cut

sub logical_length {
    my $self = shift;
    $self->warn('$hit->logical_length not implemented for Model-based searches');
    return;
}

1;
