#
# BioPerl module for Bio::Search::Hit::HMMERHit
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Hit::HMMERHit - A Hit module for HMMER hits

=head1 SYNOPSIS

    use Bio::Search::Hit::HMMERHit;
    my $hit = Bio::Search::Hit::HMMERHit->new();
    # use it in the same way as Bio::Search::Hit::GenericHit

=head1 DESCRIPTION

This is a specialization of L<Bio::Search::Hit::GenericHit>.  There
are a few news methods L<next_domain> and L<domains>.  Note that
L<bits> and L<iteration> make no sense for this object and will
return 0.

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

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::Hit::HMMERHit;
use strict;


use base qw(Bio::Search::Hit::GenericHit);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::Hit::HMMERHit->new();
 Function: Builds a new Bio::Search::Hit::HMMERHit object 
 Returns : Bio::Search::Hit::HMMERHit
 Args    : 

 Plus the Bio::Search::Hit::GenericHit inherited params
           -name         => Name of Hit (required)
           -description  => Description (optional)
           -accession    => Accession number (optional)
           -length       => Length of the Hit (optional)
           -score        => Raw Score for the Hit (optional)
           -significance => Significance value for the Hit (optional)
           -algorithm    => Algorithm used (BLASTP, FASTX, etc...)
           -hsps         => Array ref of HSPs for this Hit. 


=cut


=head2 next_domain

 Title   : next_domain 
 Usage   : my $domain = $hit->next_domain();
 Function: An alias for L<next_hsp()>, this will return the next HSP
 Returns : L<Bio::Search::HSP::HSPI> object
 Args    : none


=cut

sub next_domain{ shift->next_hsp }

=head2 domains

 Title   : domains
 Usage   : my @domains = $hit->domains();
 Function: An alias for L<hsps()>, this will return the full list of hsps
 Returns : array of L<Bio::Search::HSP::HSPI> objects
 Args    : none


=cut

sub domains{ shift->hsps() }


=head2 inherited Bio::Search::Hit::GenericHit methods

=cut

=head2 add_hsp

 Title   : add_hsp
 Usage   : $hit->add_hsp($hsp)
 Function: Add a HSP to the collection of HSPs for a Hit
 Returns : number of HSPs in the Hit
 Args    : Bio::Search::HSP::HSPI object


=cut

=head2 Bio::Search::Hit::HitI methods

=cut

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

See Also   : L<score()|score>

=cut

sub bits { return 0 }

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
 Purpose   : Get the number of HSPs for the present Blast hit.
 Example   : $nhsps = $hit_object->num_hsps();
 Returns   : Integer
 Argument  : n/a
 Throws    : Exception if the HSPs have not been collected.

See Also   : L<hsps()|hsps>

=cut

=head2 rewind

 Title   : rewind
 Usage   : $hit->rewind;
 Function: Allow one to reset the HSP iteration to the beginning
           Since this is an in-memory implementation
 Returns : none
 Args    : none

=cut

=head2 iteration

 Title   : iteration
 Usage   : $obj->iteration($newval)
 Function: PSI-BLAST iteration
 Returns : value of iteration
 Args    : newvalue (optional)


=cut


sub iteration { return 0 }

1;
