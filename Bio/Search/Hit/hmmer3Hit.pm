# $Id: bioperl.lisp 15559 2009-02-23 12:11:20Z maj $
#
# BioPerl module for Bio::Search::Hit::hmmer3Hit
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Thomas Sharpton <thomas.sharpton@gmail.com>
#
# Copyright Thomas Sharpton
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Hit::hmmer3Hit - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support

Please direct usage questions or support issues to the mailing list:

L<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Thomas Sharpton

Email thomas.sharpton@gmail.com

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::Hit::hmmer3Hit;
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

=head2 bits

 Usage     : $hit_object->bits();
 Purpose   : Gets the bit score of the best HSP for the current hit.
 Example   : $bits = $hit_object->bits();
 Returns   : Integer or undef if bit score is not set
 Argument  : n/a

See Also   : L<score()|score>

=cut

sub bits { return 0 }

=head2 iteration

 Title   : iteration
 Usage   : $obj->iteration($newval)
 Function: PSI-BLAST iteration
 Returns : value of iteration
 Args    : newvalue (optional)


=cut

sub iteration { return 0 }

1;
