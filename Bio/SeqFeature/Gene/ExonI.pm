#
# BioPerl module for Bio::SeqFeature::Gene::ExonI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Gene::ExonI - Interface for a feature representing an exon

=head1 SYNOPSIS

See documentation of methods.

=head1 DESCRIPTION

A feature representing an exon. An exon in this definition is
transcribed and at least for one particular transcript not spliced out
of the pre-mRNA. However, it does not necessarily code for amino acid.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

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

=head1 AUTHOR - Hilmar Lapp

Email hlapp@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Gene::ExonI;
use strict;

use base qw(Bio::SeqFeatureI);


=head2 is_coding

 Title   : is_coding
 Usage   : if($exon->is_coding()) {
                   # do something
           }
 Function: Whether or not the exon codes for amino acid.
 Returns : TRUE if the object represents a feature translated into protein,
           and FALSE otherwise.
 Args    : 


=cut

sub is_coding {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 cds

 Title   : cds()
 Usage   : $cds = $exon->cds();
 Function: Get the coding sequence of the exon as a sequence object.

           The returned sequence object must be in frame 0, i.e., the first
           base starts a codon.

           An implementation may return undef, indicating that a coding
           sequence does not exist, e.g. for a UTR (untranslated region).

 Returns : A L<Bio::PrimarySeqI> implementing object.
 Args    : 


=cut

sub cds {
    my ($self) = @_;
    $self->throw_not_implemented();
}

1;
