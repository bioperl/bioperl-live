
#
# BioPerl module for Bio::AnnSeqI
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::AnnSeqI - Abstract Interface of Annotated Sequence

=head1 SYNOPSIS

    # get a Bio::AnnSeqI somehow, eg, via AnnSeqIO

    # features must implement Bio::SeqFeatureI

    @features = $annseq->top_SeqFeatures(); # just top level
    @features = $annseq->all_SeqFeatures(); # descend into sub features

    # get out a specific subset
    sub test_feat {
	my $feat = shift;
        if( $feat->start > 50 && $feat->score < 20 ) {
	    return 1;
	}

	return 0;
    }

    @features = $annseq->fetch_SeqFeatures(\&test_feat);
    @features = $annseq->fetch_SeqFeatures("primary eq exon");

    $seq      = $annseq->seq(); # actual sequence

    $ann      = $annseq->annotation(); # annotation object

=head1 DESCRIPTION

AnnSeqI is the abstract interface of annotated Sequence. These methods
are those which you can be guarenteed to get for any annseq. There aren't
many here, because too many complicated functions here prevent implementations
which are just wrappers around a database or similar delayed mechanisms.

Most of the clever stuff happens inside the SeqFeatureI system.

A good reference implementation is Bio::AnnSeq which is a pure perl 
implementation of this class with alot of extra pieces for extra manipulation.
However, if you want to be able to use any annseq in your analysis, if
you can do it just using these methods, then you know you will be future
proof and compatible with other implementations of AnnSeq.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::AnnSeqI;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

@ISA = qw(Exporter);

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  
  confess "Abstract method '$caller' defined in interface Bio::AnnSeqI not implemented by pacakge $package. Not your fault - author of $package should be blamed!";
}

=head2 top_SeqFeatures

 Title   : top_SeqFeatures
 Usage   : 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub top_SeqFeatures{
   my ($self) = @_;


}


=head2 all_SeqFeatures

 Title   : all_SeqFeatures
 Usage   : @features = $annseq->all_SeqFeatures()
 Function: returns all SeqFeatures, included sub SeqFeatures
 Returns : an array
 Args    : none


=cut

sub all_SeqFeatures{
   my ($self) = @_;
   
   $self->_abstractDeath();

}




