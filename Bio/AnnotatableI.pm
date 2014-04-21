#
# BioPerl module for Bio::AnnotatableI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::AnnotatableI - the base interface an annotatable object must implement

=head1 SYNOPSIS

    use Bio::SeqIO;
    # get an annotatable object somehow: for example, Bio::SeqI objects
    # are annotatable
    my $seqio = Bio::SeqIO->new(-fh => \*STDIN, -format => 'genbank');
    while (my $seq = $seqio->next_seq()) {
        # $seq is-a Bio::AnnotatableI, hence:
        my $ann_coll = $seq->annotation();
        # $ann_coll is-a Bio::AnnotationCollectionI, hence:
        my @all_anns = $ann_coll->get_Annotations();
        # do something with the annotation objects
    }

=head1 DESCRIPTION

This is the base interface that all annotatable objects must implement. A 
good example is Bio::Seq which is an AnnotableI object.

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

=head1 AUTHOR

 Hilmar Lapp E<lt>hlapp@gmx.netE<gt>
 Allen Day E<lt>allenday@ucla.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::AnnotatableI;
use strict;

use base qw(Bio::Root::RootI);

=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($newval)
 Function: Get the annotation collection for this annotatable object.
 Example : 
 Returns : a Bio::AnnotationCollectionI implementing object, or undef
 Args    : on set, new value (a Bio::AnnotationCollectionI
           implementing object, optional) (an implementation may not
           support changing the annotation collection)

See L<Bio::AnnotationCollectionI>

=cut

sub annotation{
  shift->throw_not_implemented();
}

1;
