#
# BioPerl module for Bio::SeqEvolution::EvolutionI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki at bioperl dot org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqEvolution::EvolutionI - the interface for evolving sequences

=head1 SYNOPSIS

    # not an instantiable class

=head1 DESCRIPTION

This is the interface that all classes that mutate sequence objects in
constant fashion must implement. A good example is
Bio::SeqEvolution::DNAPoint.

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

  Heikki Lehvaslaiho E<lt>heikki at bioperl dot orgE<gt>

=head1 CONTRIBUTORS

Additional contributor's names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqEvolution::EvolutionI;
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


=head2 seq

 Title   : seq
 Usage   : $obj->seq($newval)
 Function: Set the sequence object for the original sequence
 Returns : The sequence object
 Args    : newvalue (optional)

Setting this will reset mutation and generated mutation counters.

=cut

sub seq { shift->throw_not_implemented(); }

=head2 next_seq

  Title   : next_seq
  Usage   : $obj->next_seq
  Function: Evolve the reference sequence to desired level
  Returns : A new sequence object mutated from the reference sequence
  Args    : -

=cut

sub next_seq{ shift->throw_not_implemented(); }


=head2 mutate

  Title   : mutate
  Usage   : $obj->mutate
  Function: mutate the sequence at the given location according to the model
  Returns : true
  Args    : integer, start location of the mutation, required argument

Called from next_seq().

=cut

sub mutate { shift->throw_not_implemented(); }


1;
