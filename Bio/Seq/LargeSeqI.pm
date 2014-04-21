#
# BioPerl module for Bio::Seq::LargeSeqI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Albert Vilella
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::LargeSeqI - Interface class for sequences that cache their
residues in a temporary file

=head1 SYNOPSIS

 #

=head1 DESCRIPTION

The interface class defines a group of sequence classes that do not
keep their sequence information in memory but store it in a file. This
makes it possible to work with very large files even with limited RAM.

The most important consequence of file caching for sequences is that
you do not want to inspect the sequence unless absolutely
necessary. These sequences typically override the length() method not
to check the sequence.

The seq() method is not resetable, if you want to add to the end of the
sequence you have to use add_sequence_as_string(), for any other sequence changes you'll
have to create a new object.

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
email or the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Albert Vilella

Email avilella-AT-gmail-DOT-com

=head1 CONTRIBUTORS

Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::LargeSeqI;
use strict;


use base qw(Bio::Root::RootI);


=head2 add_sequence_as_string

 Title   : add_sequence_as_string
 Usage   : $seq->add_sequence_as_string("CATGAT");
 Function: Appends additional residues to an existing  object.
           This allows one to build up a large sequence without
           storing entire object in memory.
 Returns : Current length of sequence
 Args    : string to append

=cut

sub add_sequence_as_string {
    my ($self) = @_;
    $self->throw_not_implemented();
}


1;
