# $Id$
#
# BioPerl module for Bio::Factory::SequenceBuilderI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Factory::SequenceBuilderI - This interface allows for generic building of sequences in factories which create sequences (like SeqIO) 

=head1 SYNOPSIS

# do not use this object directly it is an interface
# get a Bio::Factory::SequenceBuilderI object like

    use Bio::Seq::PrimarySequenceBuilder;
    my $seqbuilder = new Bio::Seq:PrimarySequenceBuilder();

    my $seq = $seqbuilder->new_sequence(-seq => 'ACTGAT',
					-display_id => 'exampleseq');

    print "seq is a ", ref($seq), "\n";

=head1 DESCRIPTION

Describe the interface here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Factory::SequenceBuilderI;
use vars qw(@ISA);
use strict;
@ISA = qw(Bio::Root::RootI);

=head2 new_sequence

 Title   : new_sequence
 Usage   : my $seq = $seqbuilder->new_sequence(-seq => 'CAGT', -id => 'name');
 Function: Instantiates new Bio::PrimarySeqI (or one of its child classes)
           This object allows us to genericize the instantiation of sequence
           objects.
 Returns : Bio::PrimarySeqI
 Args    : initialization parameters specific to the type of sequence
           object we want.  Typically 
           -seq        => $str,
           -display_id => $name

=cut

sub new_sequence{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

1;
