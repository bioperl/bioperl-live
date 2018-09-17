#
# BioPerl module for Bio::Factory::SequenceFactoryI
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

Bio::Factory::SequenceFactoryI - This interface allows for generic building of sequences in factories which create sequences (like SeqIO) 

=head1 SYNOPSIS

# do not use this object directly it is an interface
# get a Bio::Factory::SequenceFactoryI object like

    use Bio::Seq::SeqFactory;
    my $seqbuilder = Bio::Seq::SeqFactory->new('-type' => 'Bio::PrimarySeq');

    my $seq = $seqbuilder->create(-seq => 'ACTGAT',
				  -display_id => 'exampleseq');

    print "seq is a ", ref($seq), "\n";

=head1 DESCRIPTION

A generic way to build Sequence objects via a pluggable factory.  This
reduces the amount of code that looks like

  if( $type  eq 'Bio::PrimarySeq' ) { ... } 
  elsif( $type eq 'Bio::Seq::RichSeq' ) { ... }

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

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Factory::SequenceFactoryI;

use strict;

use base qw(Bio::Factory::ObjectFactoryI);

=head2 create

 Title   : create
 Usage   : my $seq = $seqbuilder->create(-seq => 'CAGT', 
					 -id => 'name');
 Function: Instantiates new Bio::PrimarySeqI (or one of its child classes)
           This object allows us to genericize the instantiation of sequence
           objects.
 Returns : Bio::PrimarySeqI
 Args    : initialization parameters specific to the type of sequence
           object we want.  Typically 
           -seq        => $str,
           -display_id => $name

=cut

1;
