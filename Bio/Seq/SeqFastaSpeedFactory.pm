# $Id$
#
# BioPerl module for Bio::Seq::SeqFactory
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::SeqFactory - Instantiates a new Bio::PrimarySeqI (or derived class) through a factory

=head1 SYNOPSIS

    use Bio::Seq::SeqFactory;
    my $factory = new Bio::Seq::SeqFactory;
    my $seq = $factory->create(-seq => 'WYRAVLC',
			       -id  => 'name');

    # If you want the factory to create Bio::Seq objects instead
    # of the default Bio::PrimarySeq objects, use the -type parameter:

    my $factory = new Bio::Seq::SeqFactory(-type => 'Bio::Seq');


=head1 DESCRIPTION

This object will build Bio::Seq objects generically.

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


package Bio::Seq::SeqFastaSpeedFactory;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::Factory::SequenceFactoryI;
use Bio::Seq;
use Bio::PrimarySeq;

@ISA = qw(Bio::Root::Root Bio::Factory::SequenceFactoryI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Seq::SeqFactory();
 Function: Builds a new Bio::Seq::SeqFactory object 
 Returns : Bio::Seq::SeqFactory
 Args    : -type => string, name of a PrimarySeqI derived class
                    This is optional. Default=Bio::PrimarySeq.

=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  return $self;
}


=head2 create

 Title   : create
 Usage   : my $seq = $seqbuilder->create(-seq => 'CAGT', -id => 'name');
 Function: Instantiates a new Bio::Seq object, correctly built but very
           fast, knowing stuff about Bio::PrimarySeq and Bio::Seq
 Returns : Bio::Seq

 Args    : initialization parameters specific to the type of sequence
           object we want.  Typically 
           -seq        => $str,
           -id         => $name

=cut

sub create {
    my ($self,%param) = @_;

    my $sequence = $param{'-seq'}  || $param{'-SEQ'};
    my $fulldesc = $param{'-desc'} || $param{'-DESC'};
    my $id       = $param{'-id'}   || $param{'-ID'};



    my $seq = {};
    bless $seq,"Bio::Seq";
    my $t_pseq = $seq->{'primary_seq'} = {};
    bless $t_pseq,"Bio::PrimarySeq";
    $t_pseq->{'seq'}  = $sequence;
    $t_pseq->{'desc'} = $fulldesc;
    $t_pseq->{'display_id'} = $id;
    $t_pseq->{'primary_id'} = $id;
    if( $sequence ) {
	$t_pseq->_guess_alphabet();
    }

    return $seq;
}

1;

