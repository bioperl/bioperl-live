#
# BioPerl module for Bio::Seq::SeqFastaSpeedFactory
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

Bio::Seq::SeqFastaSpeedFactory - Rapid creation of Bio::Seq objects through a factory

=head1 SYNOPSIS

    use Bio::Seq::SeqFastaSpeedFactory;
    my $factory = Bio::Seq::SeqFastaSpeedFactory->new();
    my $seq = $factory->create( -seq => 'WYRAVLC',
                                -id  => 'name'     );

=head1 DESCRIPTION

This factory was designed to build Bio::Seq objects as quickly as possible, but
is not as generic as L<Bio::Seq::SeqFactory>. It can be used to create sequences
from non-rich file formats. The L<Bio::SeqIO::fasta> sequence parser uses this
factory.

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


package Bio::Seq::SeqFastaSpeedFactory;
use strict;

use Bio::Seq;
use Bio::PrimarySeq;

use base qw(Bio::Root::Root Bio::Factory::SequenceFactoryI);


=head2 new

 Title   : new
 Usage   : my $obj = Bio::Seq::SeqFastaSpeedFactory->new();
 Function: Builds a new Bio::Seq::SeqFastaSpeedFactory object 
 Returns : Bio::Seq::SeqFastaSpeedFactory
 Args    : None

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
 Returns : A Bio::Seq object
 Args    : Initialization parameters for the sequence object we want:
              -id
              -primary_id
              -display_id
              -desc
              -seq
              -alphabet

=cut

sub create {
    my ($self,@args) = @_;
    
    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
    
    my $sequence = $param{'-seq'};
    my $fulldesc = $param{'-desc'};
    my $id       = defined $param{'-id'} ? $param{'-id'} : $param{'-primary_id'};
    my $alphabet = $param{'-alphabet'};

    my $seq = bless {}, 'Bio::Seq';
    my $t_pseq = $seq->{'primary_seq'} = bless {}, 'Bio::PrimarySeq';
    $t_pseq->{'seq'}        = $sequence;
    $t_pseq->{'length'}     = CORE::length($sequence);
    $t_pseq->{'desc'}       = $fulldesc;
    $t_pseq->{'display_id'} = $id;
    $t_pseq->{'primary_id'} = $id;
    $seq->{'primary_id'}    = $id; # currently Bio::Seq does not delegate this
    if( $sequence and !$alphabet ) {
        $t_pseq->_guess_alphabet();
    } elsif ( $sequence and $alphabet ) {
        $t_pseq->{'alphabet'} = $alphabet;
    }

    return $seq;
}

1;

