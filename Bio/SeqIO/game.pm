# $Id$
#
# BioPerl module for Bio::SeqIO::game
#
# Cared for by Sheldon McKay <smckay@bcgsc.bc.ca>
#
# Copyright Sheldon McKay
#
# You may distribute this module under the same terms as perl itself
#

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::game -- a class for parsing and writing game-XML

=head1 SYNOPSIS

 use Bio::SeqIO;

 my $in = Bio::SeqIO->new ( -file => 'file.xml', -format => 'game' );

 $in->verbose(1) # switch on non-fatal error messages

 my $seq = $in->next_seq;

=head1 DESCRIPTION

Bio::SeqIO::game will parse game XML (version 1.2) or write game XML from 
another  Bio::SeqI format.  The game-XML format current;y used by apollo contains
a single 'main' annotated sequence, so we do not really get a sequence stream.

This modules is not used directly

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
to one of the Bioperl mailing lists.

Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/MailList.shtml      - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution.

Bug reports can be submitted via email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Sheldon McKay

Email smckay@bcgsc.bc.ca

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::SeqIO::game;

use Data::Dumper;
use Bio::SeqIO;
use Bio::SeqIO::game::gameHandler;
use Bio::SeqIO::game::gameWriter;
use strict;

use vars qw{ @ISA };
@ISA = qw { Bio::SeqIO };

=head2 _initialize

 Title   : _initialize
 Function: passed arguments to the XML handler

=cut

sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);
}

=head2 next_seq

 Title   : next_seq
 Usage   : my $seq = $seqio->next_seq;
 Function: get the main sequence object
 Returns : a Bio::Seq::RichSeq object
 Args    : none

=cut

sub next_seq {
    my $self   = shift;
    
    my $seq_l  = $self->_getseqs;
    my $annseq = shift @{$seq_l};
    my $seq    = $annseq->[0];
    my $feats  = $annseq->[1];
    
    for ( @{$feats} ) {
	$seq->add_SeqFeature( $_ );
    }

    return $seq;
}   

=head2 write_seq

 Title   : write_seq
 Usage   : $seqio->write_seq($seq)
 Function: writes a sequence object as game XML
 Returns : nothing
 Args    : a Bio::SeqI compliant object

=cut

sub write_seq {
    my ($self, $seq) = @_;
    my $writer = Bio::SeqIO::game::gameWriter->new($seq);
    $writer->write_to_game;
}

=head2 _getseqs

 Title   : _getseqs
 Usage   : $self->_getseqs
 Function: An internal method to invoke the PerlSAX XML handler and retrieve
           the sequence objects
 Returns : an reference to an array with sequence object and annotations
 Args    : none

=cut

sub _getseqs {
    my $self = shift;
    my $file = $self->{_file}
      || $self->throw("No input file specified");

    # check to see if there is a source feature
    open XML, $file;
    my $text = join '', <XML>;
    close XML;
    my $source = $text =~ /type>(source|origin|region)<\/type/gm ? 1 : 0;

    # check for verbose reporting
    my $verbose = $self->{verbose};

    if ( defined $self->{seq_l} ) {
        return $self->{seq_l};
    }
    else {
        my $handler = Bio::SeqIO::game::gameHandler->new;
	$handler->{has_source} = $source if $source;
	$handler->{verbose} = 1 if $verbose;
        my $parser  = XML::Parser::PerlSAX->new( Handler => $handler );
        my $game    = $parser->parse(Source => { SystemId => $file });
	$self->{seq_l} = $game->load;
    }
}

=head2 _hide_dna

 Title   : _hide_dna
 Usage   : $seqio->_hide_dna
 Function: Hide the DNA for really huge sequences
 Returns : nothing 
 Args    : none

=cut

sub _hide_dna {
    my $self = shift;
    
    my $annseqs = $self->_getseqs;

    for ( @{$annseqs} ) {
        my $seq = $_->[0];
        $seq->seq('');
    }
    return 0;
}


=head2 verbose

 Title   : verbose
 Usage   : $seqio->verbose(1)
 Function: turn on non-fatal warnings
 Returns : true if verbosity it set to 'on'
 Args    : a true value (optional)

=cut

sub verbose {
    my ($self, $v) = @_;

    if ($v) {
	$self->{verbose} = $v;
    }
    
    $self->{verbose};
}

1;
