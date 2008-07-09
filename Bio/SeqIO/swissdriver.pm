# make as generic as possible (along with gbhandler, emblhandler)

# Let the code begin...

package Bio::SeqIO::swissdriver;
use vars qw(%FTQUAL_NO_QUOTE);
use strict;
use Bio::SeqIO::Handler::GenericRichSeqHandler;
use Data::Dumper;

use base qw(Bio::SeqIO);

# signals to process what's in the hash prior to next round, maps ann => names 
my %SEC = (
    OC      => 'CLASSIFICATION',
    OH      => 'HOST', # not currently handled, bundled with organism data for now
    OG      => 'ORGANELLE',
    OX      => 'CROSSREF',
    RA      => 'AUTHORS',
    RC      => 'COMMENT',
    RG      => 'CONSRTM',
    RP      => 'POSITION',
    RX      => 'CROSSREF',
    RT      => 'TITLE',
    RL      => 'JOURNAL',
    AS      => 'ASSEMBLYINFO',  # Third party annotation
    '//'    => 'RECORDEND'
    );

# add specialized delimiters here for easier postprocessing
my %DELIM = (
    CC      => "\n",
    DR      => "\n",
    DT      => "\n",
            );

sub _initialize {
    my($self,@args) = @_;

    $self->SUPER::_initialize(@args);
    my $handler = $self->_rearrange([qw(HANDLER)],@args);
    # hash for functions for decoding keys.
    $handler ? $self->seqhandler($handler) :
    $self->seqhandler(Bio::SeqIO::Handler::GenericRichSeqHandler->new(
                    -format => 'swiss',
                    -verbose => $self->verbose,
                    -builder => $self->sequence_builder
                    ));
    if( ! defined $self->sequence_factory ) {
        $self->sequence_factory(Bio::Seq::SeqFactory->new
                (-verbose => $self->verbose(),
                 -type => 'Bio::Seq::RichSeq'));
    }
}

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    : none

=cut

sub next_seq {
    my $self = shift;
    my $hobj = $self->seqhandler;
    local($/) = "\n";
    # these contain values that need to carry over each round
    my ($featkey, $qual, $annkey, $seqdata, $location);
    my $lastann = '';
    my $ct = 0;
    # main parser
    PARSER:
    while(defined(my $line = $self->_readline)) {
        chomp $line;
        my ($ann, $data) = split(m{\s+}, $line, 2);
        if ($ann) {
            if ($ann eq 'FT') {
                # sequence features
                if ($data =~ m{^(\w+)\s+([\d\?\<]+)\s+([\d\?\>]+)(?:\s+?(\S.*))?}ox) {
                    # has location data and desc
                    if ($seqdata) {
                        $hobj->data_handler($seqdata);
                        $seqdata = ();
                    }
                    ($seqdata->{FEATURE_KEY}, my $loc1, my $loc2, $data) = ($1, $2, $3, $4);
                    $qual = 'description';
                    $seqdata->{$qual} = $data;
                    $seqdata->{NAME} = $ann;
                    $seqdata->{LOCATION} = "$loc1..$loc2" if defined $loc1;
                    next PARSER;
                } elsif ($data =~ m{^\s+/([^=]+)(?:=(.+))?}ox) {
                    # has qualifer
                    ($qual, $data) = ($1, $2 || '');
                    $ct = ($seqdata->{$qual}) ? 
                        ((ref($seqdata->{$qual}))  ? scalar(@{ $seqdata->{$qual} }) : 1)
                        : 0 ;
                }
                $data =~ s{\.$}{};
                if ($ct == 0) {
                    $seqdata->{$qual} .= ($seqdata->{$qual}) ?
                        ' '.$data : $data;                    
                } else {
                    if (!ref($seqdata->{$qual})) {
                        $seqdata->{$qual} = [$seqdata->{$qual}];
                    }
                    ($seqdata->{$qual}->[$ct]) ?
                        ($seqdata->{$qual}->[$ct] .= ' '.$data) :
                        ($seqdata->{$qual}->[$ct] .= $data);
                }
            } else {
                # simple annotations
                if ($ann ne $lastann) {
                    if (!$SEC{$ann} && $seqdata) {
                        $hobj->data_handler($seqdata);
                        # can't use undef here; it can lead to subtle mem leaks
                        $seqdata = ();
                    }
                    $annkey = (!$SEC{$ann})    ? 'DATA'     : # primary data
                              $SEC{$ann};
                    $seqdata->{'NAME'} = $ann if !$SEC{$ann};
                }
                last PARSER if $ann eq '//';
                next PARSER if $ann eq 'SQ';
                my $delim = $DELIM{$ann} || ' ';
                $seqdata->{$annkey} .= ($seqdata->{$annkey}) ?
                    $delim.$data : $data;
                $lastann = $ann;
            } 
        } else {
            # this should only be sequence (fingers crossed!)
            SEQUENCE:
            while (defined ($line = $self->_readline)) {
                if (index($line, '//') == 0) {
                    $data =~ tr{0-9 \n}{}d;
                    $seqdata->{DATA} = $data;
                    #$self->debug(Dumper($seqdata));
                    $hobj->data_handler($seqdata);
                    $seqdata = ();
                    last PARSER;
                } else {                        
                    $data .= $line;
                    $line = undef;
                }
            }
        }
    }
    # some files have no // for the last file; this catches the last bit o' data
    $hobj->data_handler($seqdata) if $seqdata;
    return $hobj->build_sequence;
}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq($seq)
 Function: writes the $seq object (must be seq) to the stream
 Returns : 1 for success and 0 for error
 Args    : array of 1 to n Bio::SeqI objects

=cut

sub write_seq {
    shift->throw("Use Bio::SeqIO::swiss write_seq() for output");
    # maybe make a Writer class as well????
}

=head2 seqhandler

 Title   : seqhandler
 Usage   : $stream->seqhandler($handler)
 Function: Get/Set teh Bio::Seq::HandlerBaseI object
 Returns : Bio::Seq::HandlerBaseI 
 Args    : Bio::Seq::HandlerBaseI 

=cut

sub seqhandler {
    my ($self, $handler) = @_;
    if ($handler) {
        $self->throw("Not a Bio::HandlerBaseI") unless
        ref($handler) && $handler->isa("Bio::HandlerBaseI");
        $self->{'_seqhandler'} = $handler;
    }
    return $self->{'_seqhandler'};
}

1;

__END__

