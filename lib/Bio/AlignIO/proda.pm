#
# BioPerl module for Bio::AlignIO::proda
#
#   based on the Bio::SeqIO modules
#       by Ewan Birney <birney@ebi.ac.uk>
#       and Lincoln Stein  <lstein@cshl.org>
#       and the Bio::SimpleAlign module of Ewan Birney
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::proda - proda sequence input/output stream

This provides the basic capabilities to parse the output alignments
from the ProDA multiple sequence alignment program
(http://proda.stanford.edu)

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::AlignIO class.

=head1 DESCRIPTION

This object can transform Bio::Align::AlignI objects to and from proda
files.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHORS - Albert Vilella

Email: avilella-at-gmail-dot-com


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::proda;
use vars qw($LINELENGTH);
use strict;


$LINELENGTH          = 60;
use base qw(Bio::AlignIO);

=head2 new

 Title   : new
 Usage   : $alignio = Bio::AlignIO->new(-format => 'proda',
                       -file => 'filename');
 Function: returns a new Bio::AlignIO object to handle proda files
 Returns : Bio::AlignIO::proda object
 Args    : -verbose => verbosity setting (-1, 0, 1, 2)
           -file    => name of file to read in or to write, with ">"
           -fh      => alternative to -file param - provide a filehandle
                       to read from or write to
           -format  => alignment format to process or produce
           -percentages => display a percentage of identity
                           in each line of the alignment (proda only)
           -linelength=> alignment output line length (default 60)

=cut

sub _initialize {
    my ( $self, @args ) = @_;
    $self->SUPER::_initialize(@args);
    my ( $percentages, $ll ) =
      $self->_rearrange( [qw(PERCENTAGES LINELENGTH)], @args );
    defined $percentages && $self->percentages($percentages);
    $self->line_length( $ll || $LINELENGTH );
}

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream
 Returns : Bio::Align::AlignI object
 Args    : NONE

See L<Bio::Align::AlignI> for details

=cut

sub next_aln {
    my ($self) = @_;
    my $first_line;

    while ( $first_line = $self->_readline ) {
        last if $first_line !~ /^$/;
    }
    $self->_pushback($first_line);
    if ( defined( $first_line = $self->_readline )
        && $first_line !~ /\(/ )
    {
        $self->throw(
            "trying to parse a file which does not start with proda headers"
        );
    } else {
        # use it inside the loop
        $self->_pushback($first_line);
    }
    my %alignments;
    my $aln = Bio::SimpleAlign->new(
        -source  => 'proda',
        -verbose => $self->verbose
    );
    my $order = 0;
    my %order;
    $self->{_lastline} = '';
    my ($first_block, $seen_block, $seen_header) = (0,0,0);
    my @ids; my @ids_copy;
    while ( defined( $_ = $self->_readline ) ) {
        next if (/^\s+$/ && !$first_block);
        if (/^\s$/) {  # line contains no description
            $seen_block = 1;
            next;
        }
        $first_block = 1;

        # break the loop if we come to the end of the current alignment
        # and push back the proda header
        if (/\(/ && $seen_header) {
            $self->_pushback($_);
            last;
        }

        if (/\(/ && !$seen_header) {
            @ids = split(' ', $_);
            $seen_header = 1;
            next;
        }

        my ( $seqname, $aln_line ) = ( '', '' );
        if (/^\s*(\S+)\s*\/\s*(\d+)-(\d+)\s+(\S+)\s*$/ox) {

            # clustal 1.4 format
            ( $seqname, $aln_line ) = ( "$1:$2-$3", $4 );

            # } elsif( /^\s*(\S+)\s+(\S+)\s*$/ox ) { without trailing numbers
        }
        elsif (/^\s*(\S+)\s+(\S+)\s*\d*\s*$/ox) {    # with numbers
            ( $seqname, $aln_line ) = ( $1, $2 );
            if ( $seqname =~ /^[\*\.\+\:]+$/ ) {
                $self->{_lastline} = $_;
                next;
            }
        }
        else {
            $self->{_lastline} = $_;
            next;
        }

        # we ended up the first block and now are on the second
        @ids_copy = @ids unless(defined($ids_copy[0])); #FIXME - hacky
        my $seqname_with_coords = shift(@ids_copy);
        if ($seqname_with_coords !~ /$seqname/) {
            {
                $self->throw("header and body of the alignment dont match");
            }
        }
        $alignments{$seqname_with_coords} .= $aln_line;

        if ( !$seen_block ) {
            if (exists $order{$seqname_with_coords}) {
                $self->warn("Duplicate sequence : $seqname\n".
                            "Can't guarantee alignment quality");
            }
            else {
                $order{$seqname_with_coords} = $order++;
            }
        }

    }

    my ( $sname, $start, $end );
    foreach my $name ( sort { $order{$a} <=> $order{$b} } keys %alignments ) {
        if ( $name =~ /(\S+):(\d+)-(\d+)/ ) {
            ( $sname, $start, $end ) = ( $1, $2, $3 );
        }
        else {
            ( $sname, $start ) = ( $name, 1 );
            my $str = $alignments{$name};
            $str =~ s/[^A-Za-z]//g;
            $end = length($str);
        }
        my $seq = Bio::LocatableSeq->new(
					 -seq      => $alignments{$name},
					 -id       => $sname,
					 -start    => $start,
					 -end      => $end,
					 -alphabet => $self->alphabet,
					 );
        $aln->add_seq($seq);
    }
    
    return $aln if $aln->num_sequences;
    return;
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the proda-format object (.aln) into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Align::AlignI object

=cut

sub write_aln {
    my ($self,@aln) = @_;
    $self->throw_not_implemented();
}

=head2 percentages

 Title   : percentages
 Usage   : $obj->percentages($newval)
 Function: Set the percentages flag - whether or not to show percentages in
           each output line
 Returns : value of percentages
 Args    : newvalue (optional)


=cut

sub percentages {
    my ( $self, $value ) = @_;
    if ( defined $value ) {
        $self->{'_percentages'} = $value;
    }
    return $self->{'_percentages'};
}

=head2 line_length

 Title   : line_length
 Usage   : $obj->line_length($newval)
 Function: Set the alignment output line length
 Returns : value of line_length
 Args    : newvalue (optional)


=cut

sub line_length {
    my ( $self, $value ) = @_;
    if ( defined $value ) {
        $self->{'_line_length'} = $value;
    }
    return $self->{'_line_length'};
}

=head2 no_header

 Title   : no_header
 Usage   : $obj->no_header($newval)
 Function: Set if the alignment input has a CLUSTAL header or not
 Returns : value of no_header
 Args    : newvalue (optional)


=cut

sub no_header {
    my ( $self, $value ) = @_;
    if ( defined $value ) {
        $self->{'_no_header'} = $value;
    }
    return $self->{'_no_header'};
}


1;
