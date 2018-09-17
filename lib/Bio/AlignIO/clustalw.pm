#
# BioPerl module for Bio::AlignIO::clustalw
#
#   based on the Bio::SeqIO modules
#       by Ewan Birney <birney@ebi.ac.uk>
#       and Lincoln Stein  <lstein@cshl.org>
#       and the Bio::SimpleAlign module of Ewan Birney
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself
# History
# September 5, 2000
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::clustalw - clustalw sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::AlignIO class.

=head1 DESCRIPTION

This object can transform Bio::Align::AlignI objects to and from clustalw
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

=head1 AUTHORS - Peter Schattner

Email: schattner@alum.mit.edu


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::clustalw;
use vars qw($LINELENGTH $CLUSTALPRINTVERSION);
use strict;


$LINELENGTH          = 60;
$CLUSTALPRINTVERSION = '1.81';
use base qw(Bio::AlignIO);

=head2 new

 Title   : new
 Usage   : $alignio = Bio::AlignIO->new(-format => 'clustalw',
                       -file => 'filename');
 Function: returns a new Bio::AlignIO object to handle clustalw files
 Returns : Bio::AlignIO::clustalw object
 Args    : -verbose => verbosity setting (-1, 0, 1, 2)
           -file    => name of file to read in or to write, with ">"
           -fh      => alternative to -file param - provide a filehandle
                       to read from or write to
           -format  => alignment format to process or produce
           -percentages => display a percentage of identity
                           in each line of the alignment (clustalw only)
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
        && $first_line !~ /CLUSTAL/ )
    {
        $self->throw(
            "trying to parse a file which does not start with a CLUSTAL header"
        );
    }
    my %alignments;
    my $aln = Bio::SimpleAlign->new(
        -source  => 'clustalw',
        -verbose => $self->verbose
    );
    my $order = 0;
    my %order;
    $self->{_lastline} = '';
    my ($first_block, $seen_block) = (0,0);
    while ( defined( $_ = $self->_readline ) ) {
        next if (/^\s+$/ && !$first_block);
        if (/^\s$/) {  # line contains no description
            $seen_block = 1;
            next;
        }
        $first_block = 1;
        # break the loop if we come to the end of the current alignment
        # and push back the CLUSTAL header
        if (/CLUSTAL/) {
            $self->_pushback($_);
            last;
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

        if ( !$seen_block ) {
            if (exists $order{$seqname}) {
                $self->warn("Duplicate sequence : $seqname\n".
                            "Can't guarantee alignment quality");
            }
            else {
                $order{$seqname} = $order++;
            }
        }

        $alignments{$seqname} .= $aln_line;
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
        my $seq = Bio::LocatableSeq->new
	    (
	     '-seq'         => $alignments{$name},
	     '-display_id'  => $sname,
	     '-start'       => $start,
	     '-end'         => $end,
	    '-alphabet'     => $self->alphabet,
        );
        $aln->add_seq($seq);
    }

    # not sure if this should be a default option - or we can pass in
    # an option to do this in the future? --jason stajich
    # $aln->map_chars('\.','-');
    
    # no sequences added, so just return
    return $aln if $aln->num_sequences;
    return;
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the clustalw-format object (.aln) into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Align::AlignI object

=cut

sub write_aln {
    my ( $self, @aln ) = @_;
    my ( $count, $length, $seq, @seq, $tempcount, $line_len );
    $line_len = $self->line_length || $LINELENGTH;
    foreach my $aln (@aln) {
        if ( !$aln || !$aln->isa('Bio::Align::AlignI') ) {
            $self->warn(
"Must provide a Bio::Align::AlignI object when calling write_aln"
            );
            next;
        }
        my $matchline = $aln->match_line;
        if ( $self->force_displayname_flat ) {
            $aln->set_displayname_flat(1);
        }
        $self->_print(
            sprintf( "CLUSTAL W (%s) multiple sequence alignment\n\n\n",
                $CLUSTALPRINTVERSION )
        ) or return;
        $length = $aln->length();
        $count  = $tempcount = 0;
        @seq    = $aln->each_seq();
        my $max = 22;
        foreach $seq (@seq) {
            $max = length( $aln->displayname( $seq->get_nse() ) )
              if ( length( $aln->displayname( $seq->get_nse() ) ) > $max );
        }

        while ( $count < $length ) {
            my ( $linesubstr, $first ) = ( '', 1 );
            foreach $seq (@seq) {

              #
              #  Following lines are to suppress warnings
              #  if some sequences in the alignment are much longer than others.

                my ($substring);
                my $seqchars = $seq->seq();
              SWITCH: {
                    if ( length($seqchars) >= ( $count + $line_len ) ) {
                        $substring = substr( $seqchars, $count, $line_len );
                        if ($first) {
                            $linesubstr =
                              substr( $matchline, $count, $line_len );
                            $first = 0;
                        }
                        last SWITCH;
                    }
                    elsif ( length($seqchars) >= $count ) {
                        $substring = substr( $seqchars, $count );
                        if ($first) {
                            $linesubstr = substr( $matchline, $count );
                            $first = 0;
                        }
                        last SWITCH;
                    }
                    $substring = "";
                }
                $self->_print(
                    sprintf(
                        "%-" . $max . "s %s\n",
                        $aln->displayname( $seq->get_nse() ), $substring
                    )
                ) or return;
            }

            my $percentages = '';
            if ( $self->percentages ) {
                my ($strcpy) = ($linesubstr);
                my $count = ( $strcpy =~ tr/\*// );
                $percentages =
                  sprintf( "\t%d%%", 100 * ( $count / length($linesubstr) ) );
            }
            $self->_print(
                sprintf(
                    "%-" . $max . "s %s%s\n",
                    '', $linesubstr, $percentages
                )
            );
            $self->_print( sprintf("\n\n") ) or return;
            $count += $line_len;
        }
    }
    $self->flush if $self->_flush_on_write && defined $self->_fh;
    return 1;
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

1;
