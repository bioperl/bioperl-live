#
# BioPerl module for Bio::AlignIO::metafasta
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::metafasta - Metafasta MSA Sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::AlignIO> class.

=head1 DESCRIPTION

This object can transform L<Bio::SimpleAlign> objects to and from
metafasta flat file databases.

The format of a metafasta file is

  >test/1-25
  ABCDEFHIJKLMNOPQRSTUVWXYZ
  &charge
  NBNAANCNJCNNNONNCNNUNNXNZ
  &chemical
  LBSAARCLJCLSMOIMCHHULRXRZ

where the sequence block is followed by one or several meta blocks.
Each meta block starts with the ampersand character '&' in the first
column and is immediately followed by the name of the meta data which
continues until the new line. The meta data follows it. All
characters, except new line, are important in meta data.

=head1 SEE ALSO

L<Bio::SeqIO::metafasta>

=head1 FEEDBACK

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

=head1 AUTHOR - Heikki Lehvaslaiho

Email heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::metafasta;
use vars qw($WIDTH);
use strict;

use Bio::SimpleAlign;
use Bio::Seq::Meta;
use Bio::Seq::SeqFactory;
use Bio::Seq::SeqFastaSpeedFactory;

use base qw(Bio::AlignIO);

BEGIN { $WIDTH = 60}

sub _initialize {
  my($self,@args) = @_;
  $self->SUPER::_initialize(@args);
  my ($width) = $self->_rearrange([qw(WIDTH)], @args);
  $width && $self->width($width);
}

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
 Returns : L<Bio::Align::AlignI> object - returns 0 on end of file
	    or on error
 Args    : NONE

=cut

sub next_aln {
    my( $self ) = @_;
    my $seq;
    my $alphabet;
    local $/ = "\n>";

    my $aln =  Bio::SimpleAlign->new();

    while(defined (my $entry = $self->_readline)) {
        chomp($entry);
        if ($entry =~ m/\A\s*\Z/s) { # very first one
            return unless $entry = $self->_readline;
            chomp($entry);
        }
        $entry =~ s/^>//;

        my ($top,$sequence) = split(/\n/,$entry,2);
        defined $sequence && $sequence =~ s/>//g;

        my @metas;
        ($sequence, @metas) = split /\n&/, $sequence;

        my ($id, $start, $end);
        if ( $top =~ /(\S+)\/(\d+)-(\d+)/ ) {
            $id = $1;
            $start = $2;
            $end = $3;
        }
        elsif ($top =~ /(\S+)/) {
            $id = $1;
            $start = 1;
            $end = length($sequence);
        }

        defined $sequence && $sequence =~ s/\s//g; # Remove whitespace

        $seq = Bio::Seq::Meta->new('-seq'        => $sequence,
				   '-display_id' => $id,
				   '-start'      => $start,
				   '-end'        => $end,
				   '-alphabet'   => $self->alphabet,
				   );

        foreach my $meta (@metas) {
            my ($name,$string) = split /\n/, $meta;
            $string =~ s/\n//g;	# Remove newlines, spaces are important
            $seq->named_meta($name, $string);
        }

	$aln->add_seq($seq);
	
	# alignment needs seqs all the same length, pad with gaps
	my $alnlen = $aln->length;
	foreach my $seq ( $aln->each_seq ) {
		if ( $seq->length < $alnlen ) {
			my ($diff) = ($alnlen - $seq->length);
			$seq->seq( $seq->seq() . "-" x $diff);
		}
	}
    }
    return $aln if $aln->num_sequences;
	return;
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in fasta format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object

=cut

sub write_aln {
    my ($self,@aln) = @_;
    my $width = $self->width;

    foreach my $aln (@aln) {
	if( ! $aln || ! $aln->isa('Bio::Align::AlignI')  ) {
	    $self->warn("Must provide a Bio::Align::AlignI object when calling write_aln");
	    next;
	}
	foreach my $seq ( $aln->each_seq() ) {
	    my $name = $aln->displayname($seq->get_nse);

	    my $str  = $seq->seq();
            if(length($str) > 0) {
                $str =~ s/(.{1,$width})/$1\n/g;
            } else {
                $str = "\n";
            }
            $self->_print (">",$name,"\n",$str) or return;
            if ($seq->isa('Bio::Seq::MetaI')) {
                foreach my $meta ($seq->meta_names) {
                    my $str = $seq->named_meta($meta);
                    $str =~ s/(.{1,$width})/$1\n/g;
                    $self->_print ("&",$meta,"\n",$str);
                }
            }
	}
    }
    $self->flush if $self->_flush_on_write && defined $self->_fh;
    return 1;
}


=head2 width

 Title   : width
 Usage   : $obj->width($newval)
 Function: Get/Set the line width for METAFASTA output
 Returns : value of width
 Args    : newvalue (optional)


=cut

sub width{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'width'} = $value;
    }
    return $self->{'width'} || $WIDTH;
}

1;
