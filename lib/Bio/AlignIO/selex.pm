#
# BioPerl module for Bio::AlignIO::selex

#   based on the Bio::SeqIO::selex module
#       by Ewan Birney <birney@ebi.ac.uk>
#       and Lincoln Stein  <lstein@cshl.org>
#
#       and the SimpleAlign.pm module of Ewan Birney
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself
# _history
# September 5, 2000
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::selex - selex sequence input/output stream

=head1 SYNOPSIS

  # Do not use this module directly.  Use it via the L<Bio::AlignIO> class.

  use Bio::AlignIO;
  use strict;

  my $in = Bio::AlignIO->new(-format => 'selex',
                             -file   => 't/data/testaln.selex');
  while( my $aln = $in->next_aln ) {

  }

=head1 DESCRIPTION

This object can transform L<Bio::Align::AlignI> objects to and from selex flat
file databases.

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
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHORS - Peter Schattner

Email: schattner@alum.mit.edu

=head1 CONTRIBUTORS

Jason Stajich, jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::selex;
use strict;

use base qw(Bio::AlignIO);

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream. Tries to read *all* selex
          It reads all non whitespace characters in the alignment
          area. For selexs with weird gaps (eg ~~~) map them by using
          $al->map_chars('~','-')
 Returns : L<Bio::Align::AlignI> object
 Args    : NONE

=cut

sub next_aln {
    my $self = shift;
    my $entry;
    my ($start,$end,%align,$name,$seqname,%hash,@c2name, %accession,%desc);
    my $aln =  Bio::SimpleAlign->new(-source => 'selex');

    # in selex format, every non-blank line that does not start
    # with '#=' is an alignment segment; the '#=' lines are mark up lines.
    # Of particular interest are the '#=GF <name/st-ed> AC <accession>'
    # lines, which give accession numbers for each segment
    while( $entry = $self->_readline) {
        if( $entry =~ /^\#=GS\s+(\S+)\s+AC\s+(\S+)/ ) {
	    $accession{ $1 } = $2;
	    next;
	} elsif( $entry =~ /^\#=GS\s+(\S+)\s+DE\s+(.+)\s*$/ ) {
	    $desc{$1} .= $2;
	} elsif ( $entry =~ /^([^\#]\S+)\s+([A-Za-z\.\-\*]+)\s*/ ) {
	    my ($name,$seq) = ($1,$2);

	    if( ! defined $align{$name}  ) {
		push @c2name, $name;
	    }
	    $align{$name} .= $seq;
	}
    }
    # ok... now we can make the sequences

    foreach my $name ( @c2name ) {

	if( $name =~ /(\S+)\/(\d+)-(\d+)/ ) {
	    $seqname = $1;
	    $start = $2;
	    $end = $3;
	} else {
	    $seqname=$name;
	    $start = 1;
	    $end = length($align{$name});
	}
	my $seq = Bio::LocatableSeq->new
	    ('-seq'              => $align{$name},
	     '-display_id'       => $seqname,
	     '-start'            => $start,
	     '-end'              => $end,
	     '-description'      => $desc{$name},
	     '-accession_number' => $accession{$name},
 	     '-alphabet'         => $self->alphabet,
	     );

	$aln->add_seq($seq);
    }

    return $aln if $aln->num_sequences;
    return;
}


=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in selex format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object


=cut

sub write_aln {
    my ($self,@aln) = @_;
    my ($namestr,$seq,$add);
    my ($maxn);
    foreach my $aln (@aln) {
	$maxn = $aln->maxdisplayname_length();
	foreach $seq ( $aln->each_seq() ) {
	    $namestr = $aln->displayname($seq->get_nse());
	    $add = $maxn - length($namestr) + 2;
	    $namestr .= " " x $add;
	    $self->_print (sprintf("%s  %s\n",$namestr,$seq->seq())) or return;
	}
    }
    $self->flush if $self->_flush_on_write && defined $self->_fh;
    return 1;
}

1;
