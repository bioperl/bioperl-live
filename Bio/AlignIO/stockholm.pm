# $Id$
#
# BioPerl module for Bio::AlignIO::stockholm

#   based on the Bio::SeqIO::stockholm module
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

Bio::AlignIO::stockholm - stockholm sequence input/output stream

=head1 SYNOPSIS

  # Do not use this module directly.  Use it via the L<Bio::AlignIO> class.

  use Bio::AlignIO;
  use strict;

  my $in = Bio::AlignIO->new(-format => 'stockholm',
                             -file   => 't/data/testaln.stockholm');
  while( my $aln = $in->next_aln ) {

  }

=head1 DESCRIPTION

This object can transform L<Bio::Align::AlignI> objects to and from
stockholm flat file databases.

=head1 FEEDBACK

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHORS - Peter Schattner

Email: schattner@alum.mit.edu

=head1 CONTRIBUTORS

Andreas Kahari, ak@ebi.ac.uk
Jason Stajich, jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::stockholm;
use strict;


use base qw(Bio::AlignIO);

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
 Returns : L<Bio::Align::AlignI> object
 Args    : NONE

=cut

sub next_aln {
    my $self = shift;
    my $entry;

    my ($start,$end,%align,$name,$seqname,$seq,$count,%desc,
	%hash,@c2name, %accession);

    # in stockholm format, every non-blank line that does not start
    # with '#=' is an alignment segment; the '#=' lines are mark up lines.
    # Of particular interest are the '#=GF <name/st-ed> AC <accession>'
    # lines, which give accession numbers for each segment

    my $aln =  Bio::SimpleAlign->new(-source => 'stockholm');

    while( defined($entry = $self->_readline) ) {
        next unless $entry =~ /\w+/;

	if ($entry =~ /^#\s*STOCKHOLM\s+/) {
	    last;
	} else {
	    $self->throw("Not Stockholm format: Expecting \"# STOCKHOLM 1.0\"; Found \"$_\"");
	}
    }
#
#    Next section is same as for selex format
#
    while( defined($entry = $self->_readline) ) {
	# Double slash (//) signals end of file.  The flat Pfam-A data from
	# ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Pfam-A.full.gz consists
	# of several concatenated Stockholm-formatted files.  The following
	# line makes it possible to parse it without this module trying to
	# read the whole file into memory.  Andreas Kähäri 10/3/2003.
	last if $entry =~ m{^//};

	# Extra bonus:  Get the name of the alignment.
	# Andreas Kähäri 10/3/2003.
	if ($entry =~ /^\#=GF\s+AC\s+(\S+)/) {
	    $aln->id($1);
	} elsif( $entry =~ /^\#=GS\s+(\S+)\s+AC\s+(\S+)/ ) {
	    $accession{$1} = $2;
	} elsif( $entry =~ /^\#=GS\s+(\S+)\s+DE\s+(.+)\s*$/ ) {
	    $desc{$1} .= $2;
	} elsif( $entry =~ /^([A-Za-z.\-\*]+)$/ ) {
	    $align{$name} .= $1;
	} elsif( $entry =~ /^([^\#]\S+)\s+([A-Za-z.\-\*]+)\s*/ ) {
	    ($name,$seq) = ($1,$2);
	    if( ! defined $align{$name}  ) {
		push @c2name, $name;
	    }
	    $align{$name} .= $seq;
	}
    }

    # ok... now we can make the sequences
    for my $name ( @c2name ) {
	if( $name =~ m{(\S+)/(\d+)-(\d+)} ) {
	    $seqname = $1;
	    $start = $2;
	    $end = $3;
	} else {
	    $seqname=$name;
	    $start = 1;
	    $end = length($align{$name});
	}
	$seq = new Bio::LocatableSeq
	    ('-seq'              => $align{$name},
	     '-display_id'       => $seqname,
	     '-start'            => $start,
	     '-end'              => $end,
	     '-description'      => $desc{$name},
	     '-accession_number' => $accession{$name},
	     );

	$aln->add_seq($seq);
    }

#  If $end <= 0, we have either reached the end of
#  file in <fh> or we have encountered some other error
#
    return if ($end <= 0);
    return $aln;
}


=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in stockholm format  ###Not yet implemented!###
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object


=cut

sub write_aln {
    my ($self,@aln) = @_;
    $self->throw_not_implemented();
}

1;
