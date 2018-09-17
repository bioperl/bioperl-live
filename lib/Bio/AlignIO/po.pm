# $Id: po.pm
#
# BioPerl module for Bio::AlignIO::po

#   based on the Bio::AlignIO::fasta module
#       by Peter Schattner (and others?)
#
#       and the SimpleAlign.pm module of Ewan Birney
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::po - po MSA Sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::AlignIO> class.

=head1 DESCRIPTION

This object can transform L<Bio::SimpleAlign> objects to and from
'po' format flat file databases. 'po' format is the native format of
the POA alignment program (Lee C, Grasso C, Sharlow MF, 'Multiple
sequence alignment using partial order graphs', Bioinformatics (2002),
18(3):452-64).

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

=head1 AUTHORS - Matthew Betts

Email: matthew.betts@ii.uib.no

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::po;
use strict;

use Bio::SimpleAlign;

use base qw(Bio::AlignIO);


=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
 Returns : L<Bio::Align::AlignI> object - returns undef on end of file
	    or on error
 Args    : NONE

=cut

sub next_aln {
    my $self = shift;

    my $aln;
    my $entry;
    my $name;
    my $seqs;
    my $seq;
    my $nodes;
    my $list;
    my $node;
    my @chars;
    my $s;
    my $a;

    $aln =  Bio::SimpleAlign->new();

    # get to the first 'VERSION' line
    while(defined($entry = $self->_readline)) {
	if($entry =~ /^VERSION=(\S+)/) {
	    $aln->source($1);

	    if(defined($entry = $self->_readline) and $entry =~ /^NAME=(\S+)/) {
		$aln->id($1);
	    }

	    last;
	}
    }

    # read in the sequence names and node data, up to the end of
    # the file or the next 'VERSION' line, whichever comes first
    $seqs = [];
    $nodes = [];
    while(defined($entry = $self->_readline)) {
	if($entry =~ /^VERSION/) {
	    # start of a new alignment, so...
	    $self->_pushback($entry);
	    last;
	}
	elsif($entry =~ /^SOURCENAME=(\S+)/) {
	    $name = $1;

	    if($name =~ /(\S+)\/(\d+)-(\d+)/) {
		$seq = Bio::LocatableSeq->new(
					      '-display_id' => $1,
					      '-start'      => $2,
					      '-end'        => $3,
					      '-alphabet'   => $self->alphabet,
					    );

	    } else {
		$seq = Bio::LocatableSeq->new('-display_id'=> $name,
					      '-alphabet'  => $self->alphabet);
	    }

	    # store sequences in a list initially, because can't guarantee
	    # that will get them back from SimpleAlign object in the order
	    # they were read, and also can't add them to the SimpleAlign
	    # object here because their sequences are currently unknown
	    push @{$seqs}, {
			    'seq' => $seq,
			    'str' => '',
			};
	}
	elsif($entry =~ /^SOURCEINFO=(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.*)/) {
	    $seq->desc($5);
	}
	elsif($entry =~ /^(\S):(\S+)/) {
	    $node = {
		     'aa'     => $1,
		     'L'      => [],
		     'S'      => [],
		     'A'      => [],
		     'status' => 'unvisited',
		    };

	    $list = $2;
	    if($list =~ /^([L\d]*)([S\d]*)([A\d]*)/) {
		push(@{$node->{'L'}}, split(/L/, $1));
		push(@{$node->{'S'}}, split(/S/, $2));
		push(@{$node->{'A'}}, split(/A/, $3));

		(@{$node->{'L'}} > 0) and shift @{$node->{'L'}};
		(@{$node->{'S'}} > 0) and shift @{$node->{'S'}};
		(@{$node->{'A'}} > 0) and shift @{$node->{'A'}};
	    }

	    push @{$nodes}, $node;
	}
    }

    # process the nodes
    foreach $node (@{$nodes}) {
	($node->{'status'} ne 'unvisited') and next;

	@chars = ($aln->gap_char) x @{$seqs}; # char for each seq defaults to a gap

	# set the character for each sequence represented by this node
	foreach $s (@{$node->{'S'}}) {
	    $chars[$s] = $node->{'aa'};
	}
	$node->{'status'} = 'visited';

	# do the same for each node in the same align ring
	while(defined($a = $node->{'A'}->[0])) {
	    $node = $nodes->[$a];
	    ($node->{'status'} ne 'unvisited') and last;

	    foreach $s (@{$node->{'S'}}) {
		$chars[$s] = $node->{'aa'};
	    }

	    $node->{'status'} = 'visited';
	}

	# update the sequences
	foreach $seq (@{$seqs}) {
	    $seq->{'str'} .= shift @chars;
	}
    }

    # set the sequences of the bioperl objects
    # and add them to the alignment
    foreach $seq (@{$seqs}) {
	$seq->{'seq'}->seq($seq->{'str'});
	$aln->add_seq($seq->{'seq'});
    }

    # has an alignment been read?...

    return $aln if $aln->num_sequences;
	return;
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in po format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object

=cut

sub write_aln {
    my $self = shift;
    my @alns = @_;

    my $aln;
    my $seqs;
    my $nodes;
    my $seq;
    my $node;
    my $col;
    my $ring;
    my $i;
    my $char;

    foreach $aln (@alns) {
	if(!$aln or !$aln->isa('Bio::Align::AlignI')) {
	    $self->warn("Must provide a Bio::Align::AlignI object when calling write_aln");
	    next;
	}

	# store the seqs on a list, because po format
	# refers to them by position on this list
	$seqs  = [];
	foreach $seq ($aln->each_seq()) {
	    push @{$seqs}, {
			    'seq'      => $seq,
			    'n_nodes'  => 0,
			    'first'    => undef,
			    'previous' => undef,
			   };
	}

	# go through each column in the alignment and construct
	# the nodes for the equivalent poa alignment ring
	$nodes = [];
	for($col = 0; $col < $aln->length; $col++) {
	    $ring = {
		     'nodes' => {},
		     'first' => scalar @{$nodes},
		     'last'  => scalar @{$nodes},
		    };

	    for($i = 0; $i < @{$seqs}; $i++) {
		$seq = $seqs->[$i];

		$char = $seq->{'seq'}->subseq($col + 1, $col + 1);

		($char eq $aln->gap_char) and next;

		if(!defined($node = $ring->{'nodes'}->{$char})) {
		    $node = {
			     'n'  => scalar @{$nodes},
			     'aa' => $char,
			     'L'  => {},
			     'S'  => [],
			     'A'  => [],
			    };

		    # update the ring
		    $ring->{'nodes'}->{$char} = $node;
		    $ring->{'last'} = $node->{'n'};

		    # add the node to the node list
		    push @{$nodes}, $node;
		}

		# add the sequence to the node
		push @{$node->{'S'}}, $i;

		# add the node to the sequence
		defined($seq->{'first'}) or ($seq->{'first'} = $node);
		$seq->{'n_nodes'}++;

		# add an edge from the previous node in the sequence to this one.
		# Then set the previous node to the current one, ready for the next
		# residue in this sequence
		defined($seq->{'previous'}) and ($node->{'L'}->{$seq->{'previous'}->{'n'}} = $seq->{'previous'});
		$seq->{'previous'} = $node;
	    }

	    # set the 'next node in ring' field for each node in the ring
	    if($ring->{'first'} < $ring->{'last'}) {
		for($i = $ring->{'first'}; $i < $ring->{'last'}; $i++) {
		    push @{$nodes->[$i]->{'A'}}, $i + 1;
		}
		push @{$nodes->[$ring->{'last'}]->{'A'}}, $ring->{'first'};
	    }
	}

	# print header information
	$self->_print(
		      'VERSION=', ($aln->source and ($aln->source !~ /\A\s*\Z/)) ? $aln->source : 'bioperl', "\n",
		      'NAME=', $aln->id, "\n",
		      'TITLE=', ($seqs->[0]->{'seq'}->description or $aln->id), "\n",
		      'LENGTH=', scalar @{$nodes}, "\n",
		      'SOURCECOUNT=', scalar @{$seqs}, "\n",
		     );

	# print sequence information
	foreach $seq (@{$seqs}) {
	    $self->_print(
			  'SOURCENAME=', $seq->{'seq'}->display_id, "\n",
			  'SOURCEINFO=',
			  $seq->{'n_nodes'},      ' ', # number of nodes in the sequence
			  $seq->{'first'}->{'n'}, ' ', # index of first node containing the sequence
			  0,                      ' ', # FIXME - sequence weight?
			  -1,                     ' ', # FIXME - index of bundle containing sequence?
			  ($seq->{'seq'}->description or 'untitled'), "\n",
			 );
	}

	# print node information
	foreach $node (@{$nodes}) {
	    $self->_print($node->{'aa'}, ':');
	    (keys %{$node->{'L'}} > 0) and $self->_print('L', join('L', sort {$a <=> $b} keys %{$node->{'L'}}));
	    (@{$node->{'S'}} > 0) and $self->_print('S', join('S', @{$node->{'S'}}));
	    (@{$node->{'A'}} > 0) and $self->_print('A', join('A', @{$node->{'A'}}));
	    $self->_print("\n");
	}
    }
    $self->flush if $self->_flush_on_write && defined $self->_fh;

    return 1;
}

1;
