# $Id$
#
# BioPerl module for Bio::AlignIO::stockholm

#	based on the Bio::SeqIO::stockholm module
#       by Ewan Birney <birney@sanger.ac.uk>
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

Do not use this module directly.  Use it via the Bio::AlignIO class.

=head1 DESCRIPTION

This object can transform Bio::SimpleAlign objects to and from stockholm flat
file databases.

=head1 FEEDBACK

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS - Peter Schattner

Email: schattner@alum.mit.edu


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::stockholm;
use vars qw(@ISA);
use strict;

use Bio::AlignIO;

@ISA = qw(Bio::AlignIO);

# AlignIO is special new must be explict 

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;
    $self->_initialize(@args);
    return $self;
}
# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  return unless my $make = $self->SUPER::_initialize(@args);
}

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
 Returns : SimpleAlign object
 Args    : NONE

=cut

sub next_aln {
    my $self = shift;
    my $entry;

    my ($start,$end,%align,$name,$seqname,$seq,$count,
	%hash,%c2name, %accession, $no);

    # in stockholm format, every non-blank line that does not start
    # with '#=' is an alignment segment; the '#=' lines are mark up lines.
    # Of particular interest are the '#=GF <name/st-ed> AC <accession>'
    # lines, which give accession numbers for each segment



    my $aln =  Bio::SimpleAlign->new();


    while( $entry = $self->_readline) {
        $entry !~ /\w+/ && next;

	if ($entry =~ /^\#\s*STOCKHOLM\s+/) {
	    last;
	}
	else {
	    $self->throw("Not Stockholm format: Expecting \"# STOCKHOLM 1.0\"; Found \"$_\"");
	}
    }
#
#    Next section is same as for selex format
#
    while( $entry = $self->_readline) {
        $entry =~ /^\#=GS\s+(\S+)\s+AC\s+(\S+)/ && do {
	    				$accession{ $1 } = $2;
	    				next;
					};
	$entry !~ /^([^\#]\S+)\s+([A-Za-z\.\-]+)\s*/ && next;
	
	$name = $1;
	$seq = $2;

	if( ! defined $align{$name}  ) {
	    $count++;
	    $c2name{$count} = $name;
	}
	$align{$name} .= $seq;
    }

    # ok... now we can make the sequences

     foreach $no ( sort { $a <=> $b } keys %c2name ) {
	$name = $c2name{$no};

	if( $name =~ /(\S+)\/(\d+)-(\d+)/ ) {
	    $seqname = $1;
	    $start = $2;
	    $end = $3;
	} else {
	    $seqname=$name;
	    $start = 1;
	    $end = length($align{$name});
	}
	$seq = new Bio::LocatableSeq('-seq'=>$align{$name},
			    '-id'=>$seqname,
			    '-start'=>$start,
			    '-end'=>$end,
			    '-type'=>'aligned',
				     '-accession_number' => $accession{$name},

			    );

       $aln->addSeq($seq);

   }

#  If $end <= 0, we have either reached the end of
#  file in <fh> or we have encountered some other error
#
   if ($end <= 0) { undef $aln;}

   return $aln;
}


=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in stockholm format  ###Not yet implemented!###
 Returns : 1 for success and 0 for error
 Args    : Bio::SimpleAlign object


=cut

sub write_aln {
    my ($self,@aln) = @_;

    $self->throw("Sorry: stockholm-format output, not yet implemented! /n");
}

1;
