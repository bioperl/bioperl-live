#
# BioPerl module for Bio::AlignIO::maf
#
# Copyright Allen Day
#

=head1 NAME

Bio::AlignIO::maf - Multiple Alignment Format sequence input stream

=head1 SYNOPSIS

 Do not use this module directly.  Use it via the Bio::AlignIO class.

 use Bio::AlignIO;

 my $alignio = Bio::AlignIO->new(-fh => \*STDIN, -format => 'maf');

 while(my $aln = $alignio->next_aln()){
   my $match_line = $aln->match_line;

   print $aln, "\n";

   print $aln->length, "\n";
   print $aln->num_residues, "\n";
   print $aln->is_flush, "\n";
   print $aln->num_sequences, "\n";

   $aln->splice_by_seq_pos(1);

   print $aln->consensus_string(60), "\n";
   print $aln->get_seq_by_pos(1)->seq, "\n";
   print $aln->match_line(), "\n";

   print "\n";
 }

=head1 DESCRIPTION

This class constructs Bio::SimpleAlign objects from an MAF-format
multiple alignment file.

Writing in MAF format is currently unimplemented.

Spec of MAF format is here:
  http://genome.ucsc.edu/FAQ/FAQformat

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

=head1 AUTHORS - Allen Day

Email: allenday@ucla.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::maf;
use strict;

use Bio::SimpleAlign;

use base qw(Bio::AlignIO);

=head2 new

 Title   : new
 Usage   : my $alignio = Bio::AlignIO->new(-format => 'maf'
					  -file   => '>file',
					  -idlength => 10,
					  -idlinebreak => 1);
 Function: Initialize a new L<Bio::AlignIO::maf> reader
 Returns : L<Bio::AlignIO> object
 Args    :

=cut

sub _initialize {
  my($self,@args) = @_;
  $self->SUPER::_initialize(@args);

  1;
}

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
           Throws an exception if trying to read in PHYLIP
           sequential format.
 Returns : L<Bio::SimpleAlign> object
 Args    : 

=cut

sub next_aln {
    my $self = shift;

	# check beginning of file for proper header
    if(!$self->{seen_header}){
	  my $line = $self->_readline;
	  $self->throw("This doesn't look like a MAF file.  First line should start with ##maf, but it was: ".$line)
		  unless $line =~ /^##maf/;
	  $self->{seen_header} = 1;
	  # keep in case we parse this later
	  $self->_pushback($line);
    }
	
    my $aln =  Bio::SimpleAlign->new(-source => 'maf');

    my($aline, @slines, $seen_aline);
    while(my $line = $self->_readline()){
	  if ($line =~ /^a\s/xms) {
		# next block?
		if ($seen_aline) {
		  $self->_pushback($line);
		  last;
		}
		$aline = $line;
		$seen_aline++;
	  } elsif ($line =~ /^s\s/xms) {
		push @slines, $line;
	  } else {
		# missed lines
		$self->debug($line);
	  }
    }
	
	# all MAF starts with 'a' line
    return unless $aline;

    my($kvs) = $aline =~ /^a\s+(.+)$/;
    my @kvs  = split /\s+/, $kvs if $kvs;
    my %kv;
    foreach my $kv (@kvs){
	my($k,$v) = $kv =~ /(.+)=(.+)/;
	$kv{$k} = $v;
    }

    $aln->score($kv{score});

    foreach my $sline (@slines){
	my($s,$src,$start,$size,$strand,$srcsize,$text) =
	    split /\s+/, $sline;
	# adjust coordinates to be one-based inclusive
        $start = $start + 1;
    $strand = $strand eq '+' ? 1 : $strand eq '-' ? -1 : 0;
	my $seq = Bio::LocatableSeq->new('-seq'          => $text,
					 '-display_id'   => $src,
					 '-start'        => $start,
					 '-end'          => $start + $size - 1,
					 '-strand'       => $strand,
					 '-alphabet'     => $self->alphabet,
					);
	$aln->add_seq($seq);
    }

    return $aln if $aln->num_sequences;
	return;
}

sub write_aln {
  shift->throw_not_implemented
}

1;
