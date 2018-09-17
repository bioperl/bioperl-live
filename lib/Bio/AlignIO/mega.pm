#
# BioPerl module for Bio::AlignIO::mega
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::mega - Parse and Create MEGA format data files

=head1 SYNOPSIS

    use Bio::AlignIO;
    my $alignio = Bio::AlignIO->new(-format => 'mega',
                                   -file   => 't/data/hemoglobinA.meg');

    while( my $aln = $alignio->next_aln ) {
    # process each alignment or convert to another format like NEXUS
    }

=head1 DESCRIPTION

This object handles reading and writing data streams in the MEGA
format (Kumar and Nei).


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::AlignIO::mega;
use vars qw($MEGANAMELEN %VALID_TYPES $LINELEN $BLOCKLEN);
use strict;

use Bio::SimpleAlign;
use Bio::LocatableSeq;

# symbols are changed due to MEGA's use of '.' for redundant sequences

BEGIN {
  $MEGANAMELEN = 10;
  $LINELEN = 60;
  $BLOCKLEN = 10;
  %VALID_TYPES =  map {$_, 1} qw( dna rna protein standard);
}
use base qw(Bio::AlignIO);


=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
           Supports the following MEGA format features:
           - The file has to start with '#mega'
           - Reads in the name of the alignment from a comment
             (anything after '!TITLE: ') .
           - Reads in the format parameters datatype

 Returns : L<Bio::Align::AlignI> object - returns 0 on end of file
	    or on error
 Args    : NONE


=cut

sub next_aln{
   my ($self) = @_;
   my $entry;
   my ($alphabet,%seqs);
   local $Bio::LocatableSeq::OTHER_SYMBOLS = '\*\?\.';
   local $Bio::LocatableSeq::GAP_SYMBOLS = '\-';
   my $aln = Bio::SimpleAlign->new(-source => 'mega');

   while( defined($entry = $self->_readline()) && ($entry =~ /^\s+$/) ) {}

   $self->throw("Not a valid MEGA file! [#mega] not starting the file!")
       unless $entry =~ /^#mega/i;

   while( defined($entry = $self->_readline() ) ) {
       local($_) = $entry;
       if(/\!Title:\s*([^\;]+)\s*/i) { $aln->id($1)}
       elsif( s/\!Format\s+([^\;]+)\s*/$1/ ) {
	   my (@fields) = split(/\s+/,$1);
	   foreach my $f ( @fields ) {
	       my ($name,$value) = split(/\=/,$f);
	       if( $name eq 'datatype' ) {
		   $alphabet = $value;
	       } elsif( $name eq 'identical' ) {
		   $aln->match_char($value);
	       } elsif( $name eq 'indel' ) {
		   $aln->gap_char($value);
	       }
	   }
       } elsif( /^\#/ ) {
	   last;
       }
   }
   my @order;
   while( defined($entry) ) {
       if( $entry !~ /^\s+$/ ) {
	   # this is to skip the leading '#'
	   my $seqname = substr($entry,1,$MEGANAMELEN-1);
	   $seqname =~ s/(\S+)\s+$/$1/g;
	   my $line = substr($entry,$MEGANAMELEN);
	   $line =~ s/\s+//g;
	   if( ! defined $seqs{$seqname} ) {push @order, $seqname; }
	   $seqs{$seqname} .= $line;
       }
       $entry = $self->_readline();
   }

   foreach my $seqname ( @order ) {
       my $s = $seqs{$seqname};
       $s =~ s/[$Bio::LocatableSeq::GAP_SYMBOLS]+//g;
       my $end = length($s);
       my $seq = Bio::LocatableSeq->new('-alphabet'   => $alphabet,
					'-display_id' => $seqname,
					'-seq'        => $seqs{$seqname},
					'-start'      => 1,
					'-end'        => $end);

       $aln->add_seq($seq);
   }
   $aln->unmatch;
   return $aln if $aln->num_sequences;
   return;
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in MEGA format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object

=cut

sub write_aln{
   my ($self,@aln) = @_;
   my $count = 0;
   my $wrapped = 0;
   my $maxname;

   foreach my $aln ( @aln ) {
       if( ! $aln || ! $aln->isa('Bio::Align::AlignI')  ) {
	   $self->warn("Must provide a Bio::Align::AlignI object when calling write_aln");
	   return 0;
       } elsif( ! $aln->is_flush($self->verbose) ) {
	   $self->warn("All Sequences in the alignment must be the same length");
	   return 0;
       }
       $aln->match();
       my $len = $aln->length();
       my $format = sprintf('datatype=%s identical=%s indel=%s;',
			    $aln->get_seq_by_pos(1)->alphabet(),
			    $aln->match_char, $aln->gap_char);

       $self->_print(sprintf("#mega\n!Title: %s;\n!Format %s\n\n\n",
			     $aln->id, $format));

       my ($count, $blockcount,$length) = ( 0,0,$aln->length());
       $aln->set_displayname_flat();
       while( $count < $length ) {
	   foreach my $seq ( $aln->each_seq ) {
	       my $seqchars = $seq->seq();
	       $blockcount = 0;
	       my $substring = substr($seqchars, $count, $LINELEN);
	       my @blocks;
	       while( $blockcount < length($substring) ) {
		   push @blocks, substr($substring, $blockcount,$BLOCKLEN);
		   $blockcount += $BLOCKLEN;
	       }
	       $self->_print(sprintf("#%-".($MEGANAMELEN-1)."s%s\n",
				     substr($aln->displayname($seq->get_nse()),
					    0,$MEGANAMELEN-2),
				     join(' ', @blocks)));
	   }
	   $self->_print("\n");
	   $count += $LINELEN;
       }
   }
   $self->flush if $self->_flush_on_write && defined $self->_fh;
   return 1;
}


1;
