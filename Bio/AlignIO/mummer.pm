#
# BioPerl module for Bio::AlignIO::mummer
#
# You may distribute this module under the same terms as perl itself.
#
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::mummer - MUMmer multiple whole sequence alignment input/output.

=head1 SYNOPSIS

Do not use this module directly. Use it via the L<Bio::AlignIO|Bio::AlignIO>
class.

=head1 DESCRIPTION

This object can read and write MUMmer command output. MUMmer can quickly
compute multiple unique (exact) matches via suffix trees. Only the mummer
tool is supported by this module due to the complexity and the many different
output formats of the various MUMmer tools.

Due to the various output formats of the mummer tool, it may not be possible
to always convert from one Bio::AlignIO::mummer stream to another, let alone
from one of the other AlignIO objects. The output of the mummer tool is
variable and if one run produces limited output, it is not possible to convert
to an output with verbose output.

For more information:

 http://mummer.sourceforge.net

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an intergral part of the evolution of this and other Bioperl
modules. Send your comments and suggestions preferable to one of the Bioperl
mailing lists. Your participation is much appreciated.

 bioperl-l@bioperl.org                  - General discussion
 http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the
bugs and their resolution. Bug reports can be submitted via the web:

 http://bugzilla.open-bio.org

=head1 AUTHORS

 Jason Switzer  - jswitzer@gmail.com
 Joshua Wu      - nike284@gmail.com
 Aimee Seuffer  - aimeeseuf@msn.com

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a '_'.

=cut

package Bio::AlignIO::mummer;

use base qw(Bio::AlignIO);
use Bio::SimpleAlign;

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream
 Returns : Bio::Align::AlignI object
 Args    : NONE

 See L<Bio::Align::AlignI> for details

=cut

sub next_aln {
   my $self = shift;
   my ($startg, $startl, $ref_seqname, $qry_seqname);
   my ($seq, @lines, $tempdesc, $querylen, $reversed);
   my $aln = Bio::SimpleAlign->new();

   my @lines = <$mumfh>;
   chomp @lines;
   #only read one alignment, denoted by >ID
   while(defined($line = $self->_readline)) {
      chomp $line;
      if($line =~ /^\s*\>\s*(.*)/) {
         #parse the header
         split /\s+/, $1;
         if($seqname = shift) {
         }
         if($tempdesc = shift) {
            if($tempdesc =~ /reverse/i) {
               $reversed = -1;
               $tempdesc = shift;
            }
            if($tempdesc =~ /len/i) {
               shift;
               if($querylen = shift) {
               }
            }
         }

         #FIXME kludge to get the query string info in the alignment ID
         my $alnname = $seqname;
         $alnname .= " Reverse" if $reversed;
         $alnname .= " Len = " . $querylen if $querylen;
         $aln->id($alnname);

         my $count = 0;
         while(defined($line = $self->_readline)) {
            #lookahead, if it's another alignment, put it back
            if(!$line || $line =~ /^\s*\>\s*(.*)/) {
               $self->_pushback($line);
               return $aln;
            }

            #clean up the line, remove redundancy
            chomp $line;
            $line =~ s/(\s)+/$1/g;
            $line =~ s/^\s+(.*)/$1/g;
            my @tokens = split / /, $line;
            if(3 == @tokens) {
               #three column output
               $startg = shift @tokens;
               $startl = shift @tokens;
               $matchlen = shift @tokens;
            } elsif(4 == @tokens) {
               #four column output
               $ref_seqname = shift @tokens;
               $startg = shift @tokens;
               $startl = shift @tokens;
               $matchlen = shift @tokens;
            }

            #lookahead, if it's another alignment, put it back
            $line = $self->_readline;
            $seq = "";
            if(!$line || $line =~ /^\s*>\s*(.*)/) {
               $self->_pushback($line);
            } else {
               #we may or may not have a sequence string
               my @tokens = split / /, $line;
               $line =~ s/^\s+(.*)/$1/g;
               $line =~ s/^\s+(.*)/$1/g;
               if(1 == @tokens) {
                  chomp $line;
                  $seq = $line;
               } else {
                  $self->_pushback($line);
               }
            }

            #FIXME kludge for SimpleAlign using a crappy hash key name
            $count++;
            my $id = $seqname . ($reversed ? "#$count" : '');

            #completed locatableseq, add to the alignment
            my $t = Bio::Seq::RefLocatableSeq->new('-seq' => $seq,
                                                   '-id' => $id,
                                                   '-strand' => $reversed,
                                                   '-start' => $startl,
                                                   '-end' => $startl + $matchlen,
                                                   '-ref_id' => $ref_seqname,
                                                   '-ref_start' => $startg,
                                                   '-ref_end' => $startg + $querylen);
            #print "READ: " . Dumper($t) . "\n";
            $aln->add_seq($t);

         }
         return $aln;
      }
   }
   return;
}

=head2 _pop_stackframe

 Title   : _pop_stackframe
 Usage   : print "WhoCalledMe: " . &_pop_stackframe . "\n";
 Function: Prints the calling function's name for debugging.
 Returns : Stringified name of the calling function.
 Args    : NONE

=cut

sub _pop_stackframe {
   (caller(2))[3];
}

=head2 _check_subtype

 Title   : _check_subtype
 Usage   : $self->_check_subtype($aln, "Bio::Align::AlignI", 1)
 Function: Checks if an object is of the proper subtype.
 Returns : 1 if UNIVERSAL::isa match, else 0.
 Args    : self - reference to self
           seq - object to check the isa relationship
           subclass - isa relationship to check (string version)
           warn - whether to issue a warning or not

=cut

sub _check_subtype {
   my ($self, $seq, $subclass, $warn) = @_;
   if(!$seq || !$seq->isa($subclass)) {
      $self->warn("Must provide a $subclass object to " . &_pop_stackframe);
      return 0;
   }
   return 1;
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the mummer-format object (.aln) into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Align::AlignI object

=cut

sub write_aln {
   my ($self, @alignments) = @_;
   my $subref = "Bio::Align::AlignI";
   my ($seq, $header, $iseq, $name);

   for my $aln(@alignments) {
      next if(!$self->_check_subtype($aln, $subref, 1));

      $aln->set_displayname_flat($self->force_displayname_flat);

      my @eachseq = $aln->each_seq();
      #pick off the header information without touching anything
      #each seq has the same header information
      $header = "> " . $aln->id;
      $self->_print("$header\n");
      next if(!@eachseq);

      for my $iseq(@eachseq) {
         next if(!$self->_check_subtype($iseq, "Bio::Seq::RefLocatableSeq", 1));
         my $nextline = "";
         #print "WRITE: " . Dumper($iseq) . "\n";
         $nextline .= $iseq->ref_id . "    " if $iseq->ref_id;
         $nextline .= $iseq->ref_start . "    ";
         $nextline .= $iseq->start . "    ";
         $nextline .= ($iseq->end - $iseq->start) . "\n";
         $self->_print($nextline);
         $nextline = $iseq->seq . "\n" if $iseq->seq;
         $self->_print($nextline) if $iseq->seq;
      }
   }
   $self->flush if $self->_flush_on_write && $self->_fh;
   return 1;
}

1;

