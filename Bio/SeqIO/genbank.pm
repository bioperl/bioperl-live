#
# BioPerl module for Bio::SeqIO::fasta
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#          and Lincoln Stein <lstein@cshl.org>
#
# Copyright Ewan Birney, Lincoln Stein, Chris Dagdigian and Kai Kumpf
#
# You may distribute this module under the same terms as perl itself
#
# _history
# October 18, 1999  Largely rewritten by Lincoln Stein

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::gcg - GCG sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::SeqIO class.

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from GCG flat
file databases.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS - Ewan Birney, Lincoln Stein, Chris Dagdigian and Kai Kumpf

Email: birney@sanger.ac.uk
       lstein@cshl.org
       dag@sonsorol.org
       kku@imb-jena.de

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::SeqIO::genbank;
use vars qw(@ISA);
use strict;
use Bio::SeqIO;

@ISA = qw(Bio::SeqIO);

# new() is inherited from Bio::Root::Object

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    :


=cut

sub next_seq{
   my ($self,@args)    = @_;

   local $/ = "//\n";  # read a record at a time
   return unless $_ =  $self->_readline;
   my ($type)       =  /^LOCUS\s+\w+\s+\d+\s+bp\s+(\w+)/m;
   my ($definition) =  /^DEFINITION\s+((.+)(\n\s+.+)*)/m;  $definition =~ s/\n\s+/ /g;
   my ($id)         =  /^ACCESSION\s+(\S+)/m;
   my ($sequence)   =  /^ORIGIN\s+\n(.+)\Z/sm;             $sequence =~ s/[^a-zA-Z]//g;

   if( length($sequence) == 0 ) { 
    $self->throw("Looks like sequence is empty."); 
   }

   ## Turn our parsed type into the appropriate
   ## keyword that the constructor expects...
   if(defined $type) {
       if($type eq "DNA") { $type = "dna";      }
       if($type =~ /.?RNA/) { $type = "rna";      }
   }

   return Bio::Seq->new(-seq  => $sequence, 
                        -id   => $id, 
                        -desc => $definition, 
                        -type => $type );
}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq(@seq)
 Function: writes the fromatted $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object, sequence string

=cut

sub write_seq {
   my ($self,@seq) = @_;

   require POSIX;
   POSIX->import('strftime');  # for damn date format

   for my $seq (@seq) {

     my $str         = lc $seq->seq;
     my $comment     = $seq->desc; 
     my $id          = $seq->id;
     my $type        = $seq->type();
     my $timestamp   = strftime('%d-%b-%Y',localtime);
     my($sum,$offset,$len,$i,$j,$cnt,@out);

     $len = length($str);
     $offset=0;
     ## Set the offset if we have any non-standard numbering going on
     if($seq->start < 0)   { $offset = ( 0 + $seq->start); }
     if($seq->start >= 1)  { $offset = $seq->start;        }  
     if($seq->start == 0)  { $offset = -1;                 }

     #Output the sequence header info
     my $definition = $seq->desc;
     # add as many linefeeds as necessary to flesh out definition line
     $definition =~ s/(.{50,60})\s+/$1\n            /g;
     push(@out,sprintf("%-12s%-16s%-9s%-27s%-11s\n",'LOCUS',$id,"$len bp",'DNA',"\U$timestamp"));
     push(@out,"DEFINITION  $definition\n");
     push(@out,"ORIGIN\n");

     #Format the sequence
     $i = $#out + 1;
     for($j = 0 ; $j < $len ; ) {
       if( $j % 60 == 0) {
	 $out[$i] = sprintf("%8d  ",($j+$offset)); #numbering 
       }
       $out[$i] .= sprintf("%s",substr($str,$j,10));
       $j += 10;
       if( $j < $len && $j % 60 != 0 ) {
	 $out[$i] .= " ";
       }elsif($j % 60 == 0 ) {
	 $out[$i++] .= "\n";
       }                           
     }
     local($^W) = 0;
     if($j % 60 != 0 ) {
       $out[$i] .= "\n";
     }
     $out[$i] .= "//\n";
     return unless $self->_print (@out);
   }
   return 1;
}
