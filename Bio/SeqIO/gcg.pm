#
# BioPerl module for Bio::SeqIO::fasta
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#          and Lincoln Stein <lstein@cshl.org>
#
# Copyright Ewan Birney & Lincoln Stein
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

=head1 AUTHORS - Ewan Birney & Lincoln Stein

Email: birney@sanger.ac.uk
       lstein@cshl.org


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::gcg;
use vars '@ISA';
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
   my($id,$type,$desc,$line,$chksum,$sequence);


   while( defined($_ = $self->_readline()) ) {

       ## Get the descriptive info (anything before the line with '..')
       unless( /\.\.$/ ) { $desc.= $_; }

       ## Pull ID, Checksum & Type from the line containing '..'
       /\.\.$/ && do     { $line = $_; chomp; 
                           if(/Check\:\s(\d+)\s/) { $chksum = $1; }
                           if(/Type:\s(\w)\s/)    { $type   = $1; }
                           if(/(.*)\s+Length/)    { $id     = $1; }
                           last; 
                         }
   }


   while( defined($_ = $self->_readline()) ) {

       ## This is where we grab the sequence info.

       if( /\.\.$/ ) { 
        $self->throw("Looks like start of another sequence. See documentation. "); 
       }

       next if($_ eq "\n");       ## skip whitespace lines in formatted seq
       s/[^a-zA-Z]//g;            ## remove anything that is not alphabet char
       $_ = uc($_);               ## uppercase sequence
       $sequence .= $_;
   }


   ##If we parsed out a checksum, we might as well test it

   if(defined $chksum) { 
       unless(_validate_checksum($sequence,$chksum)) {
	   $self->throw("Checksum failure on parsed sequence.");
       }
   }

   ## Remove whitespace from identifier because the constructor
   ## will throw a warning otherwise...
   if(defined $id) { $id =~ s/\s+//g;}

   ## Turn our parsed "Type: N" or "Type: P" (if found) into the appropriate
   ## keyword that the constructor expects...
   if(defined $type) {
       if($type eq "N") { $type = "dna";      }
       if($type eq "P") { $type = "amino";    }
   }


   return Bio::Seq->new(-seq  => $sequence, 
                        -id   => $id, 
                        -desc => $desc, 
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
   for my $seq (@seq) {
     my $str         = $seq->seq;
     my $comment     = $seq->desc; 
     my $id          = $seq->id;
     my $type        = $seq->moltype();
     my($timestamp)  = localtime;
     my($sum,$offset,$len,$i,$j,$cnt,@out);

     $len = length($str);
     $offset=1;
     ## Set the offset if we have any non-standard numbering going on
     $sum=0;
     
     #Generate the GCG Checksum value
     for($i=0; $i<$len ;$i++) {             
       $cnt++;
       $sum += $cnt * ord(substr($str,$i,1));
       ($cnt == 57) && ($cnt=0);
     }
     $sum %= 10000;

     #Output the sequence header info
     push(@out,"$comment\n");                        
     push(@out," $id Length: $len  $timestamp  $type Check: $sum  ..\n\n");
     
     #Format the sequence
     $i = $#out + 1;
     for($j = 0 ; $j < $len ; ) {
       if( $j % 50 == 0) {
	 $out[$i] = sprintf("%8d  ",($j+$offset)); #numbering 
       }
       $out[$i] .= sprintf("%s",substr($str,$j,10));
       $j += 10;
       if( $j < $len && $j % 50 != 0 ) {
	 $out[$i] .= " ";
       }elsif($j % 50 == 0 ) {
	 $out[$i++] .= "\n";
       }                           
     }
     local($^W) = 0;
     if($j % 50 != 0 ) {
       $out[$i] .= "\n";
     }
     $out[$i] .= "\n";


     return unless $self->_print(@out);
   }
   return 1;
}

=head2 _validate_checksum

 Title   : _validate_checksum
 Usage   : n/a - internal method
 Function: if parsed gcg sequence contains a checksum field
         : we compare it to a value computed here on the parsed
         : sequence. A checksum mismatch would indicate some
         : type of parsing failure occured.
         :
 Returns : 1 for success, 0 for failure
 Args    : string containing parsed seq, value of parsed cheksum


=cut

sub _validate_checksum {
    my($seq,$parsed_sum) = @_;
    my($i,$len,$computed_sum,$cnt);

    $len = length($seq);

    #Generate the GCG Checksum value

    for($i=0; $i<$len ;$i++) {             
	$cnt++;
	$computed_sum += $cnt * ord(substr($seq,$i,1));
	($cnt == 57) && ($cnt=0);
    }
    $computed_sum %= 10000;

    ## Compare and decide if success or failure

    if($parsed_sum = $computed_sum) {
	return 1;
    } else { return 0; }


}

1;
