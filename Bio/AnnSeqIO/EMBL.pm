

#
# BioPerl module for Bio::SeqIO::EMBL
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::AnnSeqIO::EMBL - EMBL sequence input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the AnnSeqIO handler system. Go:

    $stream = Bio::AnnSeqIO->new(-file => $filename, -format => 'EMBL');

    while my $annseq ( $stream->next_annseq() ) {
	# do something with $annseq
    }

=head1 DESCRIPTION

This object can transform Bio::AnnSeq objects to and from EMBL flat
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

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::AnnSeqIO::EMBL;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::AnnSeq;
use Bio::AnnSeqIO::FTHelper;
use Bio::SeqFeature::Generic;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use FileHandle;

@ISA = qw(Bio::Root::Object Exporter);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;
   
  my ($file,$fh) = $self->_rearrange([qw(
					 FILE
					 FH
					 )],
				     @args,
				     );
  if( $file && $fh ) {
      $self->throw("Providing both a file and a filehandle for reading from - oly one please!");
  }

  if( !$file && !$fh ) {
      $self->throw("Neither a file (-file) nor a filehandle (-fh) provided to EMBL opening");
  }

  if( $file ) {
      $fh = new FileHandle;
      $fh->open($file) || $self->throw("Could not open $file for EMBL stream reading $!");
  }
  
  # hash for functions for decoding keys.
  $self->{'_func_ftunit_hash'} = {}; 
  $self->_filehandle($fh);
      

# set stuff in self from @args
 return $make; # success - we hope!
}


=head2 next_annseq

 Title   : next_annseq
 Usage   : $seq = $stream->next_annseq()
 Function: returns the next sequence in the stream
 Returns : Bio::AnnSeq object
 Args    :


=cut

sub next_annseq{
   my ($self,@args) = @_;
   my ($seq,$fh,$c,$line,$name,$desc,$seqc);

   $fh = $self->_filehandle();

   if( eof $fh ) {
       return undef; # no throws - end of file
   }

   $line = <$fh>;
   $line =~ /^ID\s+(\S+)/ || $self->throw("EMBL stream with no ID. Not embl in my book");
   $name = $1;

#   $self->warn("not parsing upper annotation in EMBL file yet!");

   my $buffer = $_;

   my $annseq = Bio::AnnSeq->new();
   
   BEFORE_FEATURE_TABLE :
   while( !eof($fh) ) {
       $_ = $buffer;

       /^DE\s+(\S.*\S)/ && do { $desc .= $1; $buffer = <$fh>; next;};
       # accession numbers...
       /^R/ && do {
	   my @refs = &_read_EMBL_References(\$buffer,$fh);
	   $annseq->annotation->add_Reference(@refs);
       };

       # we need to do the upper level annotation
       /^FH/ && last;
       # next line.

       $buffer = <$fh>;
   }

   while( <$fh> ) {
       /FT   \w/ && last;
   }
   $buffer = $_;
   
   FEATURE_TABLE :   
   while( !eof($fh) ) {
       my $ftunit = &_read_FTHelper_EMBL($fh,\$buffer);
       # process ftunit
       &_generic_seqfeature($annseq,$ftunit);

       if( $buffer !~ /^FT/ ) {
	   last;
       }
	   
   }
	
   if( $buffer !~ /^SQ/  ) {
       while( <$fh> ) {
	   /^SQ/ && last;
       }
   }

   $seqc = "";	       
   while( <$fh> ) {
       /^\/\// && last;
       $_ = uc($_);
       s/[^A-Za-z]//g;
       $seqc .= $_;
       eof $fh && last;
   }

   $seq = Bio::Seq->new(-seq => $seqc , -id => $name, -desc => $desc);
   $annseq->seq($seq);
   return $annseq;

}

=head2 write_annseq

 Title   : write_annseq
 Usage   : $stream->write_annseq($seq)
 Function: writes the $seq object (must be annseq) to the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::AnnSeq


=cut

sub write_annseq{
   my ($self,$annseq) = @_;
   my $fh = $self->_filehandle();
   my $seq = $annseq->seq();
   my $i;
   my $str = $seq->seq;

   print $fh "ID   ", $seq->id(), "\nXX   \nDE   ", $seq->desc(), "\nXX   \n";
   my $t = 1;
   foreach my $ref ( $annseq->annotation->each_Reference() ) {
       print $fh "RN   [$t]\n";
       &_write_line_EMBL_regex($fh,"RA   ","RA   ",$ref->authors,'\s+|$',80);       
       &_write_line_EMBL_regex($fh,"RT   ","RT   ",$ref->title,'\s+|$',80);       
       &_write_line_EMBL_regex($fh,"RL   ","RL   ",$ref->location,'\s+|$',80);
       print $fh "XX   \n";
       $t++;
   }

   print $fh "FH   Key             Location/Qualifiers\n";
   print $fh "FH   \n";

   foreach my $sf ( $annseq->top_SeqFeatures ) {
       my @fth = Bio::AnnSeqIO::FTHelper::from_SeqFeature($sf);
       foreach my $fth ( @fth ) {
	   &_print_EMBL_FTHelper($fth,$fh);
       }
   }

   print $fh "SQ   \n";
   print $fh "    ";
   for ($i = 10; $i < length($str); $i += 10) {
       print $fh substr($str,$i,10), " ";
       if( $i%50 == 0 ) {
	   print $fh "\n    ";
       }
   }
   print $fh "\n//\n";
   return 1;
}

=head2 _print_EMBL_FTHelper

 Title   : _print_EMBL_FTHelper
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _print_EMBL_FTHelper{
   my ($fth,$fh) = @_;

   #print $fh "FH   Key             Location/Qualifiers\n";
   #print $fh  sprintf("FT   %-15s  %s\n",$fth->key,$fth->loc);
   &_write_line_EMBL_regex($fh,sprintf("FT   %-15s ",$fth->key),"FT                   ",$fth->loc,',|$',80);
   foreach my $tag ( keys %{$fth->field} ) {
       foreach my $value ( @{$fth->field->{$tag}} ) {
	   &_write_line_EMBL_regex($fh,"FT                   ","FT                   ","/$tag=\"$value\"",'.|$',80);       
	  # print $fh "FT                   /", $tag, "=\"", $value, "\"\n";
       }
   }

}


=head2 _read_EMBL_References

 Title   : _read_EMBL_References
 Usage   :
 Function: Reads references from EMBL format. Internal function really
 Example :
 Returns : 
 Args    :


=cut

sub _read_EMBL_References{
   my ($buffer,$fh) = @_;
   my (@refs);
   
   # assumme things are starting with RN

   if( $$buffer !~ /^RN/ ) {
       warn("Not parsing line $$buffer which maybe important");
   }
   my $title;
   my $loc;
   my $au;
   while( <$fh> ) {
       /^R/ || last;
       /^RA   (.*)/ && do { $au .= $1;   next;};
       /^RT   (.*)/ && do { $title .= $1; next;};
       /^RL   (.*)/ && do { $loc .= $1; next;};
   }

   my $ref = new Bio::Annotation::Reference;
   $au =~ s/;\s*$//g;
   $title =~ s/;\s*$//g;

   $ref->authors($au);
   $ref->title($title);
   $ref->location($loc);
   push(@refs,$ref);
   $$buffer = $_;
   
   return @refs;
}


=head2 _filehandle

 Title   : _filehandle
 Usage   : $obj->_filehandle($newval)
 Function: 
 Example : 
 Returns : value of _filehandle
 Args    : newvalue (optional)


=cut

sub _filehandle{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_filehandle'} = $value;
    }
    return $obj->{'_filehandle'};

}

=head2 _read_FTHelper_EMBL

 Title   : _read_FTHelper_EMBL
 Usage   : &_read_FTHelper_EMBL($fh,$buffer)
 Function: reads the next FT key line
 Example :
 Returns : Bio::AnnSeqIO::FTHelper object 
 Args    : filehandle and reference to a scalar


=cut

sub _read_FTHelper_EMBL{
   my ($fh,$buffer) = @_;
   my ($key,$loc,$out);

   if( $$buffer =~ /^FT\s+(\S+)\s+(\S+)/ ) {
       $key = $1;
       $loc = $2;
   }

   if( $loc =~ /\(/ && $loc !~ /complement\(\d+\.\.\d+\)/ ) {
       # more location to read
       while( <$fh> ) {
	   /^FT\s+(\S+)/ || do { &Bio::Root::Object::throw("Weird location line in EMBL feature table"); };
	   $loc .= $1;
	   $loc =~ /\)/ && last;
       }

   }

   $out = new Bio::AnnSeqIO::FTHelper();
   $out->key($key);
   $out->loc($loc);

   # should read in other fields
   my $current = "";
   while( <$fh> ) {
       /^FT/  || last; # leave on non FT lines!
       /^FT   \w/ && last;
       # field on one line
       /^FT\s+\/(\S+)=\"(.*)\"/ && do {
	   my $key = $1;
	   my $value = $2;
	   if(! defined $out->field->{$key} ) {
	       $out->field->{$key} = [];
	   }
	   push(@{$out->field->{$key}},$value);
	   next;
       };
       # field on on multilines:
       /^FT\s+\/(\S+)=\"(.*)/ && do {
	   my $key = $1;
	   my $value = $2;
	   while ( <$fh> ) {
	       /FT\s+(.*)\"/ && do { $value .= $1; last; };
	       /FT\s+(.*?)\s*/ && do {$value .= $1; };
	   }
	   if(! defined $out->field->{$key} ) {
	       $out->field->{$key} = [];
	   }
	   push(@{$out->field->{$key}},$value);
	   next;
       };
       #field with no quoted value
       /^FT\s+\/(\S+)=(\S+)/ && do {
	   my $key = $1;
	   my $value = $2;
	   if(! defined $out->field->{$key} ) {
	       $out->field->{$key} = [];
	   }
	   push(@{$out->field->{$key}},$value);
	   next;
       };


   }

   $$buffer = $_;
   
   return $out;
}

=head2 _generic_seqfeature

 Title   : _generic_seqfeature
 Usage   : &_generic_seqfeature($annseq,$fthelper)
 Function: processes fthelper into a generic seqfeature
 Returns : nothing (places a new seqfeature into annseq)
 Args    : Bio::AnnSeq,Bio::AnnSeqIO::FTHelper


=cut

sub _generic_seqfeature{
   my ($annseq,$fth) = @_;
   my ($sf);

  # print "Converting from", $fth->key, "\n";

   $sf = new Bio::SeqFeature::Generic;
   if( $fth->loc =~ /join/ ) {
       my $strand;
       if ( $fth->loc =~ /complement/ ) {
	   $strand = -1;
       } else {
	   $strand = 1;
       }

       $sf->strand($strand);
       $sf->primary_tag($fth->key . "_span");
       $sf->source_tag('EMBL');
       $sf->has_tag("parent",1);
       $sf->_parse->{'parent_homogenous'} = 1;

       # we need to make sub features
       my $loc = $fth->loc;
       while( $loc =~ /(\d+)\.\.(\d+)/g ) {
	   my $start = $1;
	   my $end   = $2;
	   #print "Processing $start-$end\n";
	   my $sub = new Bio::SeqFeature::Generic;
	   $sub->primary_tag($fth->key);
	   $sub->start($start);
	   $sub->end($end);
	   $sub->strand($strand);
	   $sub->source_tag('EMBL');
	   $sf->add_sub_SeqFeature($sub,'EXPAND');
       }

   } else {
       my $lst;
       my $len;

       if( $fth->loc =~ /^(\d+)$/ ) {
	   $lst = $len = $1;
       } else {
	   $fth->loc =~ /(\d+)\.\.(\d+)/ || do {
	       $annseq->throw("Weird location line [" . $fth->loc . "] in reading EMBL");
	       last;
	   };
	   $lst = $1;
	   $len = $2;
       }

       $sf->start($lst);
       $sf->end($len);
       $sf->source_tag('EMBL');
       $sf->primary_tag($fth->key);
       if( $fth->loc =~ /complement/ ) {
	   $sf->strand(-1);
       } else {
	   $sf->strand(1);
       }
   }

   #print "Adding B4 ", $sf->primary_tag , "\n";

   foreach my $key ( keys %{$fth->field} ){
       foreach my $value ( @{$fth->field->{$key}} ) {
	   $sf->add_tag_value($key,$value);
       }
   }



   $annseq->add_SeqFeature($sf);
}

=head2 _write_line_EMBL

 Title   : _write_line_EMBL
 Usage   :
 Function: internal function
 Example :
 Returns : 
 Args    :


=cut

sub _write_line_EMBL{
   my ($fh,$pre1,$pre2,$line,$length) = @_;

   $length || die "Miscalled write_line_EMBL without length. Programming error!";
   my $subl = $length - length $pre2;
   my $linel = length $line;
   my $i;

   my $sub = substr($line,0,$length - length $pre1);

   print $fh "$pre1$sub\n";
   
   for($i= ($length - length $pre1);$i < $linel;) {
       $sub = substr($line,$i,($subl));
       print $fh "$pre2$sub\n";
       $i += $subl;
   }

}

=head2 _write_line_EMBL_regex

 Title   : _write_line_EMBL_regex
 Usage   :
 Function: internal function for writing lines of specified
           length, with different first and the next line 
           left hand headers and split at specific points in the
           text
 Example :
 Returns : nothing
 Args    : file handle, first header, second header, text-line, regex for line breaks, total line length


=cut

sub _write_line_EMBL_regex {
   my ($fh,$pre1,$pre2,$line,$regex,$length) = @_;

   
   #print STDOUT "Going to print with $line!\n";

   $length || die "Miscalled write_line_EMBL without length. Programming error!";

   if( length $pre1 != length $pre2 ) {
       die "Programming error - cannot called write_line_EMBL_regex with different length pre1 and pre2 tags!";
   }

   my $subl = $length - length $pre1;
   my @lines;

   while($line =~ m/(.{1,$subl})($regex)/g) {
       push(@lines, $1.$2);
   }
   
   my $s = shift @lines;
   print $fh "$pre1$s\n";
   foreach my $s ( @lines ) {
       print $fh "$pre2$s\n";
   }
}

    





