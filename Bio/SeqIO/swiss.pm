

#
# BioPerl module for Bio::SeqIO::swissprot
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::swissprot - Swissprot sequence input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the SeqIO handler system. Go:

    $stream = Bio::SeqIO->new(-file => $filename, -format => 'swiss');

    while ( my $seq = $stream->next_seq() ) {
	# do something with $seq
    }

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from swissprot flat
file databases.

There is alot of flexibility here about how to dump things which I need
to document fully.


=head2 Optional functions

=over

=item _show_dna()

(output only) shows the dna or not

=item _post_sort

(output only) provides a sorting func which is applied to the FTHelpers
before printing

=item _id_generation_func

This is function which is called as 

   print "ID   ", $func($seq), "\n";

To generate the ID line. If it is not there, it generates a sensible ID
line using a number of tools.

 
=end

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

=head1 AUTHOR - Elia Stupka

Email elia@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqIO::swiss;
use vars qw(@ISA);
use strict;
use Bio::Seq;
use Bio::SeqIO::FTHelper;
use Bio::SeqFeature::Generic;
use Bio::Species;
#use Bio::SeqIO::StreamI;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;
use FileHandle;

@ISA = qw(Bio::SeqIO);
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
      $self->throw("Neither a file (-file) nor a filehandle (-fh) provided to swissprot opening");
  }

  if( $file ) {
      $fh = new FileHandle;
      $fh->open($file) || $self->throw("Could not open $file for swissprot stream reading $!");
  }
  
  # hash for functions for decoding keys.
  $self->{'_func_ftunit_hash'} = {}; 
  $self->_filehandle($fh);
  $self->_show_dna(1); # sets this to one by default. People can change it

# set stuff in self from @args
 return $make; # success - we hope!
}


=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    :


=cut

sub next_seq{
   my ($self,@args) = @_;
   my ($pseq,$fh,$c,$line,$name,$desc,$acc,$seqc,$mol,$div, $date, $comment, @date_arr);
   my $seq = Bio::Seq->new();

   $fh = $self->_filehandle();

   if( eof $fh ) {
       return undef; # no throws - end of file
   }

   $line = <$fh>;

   if( $line =~ /^\s+$/ ) {
       while( <$fh> ) {
   	   /\S/ && last;
       }
       $line = $_;
   }
   
   if( !defined $line ) {
       return undef; # end of file
   }

   $line =~ /^ID\s+(\S+)/ || $self->throw("swissprot stream with no ID. Not swissprot in my book");
   $name = $1;
   
   my $buffer = $line;

   BEFORE_FEATURE_TABLE :
   until( eof($fh) ) {
       $_ = $buffer;
       
       # Exit at start of Feature table
       last if /^FT/;

       # Description line(s)
       if (/^DE\s+(\S.*\S)/) {
           $desc .= $desc ? " $1" : $1;
       }

       #Gene name
       if (/^GN\s+(\S+)/) {
	   $seq->annotation->gene_name($1);
       }  

       #accession number
       if( /^AC\s+(\S+); (\S+)?/) {
	   $acc = $1;
	   my $acc2 = $2 if $2;
	   $acc2 =~ s/\;//;
	   $seq->accession($acc);
	   $seq->add_secondary_accession($acc2);
       }
       
       #version number
       if( /^SV\s+(\S+);?/ ) {
	   my $sv = $1;
	   $sv =~ s/\;//;
	   $seq->sv($sv);
       }

       #date (NOTE: takes last date line)
       if( /^DT\s+(\S+)/ ) {
	   my $date = $1;
	   $date =~ s/\;//;
	   $seq->add_date($date);
       }
       
     

       # Organism name and phylogenetic information
       if (/^O[SC]/) {
           my $species = _read_swissprot_Species(\$buffer, $fh);
           $seq->species( $species );
       }

       # References
       if (/^R/) {
	   my @refs = &_read_swissprot_References(\$buffer,$fh);
	   $seq->annotation->add_Reference(@refs);
       }

     
       #Comments
       if (/^CC\s+(.*)/) {
	   $comment .= $1;
	   $comment .= " ";
	   while (<$fh>) {
	       if (/^CC\s+(.*)/) {
		   $comment .= $1;
		   $comment .= " ";
	       }
	       else { 
		   last;
	       }
	   }
	   my $commobj = Bio::Annotation::Comment->new();
	   $comment =~ s/-\!- //g;
	   $comment =~ s/    //g;
	   $commobj->text($comment);
	   $seq->annotation->add_Comment($commobj);
	   $comment = "";
       }

       #DBLinks
       if (/^DR\s+(\S+)\; (\S+)\; (\S+)\; (\S+)/) {
	   my $dblinkobj =  Bio::Annotation::DBLink->new();
	   $dblinkobj->database($1);
	   $dblinkobj->primary_id($2);
	   $dblinkobj->optional_id($3);
	   $dblinkobj->comment($4);
	   $seq->annotation->add_DBLink($dblinkobj);
       }

       #keywords
       if( /^KW   (.*)\S*$/ ) {
	   my $keywords = $1;
	   $seq->keywords($keywords);
       }

       # Get next line.
       $buffer = <$fh>;
   }
   
   $buffer = $_;
      
   FEATURE_TABLE :   
   while (<$fh>) {

       my $ftunit = &_read_FTHelper_swissprot($fh,\$buffer);
       
       # process ftunit
       $ftunit->_generic_seqfeature($seq);
       
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
   }

   $pseq = Bio::PrimarySeq->new(-seq => $seqc , -id => $name, -desc => $desc);
   $seq->primary_seq($pseq);
   return $seq;

}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq($seq)
 Function: writes the $seq object (must be seq) to the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq


=cut

sub write_seq {
   my ($self,$seq) = @_;

   if( !defined $seq ) {
       $self->throw("Attempting to write with no seq!");
   }

   if( ! ref $seq || ! $seq->isa('Bio::SeqI') ) {
       $self->warn(" $seq is not a SeqI compliant module. Attempting to dump, but may fail!");
   }

   my $fh = $self->_filehandle();
   my $i;
   my $str = $seq->seq;
   
   my $mol;
   my $div;
   my $len = $seq->length();

   if ( !$seq->can('division') || ! defined ($div = $seq->division()) ) {
       $div = 'UNK';
   }
   else {
       $div=$seq->division;
   }
   $mol = $seq->molecule;
  
   if(!$seq->can('molecule') || (!defined($mol))) {
       $mol = 'XXX';
   }
   else {
       $mol = $seq->molecule;
   }
   
   my $temp_line;
   if( $self->_id_generation_func ) {
       $temp_line = &{$self->_id_generation_func}($seq);
   } else {
       $temp_line = sprintf ("%s_$div    STANDARD;       $mol;   %d AA.",$seq->id(),$len);
   } 

   print $fh "ID   $temp_line\n";   

   # if there, write the accession line
   local($^W) = 0;   # supressing warnings about uninitialized fields

   if( $self->_ac_generation_func ) {
       $temp_line = &{$self->_ac_generation_func}($seq);
       print $fh "AC   $temp_line\n";   
   } else {
       if ($seq->can('accession') ) {
	   print "AC   ",$seq->accession,";";
	   if ($seq->can('each_secondary_accession') ) {
	       print " ",$seq->each_secondary_accession,";\n";
	   }
	   else {
	       print "\n";
	   }
       }
       # otherwise - cannot print <sigh>
   } 

   #Date lines
   foreach my $dt ( $seq->each_date() ) {
       _write_line_swissprot_regex($fh,"DT   ","DT   ",$dt,"\\s\+\|\$",80);
   }
   
   #Definition lines
   _write_line_swissprot_regex($fh,"DE   ","DE   ",$seq->desc(),"\\s\+\|\$",80);

   #Gene name
   if ($seq->annotation->can('gene_name')) {
       print "GN   ",$seq->annotation->gene_name,"\n";
   }
   
   # Organism lines
   if (my $spec = $seq->species) {
       my($sub_species, $species, $genus, @class) = $spec->classification();
       my $OS = "$genus $species $sub_species";
       if (my $common = $spec->common_name) {
	   $OS .= " ($common)";
       }
       print $fh "OS   $OS\n";
       my $OC = join('; ', reverse(@class));
       $OC =~ s/\;\s+$//;
       $OC .= ".";
       _write_line_swissprot_regex($fh,"OC   ","OC   ",$OC,"\; \|\$",80);
	if ($spec->organelle) {
	    _write_line_swissprot_regex($fh,"OG   ","OG   ",$spec->organelle,"\; \|\$",80);
	}
   }
   
   # Reference lines
   my $t = 1;
   foreach my $ref ( $seq->annotation->each_Reference() ) {
       print $fh "RN   [$t]\n";
       if ($ref->start && $ref->end) {
	   print "RP   SEQUENCE OF ",$ref->start,"-",$ref->end," FROM N.A.\n";
       }
       elsif ($ref->rp) {
	   print "RP   ",$ref->rp,"\n";
       } 

       &_write_line_swissprot_regex($fh,"RA   ","RA   ",$ref->authors,"\\s\+\|\$",80);       
       &_write_line_swissprot_regex($fh,"RT   ","RT   ",$ref->title,"\\s\+\|\$",80);       
       &_write_line_swissprot_regex($fh,"RL   ","RL   ",$ref->location,"\\s\+\|\$",80);
       if ($ref->comment) {
	   &_write_line_swissprot_regex($fh,"RC   ","RC   ",$ref->comment,"\\s\+\|\$",80); 
       }
       $t++;
   }
   
   # Comment lines

   foreach my $comment ( $seq->annotation->each_Comment() ) {
       _write_line_swissprot_regex($fh,"CC   ","CC   ",$comment->text,"\\s\+\|\$",80);
   }

   foreach my $dblink ( $seq->annotation->each_DBLink() ) {
       print $fh "DR   ",$dblink->database,"; ",$dblink->primary_id,"; ",$dblink->optional_id,"; ",$dblink->comment,"\n";
   }
   

   # if there, write the kw line
   
   if( $self->_kw_generation_func ) {
       $temp_line = &{$self->_kw_generation_func}($seq);
       print $fh "KW   $temp_line\n";   
   } else {
       if( $seq->can('keywords') ) {
	   print $fh "KW   ",$seq->keywords,"\n";
       }
   } 
   if( defined $self->_post_sort ) {
       # we need to read things into an array. Process. Sort them. Print 'em

       my $post_sort_func = $self->_post_sort();
       my @fth;

       foreach my $sf ( $seq->top_SeqFeatures ) {
	   push(@fth,Bio::SeqIO::FTHelper::from_SeqFeature($sf,$seq));
       }

       @fth = sort { &$post_sort_func($a,$b) } @fth;
       
       foreach my $fth ( @fth ) {
	   &_print_swissprot_FTHelper($fth,$fh);
       }
   } else {
       # not post sorted. And so we can print as we get them.
       # lower memory load...
       
       foreach my $sf ( $seq->top_SeqFeatures ) {
	   my @fth = Bio::SeqIO::FTHelper::from_SeqFeature($sf,$seq);
	   foreach my $fth ( @fth ) {
	       if( ! $fth->isa('Bio::SeqIO::FTHelper') ) {
		   $sf->throw("Cannot process FTHelper... $fth");
	       }

	       &_print_swissprot_FTHelper($fth,$fh);
	   }
       }
   }

   if( $self->_show_dna() == 0 ) {
       return;
   }

   # finished printing features.

   print $fh "SQ   SEQUENCE   $len AA;\n";
   print $fh "     ";
   my $linepos;
   for ($i = 0; $i < length($str); $i += 10) {
       print $fh substr($str,$i,10), " ";
       $linepos += 11;
       if( ($i+10)%60 == 0 ) {
	   print $fh "\n     ";
      }
   }
   print $fh "\n//\n";
   return 1;
}

=head2 _print_swissprot_FTHelper

 Title   : _print_swissprot_FTHelper
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _print_swissprot_FTHelper {
   my ($fth,$fh,$always_quote) = @_;
   
   if( ! ref $fth || ! $fth->isa('Bio::SeqIO::FTHelper') ) {
       $fth->warn("$fth is not a FTHelper class. Attempting to print, but there could be tears!");
   }
   my $loc= $fth->loc;
   $loc =~ /(\d+)..(\d+)/;
   my $loc1 = $1;
   my $loc2 = $2;
#   my $line = sprintf ("FT   %-12s %-5s %s",$fth->key,$loc1,$loc2);
#   print $fh "$line";
   my $switch=0; 
   foreach my $tag ( keys %{$fth->field} ) {
       foreach my $value ( @{$fth->field->{$tag}} ) {
	   $value =~ s/\"/\"\"/g;
	   if ($switch == 0) {
	       my $line = sprintf ("FT   %-12s %-5s %-10s%s",$fth->key,$loc1,$loc2,$value);
	       print $fh "$line\n";
	   }
           else {
	       print $fh "FT                                /$tag\n";
           } 
	   $switch = 1;
	  # print $fh "FT                   /", $tag, "=\"", $value, "\"\n";
       }
   }

}


=head2 _read_swissprot_References

 Title   : _read_swissprot_References
 Usage   :
 Function: Reads references from swissprot format. Internal function really
 Example :
 Returns : 
 Args    :


=cut

sub _read_swissprot_References{
   my ($buffer,$fh) = @_;
   my (@refs);
   my ($b1, $b2, $rp, $title, $loc, $au, $med, $com);
   
   if ($$buffer !~ /^RP/) {
       $$buffer = <$fh>;
   }
   if( $$buffer =~ /^RP/ ) {
       if ($$buffer =~ /^RP   SEQUENCE OF (\d+)-(\d+)/) { 
	   $b1=$1;
	   $b2=$2; 
       }
       elsif ($$buffer =~ /^RP   (.*)/) {
	   $rp=$1;
       }
       
   }
   while( <$fh> ) {
       /^CC/ && goto OUT;
       /^RN/ && last;
       /^RX   MEDLINE;\s+(\d+)/ && do {$med=$1};
       /^RA   (.*)/ && do { $au .= $1;   next;};
       /^RT   (.*)/ && do { $title .= $1; next;};
       /^RL   (.*)/ && do { $loc .= $1; next;};
       /^RC   (.*)/ && do { $com .= $1; next;};
   }
   
   OUT: my $ref = new Bio::Annotation::Reference;
   $au =~ s/;\s*$//g;
   if( defined $title ) {
       $title =~ s/;\s*$//g;
   }
   
   $ref->start($b1);
   $ref->end($b2);
   $ref->authors($au);
   $ref->title($title);
   $ref->location($loc);
   $ref->medline($med);
   $ref->comment($com);
   $ref->rp($rp);

   push(@refs,$ref);
   $$buffer = $_;
   return @refs;
}


=head2 _read_swissprot_Species

 Title   : _read_swissprot_Species
 Usage   :
 Function: Reads the swissprot Organism species and classification
           lines.
 Example :
 Returns : A Bio::Species object
 Args    :

=cut

sub _read_swissprot_Species {
    my( $buffer, $fh ) = @_;
    my $org;

    $_ = $$buffer;
    my( $sub_species, $species, $genus, $common, @class );
    while (defined( $_ ||= <$fh> )) {
        
        if (/^OS\s+(\S+)\s+(\S+)\s+(\S+)?(?:\s+\((.*)\))?/) {
            $genus   = $1;
	    if ($2) {
		$species = $2
		}
	    else {
		$species = "sp.";
	    }
	    $sub_species = $3 if $3;
            $common  = $4 if $4;
        }
        elsif (s/^OC\s+//) {
            push(@class, split /[\;\s\.]+/);
        }
	elsif (/^OG\s+(.*)/) {
	    $org = $1;
	}
        else {
            last;
        }
        
        $_ = undef; # Empty $_ to trigger read of next line
    }
    
    $$buffer = $_;
    
    # Don't make a species object if it's "Unknown" or "None"
    return if $genus =~ /^(Unknown|None)$/i;

    # Bio::Species array needs array in Species -> Kingdom direction
    if ($sub_species) {
	push( @class, $genus, $species, $sub_species);
    }
    else {
	push( @class, $genus, $species, "");
    }
    @class = reverse @class;
    
    my $make = Bio::Species->new();
    $make->classification( @class );
    $make->common_name( $common ) if $common;
    $make->organelle($org) if $org;
    return $make;
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

=head2 _read_FTHelper_swissprot

 Title   : _read_FTHelper_swissprot
 Usage   : &_read_FTHelper_swissprot($fh,$buffer)
 Function: reads the next FT key line
 Example :
 Returns : Bio::SeqIO::FTHelper object 
 Args    : filehandle and reference to a scalar


=cut

sub _read_FTHelper_swissprot {
   my ($fh,$buffer) = @_;
   my ($key,$loc,$out,$value);
   if ( $$buffer !~ /^FT\s+(\S+)/) {
       $out->throw("Weird location line in swissprot feature table: '$_'");
   }
   #Read key and location lines
   if( $$buffer =~ /^FT\s+(\S+)\s+(\d+)\s+(\d+)\s+(.*)/ ) {
       $key = $1;
       $loc = $2;
       $loc .= "..";
       $loc .= $3;
       $value = $4;
   }

   $out = new Bio::SeqIO::FTHelper();
   $loc =~ s/<//;
   $loc =~ s/>//;
   $out->key($key);
   $out->loc($loc);
   # $out->field->{$key};   # What's the purpose of this? I changed to next line. SAC 2/21/00
   $out->field($key);
   $value =~ s/\"\"/\"/g;
   push (@{$out->field->{$key}},$value);
   
   # Now read in other fields
   # Loop reads $_ when defined (i.e. only in first loop), then $fh, until end of file
   while( defined($_ ||= <$fh>) ) {
       
       # Exit loop on non FT lines!
       /^FT/  || last;
       
       # Exit loop on new primary key
       /^FT   \w/ && last;
       
       # Field on one line
       if (/^FT\s+\/(\S+)=\"(.+)\"/) {
	   my $key = $1;
	   my $value = $2;
	   if(! defined $out->field->{$key} ) {
	       $out->field->{$key} = [];
	   }
	   $value =~ s/\"\"/\"/g;
	   
	   push(@{$out->field->{$key}},$value);
       }
       # Field on on multilines:
       elsif (/^FT\s+\/(\S+)=\"(.*)/) {
	   my $key = $1;
	   my $value = $2;
	   while ( <$fh> ) {
	       s/\"\"/__DOUBLE_QUOTE_STRING__/g;
	       /FT\s+(.*)\"/ && do { $value .= $1; last; };
	       /FT\s+(.*)/ && do {$value .= $1; };
	   }
	   $value =~ s/__DOUBLE_QUOTE_STRING__/\"/g;

	   if(! defined $out->field->{$key} ) {
	       $out->field->{$key} = [];
	   }
	   push(@{$out->field->{$key}},$value);
       }
       # Field with no quoted value
       elsif (/^FT\s+\/(\S+)=?(\S+)?/) {
	   my $key = $1;
	   my $value = $2 if $2;
	   $value = "_no_value" unless $2;
	   if(! defined $out->field->{$key} ) {
	       $out->field->{$key} = [];
	   }
	   push(@{$out->field->{$key}},$value);
       }
       
       # Empty $_ to trigger read from $fh
       undef $_;
   }

   $$buffer = $_;
   return $out;
}


=head2 _write_line_swissprot

 Title   : _write_line_swissprot
 Usage   :
 Function: internal function
 Example :
 Returns : 
 Args    :


=cut

sub _write_line_swissprot{
   my ($fh,$pre1,$pre2,$line,$length) = @_;

   $length || die "Miscalled write_line_swissprot without length. Programming error!";
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

=head2 _write_line_swissprot_regex

 Title   : _write_line_swissprot_regex
 Usage   :
 Function: internal function for writing lines of specified
           length, with different first and the next line 
           left hand headers and split at specific points in the
           text
 Example :
 Returns : nothing
 Args    : file handle, first header, second header, text-line, regex for line breaks, total line length


=cut

sub _write_line_swissprot_regex {
   my ($fh,$pre1,$pre2,$line,$regex,$length) = @_;

   
   #print STDOUT "Going to print with $line!\n";

   $length || die "Miscalled write_line_swissprot without length. Programming error!";

   if( length $pre1 != length $pre2 ) {
       die "Programming error - cannot called write_line_swissprot_regex with different length pre1 and pre2 tags!";
   }

   my $subl = $length - (length $pre1) -1 ;
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

=head2 _post_sort

 Title   : _post_sort
 Usage   : $obj->_post_sort($newval)
 Function: 
 Returns : value of _post_sort
 Args    : newvalue (optional)


=cut

sub _post_sort{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_post_sort'} = $value;
    }
    return $obj->{'_post_sort'};

}

=head2 _show_dna

 Title   : _show_dna
 Usage   : $obj->_show_dna($newval)
 Function: 
 Returns : value of _show_dna
 Args    : newvalue (optional)


=cut

sub _show_dna{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_show_dna'} = $value;
    }
    return $obj->{'_show_dna'};

}

=head2 _id_generation_func

 Title   : _id_generation_func
 Usage   : $obj->_id_generation_func($newval)
 Function: 
 Returns : value of _id_generation_func
 Args    : newvalue (optional)


=cut

sub _id_generation_func{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_id_generation_func'} = $value;
    }
    return $obj->{'_id_generation_func'};

}

=head2 _ac_generation_func

 Title   : _ac_generation_func
 Usage   : $obj->_ac_generation_func($newval)
 Function: 
 Returns : value of _ac_generation_func
 Args    : newvalue (optional)


=cut

sub _ac_generation_func{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_ac_generation_func'} = $value;
    }
    return $obj->{'_ac_generation_func'};

}

=head2 _sv_generation_func

 Title   : _sv_generation_func
 Usage   : $obj->_sv_generation_func($newval)
 Function: 
 Returns : value of _sv_generation_func
 Args    : newvalue (optional)


=cut

sub _sv_generation_func{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_sv_generation_func'} = $value;
    }
    return $obj->{'_sv_generation_func'};

}

=head2 _kw_generation_func

 Title   : _kw_generation_func
 Usage   : $obj->_kw_generation_func($newval)
 Function: 
 Returns : value of _kw_generation_func
 Args    : newvalue (optional)


=cut

sub _kw_generation_func{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_kw_generation_func'} = $value;
    }
    return $obj->{'_kw_generation_func'};

}

1;
    





