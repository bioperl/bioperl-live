

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

 
=back

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

sub next_seq {
   my ($self,@args) = @_;
   my ($pseq,$c,$line,$name,$desc,$acc,$seqc,$mol,$div, $date, $comment, @date_arr);
   my $seq = Bio::Seq->new();
   $line = $self->_readline;

   if( !defined $line) {
       return undef; # no throws - end of file
   }
   
   if( $line =~ /^\s+$/ ) {
       while( defined ($line = $self->_readline) ) {
   	   $line =~ /\S/ && last;
       }
   }   
   if( !defined $line ) {
       return undef; # end of file
   }

   $line =~ /^ID\s+(\S+)/ || $self->throw("swissprot stream with no ID. Not swissprot in my book");
   $name = $1;
    # this is important to have the id for display in e.g. FTHelper, otherwise
    # you won't know which entry caused an error
   $seq->display_id($name);
   
   my $buffer = $line;

   BEFORE_FEATURE_TABLE :
   until( !defined ($buffer) ) {
       $_ = $buffer;
       
       # Exit at start of Feature table
       last if /^FT/;
       # and at the sequence at the latest HL 05/11/2000
       last if /^SQ/;

       # Description line(s)
       if (/^DE\s+(\S.*\S)/) {
           $desc .= $desc ? " $1" : $1;
       }
       #Gene name
       elsif (/^GN\s+(\S+)/) {
	   $seq->annotation->gene_name($1);
       }  
       #accession number
       elsif( /^AC\s+(.+)/) {
	   my( $acc, @sec ) = split /[\s;]+/, $1;
	   $seq->accession($acc);
           foreach my $s (@sec) {
	      $seq->add_secondary_accession($s);
           }
       }
       #version number
       elsif( /^SV\s+(\S+);?/ ) {
	   my $sv = $1;
	   $sv =~ s/\;//;
	   $seq->sv($sv);
       }
       #date (NOTE: takes last date line)
       elsif( /^DT\s+(\S+)/ ) {
	   my $date = $1;
	   $date =~ s/\;//;
	   $seq->add_date($date);
       }
       # Organism name and phylogenetic information
       elsif (/^O[SC]/) {
           my $species = $self->_read_swissprot_Species(\$buffer);
           $seq->species( $species );
	   # now we are one line ahead -- so continue without reading the next
	   # line   HL 05/11/2000
	   next;
       }
       # References
       elsif (/^R/) {
	   my @refs = $self->_read_swissprot_References(\$buffer);
	   $seq->annotation->add_Reference(@refs);
	   # now we are one line ahead -- so continue without reading the next
	   # line   HL 05/11/2000
	   next;
       }
       #Comments
       elsif (/^CC\s+(.*)/) {
	   $comment .= $1;
	   $comment .= " ";
	   while (defined ($buffer = $self->_readline)) {
	       if ($buffer =~ /^CC\s+(.*)/) {
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
	   # now we are one line ahead -- so continue without reading the next
	   # line   HL 05/11/2000
	   next;
       }
       #DBLinks
       elsif (/^DR\s+(\S+)\; (\S+)\; (\S+)\; (\S+)/) {
	   my $dblinkobj =  Bio::Annotation::DBLink->new();
	   $dblinkobj->database($1);
	   $dblinkobj->primary_id($2);
	   $dblinkobj->optional_id($3);
	   $dblinkobj->comment($4);
	   $seq->annotation->add_DBLink($dblinkobj);
       }
       elsif (/^DR\s+(\S+)\; (\S+)\;/) {
	   my $dblinkobj =  Bio::Annotation::DBLink->new();
	   $dblinkobj->database($1);
	   $dblinkobj->primary_id($2);
	   $seq->annotation->add_DBLink($dblinkobj);
       }
       #keywords
       elsif( /^KW   (.*)\S*$/ ) {
	   my $keywords = $1;
	   $seq->keywords($keywords);
       }

       # Get next line. Getting here assumes that we indeed need to read the
       # line.
       $buffer = $self->_readline;
   }
   
   $buffer = $_;
      
   FEATURE_TABLE :
   # if there is no feature table, or if we've got beyond, exit loop or don't
   # even enter    HL 05/11/2000
   while (defined ($buffer) && ($buffer =~ /^FT/)) {
       my $ftunit = $self->_read_FTHelper_swissprot(\$buffer);
       
       # process ftunit
       # when parsing of the line fails we get undef returned
       if($ftunit) {
	   $ftunit->_generic_seqfeature($seq, "SwissProt");
       } else {
	   $self->warn("failed to parse feature table line for seq " .
		       $seq->display_id());
       }
   }
   if( $buffer !~ /^SQ/  ) {
       while( defined($_ = $self->_readline) ) {
	   /^SQ/ && last;
       }
   }
   $seqc = "";	       
   while( defined ($_ = $self->_readline) ) {
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
       #$temp_line = sprintf ("%10s     STANDARD;      %3s;   %d AA.",
       #		     $seq->primary_id()."_".$div,$mol,$len);
       # Reconstructing the ID relies heavily upon the input source having
       # been in a format that is parsed as this routine expects it -- that is,
       # by this module itself. This is bad, I think, and immediately breaks
       # if e.g. the Bio::DB::GenPept module is used as input.
       # Hence, switch to display_id(); _every_ sequence is supposed to have
       # this. HL 2000/09/03
       $temp_line = sprintf ("%10s     STANDARD;      %3s;   %d AA.",
			     $seq->display_id(),$mol,$len);
   } 

   $self->_print( "ID   $temp_line\n");   

   # if there, write the accession line
   local($^W) = 0;   # supressing warnings about uninitialized fields

   if( $self->_ac_generation_func ) {
       $temp_line = &{$self->_ac_generation_func}($seq);
       $self->_print( "AC   $temp_line\n");   
   } else {
       if ($seq->can('accession') ) {
	   $self->_print("AC   ",$seq->accession,";");
	   if ($seq->can('each_secondary_accession') ) {
	       $self->_print(" ",$seq->each_secondary_accession,";\n");
	   }
	   else {
	       $self->_print("\n");
	   }
       }
       # otherwise - cannot print <sigh>
   } 

   #Date lines
   foreach my $dt ( $seq->each_date() ) {
       $self->_write_line_swissprot_regex("DT   ","DT   ",$dt,"\\s\+\|\$",80);
   }
   
   #Definition lines
   $self->_write_line_swissprot_regex("DE   ","DE   ",$seq->desc(),"\\s\+\|\$",80);

   #Gene name
   if ($seq->annotation->can('gene_name')) {
       $self->_print("GN   ",$seq->annotation->gene_name,"\n");
   }
   
   # Organism lines
   if ($seq->can('species') && (my $spec = $seq->species)) {
        my($species, @class) = $spec->classification();
        my $genus = $class[0];
        my $OS = "$genus $species";
	if (my $ssp = $spec->sub_species) {
            $OS .= " $ssp";
        }
        if (my $common = $spec->common_name) {
            $OS .= " ($common)";
        }
        $self->_print( "OS   $OS\n");
        my $OC = join('; ', reverse(@class)) .'.';
        $self->_write_line_swissprot_regex("OC   ","OC   ",$OC,"\; \|\$",80);
	if ($spec->organelle) {
	    $self->_write_line_swissprot_regex("OG   ","OG   ",$spec->organelle,"\; \|\$",80);
	}
   }
   
   # Reference lines
   my $t = 1;
   foreach my $ref ( $seq->annotation->each_Reference() ) {
       $self->_print( "RN   [$t]\n");
       if ($ref->start && $ref->end) {
	   # changed by jason <jason@chg.mc.duke.edu> on 3/16/00 
#	   print "RP   SEQUENCE OF ",$ref->start,"-",$ref->end," FROM N.A.\n";
	   # to
	   $self->_print("RP   SEQUENCE OF ",$ref->start,"-",$ref->end," FROM N.A.\n");
       }
       elsif ($ref->rp) {
	   # changed by jason <jason@chg.mc.duke.edu> on 3/16/00 
#	   print "RP   ",$ref->rp,"\n";
	   # to
	   $self->_print("RP   ",$ref->rp,"\n");
       } 

       $self->_write_line_swissprot_regex("RA   ","RA   ",$ref->authors,"\\s\+\|\$",80);       
       $self->_write_line_swissprot_regex("RT   ","RT   ",$ref->title,"\\s\+\|\$",80);       
       $self->_write_line_swissprot_regex("RL   ","RL   ",$ref->location,"\\s\+\|\$",80);
       if ($ref->comment) {
	   $self->_write_line_swissprot_regex("RC   ","RC   ",$ref->comment,"\\s\+\|\$",80); 
       }
       $t++;
   }
   
   # Comment lines

   foreach my $comment ( $seq->annotation->each_Comment() ) {
       $self->_write_line_swissprot_regex("CC   ","CC   ",$comment->text,"\\s\+\|\$",80);
   }

   foreach my $dblink ( $seq->annotation->each_DBLink() ) {
       $self->_print("DR   ",$dblink->database,"; ",$dblink->primary_id,"; ",$dblink->optional_id,"; ",$dblink->comment,"\n");
   }
   

   # if there, write the kw line
   
   if( $self->_kw_generation_func ) {
       $temp_line = &{$self->_kw_generation_func}($seq);
       $self->_print( "KW   $temp_line\n");   
   } else {
       if( $seq->can('keywords') ) {
	   $self->_print( "KW   ",$seq->keywords,"\n");
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
	   $self->_print_swissprot_FTHelper($fth);
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

	       $self->_print_swissprot_FTHelper($fth);
	   }
       }
   }

   if( $self->_show_dna() == 0 ) {
       return;
   }

   # finished printing features.

   $self->_print( "SQ   SEQUENCE   $len AA;\n");
   $self->_print( "     ");
   my $linepos;
   for ($i = 0; $i < length($str); $i += 10) {
       $self->_print( substr($str,$i,10), " ");
       $linepos += 11;
       if( ($i+10)%60 == 0 ) {
	   $self->_print( "\n     ");
      }
   }
   $self->_print( "\n//\n");
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
   my ($self,$fth,$always_quote) = @_;
   
    ### FIXME - not implemented
    warn "_print_swissprot_FTHelper NOT IMPLEMENTED";

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
   my ($self,$buffer) = @_;
   my (@refs);
   my ($b1, $b2, $rp, $title, $loc, $au, $med, $com);
   
   if ($$buffer !~ /^RP/) {
       $$buffer = $self->_readline;
   }
   if( !defined $$buffer ) { return undef; }
   if( $$buffer =~ /^RP/ ) {
       if ($$buffer =~ /^RP   SEQUENCE OF (\d+)-(\d+)/) { 
	   $b1=$1;
	   $b2=$2; 
       }
       elsif ($$buffer =~ /^RP   (.*)/) {
	   $rp=$1;
       }
       
   }
   while( defined ($_ = $self->_readline) ) {
       #/^CC/ && last;
       #/^RN/ && last;
       #/^SQ/ && last; # there may be sequences without CC lines! HL 05/11/2000
       /^[^R]/ && last; # may be the safest exit point HL 05/11/2000
       /^RX   MEDLINE;\s+(\d+)/ && do {$med=$1};
       /^RA   (.*)/ && do { $au .= $1;   next;};
       /^RT   (.*)/ && do { $title .= $1; next;};
       /^RL   (.*)/ && do { $loc .= $1; next;};
       /^RC   (.*)/ && do { $com .= $1; next;};
   }
   
   my $ref = new Bio::Annotation::Reference;
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
    my( $self, $buffer ) = @_;
    my $org;

    $_ = $$buffer;
    my( $sub_species, $species, $genus, $common, @class );
    while (defined( $_ ||= $self->_readline )) {
        if (/^OS\s+(\S+)(?:\s+([^\(]\S*))?(?:\s+([^\(]\S*))?(?:\s+\((.*)\))?/) {
            $genus   = $1;
	    if ($2) {
		$species = $2;
		# remove trailing dot -- TrEMBL has that. HL 05/11/2000
		$species =~ s/.$//;
	    } else {
		$species = "sp.";
	    }
	    $sub_species = $3 if $3;
            $common      = $4 if $4;
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
    
    # Don't make a species object if it is "Unknown" or "None"
    return if $genus =~ /^(Unknown|None)$/i;

    # Bio::Species array needs array in Species -> Kingdom direction
    if ($class[$#class] eq $genus) {
        push( @class, $species );
    } else {
        push( @class, $genus, $species );
    }
    @class = reverse @class;
    
    my $make = Bio::Species->new();
    $make->classification( @class );
    $make->common_name( $common      ) if $common;
    $make->sub_species( $sub_species ) if $sub_species;
    $make->organelle  ( $org         ) if $org;
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

# inherited from SeqIO.pm ! HL 05/11/2000

=head2 _read_FTHelper_swissprot

 Title   : _read_FTHelper_swissprot
 Usage   : _read_FTHelper_swissprot(\$buffer)
 Function: reads the next FT key line
 Example :
 Returns : Bio::SeqIO::FTHelper object 
 Args    : filehandle and reference to a scalar


=cut

sub _read_FTHelper_swissprot {
    # initial version implemented by HL 05/10/2000
    # FIXME this may not be perfect, so please review 
    my ($self,$buffer) = @_;
    my ($key,   # The key of the feature
        $loc,   # The location line from the feature
        $desc,  # The descriptive text
        );
    
    if ($$buffer =~ /^FT   (\w+)\s+([\d\?\<]+)\s+([\d\?\>]+)\s*(.*)$/) {
        $key = $1;
        my $loc1 = $2;
        my $loc2 = $3;
	$loc = "$loc1..$loc2";
	if($4 && (length($4) > 0)) {
	    $desc = 4;
	    chomp($desc);
	} else {
	    $desc = "";
	}
	# Read all the continuation lines up to the next feature
	while (defined($_ = $self->_readline) && /^FT\s{20,}(\S.*)$/) {
	    $desc .= $1;
	    chomp($desc);
	}
	$desc =~ s/\.$//;
    } else {
        # No feature key. What's this?
	$self->warn("No feature key in putative feature table line: $_");
        return;
    } 
    
    # Put the first line of the next feature into the buffer
    $$buffer = $_;

    # Make the new FTHelper object
    my $out = new Bio::SeqIO::FTHelper();
    $out->key($key);
    $out->loc($loc);
    
    # store the description if there is one
    if($desc && (length($desc) > 0)) {
	$out->field->{"description"} ||= [];
	push(@{$out->field->{"description"}}, $desc);
    }
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
   my ($self,$pre1,$pre2,$line,$length) = @_;

   $length || die "Miscalled write_line_swissprot without length. Programming error!";
   my $subl = $length - length $pre2;
   my $linel = length $line;
   my $i;

   my $sub = substr($line,0,$length - length $pre1);

   $self->_print( "$pre1$sub\n");
   
   for($i= ($length - length $pre1);$i < $linel;) {
       $sub = substr($line,$i,($subl));
       $self->_print( "$pre2$sub\n");
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
   my ($self,$pre1,$pre2,$line,$regex,$length) = @_;

   
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
   $self->_print( "$pre1$s\n");
   foreach my $s ( @lines ) {
       $self->_print( "$pre2$s\n");
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
    





