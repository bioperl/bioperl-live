

# BioPerl module for Bio::SeqIO::GenBank
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::AnnSeqIO::GenBank - GenBank sequence input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the AnnSeqIO handler system. Go:

    $stream = Bio::AnnSeqIO->new(-file => $filename, -format => 'GenBank');

    while ( my $annseq = $stream->next_annseq() ) {
	# do something with $annseq
    }

=head1 DESCRIPTION

This object can transform Bio::AnnSeq objects to and from GenBank flat
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

   print "ID   ", $func($annseq), "\n";

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


package Bio::AnnSeqIO::GenBank;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::AnnSeq;
use Bio::AnnSeqIO::FTHelper;
use Bio::SeqFeature::Generic;
use Bio::Species;
use Bio::AnnSeqIO::StreamI;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use FileHandle;

@ISA = qw(Bio::Root::Object Bio::AnnSeqIO::StreamI);
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
	$self->throw("Neither a file (-file) nor a filehandle (-fh) provided to GenBank opening");
    }
    
    if( $file ) {
	$fh = new FileHandle;
	$fh->open($file) || $self->throw("Could not open $file for GenBank stream reading $!");
    }
    
    # hash for functions for decoding keys.
    $self->{'_func_ftunit_hash'} = {}; 
    $self->_filehandle($fh);
    $self->_show_dna(1); # sets this to one by default. People can change it
    
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
    my ($seq,$fh,$c,$line,$name,$desc,$acc,$seqc);
    
    $fh = $self->_filehandle();
    $self->throw("This function has not been implemented yet!");
    if( eof $fh ) {
	return undef; # no throws - end of file
    }
    
    $line = <$fh>;
    $line =~ /^ID\s+(\S+)/ || $self->throw("GenBank stream with no ID. Not GenBank in my book");
    $name = $1;
    
#   $self->warn("not parsing upper annotation in GenBank file yet!");
	
    my $buffer = $line;
    
    my $annseq = Bio::AnnSeq->new();
    
    BEFORE_FEATURE_TABLE :
	until( eof($fh) ) {
	    $_ = $buffer;
	    

	    # Exit at start of Feature table
	    last if /^FEATURES/;
	    
	    # Description line(s)
	    if (/^DEFINITION\s+(\S.*\S)/) {
		$desc .= $desc ? " $1" : $1;
	    }
	    if( /^ACCESSION\s+(\S+);?/ ) {
		$acc = $1;
		$annseq->accession($acc);
	    }
	    
	    if( /^VERSION\s+(\S+);?/ ) {
		my $sv = $1;
		$annseq->sv($sv);
	    }
	    
	    if( /^KEYWORDS    (.*)\S*$/ ) {
		my $keywords = $1;
		$annseq->keywords($keywords);
	    }
	    
	    # accession numbers...
       
	    # References
	    elsif (/^R/) {
		my @refs = &_read_GenBank_References(\$buffer,$fh);
		$annseq->annotation->add_Reference(@refs);
	    }
	    
	    # Organism name and phylogenetic information
	    elsif (/^O[SC]/) {
		my $species = _read_GenBank_Species(\$buffer, $fh);
		$annseq->species( $species );
	    }
	    
	    # Get next line.
	    $buffer = <$fh>;
	}

    while( <$fh> ) {
	/FT   \w/ && last;
    }
    $buffer = $_;
    
    FEATURE_TABLE :   
	until( eof($fh) ) {
	    my $ftunit = &_read_FTHelper_GenBank($fh,\$buffer);
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

sub write_annseq {
    my ($self,$annseq) = @_;
    
    if( !defined $annseq ) {
	$self->throw("Attempting to write with no annseq!");
    }
    
    if( ! ref $annseq || ! $annseq->isa('Bio::AnnSeqI') ) {
	$self->warn(" $annseq is not a AnnSeqI compliant module. Attempting to dump, but may fail!");
    }
    
    my $fh = $self->_filehandle();
    my $seq = $annseq->seq();
    my $i;
    my $str = $seq->seq;
    
    my $div;
    my $len = $seq->seq_len();
    
    if( !$annseq->can('division') || ($div = $annseq->division()) == undef ) {
	$div = 'UNK';
    }
    
    my $temp_line;
    if( $self->_id_generation_func ) {
	$temp_line = &{$self->_id_generation_func}($annseq);
    } else {
	#Note: no date field at the end of header
	$temp_line = sprintf ("%-12s%-10s%10s%8s%15s\n", 'LOCUS',$seq->id(),$len,'DNA',$div);
    } 
    
    print $fh "$temp_line\n";   
    
# this next line screws up perl mode parsing. Sorry. It is a pain!
    _write_line_GenBank_regex($fh,"DEFINITION  ","            ",$seq->desc(),"\\s\+\|\$",80);
    
    # if there, write the accession line
    
    if( $self->_ac_generation_func ) {
	$temp_line = &{$self->_ac_generation_func}($annseq);
	print $fh "ACCESSION    $temp_line\n";   
    } else {
	# nothing at the moment
    } 
    
    # if there, write the version line
    
    if( $self->_sv_generation_func ) {
	$temp_line = &{$self->_sv_generation_func}($annseq);
	if( $temp_line ) {
	    print $fh "VERSION     $temp_line\n";   
	}
    } else {
	# nothing at the moment
    } 
    
    # if there, write the keywords line
    
    if( $self->_kw_generation_func ) {
	$temp_line = &{$self->_kw_generation_func}($annseq);
	print $fh "KEYWORDS    $temp_line\n";   
    } else {
	# nothing at the moment
    } 
    
    
    # Organism lines
    if (my $spec = $annseq->species) {
        my($species, $genus, @class) = $spec->classification();
        my $OS = "$genus $species";
        if (my $common = $spec->common_name) {
	    print $fh "SOURCE      $common\n";
        }
	else {
	    print $fh "SOURCE      $OS\n";
	}
	print $fh "  ORGANISM  $OS\n";
        my $OC = join('; ', reverse(@class)). '.';
        _write_line_GenBank_regex($fh,"            ","            ",$OC,"\; \|\$",80);
    }
    
    # Reference lines
    my $t = 1;
    foreach my $ref ( $annseq->annotation->each_Reference() ) {
	print $fh "REFERENCE   $t\n";
	&_write_line_GenBank_regex($fh,"  AUTHORS   ","            ",$ref->authors,"\\s\+\|\$",80);
	&_write_line_GenBank_regex($fh,"  TITLE     ","            ",$ref->title,"\\s\+\|\$",80);
	&_write_line_GenBank_regex($fh,"  JOURNAL   ","            ",$ref->location,"\\s\+\|\$",80);
	$t++;
    }
    # Comment lines
    
    foreach my $comment ( $annseq->annotation->each_Comment() ) {
	_write_line_GenBank_regex($fh,"COMMENT     ","            ",$comment->text(),"\\s\+\|\$",80);
    }
    print $fh "FEATURES             Location/Qualifiers\n";
    
    if( defined $self->_post_sort ) {
	# we need to read things into an array. Process. Sort them. Print 'em
	
	my $post_sort_func = $self->_post_sort();
	my @fth;
	
	foreach my $sf ( $annseq->top_SeqFeatures ) {
	    push(@fth,Bio::AnnSeqIO::FTHelper::from_SeqFeature($sf,$annseq));
	}
	
	@fth = sort { &$post_sort_func($a,$b) } @fth;
	
	foreach my $fth ( @fth ) {
	    &_print_GenBank_FTHelper($fth,$fh);
	}
    } else {
	# not post sorted. And so we can print as we get them.
	# lower memory load...
	
	foreach my $sf ( $annseq->top_SeqFeatures ) {
	    my @fth = Bio::AnnSeqIO::FTHelper::from_SeqFeature($sf,$annseq);
	    foreach my $fth ( @fth ) {
		if( ! $fth->isa('Bio::AnnSeqIO::FTHelper') ) {
		    $sf->throw("Cannot process FTHelper... $fth");
		}
		
		&_print_GenBank_FTHelper($fth,$fh);
	    }
	}
    }
    
    if( $self->_show_dna() == 0 ) {
       return;
   }
    
# finished printing features.
    
    $str =~ tr/A-Z/a-z/;
    my $a = $str; 
    $a =~ s/[^a]//g;
    my $alen = length $a;
    
    my $t = $str; 
    $t =~ s/[^t]//g;
    my $tlen = length $t;
    
    my $g = $str; 
    $g =~ s/[^g]//g;
    my $glen = length $g;
    
    my $c = $str; 
    $c =~ s/[^c]//g;
    my $clen = length $c;
    
    my $olen = $len - ($alen + $tlen + $clen + $glen);
    if( $olen < 0 ) {
	$self->warn("Weird. More atgc than bases. Problem!");
    }
    printf $fh ("BASE COUNT %8s a %6s c %6s g %6s t\n",$alen,$clen,$glen,$tlen); 
    printf $fh "ORIGIN\n";
    my $di;
    for ($i = 0; $i < length($str); $i += 10) {
	
	$di=$i+11;

	#first line
	if ($i==0) {
	    print $fh sprintf("%9d ",1);
	}
	#print sequence, spaced by 10
	print $fh substr($str,$i,10), " ";
	
	#break line and print number at beginning of next line
	if(($i+10)%60 == 0) {
	    print $fh "\n";
	    print $fh sprintf("%9d ",$di);
	}
    }
    
    
    print $fh "//\n";
    return 1;
}

=head2 _print_GenBank_FTHelper

 Title   : _print_GenBank_FTHelper
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _print_GenBank_FTHelper {
   my ($fth,$fh,$always_quote) = @_;
   
   if( ! ref $fth || ! $fth->isa('Bio::AnnSeqIO::FTHelper') ) {
       $fth->warn("$fth is not a FTHelper class. Attempting to print, but there could be tears!");
   }
    &_write_line_GenBank_regex($fh,sprintf("     %-16s",$fth->key),"                     ",$fth->loc,"\,\|\$",80);

   foreach my $tag ( keys %{$fth->field} ) {
       foreach my $value ( @{$fth->field->{$tag}} ) {
           if( $always_quote == 1 || $value !~ /^\d+$/ ) {
	      &_write_line_GenBank_regex($fh,"                     ","                     ","/$tag=\"$value\"","\.\|\$",80);
           } else {
              &_write_line_GenBank_regex($fh,"                     ","                     ","/$tag=$value","\.\|\$",80);
           }
       }
   }

}


=head2 _read_GenBank_References

 Title   : _read_GenBank_References
 Usage   :
 Function: Reads references from GenBank format. Internal function really
 Example :
 Returns : 
 Args    :


=cut

sub _read_GenBank_References{
   my ($buffer,$fh) = @_;
   my (@refs);
   
   # assumme things are starting with RN

   if( $$buffer !~ /^RN/ ) {
       warn("Not parsing line '$$buffer' which maybe important");
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


=head2 _read_GenBank_Species

 Title   : _read_GenBank_Species
 Usage   :
 Function: Reads the GenBank Organism species and classification
           lines.
 Example :
 Returns : A Bio::Species object
 Args    :

=cut

sub _read_GenBank_Species {
    my( $buffer, $fh ) = @_;
    
    $_ = $$buffer;
    
    my( $species, $genus, $common, @class );
    while (defined( $_ ||= <$fh> )) {
        
        if (/^OS\s+(\S+)\s+(\S+)(?:\s+\((\s+)\))?/) {
            $genus   = $1;
            $species = $2;
            $common  = $3 if $3;
        }
        elsif (s/^OC\s+//) {
            push(@class, split /[\s\.;]+/, $_);
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
    push( @class, $genus, $species );
    @class = reverse @class;
    
    my $make = Bio::Species->new();
    $make->classification( @class );
    $make->common_name( $common ) if $common;
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

=head2 _read_FTHelper_GenBank

 Title   : _read_FTHelper_GenBank
 Usage   : &_read_FTHelper_GenBank($fh,$buffer)
 Function: reads the next FT key line
 Example :
 Returns : Bio::AnnSeqIO::FTHelper object 
 Args    : filehandle and reference to a scalar


=cut

sub _read_FTHelper_GenBank {
   my ($fh,$buffer) = @_;
   my ($key,$loc,$out);

   if( $$buffer =~ /^FT\s+(\S+)\s+(\S+)/ ) {
       $key = $1;
       $loc = $2;
   }

   $out = new Bio::AnnSeqIO::FTHelper();

   # Read the rest of the location
   while( <$fh> ) {
       # Exit loop on first qualifier - line is left in $_
       last if /^XX/; # sometimes a trailing XX comment left in!
       last if /^FT\s+\//;
       # Get location line
       /^FT\s+(\S+)/ or $out->throw("Weird location line in GenBank feature table: '$_'");
       $loc .= $1;
   }

   $out->key($key);
   $out->loc($loc);

   # Now read in other fields
   while( defined($_ ||= <$fh>) ) {
      
       # Exit loop on non FT lines!
       /^FT/  || last;
       
       # Exit loop on new primary key
       /^FT   \w/ && last;
       
       # Field on one line
       if (/^FT\s+\/(\S+)=\"(.*)\"/) {
	   my $key = $1;
	   my $value = $2;
	   if(! defined $out->field->{$key} ) {
	       $out->field->{$key} = [];
	   }
	   push(@{$out->field->{$key}},$value);
       }
       # Field on on multilines:
       elsif (/^FT\s+\/(\S+)=\"(.*)/) {
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
       }
       # Field with no quoted value
       elsif (/^FT\s+\/(\S+)=(\S+)/) {
	   my $key = $1;
	   my $value = $2;
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

=head2 _generic_seqfeature

 Title   : _generic_seqfeature
 Usage   : &_generic_seqfeature($annseq,$fthelper)
 Function: processes fthelper into a generic seqfeature
 Returns : nothing (places a new seqfeature into annseq)
 Args    : Bio::AnnSeq,Bio::AnnSeqIO::FTHelper


=cut

sub _generic_seqfeature {
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
       $sf->source_tag('GenBank');
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
	   $sub->source_tag('GenBank');
	   $sf->add_sub_SeqFeature($sub,'EXPAND');
       }

   } else {
       my $lst;
       my $len;

       if( $fth->loc =~ /^(\d+)$/ ) {
	   $lst = $len = $1;
       } else {
	   $fth->loc =~ /(\d+)\.\.(\d+)/ || do {
	       $annseq->throw("Weird location line [" . $fth->loc . "] in reading GenBank");
	       last;
	   };
	   $lst = $1;
	   $len = $2;
       }

       $sf->start($lst);
       $sf->end($len);
       $sf->source_tag('GenBank');
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

=head2 _write_line_GenBank

 Title   : _write_line_GenBank
 Usage   :
 Function: internal function
 Example :
 Returns : 
 Args    :


=cut

sub _write_line_GenBank{
   my ($fh,$pre1,$pre2,$line,$length) = @_;

   $length || die "Miscalled write_line_GenBank without length. Programming error!";
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

=head2 _write_line_GenBank_regex

 Title   : _write_line_GenBank_regex
 Usage   :
 Function: internal function for writing lines of specified
           length, with different first and the next line 
           left hand headers and split at specific points in the
           text
 Example :
 Returns : nothing
 Args    : file handle, first header, second header, text-line, regex for line breaks, total line length


=cut

sub _write_line_GenBank_regex {
   my ($fh,$pre1,$pre2,$line,$regex,$length) = @_;

   
   #print STDOUT "Going to print with $line!\n";

   $length || die "Miscalled write_line_GenBank without length. Programming error!";

   if( length $pre1 != length $pre2 ) {
       die "Programming error - cannot called write_line_GenBank_regex with different length pre1 and pre2 tags!";
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
    





