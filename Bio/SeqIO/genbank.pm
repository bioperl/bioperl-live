

# BioPerl module for Bio::SeqIO::GenBank
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::GenBank - GenBank sequence input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the SeqIO handler system. Go:

    $stream = Bio::SeqIO->new(-file => $filename, -format => 'GenBank');

    while ( my $seq = $stream->next_seq() ) {
	# do something with $seq
    }

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from GenBank flat
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


package Bio::SeqIO::genbank;
use vars qw(@ISA);
use strict;

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::FTHelper;
use Bio::SeqFeature::Generic;
use Bio::Species;


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


=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    :


=cut

sub next_seq{
    my ($self,@args) = @_;
    my ($pseq,$c,$line,$name,$desc,$acc,$seqc,$mol,$div,$date);
    my $seq = Bio::Seq->new();

    $line = $self->_readline;
    
    if( !defined $line) {
	return undef; # no throws - end of file
    }
    
    if( $line =~ /^\s+$/ ) {
	while( defined($line = $self->_readline) ) {
	    $line =~ /\S/ && last;
	}
    }
    if( !defined $line ) {
	return undef; # end of file
    }
    
    $line =~ /^LOCUS\s+\S+/ || $self->throw("GenBank stream with no LOCUS. Not GenBank in my book. Got $line");
    $line =~ /^LOCUS\s+(\S+)\s+\S+\s+bp\s+(\S+)\s+(\S+)\s+(\S+)/;

    $name = $1;

    $mol=$2; 
    if ($mol) {
	$seq->molecule($mol);
    }
    
    $div=$3;
    if ($div) {
	$seq->division($div);
    }
    
    $date=$4;
    if ($date) {
	$seq->add_date($date);
    }

    my $buffer = $line;
        
    BEFORE_FEATURE_TABLE :
	until( !defined ( $buffer )  ) {
	    $_ = $buffer;
	    
	    # Exit at start of Feature table
	    last if /^FEATURES/;
	    
	    # Description line(s)
	    if (/^DEFINITION\s+(\S.*\S)/) {
		$desc .= $desc ? " $1" : $1;
		$desc .= " ";
		while ( defined($_ = $self->_readline) ) { 
		    /^\s+(.*)/ && do { $desc .= $1; next;};
		    last;
		}
	    }		 
	    
	    if( /^ACCESSION\s+(\S+)/ ) {
		$acc = $1;
		$seq->accession($acc);
	    }
	    
	    #Version number
	    if( /^VERSION\s+(\S+)/ ) {
		my $sv = $1;
		$seq->sv($sv);
	    }

	    #Keywords
	    if( /^KEYWORDS\s+(.*)/ ) {
		my $keywords = $1;
		$keywords =~ s/\;//g;
		$seq->keywords($keywords);
	    }

	    # Organism name and phylogenetic information
	    if (/^SOURCE/) {
		my $species = $self->_read_GenBank_Species(\$buffer);
		$seq->species( $species );
		next;
	    }

	    #References
	    if (/^REFERENCE/) {
		my @refs = $self->_read_GenBank_References(\$buffer);
		$seq->annotation->add_Reference(@refs);
		next;
	    }
	    
	    #Comments
	    if (/^COMMENT\s+(.*)/) {
		my $comment = $1;
		while (defined($_ = $self->_readline)) {
		    if (/^FEATURES/) {
			$comment =~ s/\n/ /g;
			$comment =~ s/  +/ /g;
			my $commobj = Bio::Annotation::Comment->new();
			$commobj->text($comment);
			$seq->annotation->add_Comment($commobj);
			last BEFORE_FEATURE_TABLE;
		    }
		    $comment .= $_; 
		}
	    }
	    
	    # Get next line.
	    $buffer = $self->_readline;
	}

    # need to read the first line of the feature table
    FEATURE_TABLE :   
	while( defined ( $buffer = $self->_readline) ) {	    
	    my $ftunit = $self->_read_FTHelper_GenBank(\$buffer);
	    # process ftunit
	    $ftunit->_generic_seqfeature($seq);
	    last if /^BASE/;
	}
    $seqc = "";	
    while (defined( $_ = $self->_readline)) {
	last if /^ORIGIN/;
    }

    while( defined($_ = $self->_readline) ) {
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
    
    my ($div, $mol);
    my $len = $seq->length();
    

    if ( $seq->can('division') ) {
	$div=$seq->division;
    } 
    if( !defined $div || ! $div ) { $div = 'UNK'; }

    if( !$seq->can('molecule') || ! defined ($mol = $seq->molecule()) ) {
	$mol = 'DNA';
    }
    else {
	$mol = $seq->molecule;
    }
    
    local($^W) = 0;   # supressing warnings about uninitialized fields.
    
    my $temp_line;
    if( $self->_id_generation_func ) {
	$temp_line = &{$self->_id_generation_func}($seq);
    } else {
	my @dates = $seq->each_date();
	my $date = shift @dates;
	$temp_line = sprintf ("%-12s%-10s%10s bp%8s%5s %3s ", 'LOCUS',$seq->id(),$len,$mol,$div,$date);
    } 
    
    $self->_print("$temp_line\n");   
    $self->_write_line_GenBank_regex("DEFINITION  ","            ",$seq->desc(),"\\s\+\|\$",80);
    
    # if there, write the accession line

    if( $self->_ac_generation_func ) {
	$temp_line = &{$self->_ac_generation_func}($seq);
	$self->_print("ACCESSION   $temp_line\n");   
    } else {
	if( $seq->can('accession') ) {
	    $self->_print("ACCESSION   ",$seq->accession,"\n");
	}
	# otherwise - cannot print <sigh>
    } 
    
    # if there, write the version line
    
    if( defined $self->_sv_generation_func() ) {
	$temp_line = &{$self->_sv_generation_func}($seq);
	if( $temp_line ) {
	    $self->_print("VERSION     $temp_line\n");   
	}
    } else {
	if( $seq->can('sv') ) {
	    $self->_print("VERSION     ",$seq->sv,"\n");
       }
    } 
    
    # if there, write the keywords line
    
    if( defined $self->_kw_generation_func() ) {
	$temp_line = &{$self->_kw_generation_func}($seq);
	$self->_print("KEYWORDS    $temp_line\n");   
    } else {
	if( $seq->can('keywords') ) {
	    $self->_print("KEYWORDS    ",$seq->keywords,"\n");
	}
    } 
    
    
    # Organism lines
    if (my $spec = $seq->species) {
        my($species, $genus, @class) = $spec->classification();
	my $sub_species = $spec->sub_species();
	if( !defined $sub_species ) { $sub_species = ""; }

        my $OS = "$genus $species $sub_species";
	$self->_print("SOURCE      $OS\n");
	$self->_print("  ORGANISM  $OS\n");
        my $OC = join (';', reverse(@class));
	$OC =~ s/\n//g;
	$OC =~ s/\;/\; /g;
	$OC = "$OC; $genus.";
        $self->_write_line_GenBank_regex("            ","            ",
					 $OC,"\\s\+\|\$",80);
    }
    
    # Reference lines
    my $count = 1;
    foreach my $ref ( $seq->annotation->each_Reference() ) {
	$temp_line = sprintf ("REFERENCE    $count  (bases %d to %d)",$ref->start,$ref->end);
	$self->_print("$temp_line\n");
	$self->_write_line_GenBank_regex("  AUTHORS   ","            ",
					 $ref->authors,"\\s\+\|\$",80);
	$self->_write_line_GenBank_regex("  TITLE     ","            ",
					 $ref->title,"\\s\+\|\$",80);
	$self->_write_line_GenBank_regex("  JOURNAL   ","            ",
					 $ref->location,"\\s\+\|\$",80);
	if ($ref->comment) {
	    $self->_write_line_GenBank_regex("  REMARK    ","            ",
					     $ref->comment,"\\s\+\|\$",80);
	}
	    $count++;
    }
    # Comment lines
    
    foreach my $comment ( $seq->annotation->each_Comment() ) {
	$self->_write_line_GenBank_regex("COMMENT     ","            ",
					 $comment->text,"\\s\+\|\$",80);
    }
    $self->_print("FEATURES             Location/Qualifiers\n");
    
    if( defined $self->_post_sort ) {
	# we need to read things into an array. Process. Sort them. Print 'em

	my $post_sort_func = $self->_post_sort();
	my @fth;

	foreach my $sf ( $seq->top_SeqFeatures ) {
	    push(@fth,Bio::SeqIO::FTHelper::from_SeqFeature($sf,$seq));
	}

	@fth = sort { &$post_sort_func($a,$b) } @fth;

	foreach my $fth ( @fth ) {
	    $self->_print_GenBank_FTHelper($fth);
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
		
		$self->_print_GenBank_FTHelper($fth);
	    }
	}
    }
    
    if( $self->_show_dna() == 0 ) {
       return;
   }
    
# finished printing features.
    
    $str =~ tr/A-Z/a-z/;

# Count each nucleotide
    my $alen = $str =~ tr/a/a/;
    my $clen = $str =~ tr/c/c/;
    my $glen = $str =~ tr/g/g/;
    my $tlen = $str =~ tr/t/t/;
    
    my $olen = $len - ($alen + $tlen + $clen + $glen);
    if( $olen < 0 ) {
	$self->warn("Weird. More atgc than bases. Problem!");
    }
    $self->_print("BASE COUNT %8s a %6s c %6s g %6s t\n",$alen,$clen,$glen,$tlen); 
    $self->_print("ORIGIN\n");
    my $di;
    for ($i = 0; $i < length($str); $i += 10) {
	
	$di=$i+11;

	#first line
	if ($i==0) {
	    $self->_print(sprintf("%9d ",1));
	}
	#print sequence, spaced by 10
	$self->_print(substr($str,$i,10), " ");
	
	#break line and print number at beginning of next line
	if(($i+10)%60 == 0) {
	    $self->_print("\n");
	    $self->_print(sprintf("%9d ",$di));
	}
    }
    
    
    $self->_print("\n//\n");
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
   my ($self,$fth,$always_quote) = @_;
   
   if( ! ref $fth || ! $fth->isa('Bio::SeqIO::FTHelper') ) {
       $fth->warn("$fth is not a FTHelper class. Attempting to print, but there could be tears!");
   }
   $self->_write_line_GenBank_regex(sprintf("     %-16s",$fth->key),
				    "                     ",
				    $fth->loc,"\,\|\$",80);

   if( !defined $always_quote) { $always_quote = 0; }

   foreach my $tag ( keys %{$fth->field} ) {
       foreach my $value ( @{$fth->field->{$tag}} ) {
	   $value =~ s/\"/\"\"/g;
	   if ($value eq "_no_value") {
	       $self->_write_line_GenBank_regex("                     ",
						"                     ",
						"/$tag","\.\|\$",80);
	   }
           elsif( $always_quote == 1 || $value !~ /^\d+$/ ) {
	      $self->_write_line_GenBank_regex("                     ",
					       "                     ",
					       "/$tag=\"$value\"","\.\|\$",80);
           } else {
              $self->_write_line_GenBank_regex("                     ",
					       "                     ",
					       "/$tag=$value","\.\|\$",80);
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
   my ($self,$buffer) = @_;
   my (@refs);
   my $previous;
   
   # assumme things are starting with RN

   if( $$buffer !~ /^REFERENCE/ ) {
       warn("Not parsing line '$$buffer' which maybe important");
   }
   my $title;
   my $loc;
   my $au;
   my $b1;
   my $b2;
   my $com;

   if ($$buffer =~ /^REFERENCE\s+\d+\s+\(bases (\d+) to (\d+)/){
       $b1=$1;
       $b2=$2;
   }
   while( defined($_ = $self->_readline) ) {
       if (/^  AUTHORS\s+(.*)/) { 
	   $au .= $1;   
	   while ( defined($_ = $self->_readline) ) {
	       /^  TITLE/ && last; 
	       /^\s+(.*)/ && do { $au .= $1; $au =~ s/\,(\S)/ $1/g;$au .= " ";next;};
	   }    
       }
       if (/^  TITLE\s+(.*)/)  { 
	   $title .= $1;
	   while ( defined($_ = $self->_readline) ) {
	       /^  JOURNAL/ && last; 
	       /^\s+(.*)/ && do { $title .= $1;$title .= " ";next;};
	   }
       }
       if (/^  JOURNAL\s+(.*)/) { 
	   $loc .= $1;
	   while ( defined($_ = $self->_readline) ) {
	       /^  REMARK/ && last;
	       /^\s+(.*)/ && do { $loc .= $1;$loc .= " ";next;};
	       last;
	   }
       }

       if (/^  REMARK\s+(.*)/) { 
	   $com .= $1;
	   while ( defined($_ = $self->_readline) ) {	       
	       /^\s+(.*)/ && do { $com .= $1;$com .= " ";next;};
	       last;
	   }
       }

       /^REFERENCE/ && do {
	   # put away current reference
	   my $ref = new Bio::Annotation::Reference;
	   $au =~ s/;\s*$//g;
	   $title =~ s/;\s*$//g;
	   $ref->start($b1);
	   $ref->end($b2);
	   $ref->authors($au);
	   $ref->title($title);
	   $ref->location($loc);
	   $ref->comment($com);
	   push(@refs,$ref);
	   $au = "";
	   $title = "";
	   $loc = "";
	   $com = "";

	   # return to main reference loop
	   next;
       };
	   
       /^FEATURES||^COMMENT/ && last;
   }

   # put away last reference

   my $ref = new Bio::Annotation::Reference;
   $au =~ s/;\s*$//g;
   $title =~ s/;\s*$//g;

   $ref->start($b1);
   $ref->end($b2);
   $ref->authors($au);
   $ref->title($title);
   $ref->location($loc);
   $ref->comment($com);
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
    my( $self,$buffer) = @_;
    
    $_ = $$buffer;
    
    my( $sub_species, $species, $genus, $common, @class );
    while (defined( $_ = $self->_readline )) {

        if (/^SOURCE\s+(.*)/) {
	    $common = $1;
	    $common =~ s/\.//;
	}
	if (/^  ORGANISM\s+(\S+)\s+(\S+)\s+?(\S+)?/) {
            $genus = $1 if $1;
	    $genus = "None" unless $1;
	    if ($2) {
		$species = $2;
	    }
	    else {
		$species = "sp.";
	    }
	    $sub_species = $3 if $3;
        }
        elsif ($_ !~ /^REFERENCE/) {
	    $_ =~ s/\s+//;
	    $_ =~ s/^SOURCE\S+//;
            push(@class, split /[\;]+/, $_);
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
	push( @class, $genus, $species, $sub_species );
    }
    else {
	push( @class, $genus, $species, "");
    }
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
 Usage   : _read_FTHelper_GenBank($buffer)
 Function: reads the next FT key line
 Example :
 Returns : Bio::SeqIO::FTHelper object 
 Args    : filehandle and reference to a scalar


=cut

sub _read_FTHelper_GenBank {
   my ($self,$buffer) = @_;
   my ($key,$loc,$out);

   $_ = $$buffer;
   if( $$buffer =~ /^\s+(\S+)\s+(\S+)/ ) {
       $key = $1;
       $loc = $2;
       if ($$buffer !~ /^\s+\//) {
	   while ( defined($_ = $self->_readline) ) {
	       # read location line until feature, qualifier or end of feature table
	       /^\s+\S+\s+\S+/ && last; # trigger for next feature
	       /\s+\// && last; # trigger for qualifier
	       /^BASE/ && last; # trigger for end of feature table
	       /\s+\>?(\S+)/ && do {$loc .= $1;};
	       
	   }
       }
   }

   $out = new Bio::SeqIO::FTHelper();
   $out->key($key);
   $out->loc($loc);

   # Now read in other fields
   # Loop reads $_ when defined (i.e. only in first loop), then read $self->_readline, until eof 
   while ( defined($_ ||= $self->_readline) ) {
       # Exit loop on new primary key or end of features
       (/^     \S+/||/^BASE/) && last;
       
       # Field on one line
       if (/^\s+\/(\S+)=\"(.*)\"/) {
	   my $key = $1;
	   my $value = $2;
	   if(! defined $out->field->{$key} ) {
	       $out->field->{$key} = [];
	   }
	   $value =~ s/\"\"/\"/g;

	   push(@{$out->field->{$key}},$value);
       }
       # Field on on multilines:
       elsif (/^\s+\/(\S+)=\"(.*)/) {
	   my $key = $1;
	   my $value = $2;
	   while ( defined($_ = $self->_readline) ) {
	       # first subs out the ""
	       s/\"\"/__DOUBLE_QUOTE_STRING__/g;
	       /\s+(.*)\"/ && do { $value .= $1; last; };
	       /\s+(.*)/ && do {$value .= $1; };
	   }
	   $value =~ s/__DOUBLE_QUOTE_STRING__/\"/g;

	   if(! defined $out->field->{$key} ) {
	       $out->field->{$key} = [];
	   }
	   push(@{$out->field->{$key}},$value);
       }
       # Field with no quoted value
       elsif (/^\s+\/(\S+)=?(\S+)?/) {
	   my $key = $1;
	   my $value = $2 if $2;
	   $value = "_no_value" unless $2;
	   if(! defined $out->field->{$key} ) {
	       $out->field->{$key} = [];
	   }
	   push(@{$out->field->{$key}},$value);
       }

       # Empty $_ to trigger read from _readline
       undef $_;
   }   
   $$buffer = $_;
   return $out;
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
   my ($self,$pre1,$pre2,$line,$length) = @_;

   $length || die "Miscalled write_line_GenBank without length. Programming error!";
   my $subl = $length - length $pre2;
   my $linel = length $line;
   my $i;

   my $sub = substr($line,0,$length - length $pre1);

   $self->_print("$pre1$sub\n");
   
   for($i= ($length - length $pre1);$i < $linel;) {
       $sub = substr($line,$i,($subl));
       $self->_print("$pre2$sub\n");
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
   my ($self,$pre1,$pre2,$line,$regex,$length) = @_;
   
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
   $self->_print("$pre1$s\n");
   foreach my $s ( @lines ) {
       $self->_print("$pre2$s\n");
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
    





