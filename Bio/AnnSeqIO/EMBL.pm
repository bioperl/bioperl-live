

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

    while ( my $annseq = $stream->next_annseq() ) {
	# do something with $annseq
    }

=head1 DESCRIPTION

This object can transform Bio::AnnSeq objects to and from EMBL flat
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
      $self->throw("Neither a file (-file) nor a filehandle (-fh) provided to EMBL opening");
  }

  if( $file ) {
      $fh = new FileHandle;
      $fh->open($file) || $self->throw("Could not open $file for EMBL stream reading $!");
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

sub next_annseq {
   my ($self,@args) = @_;
   my ($seq,$fh,$c,$line,$name,$desc,$acc,$seqc,$mol,$div, $date, $comment, @date_arr);
   my $annseq = Bio::AnnSeq->new();

   $fh = $self->_filehandle();

   $line = <$fh>;   # This needs to be before the first eof() test

   if( eof $fh ) {
       return undef; # no throws - end of file
   }

   if( $line =~ /^\s+$/ ) {
       while( <$fh> ) {
	   /\S/ && last;
       }
       $line = $_;
   }
   
   if( !defined $line ) {
       return undef; # end of file
   }
   $line =~ /^ID\s+\S+/ || $self->throw("EMBL stream with no ID. Not embl in my book");
   $line =~ /^ID\s+(\S+)\s+\S+\; (.+)\; (\S+)/;
   $name = $1;
   
   $mol= $2;
   $mol =~ s/\;//;
   if ($mol) {
       $annseq->molecule($mol);
   }
   
   if ($3) {
       $div = $3;
       $div =~ s/\;//;
       $annseq->division($div);
   }
   
#   $self->warn("not parsing upper annotation in EMBL file yet!");
   my $buffer = $line;

 
   
   BEFORE_FEATURE_TABLE :
   until( eof($fh) ) {
       $_ = $buffer;


       # Exit at start of Feature table
       last if /^FH/;

       # Description line(s)
       if (/^DE\s+(\S.*\S)/) {
           $desc .= $desc ? " $1" : $1;
       }

       #accession number
       if( /^AC\s+(\S+);?/ ) {
	   $acc = $1;
	   $acc =~ s/\;//;
	   $annseq->accession($acc);
       }
       
       #version number
       if( /^SV\s+(\S+);?/ ) {
	   my $sv = $1;
	   $sv =~ s/\;//;
	   $annseq->sv($sv);
       }

       #date (NOTE: takes last date line)
       if( /^DT\s+(\S+)/ ) {
	   my $date = $1;
	   $date =~ s/\;//;
	   $annseq->add_date($date);
       }
       
       #keywords
       if( /^KW   (.*)\S*$/ ) {
	   my $keywords = $1;
	   $annseq->keywords($keywords);
       }

       # Organism name and phylogenetic information
       elsif (/^O[SC]/) {
           my $species = _read_EMBL_Species(\$buffer, $fh);
           $annseq->species( $species );
       }

       # References
       elsif (/^R/) {
	   my @refs = &_read_EMBL_References(\$buffer,$fh);
	   $annseq->annotation->add_Reference(@refs);
       }
       
       # DB Xrefs
       elsif (/^DR/) {
	   my @links = &_read_EMBL_DBLink(\$buffer,$fh);
	   $annseq->annotation->add_DBLink(@links);
       }
       
       # Comments
       elsif (/^CC\s+(.*)/) {
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
	   $commobj->text($comment);
	   $annseq->annotation->add_Comment($commobj);
	   $comment = "";
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
       my $ftunit = &_read_FTHelper_EMBL($fh,\$buffer);
       # process ftunit
       $ftunit->_generic_seqfeature($annseq);

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

    my $mol;
    my $div;
    my $len = $seq->seq_len();

    if ($annseq->can('division') && defined $annseq->division) {
        $div = $annseq->division();
    }
    $div ||= 'UNK';
    
    if ($annseq->can('molecule')) {
        $mol = $annseq->molecule();
    }
    $mol ||= 'XXX';
   
    my $temp_line;
    if( $self->_id_generation_func ) {
        $temp_line = &{$self->_id_generation_func}($annseq);
    } else {
        $temp_line = sprintf("%-11sstandard; $mol; $div; %d BP.", $seq->id(), $len);
    } 

    print $fh "ID   $temp_line\n",
              "XX\n";

    # Write the accession line if present
    {
        my( $acc );
        if( my $func = $self->_ac_generation_func ) {
            $acc = &{$func}($annseq);
        } elsif( $annseq->can('accession')) {
            $acc = $annseq->accession;
        }
        if (defined $acc) {
            print $fh "AC   $acc\n",
                      "XX\n";
        }
    }

    # Write the sv line if present
    {
        my( $sv );
        if (my $func = $self->_sv_generation_func) {
            $sv = &{$func}($annseq);
        } elsif( $annseq->can('sv')) {
            $sv = $annseq->sv;
        }
        if (defined $sv) {
            print $fh "SV   $sv\n",
                      "XX\n";
        }
    }

    # Date lines
    my $switch=0;
    foreach my $dt ( $annseq->each_date() ) {
        _write_line_EMBL_regex($fh,"DT   ","DT   ",$dt,'\s+|$',80);
        $switch=1;
    }
    if ($switch == 1) {
        print $fh "XX\n";
    }

    # Description lines
    _write_line_EMBL_regex($fh,"DE   ","DE   ",$seq->desc(),'\s+|$',80);
    print $fh "XX\n";

    # if there, write the kw line
    {
        my( $kw );
        if( my $func = $self->_kw_generation_func ) {
            $kw = &{$func}($annseq);
        } elsif( $annseq->can('keywords') ) {
	    $kw = $annseq->keywords;
        }
        if (defined $kw) {
            print $fh "KW   $kw\n",
                      "XX\n";
        }
    }
   
    # Organism lines
    if (my $spec = $annseq->species) {
        my($sub_species, $species, $genus, @class) = $spec->classification();
        my $OS = "$genus $species $sub_species";
        if (my $common = $spec->common_name) {
            $OS .= " ($common)";
        }
        print $fh "OS   $OS\n";
        my $OC = join('; ', reverse(@class));
	$OC =~ s/\;\s+$//;
	$OC .= ".";
        _write_line_EMBL_regex($fh,"OC   ","OC   ",$OC,'; |$',80);
	if ($spec->organelle) {
	    _write_line_EMBL_regex($fh,"OG   ","OG   ",$spec->organelle,'; |$',80);
	}
        print $fh "XX\n";
    }
   
    # Reference lines
    my $t = 1;
    foreach my $ref ( $annseq->annotation->each_Reference() ) {
        print $fh "RN   [$t]\n";
        
        # Having no RP line is legal, but we need both
        # start and end for a valid location.
        my $start = $ref->start;
        my $end   = $ref->end;
        if ($start and $end) {
            print $fh "RP   $start-$end\n";
        } elsif ($start or $end) {
            $self->throw("Both start and end are needed for a valid RP line.  Got: start='$start' end='$end'");
        }

        if (my $med = $ref->medline) {
            print $fh "RX   MEDLINE; $med\n";
        }

        &_write_line_EMBL_regex($fh, "RA   ", "RA   ", $ref->authors,  '\s+|$', 80);       

        # If there is no title to the reference, it appears
        # as a single semi-colon.  All titles must end in
        # a semi-colon.
        my $ref_title = $ref->title || '';
        $ref_title =~ s/[\s;]*$/;/;
        &_write_line_EMBL_regex($fh, "RT   ", "RT   ", $ref_title,    '\s+|$', 80);       

        &_write_line_EMBL_regex($fh, "RL   ", "RL   ", $ref->location, '\s+|$', 80);
        if ($ref->comment) {
	    &_write_line_EMBL_regex($fh, "RC   ", "RC   ", $ref->comment, '\s+|$', 80); 
        }
        print $fh "XX\n";
        $t++;
    }

    # DB Xref lines
    if (my @db_xref = $annseq->annotation->each_DBLink) {
        foreach my $dr (@db_xref) {
            my $db_name = $dr->database;
            my $prim    = $dr->primary_id;
            my $opt     = $dr->optional_id || '';
            
            my $line = "$db_name; $prim; $opt.";
            &_write_line_EMBL_regex($fh, "DR   ", "DR   ", $line, '\s+|$', 80);
        }
        print $fh "XX\n";
    }

    # Comment lines
    foreach my $comment ( $annseq->annotation->each_Comment() ) {
        _write_line_EMBL_regex($fh, "CC   ", "CC   ", $comment->text, '\s+|$', 80);
        print $fh "XX\n";
    }
    # "\\s\+\|\$"

    ## FEATURE TABLE

    print $fh "FH   Key             Location/Qualifiers\n";
    print $fh "FH\n";

    if( defined $self->_post_sort ) {
        # we need to read things into an array. Process. Sort them. Print 'em

        my $post_sort_func = $self->_post_sort();
        my @fth;

        foreach my $sf ( $annseq->top_SeqFeatures ) {
	    push(@fth,Bio::AnnSeqIO::FTHelper::from_SeqFeature($sf,$annseq));
        }

        @fth = sort { &$post_sort_func($a,$b) } @fth;

        foreach my $fth ( @fth ) {
	    &_print_EMBL_FTHelper($fth,$fh);
        }
    } else {
        # not post sorted. And so we can print as we get them.
        # lower memory load...

        foreach my $sf ( $annseq->top_SeqFeatures ) {
	    my @fth = Bio::AnnSeqIO::FTHelper::from_SeqFeature($sf,$annseq);
	    foreach my $fth ( @fth ) {
	        &_print_EMBL_FTHelper($fth,$fh);
	    }
        }
    }

    print $fh "XX\n";

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

    print $fh "SQ   Sequence $len BP; $alen A; $clen C; $glen G; $tlen T; $olen other;\n";
    print $fh "     ";
    my $linepos;
    for ($i = 0; $i < length($str); $i += 10) {

        if( $i+10 >= length($str) ) {
	    # last line.
	    print $fh substr($str,$i);
	    $linepos += length($str)-$i;
	    print $fh ' ' x (70 - $linepos);
	    print $fh sprintf(" %-5d\n",length($str));
	    last;
        }
        print $fh substr($str,$i,10), " ";
        $linepos += 11;
        if( ($i+10)%60 == 0 ) {
	    my $end = $i+10;
	    print $fh sprintf("%-5d\n     ",$end);
	    $linepos = 5;
        }
    }


    print $fh "//\n";
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

sub _print_EMBL_FTHelper {
   my ($fth,$fh,$always_quote) = @_;
   $always_quote ||= 0;
   
   if( ! ref $fth || ! $fth->isa('Bio::AnnSeqIO::FTHelper') ) {
       $fth->warn("$fth is not a FTHelper class. Attempting to print, but there could be tears!");
   }


   #print $fh "FH   Key             Location/Qualifiers\n";
   #print $fh  sprintf("FT   %-15s  %s\n",$fth->key,$fth->loc);
   &_write_line_EMBL_regex($fh,sprintf("FT   %-15s ",$fth->key),"FT                   ",$fth->loc,'\,|$',80);
   foreach my $tag ( keys %{$fth->field} ) {
       if( ! defined $fth->field->{$tag} ) { next; } 
       foreach my $value ( @{$fth->field->{$tag}} ) {
	   $value =~ s/\"/\"\"/g;
	   if ($value eq "_no_value") {
	       &_write_line_EMBL_regex($fh,"FT                   ","FT                   ","/$tag",'.|$',80);
	   }
           elsif( $always_quote == 1 || $value !~ /^\d+$/ ) {
	      &_write_line_EMBL_regex($fh,"FT                   ","FT                   ","/$tag=\"$value\"",'.|$',80);
           }
           else {
              &_write_line_EMBL_regex($fh,"FT                   ","FT                   ","/$tag=$value",'.|$',80);
           }
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

sub _read_EMBL_References {
   my ($buffer,$fh) = @_;
   my (@refs);
   
   # assumme things are starting with RN

   if( $$buffer !~ /^RN/ ) {
       warn("Not parsing line '$$buffer' which maybe important");
   }
   my $b1;
   my $b2;
   my $title;
   my $loc;
   my $au;
   my $med;
   my $com;

   while( <$fh> ) {
       /^R/ || last;
       /^RP   (\d+)-(\d+)/ && do {$b1=$1;$b2=$2;};
       /^RX   MEDLINE;\s+(\d+)/ && do {$med=$1};
       /^RA   (.*)/ && do { $au .= $1;   next;};
       /^RT   (.*)/ && do { $title .= $1; next;};
       /^RL   (.*)/ && do { $loc .= $1; next;};
       /^RC   (.*)/ && do { $com .= $1; next;};
   }
   
   my $ref = new Bio::Annotation::Reference;
   $au =~ s/;\s*$//g;
   $title =~ s/;\s*$//g;

   $ref->start($b1);
   $ref->end($b2);
   $ref->authors($au);
   $ref->title($title);
   $ref->location($loc);
   $ref->medline($med);
   $ref->comment($com);

   push(@refs,$ref);
   $$buffer = $_;
   
   return @refs;
}


=head2 _read_EMBL_Species

 Title   : _read_EMBL_Species
 Usage   :
 Function: Reads the EMBL Organism species and classification
           lines.
 Example :
 Returns : A Bio::Species object
 Args    :

=cut

sub _read_EMBL_Species {
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

=head2 _read_EMBL_DBLink

 Title   : _read_EMBL_DBLink
 Usage   :
 Function: Reads the EMBL database cross reference ("DR") lines
 Example :
 Returns : A list of Bio::Annotation::DBLink objects
 Args    :

=cut

sub _read_EMBL_DBLink {
    my( $buffer, $fh ) = @_;
    my( @db_link );

    $_ = $$buffer;
    while (defined( $_ ||= <$fh> )) {
        
        if (my($databse, $prim_id, $sec_id)
                = /^DR   ([^\s;]+);\s*([^\s;]+);\s*([^\s;]+)?\.$/) {
            my $link = Bio::Annotation::DBLink->new();
            $link->database   ( $databse );
            $link->primary_id ( $prim_id );
            $link->optional_id( $sec_id  ) if $sec_id;
            push(@db_link, $link);
	}
        else {
            last;
        }
        
        $_ = undef; # Empty $_ to trigger read of next line
    }
    
    $$buffer = $_;
    
    return @db_link;
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

sub _read_FTHelper_EMBL {
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
       last if /^FT\s+(\S+)\s+(\S+)/;
       
       # Get location line
       /^FT\s+(\S+)/ or $out->throw("Weird location line in EMBL feature table: '$_'");
       $loc .= $1;
   }

   $loc =~ s/<//;
   $loc =~ s/>//;
   $out->key($key);
   $out->loc($loc);

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


=head2 _write_line_EMBL

 Title   : _write_line_EMBL
 Usage   :
 Function: internal function
 Example :
 Returns : 
 Args    :


=cut

sub _write_line_EMBL {
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

    my $subl = $length - (length $pre1) -1 ;

    my( @lines );
    while($line =~ m/(.{1,$subl})($regex)/g) {
        push(@lines, $1.$2);
    }
    foreach (@lines) { s/\s+$//; }
    #chomp(@lines);
    
    # Print first line
    my $s = shift(@lines);
    print $fh "$pre1$s\n";
    
    # Print the rest
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
    





