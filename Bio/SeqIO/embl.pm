# $Id$
#
# BioPerl module for Bio::SeqIO::EMBL
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::embl - EMBL sequence input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the AnnSeqIO handler system. Go:

    $stream = Bio::SeqIO->new(-file => $filename, -format => 'EMBL');

    while ( (my $seq = $stream->next_seq()) ) {
	# do something with $seq
    }

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from EMBL flat
file databases.

There is alot of flexibility here about how to dump things which I need
to document fully.

There should be a common object that this and genbank share (probably
with swissprot). Too much of the magic is identical. 

=head2 Optional functions

=over 3

=item _show_dna()

(output only) shows the dna or not

=item _post_sort()

(output only) provides a sorting func which is applied to the FTHelpers
before printing

=item _id_generation_func()

This is function which is called as 

   print "ID   ", $func($annseq), "\n";

To generate the ID line. If it is not there, it generates a sensible ID
line using a number of tools.

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://www.bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqIO::embl;
use vars qw(@ISA);
use strict;
use Bio::Seq::RichSeq;
use Bio::SeqIO::FTHelper;
use Bio::SeqFeature::Generic;
use Bio::Species;

@ISA = qw(Bio::SeqIO);

sub _initialize {
  my($self,@args) = @_;

  $self->SUPER::_initialize(@args);  
  # hash for functions for decoding keys.
  $self->{'_func_ftunit_hash'} = {}; 
  $self->_show_dna(1); # sets this to one by default. People can change it 
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
   my ($pseq,$c,$line,$name,$desc,$acc,$seqc,$mol,$div, 
       $date, $comment, @date_arr);
   my $seq = Bio::Seq::RichSeq->new(-verbose =>$self->verbose());


   $line = $self->_readline;   # This needs to be before the first eof() test

   if( !defined $line ) {
       return undef; # no throws - end of file
   }

   if( $line =~ /^\s+$/ ) {
       while( defined ($line = $self->_readline) ) {
	   $line =~/\S/ && last;
       }
   }   
   if( !defined $line ) {
       return undef; # end of file
   }
   $line =~ /^ID\s+\S+/ || $self->throw("EMBL stream with no ID. Not embl in my book");
   $line =~ /^ID\s+(\S+)\s+\S+\;\s+([^;]+)\;\s+(\S+)\;/;
   $name = $1;
   $mol = $2;
   $div = $3;
   if(! $name) {
       $name = "unknown id";
   }
    # this is important to have the id for display in e.g. FTHelper, otherwise
    # you won't know which entry caused an error
   $seq->display_id($name);
   if($mol) {
       if ( $mol =~ /circular/ ) {
	   $seq->is_circular(1);
	   $mol =~  s|circular ||;
       }
       $seq->molecule($mol);
       my $alphabet;
       if (defined $seq->molecule) {
	   my $mol =$seq->molecule;
	   if ($mol =~ /DNA/) {
	       $alphabet='dna';
	   }
	   elsif ($mol =~ /RNA/) {
	       $alphabet='dna';
	   }
	   elsif ($mol =~ /AA/) {
	       $alphabet='protein';
	   }
       }
       if ($alphabet) {
	   $seq->primary_seq->alphabet($alphabet);
       }

   }
   if ($div) {
       $seq->division($div);
   }
   
#   $self->warn("not parsing upper annotation in EMBL file yet!");
   my $buffer = $line;
 
   
   BEFORE_FEATURE_TABLE :
   until( !defined $buffer ) {
       $_ = $buffer;

       # Exit at start of Feature table
       last if /^F[HT]/;

       # Description line(s)
       if (/^DE\s+(\S.*\S)/) {
           $desc .= $desc ? " $1" : $1;
       }

       #accession number
       if( /^AC\s+(\S+);?/ ) {
	   $acc = $1;
	   $acc =~ s/\;//;
	   $seq->accession_number($acc);
       }
       
       #version number
       if( /^SV\s+\S+\.(\d+);?/ ) {
	   my $sv = $1;
	   $sv =~ s/\;//;
	   $seq->seq_version($sv);
       }

       #date (NOTE: takes last date line)
       if( /^DT\s+(.+)$/ ) {
	   my $date = $1;
	   $date =~ s/\;$//;
	   $seq->add_date($date);
       }
       
       #keywords
       if( /^KW   (.*)\S*$/ ) {
	   my $keywords = $1;
	   $seq->keywords($keywords);
       }

       # Organism name and phylogenetic information
       elsif (/^O[SC]/) {
           my $species = $self->_read_EMBL_Species(\$buffer);
           $seq->species( $species );
       }

       # References
       elsif (/^R/) {
	   my @refs = $self->_read_EMBL_References(\$buffer);
	   foreach my $ref ( @refs ) {
	       $seq->annotation->add_Annotation('reference',$ref);
	   }
       }
       
       # DB Xrefs
       elsif (/^DR/) {
	   my @links = $self->_read_EMBL_DBLink(\$buffer);
	   foreach my $dblink ( @links ) {
	       $seq->annotation->add_Annotation('dblink',$dblink);
	   }
       }
       
       # Comments
       elsif (/^CC\s+(.*)/) {
	   $comment .= $1;
	   $comment .= " ";
	   while (defined ($_ = $self->_readline) ) {
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
	   $seq->annotation->add_Annotation('comment',$commobj);
	   $comment = "";
       }

       # Get next line.
       $buffer = $self->_readline;
   }
    $seq->desc($desc);

   while( defined ($_ = $self->_readline) ) {
       /^FT   \w/ && last;
       /^SQ / && last;
   }
   $buffer = $_;
      
   if (defined($buffer) && $buffer =~ /^FT /) {
     until( !defined ($buffer) ) {
	 my $ftunit = $self->_read_FTHelper_EMBL(\$buffer);
	 # process ftunit
	 $ftunit->_generic_seqfeature($seq);

	 if( $buffer !~ /^FT/ ) {
	     last;
	 }
     }
   }
   
   if( $buffer !~ /^SQ/  ) {
       while( defined ($_ = $self->_readline) ) {
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

   $seq->seq($seqc);
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

    my $str = $seq->seq;

    my $mol;
    my $div;
    my $len = $seq->length();

    if ($seq->can('division') && defined $seq->division) {
        $div = $seq->division();
    }
    $div ||= 'UNK';
    
    if ($seq->can('molecule')) {
        $mol = $seq->molecule();
	$mol = 'RNA' if $mol =~ /RNA/; # no 'mRNA' 
    }
    elsif (defined $seq->primary_seq->alphabet) {
	my $alphabet =$seq->primary_seq->alphabet;
	if ($alphabet eq 'dna') {
	    $mol ='DNA';
	}
	elsif ($alphabet eq 'rna') {
	    $mol='RNA';
	}
	elsif ($alphabet eq 'protein') {
	    $mol='AA';
	}
    }
    $mol ||= 'XXX';
    $mol = "circular $mol" if $seq->is_circular;
   
    my $temp_line;
    if( $self->_id_generation_func ) {
        $temp_line = &{$self->_id_generation_func}($seq);
    } else {
        $temp_line = sprintf("%-11sstandard; $mol; $div; %d BP.", $seq->id(), $len);
    } 

    $self->_print( "ID   $temp_line\n","XX\n");

    # Write the accession line if present
    my( $acc );
    {
        if( my $func = $self->_ac_generation_func ) {
            $acc = &{$func}($seq);
        } elsif( $seq->isa('Bio::Seq::RichSeqI') && 
		 defined($seq->accession_number) ) {
            $acc = $seq->accession_number;
        }

        if (defined $acc) {
            $self->_print("AC   $acc;\n",
			  "XX\n");
        }
    }

    # Write the sv line if present
    {
        my( $sv );
        if (my $func = $self->_sv_generation_func) {
            $sv = &{$func}($seq);
        } elsif($seq->isa('Bio::Seq::RichSeqI') && 
		defined($seq->seq_version)) {
            $sv = "$acc.". $seq->seq_version();
        }	
        if (defined $sv) {
            $self->_print( "SV   $sv\n",
			   "XX\n");
        }
    }

    # Date lines
    my $switch=0;
    if( $seq->can('get_dates') ) {
	foreach my $dt ( $seq->get_dates() ) {
	    $self->_write_line_EMBL_regex("DT   ","DT   ",$dt,'\s+|$',80);#'
            $switch=1;
        }
        if ($switch == 1) {
            $self->_print("XX\n");
        }
    }

    # Description lines
    $self->_write_line_EMBL_regex("DE   ","DE   ",$seq->desc(),'\s+|$',80); #'
    $self->_print( "XX\n");

    # if there, write the kw line
    {
        my( $kw );
        if( my $func = $self->_kw_generation_func ) {
            $kw = &{$func}($seq);
        } elsif( $seq->can('keywords') ) {
	    $kw = $seq->keywords;
        }
        if (defined $kw) {
            $self->_print( "KW   $kw\n",
			   "XX\n");
        }
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
        $self->_print("OS   $OS\n");
        my $OC = join('; ', reverse(@class)) .'.';
        $self->_write_line_EMBL_regex("OC   ","OC   ",$OC,'; |$',80); #'
	if ($spec->organelle) {
	    $self->_write_line_EMBL_regex("OG   ","OG   ",$spec->organelle,'; |$',80); #'
	}
        $self->_print("XX\n");
    }
   
    # Reference lines
    my $t = 1;
    if ( defined $seq->annotation ) {
	foreach my $ref ( $seq->annotation->get_Annotations('reference') ) {
	    $self->_print( "RN   [$t]\n");
	    
	    # Having no RP line is legal, but we need both
	    # start and end for a valid location.
	    my $start = $ref->start;
	    my $end   = $ref->end;
	    if ($start and $end) {
		$self->_print( "RP   $start-$end\n");
	    } elsif ($start or $end) {
		$self->throw("Both start and end are needed for a valid RP line.  Got: start='$start' end='$end'");
	    }
	    
	    if (my $med = $ref->medline) {
		$self->_print( "RX   MEDLINE; $med.\n");
	    }
	    if (my $pm = $ref->pubmed) {
		$self->_print( "RX   PUBMED; $pm.\n");
	    }
	    $self->_write_line_EMBL_regex("RA   ", "RA   ", 
					  $ref->authors . ";",
					  '\s+|$', 80); #'

           # If there is no title to the reference, it appears
           # as a single semi-colon.  All titles must end in
           # a semi-colon.
           my $ref_title = $ref->title || '';
           $ref_title =~ s/[\s;]*$/;/;
           $self->_write_line_EMBL_regex("RT   ", "RT   ", $ref_title,    '\s+|$', 80); #'
	   $self->_write_line_EMBL_regex("RL   ", "RL   ", $ref->location, '\s+|$', 80); #'
           if ($ref->comment) {
	       $self->_write_line_EMBL_regex("RC   ", "RC   ", $ref->comment, '\s+|$', 80); #' 
           }
           $self->_print("XX\n");
           $t++;
       }
        

       # DB Xref lines
       if (my @db_xref = $seq->annotation->get_Annotations('dblink') ) {
           foreach my $dr (@db_xref) {
              my $db_name = $dr->database;
              my $prim    = $dr->primary_id;
              my $opt     = $dr->optional_id || '';
            
              my $line = "$db_name; $prim; $opt.";
              $self->_write_line_EMBL_regex("DR   ", "DR   ", $line, '\s+|$', 80); #'
          }
          $self->_print("XX\n");
       }

       # Comment lines
       foreach my $comment ( $seq->annotation->get_Annotations('comment') ) {
           $self->_write_line_EMBL_regex("CC   ", "CC   ", $comment->text, '\s+|$', 80); #'
           $self->_print("XX\n");
       }
    }
    # "\\s\+\|\$"

    ## FEATURE TABLE

    $self->_print("FH   Key             Location/Qualifiers\n");
    $self->_print("FH\n");

    if( defined $self->_post_sort ) {
        # we need to read things into an array. Process. Sort them. Print 'em

        my $post_sort_func = $self->_post_sort();
        my @fth;

        foreach my $sf ( $seq->top_SeqFeatures ) {
	    push(@fth,Bio::SeqIO::FTHelper::from_SeqFeature($sf,$seq));
        }

        @fth = sort { &$post_sort_func($a,$b) } @fth;

        foreach my $fth ( @fth ) {
	    $self->_print_EMBL_FTHelper($fth);
        }
    } else {
        # not post sorted. And so we can print as we get them.
        # lower memory load...

        foreach my $sf ( $seq->top_SeqFeatures ) {
	    my @fth = Bio::SeqIO::FTHelper::from_SeqFeature($sf,$seq);
	    foreach my $fth ( @fth ) {
	        $self->_print_EMBL_FTHelper($fth);
	    }
        }
    }

    $self->_print( "XX\n");

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

    $self->_print("SQ   Sequence $len BP; $alen A; $clen C; $glen G; $tlen T; $olen other;\n");
    
    my $nuc = 60;               # Number of nucleotides per line
    my $whole_pat = 'a10' x 6;  # Pattern for unpacking a whole line
    my $out_pat   = 'A11' x 6;  # Pattern for packing a line
    my $length = length($str);
    
    # Calculate the number of nucleotides which fit on whole lines
    my $whole = int($length / $nuc) * $nuc;

    # Print the whole lines
    my( $i );
    for ($i = 0; $i < $whole; $i += $nuc) {
        my $blocks = pack $out_pat,
                     unpack $whole_pat,
                     substr($str, $i, $nuc);
        $self->_print(sprintf("     $blocks%9d\n", $i + $nuc));
    }

    # Print the last line
    if (my $last = substr($str, $i)) {
        my $last_len = length($last);
        my $last_pat = 'a10' x int($last_len / 10) .'a'. $last_len % 10;
        my $blocks = pack $out_pat,
                     unpack($last_pat, $last);
        $self->_print(sprintf("     $blocks%9d\n", $length));    # Add the length to the end
    }

    $self->_print( "//\n");
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
   my ($self,$fth,$always_quote) = @_;
   $always_quote ||= 0;
   
   if( ! ref $fth || ! $fth->isa('Bio::SeqIO::FTHelper') ) {
       $fth->warn("$fth is not a FTHelper class. Attempting to print, but there could be tears!");
   }


   #$self->_print( "FH   Key             Location/Qualifiers\n");
   #$self->_print( sprintf("FT   %-15s  %s\n",$fth->key,$fth->loc));
   $self->_write_line_EMBL_regex(sprintf("FT   %-15s ",$fth->key),"FT                   ",$fth->loc,'\,|$',80); #'
   foreach my $tag ( keys %{$fth->field} ) {
       if( ! defined $fth->field->{$tag} ) { next; } 
       foreach my $value ( @{$fth->field->{$tag}} ) {
	   $value =~ s/\"/\"\"/g;
	   if ($value eq "_no_value") {
	       $self->_write_line_EMBL_regex("FT                   ","FT                   ","/$tag",'.|$',80); #'
	   }
           elsif( $always_quote == 1 || $value !~ /^\d+$/ ) {
              my $pat = $value =~ /\s/ ? '\s|$' : '.|$';
	      $self->_write_line_EMBL_regex("FT                   ","FT                   ","/$tag=\"$value\"",$pat,80);
           }
           else {
              $self->_write_line_EMBL_regex("FT                   ","FT                   ","/$tag=$value",'.|$',80); #'
           }
	  # $self->_print( "FT                   /", $tag, "=\"", $value, "\"\n");
       }
   }

}

#'
=head2 _read_EMBL_References

 Title   : _read_EMBL_References
 Usage   :
 Function: Reads references from EMBL format. Internal function really
 Example :
 Returns : 
 Args    :


=cut

sub _read_EMBL_References {
   my ($self,$buffer) = @_;
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
   my $pm;
   my $com;

   while( defined ($_ = $self->_readline) ) {
       /^R/ || last;
       /^RP   (\d+)-(\d+)/ && do {$b1=$1;$b2=$2;};
       /^RX   MEDLINE;\s+(\d+)/ && do {$med=$1};
       /^RX   PUBMED;\s+(\d+)/ && do {$pm=$1};       
       /^RA   (.*)/ && do {
	   $au = $self->_concatenate_lines($au,$1); next;
       };
       /^RT   (.*)/ && do {
	   $title = $self->_concatenate_lines($title,$1); next;
       };
       /^RL   (.*)/ && do {
	   $loc = $self->_concatenate_lines($loc,$1); next;
       };
       /^RC   (.*)/ && do {
	   $com = $self->_concatenate_lines($com,$1); next;
       };
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
   $ref->pubmed($pm);

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
    my( $self, $buffer ) = @_;
    my $org;

    $_ = $$buffer;
    my( $sub_species, $species, $genus, $common, @class );
    while (defined( $_ ||= $self->_readline )) {
        
        if (/^OS\s+(\S+)(?:\s+([^\(]\S*))?(?:\s+([^\(]\S*))?(?:\s+\((.*)\))?/) {
            $genus   = $1;
	    $species = $2 || 'sp.';
	    $sub_species = $3 if $3;
            $common      = $4 if $4;
        }
        elsif (s/^OC\s+//) {
	    # only split on ';' or '.' so that 
	    # classification that is 2 words will 
	    # still get matched
	    # use map to remove trailing/leading spaces
	    chomp;
            push(@class,  map { s/^\s+//; s/\s+$//; $_; } split /[;\.]+/);
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

=head2 _read_EMBL_DBLink

 Title   : _read_EMBL_DBLink
 Usage   :
 Function: Reads the EMBL database cross reference ("DR") lines
 Example :
 Returns : A list of Bio::Annotation::DBLink objects
 Args    :

=cut

sub _read_EMBL_DBLink {
    my( $self,$buffer ) = @_;
    my( @db_link );

    $_ = $$buffer;
    while (defined( $_ ||= $self->_readline )) {
        
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
 Usage   : _read_FTHelper_EMBL($buffer)
 Function: reads the next FT key line
 Example :
 Returns : Bio::SeqIO::FTHelper object 
 Args    : filehandle and reference to a scalar


=cut

sub _read_FTHelper_EMBL {
    my ($self,$buffer) = @_;
    
    my ($key,   # The key of the feature
        $loc,   # The location line from the feature
        @qual,  # An arrray of lines making up the qualifiers
        );
    
    if ($$buffer =~ /^FT   (\S+)\s+(\S+)/) {
        $key = $1;
        $loc = $2;
        # Read all the lines up to the next feature
        while ( defined($_ = $self->_readline) ) {
            if (/^FT(\s+)(.+?)\s*$/) {
                # Lines inside features are preceeded by 19 spaces
                # A new feature is preceeded by 3 spaces
                if (length($1) > 4) {
                    # Add to qualifiers if we're in the qualifiers
                    if (@qual) {
                        push(@qual, $2);
                    }
                    # Start the qualifier list if it's the first qualifier
                    elsif (substr($2, 0, 1) eq '/') {
                        @qual = ($2);
                    }
                    # We're still in the location line, so append to location
                    else {
                        $loc .= $2;
                    }
                } else {
                    # We've reached the start of the next feature
                    last;
                }
            } else {
                # We're at the end of the feature table
                last;
            }
        }
    } else {
        # No feature key
        return;
    } 
    
    # Put the first line of the next feature into the buffer
    $$buffer = $_;

    # Make the new FTHelper object
    my $out = new Bio::SeqIO::FTHelper(-verbose => $self->verbose());
    $out->key($key);
    $out->loc($loc);

    # Now parse and add any qualifiers.  (@qual is kept
    # intact to provide informative error messages.)
  QUAL: for (my $i = 0; $i < @qual; $i++) {
        $_ = $qual[$i];
        my( $qualifier, $value ) = m{^/([^=]+)(?:=(.+))?}
            or $self->throw("Can't see new qualifier in: $_\nfrom:\n"
                . join('', map "$_\n", @qual));
        if (defined $value) {
            # Do we have a quoted value?
            if (substr($value, 0, 1) eq '"') {
                # Keep adding to value until we find the trailing quote
                # and the quotes are balanced
                while ($value !~ /"$/ or $value =~ tr/"/"/ % 2) { #"
                    $i++;
                    my $next = $qual[$i];
                    unless (defined($next)) {
                        warn("Unbalanced quote in:\n", map("$_\n", @qual),
                            "No further qualifiers will be added for this feature");
                        last QUAL;
                    }

                    # Join to value with space if value or next line contains a space
                    $value .= (grep /\s/, ($value, $next)) ? " $next" : $next;
                }
                # Trim leading and trailing quotes
                $value =~ s/^"|"$//g;
                # Undouble internal quotes
                $value =~ s/""/"/g; #"
            }
        } else {
            $value = '_no_value';
        }

        # Store the qualifier
        $out->field->{$qualifier} ||= [];
        push(@{$out->field->{$qualifier}},$value);
    }   

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
   my ($self,$pre1,$pre2,$line,$length) = @_;

   $length || die "Miscalled write_line_EMBL without length. Programming error!";
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
    my ($self,$pre1,$pre2,$line,$regex,$length) = @_;

    #print STDOUT "Going to print with $line!\n";

    $length || die "Programming error - called write_line_EMBL_regex without length.";

    if( length $pre1 != length $pre2 ) {
        die "Programming error - called write_line_EMBL_regex with different length pre1 and pre2 tags!";
    }

    my $subl = $length - (length $pre1) -1 ;



    my( @lines );
    while(defined $line && 
	  $line =~ m/(.{1,$subl})($regex)/g) {
	push(@lines, $1.$2);
    }
    foreach (@lines) { s/\s+$//; }
    
    # Print first line
    my $s = shift(@lines) || '';    
    $self->_print( "$pre1$s\n");
    
    # Print the rest
    foreach my $s ( @lines ) {
	$s = '' if( !defined $s );
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
