# $Id$
#
# BioPerl module for Bio::SeqIO::swiss
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::swiss - Swissprot sequence input/output stream

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

=over 3

=item _show_dna()

(output only) shows the dna or not

=item _post_sort()

(output only) provides a sorting func which is applied to the FTHelpers
before printing

=item _id_generation_func()

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

  bioperl-l@bioperl.org          - General discussion
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
use Bio::SeqIO;
use Bio::SeqIO::FTHelper;
use Bio::SeqFeature::Generic;
use Bio::Species;
use Bio::Tools::SeqStats;
use Bio::Seq::SeqFactory;
use Bio::Annotation::Collection;
use Bio::Annotation::Comment;
use Bio::Annotation::Reference;
use Bio::Annotation::DBLink;
use Bio::Annotation::OntologyTerm;
use Bio::Annotation::SimpleValue;

@ISA = qw(Bio::SeqIO);


sub _initialize {
  my($self,@args) = @_;
  $self->SUPER::_initialize(@args);
   
  # hash for functions for decoding keys.
  $self->{'_func_ftunit_hash'} = {};
  $self->_show_dna(1); # sets this to one by default. People can change it
  if( ! defined $self->sequence_factory ) {
      $self->sequence_factory(new Bio::Seq::SeqFactory
			      (-verbose => $self->verbose(), 
			       -type => 'Bio::Seq::RichSeq'));      
  }
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
       $date,$comment,@date_arr, @sec);
   my ($keywords,$acc_string);
   my ($annotation, %params, @features) = ( new Bio::Annotation::Collection);
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
   
   # fixed to allow _DIVISION to be optional for bug #946
   # see bug report for more information
   $line =~ /^ID\s+([^\s_]+)(_([^\s_]+))?\s+([^\s;]+);\s+([^\s;]+);/ 
     || $self->throw("swissprot stream with no ID. Not swissprot in my book");
   
   if( $3 ) {
       $name = "$1$2";
       $params{'-division'} = $3;
   } else {
       $name = $1;
       $params{'-division'} = 'UNK';
       $params{'-primary_id'} = $1;
   }
   $params{'-alphabet'} = 'protein';
    # this is important to have the id for display in e.g. FTHelper, otherwise
    # you won't know which entry caused an error
   $params{'-display_id'} = $name;
   
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
       elsif(/^GN\s+([^\.]+)/) {
           # Drop trailing spaces and dots
           s/[\. ]*$//;  # imported from swiss knife by cjm
           if (/^GN\s+(.*)/) {
               foreach my $gn (split(/ OR /,$1)) {
                   my $sim = Bio::Annotation::SimpleValue->new();
                   $sim->value($gn);		   
                   $annotation->add_Annotation('gene_name',$sim);
               }
           }
       }
       #accession number(s)
       elsif( /^AC\s+(.+)/) {
           $acc_string .= $acc_string ? " $1" : $1;
       }
       #version number
       elsif( /^SV\s+(\S+);?/ ) {
	   my $sv = $1;
	   $sv =~ s/\;//;
	   $params{'-seq_version'} = $sv;
       }
       #date
       elsif( /^DT\s+(.*)/ ) {
	   my $date = $1;
	   $date =~ s/\;//;
	   $date =~ s/\s+$//;
	   push @{$params{'-dates'}}, $date;
       }
       # Organism name and phylogenetic information
       elsif (/^O[SCG]/) {
           my $species = $self->_read_swissprot_Species(\$buffer);
           $params{'-species'}= $species;
	   # now we are one line ahead -- so continue without reading the next
	   # line   HL 05/11/2000
	   next;
       }
       # References
       elsif (/^R/) {
	   my $refs = $self->_read_swissprot_References(\$buffer);

	   foreach my $r (@$refs) {
	       $annotation->add_Annotation('reference',$r);
	   }
	   # now we are one line ahead -- so continue without reading the next
	   # line   HL 05/11/2000
	   next;
       }
       #Comments
       elsif (/^CC\s{3}(.*)/) {
	   $comment .= $1;
	   $comment .= "\n";
	   while (defined ($buffer = $self->_readline)) {
	       if ($buffer =~ /^CC\s{3}(.*)/) {
		   $comment .= $1;
		   $comment .= "\n";
	       }
	       else {
		   last;
	       }
	   }
	   my $commobj = Bio::Annotation::Comment->new();
	   # note: don't try to process comments here -- they may contain
           # structure. LP 07/30/2000
	   $commobj->text($comment);
	   $annotation->add_Annotation('comment',$commobj);
	   $comment = "";
	   # now we are one line ahead -- so continue without reading the next
	   # line   HL 05/11/2000
	   next;
       }
       #DBLinks
       elsif (/^DR\s+(\S+)\;\s+(\S+)\;\s+(\S+)[\;\.](.*)$/) {
	   my $dblinkobj =  Bio::Annotation::DBLink->new();
	   $dblinkobj->database($1);
	   $dblinkobj->primary_id($2);
	   $dblinkobj->optional_id($3);
	   my $comment = $4;
	   if(length($comment) > 0) {
	       # edit comment to get rid of leading space and trailing dot
	       if( $comment =~ /^\s*(\S+)\./ ) {
		   $dblinkobj->comment($1);
	       } else {
		   $dblinkobj->comment($comment);
	       }
	   }
	   $annotation->add_Annotation('dblink',$dblinkobj);
       }
       #keywords
       elsif( /^KW\s+(.*)$/ ) {
           $keywords .= $keywords ? " $1" : $1;
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
	   push(@features,
		$ftunit->_generic_seqfeature($self->location_factory(),
					     $params{'-seqid'}, "SwissProt"));
       } else {
	   $self->warn("failed to parse feature table line for seq " .
		       $params{'-display_id'});
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
   $acc_string =~ s/\;\s*/ /g;
   ( $acc, @sec ) = split " ",$acc_string;

   my $seq=  $self->sequence_factory->create
       (-verbose  => $self->verbose,
	%params,
	-seq      => $seqc,
	-desc     => $desc,
	-keywords => $keywords,
	-accession_number => $acc,
	-secondardy_accessions => \@sec,
	-features => \@features,
	-annotation => $annotation,
	);

   # The annotation doesn't get added by the contructor
   $seq->annotation($annotation);

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
   
   if( ! $seq->can('alphabet') || ! defined ($mol = $seq->alphabet) ) {
       $mol = 'XXX';
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
       $mol =~ s/protein/PRT/;
       $temp_line = sprintf ("%10s     STANDARD;      %3s;   %d AA.",
			     $seq->display_id(), $mol, $len);
   }

   $self->_print( "ID   $temp_line\n");

   # if there, write the accession line
   local($^W) = 0;   # supressing warnings about uninitialized fields

   if( $self->_ac_generation_func ) {
       $temp_line = &{$self->_ac_generation_func}($seq);
       $self->_print( "AC   $temp_line\n");
   } else {
       if ($seq->can('accession_number') ) {
	   $self->_print("AC   ",$seq->accession_number,";");
	   if ($seq->can('get_secondary_accessions') ) {
	     foreach my $sacc ($seq->get_secondary_accessions) {
	       $self->_print(" ",$sacc,";");
	     }
	     $self->_print("\n");
	   }
	   else {
	       $self->_print("\n");
	   }
       }
       # otherwise - cannot print <sigh>
   }

   # Date lines

   if( $seq->can('get_dates') ) {
       foreach my $dt ( $seq->get_dates() ) {
	   $self->_write_line_swissprot_regex("DT   ","DT   ",
					      $dt,"\\s\+\|\$",80);
       }
   }

   #Definition lines
   $self->_write_line_swissprot_regex("DE   ","DE   ",$seq->desc(),"\\s\+\|\$",80);

   #Gene name
   if ((my @genes = $seq->annotation->get_Annotations('gene_name') ) ) {
       $self->_print("GN   ",join(' OR ', map { $_->value } @genes),".\n");
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
            $OS .= " ($common).";
        }
	if ($class[$#class] =~ /viruses/i) { # different OS / OC syntax
	  $OS = $spec->common_name;          # for viruses LP 09/16/2000
	  unshift @class, $species;
	}
        $self->_print( "OS   $OS\n");
        my $OC = join('; ', reverse(@class)) .'.';
        $self->_write_line_swissprot_regex("OC   ","OC   ",$OC,"\; \|\$",80);
	if ($spec->organelle) {
	    $self->_write_line_swissprot_regex("OG   ","OG   ",$spec->organelle,"\; \|\$",80);
	}
	if ($spec->ncbi_taxid) {
	    $self->_print("OX   NCBI_TaxID=".$spec->ncbi_taxid.";\n");
	}
   }
   
   # Reference lines
   my $t = 1;
   foreach my $ref ( $seq->annotation->get_Annotations('reference') ) {
       $self->_print( "RN   [$t]\n");
       # changed by lorenz 08/03/00
       # j.gilbert and h.lapp agreed that the rp line in swissprot seems 
       # more like a comment than a parseable value, so print it as is
       if ($ref->rp) {
	 $self->_write_line_swissprot_regex("RP   ","RP   ",$ref->rp,
					    "\\s\+\|\$",80);
       }
       if ($ref->comment) {
	 $self->_write_line_swissprot_regex("RC   ","RC   ",$ref->comment,
					    "\\s\+\|\$",80);
       }
       if ($ref->medline) {
	 # new RX format in swissprot LP 09/17/00
	 if ($ref->pubmed) {
	 $self->_write_line_swissprot_regex("RX   ","RX   ",
					    "MEDLINE=".$ref->medline.
					    "; PubMed=".$ref->pubmed.";",
					    "\\s\+\|\$",80);
	 } else {
	 $self->_write_line_swissprot_regex("RX   MEDLINE; ","RX   MEDLINE; ",
					    $ref->medline.".","\\s\+\|\$",80);
	 }
       }
       my $author = $ref->authors .';' if($ref->authors);
       my $title = $ref->title .';' if( $ref->title);
       
       $self->_write_line_swissprot_regex("RA   ","RA   ",$author,"\\s\+\|\$",80);
       $self->_write_line_swissprot_regex("RT   ","RT   ",$title,"\\s\+\|\$",80);
       $self->_write_line_swissprot_regex("RL   ","RL   ",$ref->location,"\\s\+\|\$",80);
       $t++;
   }
   
   # Comment lines

   foreach my $comment ( $seq->annotation->get_Annotations('comment') ) {
       foreach my $cline (split ("\n", $comment->text)) {
	 while (length $cline > 74) {
	   $self->_print("CC   ",(substr $cline,0,74),"\n");
	   $cline = substr $cline,74;
	 }
	 $self->_print("CC   ",$cline,"\n");
       }
   }

   foreach my $dblink ( $seq->annotation->get_Annotations('dblink') ) 
   {
     if (defined($dblink->comment)&&($dblink->comment)) {
	 $self->_print("DR   ",$dblink->database,"; ",$dblink->primary_id,"; ",
		       $dblink->optional_id,"; ",$dblink->comment,".\n");
     } elsif($dblink->optional_id) {
	 $self->_print("DR   ",$dblink->database,"; ",$dblink->primary_id,"; ",
		       $dblink->optional_id,".\n");
     }
     else {
	 $self->_print("DR   ",$dblink->database,"; ",$dblink->primary_id,"; ",
		       "-.\n");
     }
   }   

   # if there, write the kw line
   
   if( $self->_kw_generation_func ) {
       $temp_line = &{$self->_kw_generation_func}($seq);
       $self->_print( "KW   $temp_line\n");
   } else {
       if( $seq->can('keywords') ) {
	 $self->_write_line_swissprot_regex("KW   ","KW   ",
					    $seq->keywords,"\\s\+\|\$",80);       
       }
   }

#Check if there is seqfeatures before printing the FT line
   my @feats = $seq->top_SeqFeatures;
       if ($feats[0]) {
         if( defined $self->_post_sort ) {
	     
	     # we need to read things into an array. Process. Sort them. Print 'em
	     
	     my $post_sort_func = $self->_post_sort();
	     my @fth;
	     
	   my @feats = $seq->top_SeqFeatures;
	   
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
   }
   # finished printing features.

   # molecular weight
   my $mw = ${Bio::Tools::SeqStats->get_mol_wt($seq->primary_seq)}[0];
   # checksum
   # was crc32 checksum, changed it to crc64 
   my $crc64 = $self->_crc64(\$str); 
   $self->_print( sprintf("SQ   SEQUENCE  %4d AA;  %d MW;  %16s CRC64;\n",
			  $len,$mw,$crc64));
   $self->_print( "     ");
   my $linepos;
   for ($i = 0; $i < length($str); $i += 10) {
       $self->_print( substr($str,$i,10), " ");
       $linepos += 11;
       if( ($i+10)%60 == 0 && (($i+10) < length($str))) {
	   $self->_print( "\n     ");
      }
   }
   $self->_print( "\n//\n");

   $self->_fh->flush if $self->_flush_on_write && defined $self->_fh;
   return 1;
}

# Thanks to James Gilbert for the following two. LP 08/01/2000

=head2 _generateCRCTable

 Title   : _generateCRCTable
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _generateCRCTable {
  # 10001000001010010010001110000100
  # 32
  my $poly = 0xEDB88320;
  my ($self) = shift;
        
  $self->{'_crcTable'} = [];
  foreach my $i (0..255) {
    my $crc = $i;
    for (my $j=8; $j > 0; $j--) {
      if ($crc & 1) {
	$crc = ($crc >> 1) ^ $poly;
      }
      else {
	$crc >>= 1;
      }
    }
    ${$self->{'_crcTable'}}[$i] = $crc;
  }
}


=head2 _crc32

 Title   : _crc32
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _crc32 {
  my( $self, $str ) = @_;
  
  $self->throw("Argument to crc32() must be ref to scalar")
    unless ref($str) eq 'SCALAR';
  
  $self->_generateCRCTable() unless exists $self->{'_crcTable'};
  
  my $len = length($$str);
  
  my $crc = 0xFFFFFFFF;
  for (my $i = 0; $i < $len; $i++) {
    # Get upper case value of each letter
    my $int = ord uc substr $$str, $i, 1;
    $crc = (($crc >> 8) & 0x00FFFFFF) ^ 
      ${$self->{'_crcTable'}}[ ($crc ^ $int) & 0xFF ];
  }
  return $crc;
}

=head2 _crc64

 Title   : _crc64
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _crc64{
    my ($self, $sequence) = @_;
    my $POLY64REVh = 0xd8000000;
    my @CRCTableh = 256;
    my @CRCTablel = 256;
    my $initialized;       
    

    my $seq = $$sequence;
      
    my $crcl = 0;
    my $crch = 0;
    if (!$initialized) {
	$initialized = 1;
	for (my $i=0; $i<256; $i++) {
	    my $partl = $i;
	    my $parth = 0;
	    for (my $j=0; $j<8; $j++) {
		my $rflag = $partl & 1;
		$partl >>= 1;
		$partl |= (1 << 31) if $parth & 1;
		$parth >>= 1;
		$parth ^= $POLY64REVh if $rflag;
	    }
	    $CRCTableh[$i] = $parth;
	    $CRCTablel[$i] = $partl;
	}
    }
    
    foreach (split '', $seq) {
	my $shr = ($crch & 0xFF) << 24;
	my $temp1h = $crch >> 8;
	my $temp1l = ($crcl >> 8) | $shr;
	my $tableindex = ($crcl ^ (unpack "C", $_)) & 0xFF;
	$crch = $temp1h ^ $CRCTableh[$tableindex];
	$crcl = $temp1l ^ $CRCTablel[$tableindex];
    }
    my $crc64 = sprintf("%08X%08X", $crch, $crcl);
        
    return $crc64;
      
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
   $always_quote ||= 0;
   my ($start,$end) = ('?', '?');
   
   if( ! ref $fth || ! $fth->isa('Bio::SeqIO::FTHelper') ) {
       $fth->warn("$fth is not a FTHelper class. ".
		  "Attempting to print, but there could be tears!");
   }

   if( $fth->loc =~ /(\?|\d+)?\.\.(\?|\d+)?/ ) {
       $start = $1 if defined $1;
       $end = $2 if defined $2;

       # to_FTString only returns one value when start == end, #JB955
       # so if no match is found, assume it is both start and end #JB955
   } else {
       $start = $end = $fth->loc; 
   }
   
   my $desc = "";
   $desc = @{$fth->field->{"description"}}[0]."." 
       if exists $fth->field->{"description"};
   $self->_write_line_swissprot_regex(sprintf("FT   %-8s %6s %6s       ",
					      substr($fth->key,0,8),
					      $start,$end),
				      "FT                                ",
				      $desc,'\s+|$',80);
}
#'

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
   my ($b1, $b2, $rp, $title, $loc, $au, $med, $com, $pubmed);
   
   if ($$buffer !~ /^RP/) {
       $$buffer = $self->_readline;
   }
   if( !defined $$buffer ) { return undef; }
   if( $$buffer =~ /^RP/ ) {
       if ($$buffer =~ /^RP   (SEQUENCE OF (\d+)-(\d+).*)/) { 
	   $rp=$1;
	   $b1=$2;
	   $b2=$3; 
       }
       elsif ($$buffer =~ /^RP   (.*)/) {
	   $rp=$1;
       }
       
   }
   while( defined ($_ = $self->_readline) ) {
       #/^CC/ && last;
       /^RN/ && last; # separator between references ! LP 07/25/2000
       #/^SQ/ && last; # there may be sequences without CC lines! HL 05/11/2000
       /^[^R]/ && last; # may be the safest exit point HL 05/11/2000
       /^RX   MEDLINE;\s+(\d+)/ && do {$med=$1};
       /^RX   MEDLINE=(\d+);\s+PubMed=(\d+);/ && do {$med=$1;$pubmed=$2};
       /^RA   (.*)/ && do { $au .= $au ? " $1" : $1;   next;};
       /^RT   (.*)/ && do { $title .= $title ? " $1" : $1; next;};
       /^RL   (.*)/ && do { $loc .= $loc ? " $1" : $1; next;};
       /^RC   (.*)/ && do { $com .= $com ? " $1" : $1; next;};
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
   $ref->pubmed($pubmed) if (defined $pubmed);
   $ref->comment($com);
   $ref->rp($rp);

   push(@refs,$ref);
   $$buffer = $_;
   return \@refs;
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
    my( $sub_species, $species, $genus, $common, @class, $osline, $ncbi_taxid );
    while (defined( $_ ||= $self->_readline )) {
	last unless /^O[SCGX]/;
	# believe it or not, but OS may come multiple times -- at this time
	# we can't capture multiple species
        if((! defined($genus)) &&
	   /^OS\s+((\S+)(?:\s+([^\(]\S*))?(?:\s+([^\(]\S*))?(?:\s+\((.*)\))?.*)/) {
	    $osline = $1;
            $genus   = $2;
	    if ($3) {
		$species = $3;
		# remove trailing dot (or comma if multiple species)
		$species =~ s/[\.,]$//;
	    } else {
		$species = "sp.";
	    }
	    $sub_species = $4 if $4;
            $common      = $5 if $5;
        }
        elsif (s/^OC\s+//) {
            push(@class, split /[\;\.]\s*/);
	    if($class[0] =~ /viruses/i) { # viruses have different OS/OC syntax
		$common = $osline;        # LP 09/16/2000
	    }
        }
	elsif (/^OG\s+(.*)/) {
	    $org = $1;
	}
	elsif (/^OX\s+(.*)\;/) {
            my $taxstring = $1;
            if ($taxstring =~ /NCBI_TaxID=([\w\d]+)/) {
                $ncbi_taxid = $1;
            }
            else {
                $self->throw("$taxstring doesn't look like NCBI_TaxID");
            }
	}
        
        $_ = undef; # Empty $_ to trigger read of next line
    }
    
    $$buffer = $_;
    
    # Don't make a species object if it is "Unknown" or "None"
    return if $genus =~ /^(Unknown|None)$/i;

    if ($class[0] !~ /viruses/i) { # different OS / OC syntax for viruses
      # Bio::Species array needs array in Species -> Kingdom direction
      if ($class[$#class] eq $genus) {
        push( @class, $species );
      } else {
        push( @class, $genus, $species );
      }
    }
    @class = reverse @class;
    
    my $make = Bio::Species->new();
    $make->classification( @class );
    $make->common_name( $common      ) if $common;
    $make->sub_species( $sub_species ) if $sub_species;
    $make->organelle  ( $org         ) if $org;
    $make->ncbi_taxid ( $ncbi_taxid  ) if $ncbi_taxid;
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
	    $desc = $4;
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
    my $out = new Bio::SeqIO::FTHelper(-verbose => $self->verbose());
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
       print STDERR "len 1 is ", length $pre1, " len 2 is ", length $pre2, "\n";
       die "Programming error - cannot called write_line_swissprot_regex with different length \npre1 ($pre1) and \npre2 ($pre2) tags!";
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
