# $Id$
#
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

  bioperl-l@bioperl.org                  - General discussion
  http://www.bioperl.org/MailList.shtml  - About the mailing lists

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

use Bio::Seq::RichSeq;
use Bio::SeqIO;
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

sub next_seq{
    my ($self,@args) = @_;
    my ($pseq,$c,$name,$desc,$seqc,$mol,$div,$date);
    my $line;
    my @acc = ();
    my $seq = Bio::Seq::RichSeq->new('-verbose' =>$self->verbose());
    
    while(defined($line = $self->_readline())) {
	$line =~ /^LOCUS\s+\S+/ && last;
    }
    return undef if( !defined $line ); # end of file
    
    $line =~ /^LOCUS\s+\S+/ ||
	$self->throw("GenBank stream with no LOCUS. Not GenBank in my book. Got $line");
    $line =~ /^LOCUS\s+(\S+)\s+\S+\s+(bp|aa)\s+(\S+)\s+(\S+)\s+(\S+)?/ || do {
	$line =~ /^LOCUS\s+(\S+)/ ||
	    $self->throw("GenBank stream with no LOCUS. Not GenBank in my book. Got $line");
	# fall back to at least an ID
	$name = $1;
	$self->warn("GenBank record with LOCUS line in unexpected format. ".
		    "Attributes from this line other than ID will be missing.");
    };
    $name = $1;
    # this is important to have the id for display in e.g. FTHelper, otherwise
    # you won't know which entry caused an error
    $seq->display_id($name);
    # the moltype of the entry
    if($2 eq 'bp') {
	$seq->moltype('dna');
    } else {
	# $2 eq 'aa'
	$seq->moltype('protein');
    }
    # for aa there is usually no 'molecule' (mRNA etc)
    if (($2 eq 'bp') || defined($5)) {
	$seq->molecule($3);
	$seq->division($4);
	$date = $5;
    } else {
	$seq->molecule('PRT') if($2 eq 'aa');
	$seq->division($3);
	$date = $4;
    }
    if ($date) {
	$seq->add_date($date);
    }
    my $buffer = $self->_readline();
        
    until( !defined ($buffer) ) {
	$_ = $buffer;
	
	# Description line(s)
	if (/^DEFINITION\s+(\S.*\S)/) {
	    $desc .= $desc ? " $1" : $1;
	    $desc .= " ";
	    while ( defined($_ = $self->_readline) ) { 
		/^\s+(.*)/ && do { $desc .= $1; next;};
		last;
	    }
	}		 
	
	# accession number (there can be multiple accessions)
	if( /^ACCESSION\s+(\S.*\S)/ ) {
	    push(@acc, split(' ',$1));
	}
	
	# PID
	if( /^PID\s+(\S+)/ ) {
	    $seq->pid($1);
	}
	
	#Version number
	if( /^VERSION\s+(\S+)\.(\d+)\s*(GI:\d+)?/ ) {
	    $seq->seq_version($2);
	    $seq->primary_id(substr($3, 3)) if($3);
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
		if (/^\S/) {
		    $comment =~ s/\n/ /g;
		    $comment =~ s/  +/ /g;
		    my $commobj = Bio::Annotation::Comment->new();
		    $commobj->text($comment);
		    $seq->annotation->add_Comment($commobj);
		    last;
		}
		$comment .= $_; 
	    }
	    $buffer = $_;
	    next;
	}
	# Exit at start of Feature table, or start of sequence
	last if( /^(FEATURES)|(ORIGIN)/ );
	# Get next line and loop again
	$buffer = $self->_readline;
    }
    return undef if(! defined($buffer));
    $seq->desc($desc); # don't forget!
    $seq->accession_number(shift(@acc));
    $seq->add_secondary_accession(@acc);
    
    # some "minimal" formats may not necessarily have a feature table
    if(/^FEATURES/) {
	# need to read the first line of the feature table
	$buffer = $self->_readline;

	# DO NOT read lines in the while condition -- this is done as a side
	# effect in _read_FTHelper_GenBank!
	while( defined($buffer) ) {
	    # check immediately -- not at the end of the loop
	    # note: GenPept entries obviously do not have a BASE line
	    last if(($buffer =~ /^BASE/) || ($buffer =~ /^ORIGIN/));
	    # slurp in one feature at a time -- at return, the start of
	    # the next feature will have been read already, so we need
	    # to pass a reference, and the called method must set this
	    # to the last line read before returning 
	    my $ftunit = $self->_read_FTHelper_GenBank(\$buffer);
	    # process ftunit
	    $ftunit->_generic_seqfeature($seq);
	}
	$_ = $buffer;
    }
    # advance to the section with the sequence
    if(! /^ORIGIN/) {
	$seqc = "";	
	while (defined( $_ = $self->_readline)) {
	    last if /^ORIGIN/;
	}
    }
    while( defined($_ = $self->_readline) ) {
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
	my $date = '';
	if( $seq->can('get_dates') ) { 
	    ($date) = $seq->get_dates(); # get first one from the list
	}
	$temp_line = sprintf ("%-12s%-10s%10s bp%8s%5s %3s ", 'LOCUS',$seq->id(),$len,$mol,$div,$date);
    } 
    
    $self->_print("$temp_line\n");   
    $self->_write_line_GenBank_regex("DEFINITION  ","            ",$seq->desc(),"\\s\+\|\$",80);
    
    # if there, write the accession line

    if( $self->_ac_generation_func ) {
	$temp_line = &{$self->_ac_generation_func}($seq);
	$self->_print("ACCESSION   $temp_line\n");   
    } else {
	my @acc = ();
	push(@acc, $seq->accession_number());
	if( $seq->isa('Bio::Seq::RichSeqI') ) {
	    push(@acc, $seq->get_secondary_accessions());
	}
	$self->_print("ACCESSION   ", join(" ", @acc), "\n");
	# otherwise - cannot print <sigh>
    } 

    # if PID defined, print it
    if($seq->isa('Bio::Seq::RichSeqI') && $seq->pid()) {
	$self->_print("PID         ", $seq->pid(), "\n");
    }

    # if there, write the version line
    
    if( defined $self->_sv_generation_func() ) {
	$temp_line = &{$self->_sv_generation_func}($seq);
	if( $temp_line ) {
	    $self->_print("VERSION     $temp_line\n");   
	}
    } else {
	if($seq->isa('Bio::Seq::RichSeqI') && defined($seq->seq_version)) {
	    my $id = $seq->primary_id(); # this may be a GI number
	    $self->_print("VERSION     ",
			  $seq->accession_number(), ".", $seq->seq_version,
			  ($id && ($id =~ /^\d+$/) ? " GI:".$id : ""),
			  "\n");
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
        my ($species, $genus, @class) = $spec->classification();
        my $OS = "$genus $species";
        if (my $ssp = $spec->sub_species) {
            $OS .= " $ssp";
        }
	$self->_print("SOURCE      $OS.\n");
	$self->_print("  ORGANISM  ",
		      ($spec->organelle() ? $spec->organelle()." " : ""),
		      $OS, "\n");
        my $OC = join('; ', (reverse(@class), $genus)) .'.';
        $self->_write_line_GenBank_regex("            ","            ",
					 $OC,"\\s\+\|\$",80);
    }
    
    # Reference lines
    my $count = 1;
    foreach my $ref ( $seq->annotation->each_Reference() ) {
	$temp_line = sprintf ("REFERENCE   $count  (%s %d to %d)",
			      ($seq->moltype() eq "protein" ?
			       "residues" : "bases"),
			      $ref->start,$ref->end);
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
    my $base_count = sprintf("BASE COUNT %8s a %6s c %6s g %6s t %6s others\n",
			     $alen,$clen,$glen,$tlen,$olen);
    $self->_print($base_count); 
    $self->_print("ORIGIN      \n");
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
              my $pat = $value =~ /\s/ ? '\s|$' : '.|$';
	      $self->_write_line_GenBank_regex("                     ",
					       "                     ",
					       "/$tag=\"$value\"",$pat,80);
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
   my $ref;
   
   # assumme things are starting with RN

   if( $$buffer !~ /^REFERENCE/ ) {
       warn("Not parsing line '$$buffer' which maybe important");
   }

   $_ = $$buffer;

   my $title;
   my $loc;
   my $au;
   my $com;

   while( defined($_) || defined($_ = $self->_readline) ) {
       if (/^  AUTHORS\s+(.*)/) { 
	   $au .= $1;   
	   while ( defined($_ = $self->_readline) ) {
	       /^\s{3,}(.*)/ && do {
		   $au .= $1;
		   $au =~ s/\,(\S)/ $1/g;
		   $au .= " ";
		   next;
	       };
	       last;
	   }
	   $ref->authors($au);
       }
       if (/^  TITLE\s+(.*)/)  { 
	   $title .= $1;
	   while ( defined($_ = $self->_readline) ) {
	       /^\s{3,}(.*)/ && do { $title .= " " . $1;
				     next;};
	       last;
	   }
	   $ref->title($title);
       }
       if (/^  JOURNAL\s+(.*)/) { 
	   $loc .= $1;
	   while ( defined($_ = $self->_readline) ) {
	       /^\s{3,}(.*)/ && do { $loc .= $1;$loc .= " ";next;};
	       last;
	   }
	   $ref->location($loc);
       }
       if (/^  REMARK\s+(.*)/) { 
	   $com .= $1;
	   while ( defined($_ = $self->_readline) ) {	       
	       /^\s{3,}(.*)/ && do { $com .= $1;$com .= " ";next;};
	       last;
	   }
	   $ref->comment($com);
       }
       /^REFERENCE/ && do {
	   # store current reference
	   $self->_add_ref_to_array(\@refs,$ref) if $ref;
	   # reset
	   $au = "";
	   $title = "";
	   $loc = "";
	   $com = "";
	   # create the new reference object
	   $ref = Bio::Annotation::Reference->new();
	   # check whether start and end base is given
	   if (/^REFERENCE\s+\d+\s+\([a-z]+ (\d+) to (\d+)/){
	       $ref->start($1);
	       $ref->end($2);
	   }
       };

       /^(FEATURES)|(COMMENT)/ && last;

       $_ = undef; # Empty $_ to trigger read of next line
   }

   # store last reference
   $self->_add_ref_to_array(\@refs,$ref) if $ref;

   $$buffer = $_;

   #print "\nnumber of references found: ", $#refs+1,"\n";

   return @refs;
}

#
# This is undocumented as it shouldn't be called by anywhere else as
# read_GenBank_References. For those who still want to know:
#
# Purpose: adds a Reference object to an array of Reference objects, takes
#     care of possible cleanups to be done (currently, only author and title
#     will be chopped of trailing semicolons).
# Parameters:
#     a reference to an array of Reference objects
#     the Reference object to be added
# Returns: nothing
#
sub _add_ref_to_array {
    my ($self, $refs, $ref) = @_;

    # first, polish author and title by removing possible trailing semicolons
    my $au = $ref->authors();
    my $title = $ref->title();
    $au =~ s/;\s*$//g if $au;
    $title =~ s/;\s*$//g if $title;
    $ref->authors($au);
    $ref->title($title);
    # the rest should be clean already, so go ahead and add it
    push(@{$refs}, $ref);
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
    my @organell_names = ("chloroplast", "mitochondr"); 
             # only those carrying DNA, apart from the nucleus

    $_ = $$buffer;
    
    my( $sub_species, $species, $genus, $common, $organelle, @class );
    # upon first entering the loop, we must not read a new line -- the SOURCE
    # line is already in the buffer (HL 05/10/2000)
    while (defined($_) || defined($_ = $self->_readline())) {
	# de-HTMLify (links that may be encountered here don't contain
	# escaped '>', so a simple-minded approach suffices)
        s/<[^>]+>//g;
	if (/^SOURCE\s+(.*)/) {
	    # FIXME this is probably mostly wrong (e.g., it yields things like
	    # Homo sapiens adult placenta cDNA to mRNA
	    # which is certainly not what you want)
	    $common = $1;
	    $common =~ s/\.$//; # remove trailing dot
	} elsif (/^\s+ORGANISM/) {
	    my @spflds = split(' ', $_);
	    shift(@spflds); # ORGANISM
	    if(grep { $_ =~ /^$spflds[0]/i; } @organell_names) {
		$organelle = shift(@spflds);
	    }
            $genus = shift(@spflds);
	    if(@spflds) {
		$species = shift(@spflds);
	    } else {
		$species = "sp.";
	    }
	    $sub_species = shift(@spflds) if(@spflds);
        } elsif (/^\s+(.+)/) {
            push(@class, split /[;\s\.]+/, $1);
        } else {
            last;
        }
        
        $_ = undef; # Empty $_ to trigger read of next line
    }
    
    $$buffer = $_;
    
    # Don't make a species object if it's empty or "Unknown" or "None"
    return unless $genus and  $genus !~ /^(Unknown|None)$/i;
    
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
    $make->organelle($organelle) if $organelle;
    return $make;
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
    
    my ($key,   # The key of the feature
        $loc    # The location line from the feature
        );
    my @qual = (); # An arrray of lines making up the qualifiers
    
    if ($$buffer =~ /^     (\S+)\s+(\S+)/) {
        $key = $1;
        $loc = $2;
        # Read all the lines up to the next feature
        while ( defined($_ = $self->_readline) ) {
            if (/^(\s+)(.+?)\s*$/) {
                # Lines inside features are preceeded by 21 spaces
                # A new feature is preceeded by 5 spaces
                if (length($1) > 6) {
                    # Add to qualifiers if we're in the qualifiers, or if it's
		    # the first qualifier
                    if (($#qual >= 0) || (substr($2, 0, 1) eq '/')) {
                        push(@qual, $2);
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
    my $out = new Bio::SeqIO::FTHelper(-verbose=>$self->verbose());
    $out->key($key);
    $out->loc($loc);

    # Now parse and add any qualifiers.  (@qual is kept
    # intact to provide informative error messages.)
  QUAL: for (my $i = 0; $i < @qual; $i++) {
        $_ = $qual[$i];
        my( $qualifier, $value ) = (m{^/([^=]+)(?:=(.+))?})
	    or $self->warn("cannot see new qualifier in feature $key: ".
			   $qual[$i]);
            #or $self->throw("Can't see new qualifier in: $_\nfrom:\n"
            #    . join('', map "$_\n", @qual));
	$qualifier = '' unless( defined $qualifier);
        if (defined $value) {
            # Do we have a quoted value?
            if (substr($value, 0, 1) eq '"') {
                # Keep adding to value until we find the trailing quote
                # and the quotes are balanced
                while ($value !~ /\"$/ or $value =~ tr/"/"/ % 2) {
		    if($i >= $#qual) {
			# We haven't found the closing douple quote ...
			# Even though this should be considered as a malformed
			# entry, we allow for a single exception, namely if the
			# value ends exactly at the last char position of a
			# GenBank line.
			# At least 78 chars required, of which 21 are spaces
			#Currently disabled, let's wait for an actual sequence
			#entry that really requires this.
			#if(length($qual[$i]) >= 57) {
			#    $self->warn("unbalanced quotes in feature $key ".
			#		"(location: $loc), ".
			#		"qualifier $qualifier, ".
			#		"accepting though");
			#    last;
			#} else {
			    $self->warn("Unbalanced quote in:\n" .
					join('', map("$_\n", @qual)) .
					"No further qualifiers will " .
					"be added for this feature");
			    last QUAL;
			#}
                    }
                    $i++; # modifying a for-loop variable inside of the loop
		          # is not the best programming style ...
                    my $next = $qual[$i];

                    # Join to value with space if value or next line contains a space
                    $value .= (grep /\s/, ($value, $next)) ? " $next" : $next;
                }
                # Trim leading and trailing quotes
                $value =~ s/^"|"$//g;
                # Undouble internal quotes
                $value =~ s/""/\"/g;
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
