# $Id$
#
# BioPerl module for Bio::SeqIO::GenBank
#
# Cared for by Elia Stupka <elia@tll.org.sg>
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


If you want to output annotations in genbank format they need to be
stored in a Bio::Annotation::Collection object which is accessible
through the Bio::SeqI interface method L<annotation()|annotation>.  

The following are the names of the keys which are polled from a
L<Bio::Annotation::Collection> object.

reference       - Should contain Bio::Annotation::Reference objects
comment         - Should contain Bio::Annotation::Comment objects

segment         - Should contain a Bio::Annotation::SimpleValue object
origin          - Should contain a Bio::Annotation::SimpleValue object

=back

=head1 Where does the data go?

Data parsed in Bio::SeqIO::genbank is stored in a variety of data
fields in the sequence object that is returned.  More information in
the HOWTOs about exactly what each field means and where it goes.
Here is a partial list of fields.

Items listed as RichSeq or Seq or PrimarySeq and then NAME() tell you
the top level object which defines a function called NAME() which
stores this information.

Items listed as Annotation 'NAME' tell you the data is stored the
associated Bio::AnnotationCollectionI object which is associated with
Bio::Seq objects.  If it is explictly requested that no annotations
should be stored when parsing a record of course they won't be
available when you try and get them.  If you are having this problem
look at the type of SeqBuilder that is being used to contruct your
sequence object.

Comments             Annotation 'comment'
References           Annotation 'reference'
Segment              Annotation 'segment'
Origin               Annotation 'origin'

Accessions           PrimarySeq accession_number()
Secondary accessions RichSeq get_secondary_accessions()
GI number            PrimarySeq primary_id()
LOCUS                PrimarySeq display_id()
Keywords             RichSeq get_keywords()
Dates                RichSeq get_dates()
Molecule             RichSeq molecule()
Seq Version          RichSeq seq_version()
PID                  RichSeq pid()
Division             RichSeq division()
Features             Seq get_SeqFeatures()
Alphabet             PrimarySeq alphabet()
Definition           PrimarySeq description() or desc()
Version              PrimarySeq version()

Sequence             PrimarySeq seq()

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
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Elia Stupka

Email elia@tll.org.sg

=head1 CONTRIBUTORS

Ewan Birney birney at ebi.ac.uk
Jason Stajich jason at bioperl.org
Chris Mungall cjm at fruitfly.bdgp.berkeley.edu
Lincoln Stein lstein at cshl.org
Heikki Lehvaslaiho, heikki at ebi.ac.uk
Hilmar Lapp, hlapp at gmx.net
Donald G. Jackson, donald.jackson at bms.com

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::genbank;
use vars qw(@ISA %FTQUAL_NO_QUOTE);
use strict;

use Bio::SeqIO;
use Bio::SeqIO::FTHelper;
use Bio::SeqFeature::Generic;
use Bio::Species;
use Bio::Seq::SeqFactory;
use Bio::Annotation::Collection;
use Bio::Annotation::Comment;
use Bio::Annotation::Reference;
use Bio::Annotation::DBLink;

@ISA = qw(Bio::SeqIO);

%FTQUAL_NO_QUOTE=(
		  'anticodon'    => 1,
		  'citation'     => 1,
		  'codon'        => 1,
		  'codon_start'  => 1,
		  'cons_splice'  => 1,
		  'direction'    => 1,
		  'evidence'     => 1,
		  'label'        => 1,
		  'mod_base'     => 1,
		  'number'       => 1,
		  'rpt_type'     => 1,
		  'rpt_unit'     => 1,
		  'transl_except'=> 1,
		  'transl_table' => 1,
		  'usedin'       => 1,
		  );

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
    my $builder = $self->sequence_builder();
    my $seq;
    my %params;

  RECORDSTART: while (1) {
      my $buffer;
      my (@acc, @features);
      my ($display_id, $annotation);
      my $species;

      # initialize; we may come here because of starting over
      @features = ();
      $annotation = undef;
      @acc = ();
      $species = undef;
      %params = (-verbose => $self->verbose); # reset hash
      local($/) = "\n";
      while(defined($buffer = $self->_readline())) {
	  last if index($buffer,'LOCUS       ') == 0;
      }
      return undef if( !defined $buffer ); # end of file
      $buffer =~ /^LOCUS\s+(\S.*)$/o ||
	  $self->throw("GenBank stream with bad LOCUS line. Not GenBank in my book. Got '$buffer'");
      
      my @tokens = split(' ', $1);

      # this is important to have the id for display in e.g. FTHelper,
      # otherwise you won't know which entry caused an error
      $display_id = shift(@tokens);
      $params{'-display_id'} = $display_id;
      # may still be useful if we don't want the seq
      $params{'-length'} = shift(@tokens);
      # the alphabet of the entry
      $params{'-alphabet'} = (lc(shift @tokens) eq 'bp') ? 'dna' : 'protein';
      # for aa there is usually no 'molecule' (mRNA etc)
      if (($params{'-alphabet'} eq 'dna') || (@tokens > 2)) {	
	  $params{'-molecule'} = shift(@tokens);
	  my $circ = shift(@tokens);
	  if ($circ eq 'circular') {
	      $params{'-is_circular'} = 1;
	      $params{'-division'} = shift(@tokens);
	  } else {
	      # 'linear' or 'circular' may actually be omitted altogether
	      $params{'-division'} =
		  (CORE::length($circ) == 3 ) ? $circ : shift(@tokens);
	  }
      } else {
	  $params{'-molecule'} = 'PRT' if($params{'-alphabet'} eq 'aa');
	  $params{'-division'} = shift(@tokens);
      }
      my $date = join(' ', @tokens); # we lump together the rest

      # this is per request bug #1513
      # we can handle
      # 9-10-2003
      # 9-10-03
      #09-10-2003
      #09-10-03
      if($date =~ s/\s*((\d{1,2})-(\w{3})-(\d{2,4})).*/$1/) {
	  if( length($date) < 11 ) { # improperly formatted date
	                             # But we'll be nice and fix it for them
	      my ($d,$m,$y) = ($2,$3,$4);
	      if( length($d) == 1 ) {
		  $d = "0$d";
	      }
	      # guess the century here
	      if( length($y) == 2 ) {
		  if( $y > 60 ) {  # arbitrarily guess that '60' means 1960
		      $y = "19$y";
		  } else { 
		      $y = "20$y";
		  }
		  $self->warn("Date was malformed, guessing the century for $date to be $y\n");
	      }
	      $params{'-dates'} = [join('-',$d,$m,$y)];
	  } else { 
	      $params{'-dates'} = [$date];
	  }
      }
      # set them all at once
      $builder->add_slot_value(%params);
      %params = ();

      # parse the rest if desired, otherwise start over
      if(! $builder->want_object()) {
	  $builder->make_object();
	  next RECORDSTART;
      }
      
      # set up annotation depending on what the builder wants
      if($builder->want_slot('annotation')) {
	  $annotation = new Bio::Annotation::Collection;
      }
      $buffer = $self->_readline();
      until( !defined ($buffer) ) {
	  $_ = $buffer;
	  
	  # Description line(s)
	  if (/^DEFINITION\s+(\S.*\S)/) {
	      my @desc = ($1);
	      while ( defined($_ = $self->_readline) ) { 
		  if( /^\s+(.*)/ ) { push (@desc, $1); next };		  
		  last;
	      }
	      $builder->add_slot_value(-desc => join(' ', @desc));
	      # we'll continue right here because DEFINITION always comes
	      # at the top of the entry
	  }
	  # accession number (there can be multiple accessions)
	  if( /^ACCESSION\s+(\S.*\S)/ ) {
	      push(@acc, split(/\s+/,$1));
	      while( defined($_ = $self->_readline) ) { 
		  /^\s+(.*)/ && do { push (@acc, split(/\s+/,$1)); next };
		  last;
	      }
	      $buffer = $_;
	      next;
	  }
	  # PID
	  elsif( /^PID\s+(\S+)/ ) {
	      $params{'-pid'} = $1;
	  }
	  #Version number
	  elsif( /^VERSION\s+(.+)$/ ) {
	      my ($acc,$gi) = split(' ',$1);
	      if($acc =~ /^\w+\.(\d+)/) {
		  $params{'-version'} = $1;
		  $params{'-seq_version'} = $1;
	      }
	      if($gi && (index($gi,"GI:") == 0)) {
		  $params{'-primary_id'} = substr($gi,3);
	      }
	  }
	  #Keywords
	  elsif( /^KEYWORDS\s+(.*)/ ) {
	      my @kw = split(/\s*\;\s*/,$1);
	      while( defined($_ = $self->_readline) ) { 
		  chomp;
		  /^\s+(.*)/ && do { push (@kw, split(/\s*\;\s*/,$1)); next };
		  last;
	      }
	      
	      @kw && $kw[-1] =~ s/\.$//;
	      $params{'-keywords'} = \@kw;
	      $buffer = $_;
	      next;
	  }
	  # Organism name and phylogenetic information
	  elsif (/^SOURCE/) {
	      if($builder->want_slot('species')) {
		  $species = $self->_read_GenBank_Species(\$buffer);
		  $builder->add_slot_value(-species => $species);
	      } else {
		  while(defined($buffer = $self->_readline())) {
		      last if substr($buffer,0,1) ne ' ';
		  }
	      }
	      next;
	  }
	  #References
	  elsif (/^REFERENCE/) {
	      if($annotation) {
		  my @refs = $self->_read_GenBank_References(\$buffer);
		  foreach my $ref ( @refs ) {
		      $annotation->add_Annotation('reference',$ref);
		  }
	      } else {
		  while(defined($buffer = $self->_readline())) {
		      last if substr($buffer,0,1) ne ' ';
		  }
	      }
	      next;
	  }
	  #Comments
	  elsif (/^COMMENT\s+(.*)/) {
	      if($annotation) {
		  my $comment = $1;
		  while (defined($_ = $self->_readline)) {
		      last if (/^\S/);
		      $comment .= $_; 
		  }
		  $comment =~ s/\n/ /g;
		  $comment =~ s/  +/ /g;
		  $annotation->add_Annotation(
			    'comment',
			    Bio::Annotation::Comment->new(-text => $comment));
		  $buffer = $_;
	      } else {
		  while(defined($buffer = $self->_readline())) {
		      last if substr($buffer,0,1) ne ' ';
		  }
	      }
	      next;
	  } elsif( /^SEGMENT\s+(.+)/ ) {
	      if($annotation) {
		  my $segment = $1;
		  while (defined($_ = $self->_readline)) {
		      last if (/^\S/);
		      $segment .= $_; 
		  }
		  $segment =~ s/\n/ /og;
		  $segment =~ s/ {2,}/ /og;
		  $annotation->add_Annotation(
					      'segment',
					      Bio::Annotation::SimpleValue->new(-value => $segment));
		  $buffer = $_;
	      } else {
		  while(defined($buffer = $self->_readline())) {
		      last if substr($buffer,0,1) ne ' ';
		  }
	      }
	      next;
	  }
	  # Exit at start of Feature table, or start of sequence
	  last if( /^(FEATURES|ORIGIN)/ );
	  # Get next line and loop again
	  $buffer = $self->_readline;
      }
      return undef if(! defined($buffer));

      # add them all at once for efficiency
      $builder->add_slot_value(-accession_number => shift(@acc),
			       -secondary_accessions => \@acc,
			       %params);
      $builder->add_slot_value(-annotation => $annotation) if $annotation;
      %params = (); # reset before possible re-use to avoid setting twice

      # start over if we don't want to continue with this entry
      if(! $builder->want_object()) {
	  $builder->make_object();
	  next RECORDSTART;
      }      
      
      # some "minimal" formats may not necessarily have a feature table
      if($builder->want_slot('features') && defined($_) && /^FEATURES/o) {
	  # need to read the first line of the feature table
	  $buffer = $self->_readline;
	  
	  # DO NOT read lines in the while condition -- this is done as a side
	  # effect in _read_FTHelper_GenBank!
	  while( defined($buffer) ) {
	      # check immediately -- not at the end of the loop
	      # note: GenPept entries obviously do not have a BASE line
	      last if(($buffer =~ /^BASE/o) || ($buffer =~ /^ORIGIN/o) ||
		      ($buffer =~ /^CONTIG/o) );

	      # slurp in one feature at a time -- at return, the start of
	      # the next feature will have been read already, so we need
	      # to pass a reference, and the called method must set this
	      # to the last line read before returning 
	      
	      my $ftunit = $self->_read_FTHelper_GenBank(\$buffer);
	      
	      # fix suggested by James Diggans

	      if( !defined $ftunit ) {
		  # GRRRR. We have fallen over. Try to recover
		  $self->warn("Unexpected error in feature table for ".$params{'-display_id'}." Skipping feature, attempting to recover");
		  unless( ($buffer =~ /^\s{5,5}\S+/o) or 
			  ($buffer =~ /^\S+/o)) {
		      $buffer = $self->_readline();
		  }
		  next; # back to reading FTHelpers
	      }
		
	      # process ftunit
	      my $feat =
		  $ftunit->_generic_seqfeature($self->location_factory(),
					       $display_id);
	      # add taxon_id from source if available
	      if($species && ($feat->primary_tag eq 'source') &&
		 $feat->has_tag('db_xref') && (! $species->ncbi_taxid())) {
		  foreach my $tagval ($feat->get_tag_values('db_xref')) {
		      if(index($tagval,"taxon:") == 0) {
			  $species->ncbi_taxid(substr($tagval,6));
                          last;
		      }
		  }
	      }
	      # add feature to list of features
	      push(@features, $feat);
	  }
	  $builder->add_slot_value(-features => \@features);
	  $_ = $buffer;
      }
      if( defined ($_) ) {
	  if( /^CONTIG/o && $builder->want_slot('features')) {
	      $b = "     $_"; # need 5 spaces to treat it like a feature
	      my $ftunit = $self->_read_FTHelper_GenBank(\$b);
	      if( ! defined $ftunit ) {
		  $self->warn("unable to parse the CONTIG feature\n");
	      } else { 
		  push(@features,
		       $ftunit->_generic_seqfeature($self->location_factory(),
						    $display_id));
	      }	
	  } elsif(! /^(ORIGIN|\/\/)/ ) {    # advance to the sequence, if any
	      while (defined( $_ = $self->_readline) ) {
		  last if /^(ORIGIN|\/\/)/;
	      }
	  }
      }
      if(! $builder->want_object()) {
	  $builder->make_object(); # implicit end-of-object
	  next RECORDSTART;
      }
      if($builder->want_slot('seq')) {
	  # the fact that we want a sequence does not necessarily mean that
	  # there also is a sequence ...
	  if(defined($_) && s/^ORIGIN\s+//) {
	      chomp;
	      if( $annotation && length($_) > 0 ) {
		  $annotation->add_Annotation('origin',
					      Bio::Annotation::SimpleValue->new(-value => $_));
	      }
	      my $seqc = '';
	      while( defined($_ = $self->_readline) ) {
		  /^\/\// && last;
		  $_ = uc($_);
		  s/[^A-Za-z]//g;
		  $seqc .= $_;
	      }
	      $self->debug("sequence length is ". length($seqc) ."\n");
	      $builder->add_slot_value(-seq => $seqc);
	  }
      } elsif ( defined($_) && (substr($_,0,2) ne '//')) {
	  # advance to the end of the record
	  while( defined($_ = $self->_readline) ) {
	      last if substr($_,0,2) eq '//';
	  }
      }
      # Unlikely, but maybe the sequence is so weird that we don't want it
      # anymore. We don't want to return undef if the stream's not exhausted
      # yet.
      $seq = $builder->make_object();
      next RECORDSTART unless $seq;
      last RECORDSTART;
  } # end while RECORDSTART

    return $seq;
}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq($seq)
 Function: writes the $seq object (must be seq) to the stream
 Returns : 1 for success and 0 for error
 Args    : array of 1 to n Bio::SeqI objects


=cut

sub write_seq {
    my ($self,@seqs) = @_;

    foreach my $seq ( @seqs ) {
	$self->throw("Attempting to write with no seq!") unless defined $seq;

	if( ! ref $seq || ! $seq->isa('Bio::SeqI') ) {
	    $self->warn(" $seq is not a SeqI compliant module. Attempting to dump, but may fail!");
	}

	my $str = $seq->seq;

	my ($div, $mol);
	my $len = $seq->length();

	if ( $seq->can('division') ) {
	    $div=$seq->division;
	} 
	if( !defined $div || ! $div ) { $div = 'UNK'; }
	my $alpha = $seq->alphabet;
	if( !$seq->can('molecule') || ! defined ($mol = $seq->molecule()) ) {
	    $mol =  $alpha || 'DNA';
	}

	my $circular = 'linear  ';
	$circular = 'circular' if $seq->is_circular;

	local($^W) = 0;	# supressing warnings about uninitialized fields.

	my $temp_line;
	if( $self->_id_generation_func ) {
	    $temp_line = &{$self->_id_generation_func}($seq);
	} else {
	    my $date = '';
	    if( $seq->can('get_dates') ) { 	    
		($date) = $seq->get_dates();
	    }

            $self->warn("No whitespace allowed in GenBank display id [". $seq->display_id. "]")
                if $seq->display_id =~ /\s/;

	    $temp_line = sprintf ("%-12s%-15s%13s %s%4s%-8s%-8s %3s %-s", 
				  'LOCUS', $seq->id(),$len,
				  (lc($alpha) eq 'protein') ? ('aa','', '') : 
				  ('bp', '',$mol),$circular,
				  $div,$date);
	} 

	$self->_print("$temp_line\n");
	$self->_write_line_GenBank_regex("DEFINITION  ", "            ",
					 $seq->desc(),"\\s\+\|\$",80);

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
			      ($id && ($id =~ /^\d+$/) ? "  GI:".$id : ""),
			      "\n");
	    }
	} 

	# if there, write the keywords line

	if( defined $self->_kw_generation_func() ) {
	    $temp_line = &{$self->_kw_generation_func}($seq);
	    $self->_print("KEYWORDS    $temp_line\n");   
	} else {
	    if( $seq->can('keywords') ) {
		my $kw = $seq->keywords;
		$kw .= '.' if( $kw !~ /\.$/ );
		$self->_print("KEYWORDS    $kw\n");
	    }
	} 

	# SEGMENT if it exists
	foreach my $ref ( $seq->annotation->get_Annotations('segment') ) {
	    $self->_print(sprintf ("%-11s %s\n",'SEGMENT',
				  $ref->value));
	}

	# Organism lines
	if (my $spec = $seq->species) {
	    my ($species, $genus, @class) = $spec->classification();
	    my $OS;
	    if( $spec->common_name ) {
		$OS = $spec->common_name;
	    } else { 
		$OS = "$genus $species";
	    }
	    if (my $ssp = $spec->sub_species) {
		$OS .= " $ssp";
	    }
	    $self->_print("SOURCE      $OS\n");
	    $self->_print("  ORGANISM  ",
			  ($spec->organelle() ? $spec->organelle()." " : ""),
			  "$genus $species", "\n");
	    my $OC = join('; ', (reverse(@class), $genus)) .'.';
	    $self->_write_line_GenBank_regex(' 'x12,' 'x12,
					     $OC,"\\s\+\|\$",80);
	}

	# Reference lines
	my $count = 1;
	foreach my $ref ( $seq->annotation->get_Annotations('reference') ) {
            $temp_line = "REFERENCE   $count";
	    $temp_line .= sprintf ("  (%s %d to %d)",
				  ($seq->alphabet() eq "protein" ?
				   "residues" : "bases"),
				  $ref->start,$ref->end)
                if $ref->start;
	    $self->_print("$temp_line\n");
	    $self->_write_line_GenBank_regex("  AUTHORS   ",' 'x12,
					     $ref->authors,"\\s\+\|\$",80);
	    $self->_write_line_GenBank_regex("  TITLE     "," "x12,
					     $ref->title,"\\s\+\|\$",80);
	    $self->_write_line_GenBank_regex("  JOURNAL   "," "x12,
					     $ref->location,"\\s\+\|\$",80);
	    if ($ref->comment) {
		$self->_write_line_GenBank_regex("  REMARK    "," "x12,
						 $ref->comment,"\\s\+\|\$",80);
	    }
	    if( $ref->medline) {
		$self->_write_line_GenBank_regex("  MEDLINE   "," "x12,
						 $ref->medline, "\\s\+\|\$",80);
		# I am assuming that pubmed entries only exist when there
		# are also MEDLINE entries due to the indentation
	    }
	    # This could be a wrong assumption
	    if( $ref->pubmed ) {
		$self->_write_line_GenBank_regex("   PUBMED   "," "x12,
						 $ref->pubmed, "\\s\+\|\$",
						 80);
	    }
	    $count++;
	}
	# Comment lines

	foreach my $comment ( $seq->annotation->get_Annotations('comment') ) {
	    $self->_write_line_GenBank_regex("COMMENT     "," "x12,
					     $comment->text,"\\s\+\|\$",80);
	}
	$self->_print("FEATURES             Location/Qualifiers\n");

	my $contig;
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
	if( $seq->length == 0 ) { $self->_show_dna(0) }

	if( $self->_show_dna() == 0 ) {
	    $self->_print("\n//\n");
	    return;
	}

# finished printing features.

	$str =~ tr/A-Z/a-z/;

# Count each nucleotide
	unless(  $mol eq 'protein' ) {
	    my $alen = $str =~ tr/a/a/;
	    my $clen = $str =~ tr/c/c/;
	    my $glen = $str =~ tr/g/g/;
	    my $tlen = $str =~ tr/t/t/;

	    my $olen = $len - ($alen + $tlen + $clen + $glen);
	    if( $olen < 0 ) {
		$self->warn("Weird. More atgc than bases. Problem!");
	    }

	    my $base_count = sprintf("BASE COUNT %8s a %6s c %6s g %6s t%s\n",
				     $alen,$clen,$glen,$tlen,
				     ( $olen > 0 ) ? 
				     sprintf("%6s others",$olen) : '');
	    $self->_print($base_count); 
	}

	my ($o) = $seq->annotation->get_Annotations('origin');
	$self->_print(sprintf("%-12s%s\n",
			      'ORIGIN', $o ? $o->value : ''));
        # print out the sequence
	my $nuc = 60;		# Number of nucleotides per line
	my $whole_pat = 'a10' x 6; # Pattern for unpacking a whole line
	my $out_pat   = 'A11' x 6; # Pattern for packing a line
	my $length = length($str);

	# Calculate the number of nucleotides which fit on whole lines
	my $whole = int($length / $nuc) * $nuc;

	# Print the whole lines
	my $i;
	for ($i = 0; $i < $whole; $i += $nuc) {
	    my $blocks = pack $out_pat,
	    unpack $whole_pat,
	    substr($str, $i, $nuc);
            chop $blocks;
	    $self->_print(sprintf("%9d $blocks\n", $i + $nuc - 59));
	}

	# Print the last line
	if (my $last = substr($str, $i)) {
	    my $last_len = length($last);
	    my $last_pat = 'a10' x int($last_len / 10) .'a'. $last_len % 10;
	    my $blocks = pack $out_pat,
	    unpack($last_pat, $last);
            $blocks =~ s/ +$//;
	    $self->_print(sprintf("%9d $blocks\n", $length - $last_len + 1));
	}

	$self->_print("//\n");

	$self->flush if $self->_flush_on_write && defined $self->_fh;
	return 1;
    }
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
   my ($self,$fth) = @_;
   
   if( ! ref $fth || ! $fth->isa('Bio::SeqIO::FTHelper') ) {
       $fth->warn("$fth is not a FTHelper class. Attempting to print, but there could be tears!");   
   }
   if( defined $fth->key && 
       $fth->key eq 'CONTIG' ) {
       $self->_write_line_GenBank_regex(sprintf("%-12s",$fth->key),
					' 'x12,$fth->loc,"\,\|\$",80);
   } else {
       $self->_write_line_GenBank_regex(sprintf("     %-16s",$fth->key),
					" "x21,
					$fth->loc,"\,\|\$",80);
   }

   foreach my $tag ( keys %{$fth->field} ) {
       foreach my $value ( @{$fth->field->{$tag}} ) {
	   $value =~ s/\"/\"\"/g;
	   if ($value eq "_no_value") {
	       $self->_write_line_GenBank_regex(" "x21,
						" "x21,
						"/$tag","\.\|\$",80);
	   }
           # there are almost 3x more quoted qualifier values and they
           # are more common too so we take quoted ones first
           elsif (!$FTQUAL_NO_QUOTE{$tag}) {
              my ($pat) = ($value =~ /\s/ ? '\s|$' : '.|$');
	      $self->_write_line_GenBank_regex(" "x21,
					       " "x21,
					       "/$tag=\"$value\"",$pat,80);

           } else {
              $self->_write_line_GenBank_regex(" "x21,
					       " "x21,
					       "/$tag=$value","\.\|\$",80);
           }
       }
   }

}


=head2 _read_GenBank_References

 Title   : _read_GenBank_References
 Usage   :
 Function: Reads references from GenBank format. Internal function really
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

   my (@title,@loc,@authors,@com,@medline,@pubmed);

   REFLOOP: while( defined($_) || defined($_ = $self->_readline) ) {
       if (/^\s{2}AUTHORS\s+(.*)/o) { 
	   push (@authors, $1);   
	   while ( defined($_ = $self->_readline) ) {
	       /^\s{3,}(.*)/o && do { push (@authors, $1);next;};
	       last;
	   }
	   $ref->authors(join(' ', @authors));
       }
       if (/^\s{2}TITLE\s+(.*)/o)  { 
	   push (@title, $1);
	   while ( defined($_ = $self->_readline) ) {
	       /^\s{3,}(.*)/o && do { push (@title, $1);
				     next;
				 };
	       last;
	   }
	   $ref->title(join(' ', @title));
       }
       if (/^\s{2}JOURNAL\s+(.*)/o) { 
	   push(@loc, $1);
	   while ( defined($_ = $self->_readline) ) {
	       # we only match when there are at least 4 spaces
	       # there is probably a better way to match this
	       # as it assumes that the describing tag is short enough 
	       /^\s{4,}(.*)/o && do { push(@loc, $1);
				      next;
				 };
	       last;
	   }
	   $ref->location(join(' ', @loc));
	   redo REFLOOP;
       }
       if (/^\s{2}REMARK\s+(.*)/o) { 
	   push (@com, $1);
	   while ( defined($_ = $self->_readline) ) {	       
	       /^\s{3,}(.*)/o && do { push(@com, $1);
				     next;
				 };
	       last;
	   }
	   $ref->comment(join(' ', @com));
	   redo REFLOOP;
       }
       if( /^\s{2}MEDLINE\s+(.*)/ ) {
	   push(@medline,$1);
	   while ( defined($_ = $self->_readline) ) {	       
	       /^\s{4,}(.*)/ && do { push(@medline, $1);
				     next;
				 };
	       last;
	   }
	   $ref->medline(join(' ', @medline));
	   redo REFLOOP;
       }
       if( /^\s{3}PUBMED\s+(.*)/ ) {
	   push(@pubmed,$1);
	   while ( defined($_ = $self->_readline) ) {	       
	       /^\s{5,}(.*)/ && do { push(@pubmed, $1);
				     next;
				 };
	       last;
	   }
	   $ref->pubmed(join(' ', @pubmed));
	   redo REFLOOP;
       }
       
       /^REFERENCE/o && do {
	   # store current reference
	   $self->_add_ref_to_array(\@refs,$ref) if $ref;
	   # reset
	   @authors = ();
	   @title = ();
	   @loc = ();
	   @com = ();
	   @pubmed = ();
	   @medline = ();
	   # create the new reference object
	   $ref = Bio::Annotation::Reference->new();
	   # check whether start and end base is given
	   if (/^REFERENCE\s+\d+\s+\([a-z]+ (\d+) to (\d+)/){
	       $ref->start($1);
	       $ref->end($2);
	   }
       };

       /^(FEATURES)|(COMMENT)/o && last;

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
 Args    : a reference to the current line buffer

=cut

sub _read_GenBank_Species {
    my( $self,$buffer) = @_;
    my @organell_names = ("chloroplast", "mitochondr"); 
    # only those carrying DNA, apart from the nucleus

    $_ = $$buffer;
    
    my( $sub_species, $species, $genus, $common, $organelle, @class, $ns_name );
    # upon first entering the loop, we must not read a new line -- the SOURCE
    # line is already in the buffer (HL 05/10/2000)
    while (defined($_) || defined($_ = $self->_readline())) {
	# de-HTMLify (links that may be encountered here don't contain
	# escaped '>', so a simple-minded approach suffices)
        s/<[^>]+>//g;
	if (/^SOURCE\s+(.*)/o) {
	    # FIXME this is probably mostly wrong (e.g., it yields things like
	    # Homo sapiens adult placenta cDNA to mRNA
	    # which is certainly not what you want)
	    $common = $1;
	    $common =~ s/\.$//; # remove trailing dot
	} elsif (/^\s{2}ORGANISM/o) {
	    my @spflds = split(' ', $_);
            ($ns_name) = $_ =~ /\w+\s+(.*)/o;
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
        } elsif (/^\s+(.+)/o) {
	    # only split on ';' or '.' so that 
	    # classification that is 2 words will 
	    # still get matched
	    # use map to remove trailing/leading spaces
            push(@class, map { s/^\s+//; s/\s+$//; $_; } split /[;\.]+/, $1);
        } else {
            last;
        }
        
        $_ = undef; # Empty $_ to trigger read of next line
    }
    
    $$buffer = $_;
    
    # Don't make a species object if it's empty or "Unknown" or "None"
    return unless $genus and  $genus !~ /^(Unknown|None)$/oi;
    
    # Bio::Species array needs array in Species -> Kingdom direction
    if ($class[0] eq 'Viruses') {
        push( @class, $ns_name );
    }
    elsif ($class[$#class] eq $genus) {
        push( @class, $species );
    } else {
        push( @class, $genus, $species );
    }
    @class = reverse @class;
    
    my $make = Bio::Species->new();
    $make->classification( \@class, "FORCE" ); # no name validation please
    $make->common_name( $common      ) if $common;
    unless ($class[-1] eq 'Viruses') {
        $make->sub_species( $sub_species ) if $sub_species;
    }
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
    
    if ($$buffer =~ /^\s{5}(\S+)\s+(.+?)\s*$/o) {
        $key = $1;
        $loc = $2;
        # Read all the lines up to the next feature
        while ( defined($_ = $self->_readline) ) {
            if (/^(\s+)(.+?)\s*$/o) {
                # Lines inside features are preceded by 21 spaces
                # A new feature is preceded by 5 spaces
                if (length($1) > 6) {
                    # Add to qualifiers if we're in the qualifiers, or if it's
		    # the first qualifier
                    if (@qual || (index($2,'/') == 0)) {
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
	$self->debug("no feature key!\n");
	# change suggested by JDiggans to avoid infinite loop- 
	# see bugreport 1062.
	# reset buffer to prevent infinite loop
	$$buffer = $self->_readline();
        return;
    } 
    
    # Put the first line of the next feature into the buffer
    $$buffer = $_;

    # Make the new FTHelper object
    my $out = new Bio::SeqIO::FTHelper();
    $out->verbose($self->verbose());
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
		       $self->warn("Unbalanced quote in:\n" .
				   join("\n", @qual) .
				   "No further qualifiers will " .
				   "be added for this feature");
		       last QUAL;
                    }
                    $i++; # modifying a for-loop variable inside of the loop
		          # is not the best programming style ...
                    my $next = $qual[$i];

                    # add to value with a space unless the value appears
		    # to be a sequence (translation for example)
		    if(($value.$next) =~ /[^A-Za-z\"\-]/o) {
			$value .= " ";
		    }
                    $value .= $next;
                }
                # Trim leading and trailing quotes
                $value =~ s/^"|"$//g;
                # Undouble internal quotes
                $value =~ s/""/\"/g;
            } elsif ( $value =~ /^\(/ ) { # values quoted by ()s
		# Keep addingto value until we find the trailing bracket
		# and the ()s are balanced
		my $left = ($value =~ tr/\(/\(/); # count left parens
		my $right = ($value =~ tr/\)/\)/); # count right parens
		while( $value !~ /\)$/ or $left != $right ) {
		    if( $i >= $#qual) {
			$self->warn("Unbalanced parens in:\n".
				    join("\n", @qual).
				    "No further qualifiers will ".
				    "be added for this feature");
			last QUAL;
		    }
		    $i++;
		    my $next = $qual[$i];		    
		    $value .= $next;
		    $left += ($next =~ tr/\(/\(/);
		    $right += ($next =~ tr/\)/\)/);
		}
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

   $length || $self->throw("Miscalled write_line_GenBank without length. Programming error!");
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

   $length || $self->throw( "Miscalled write_line_GenBank without length. Programming error!");

#   if( length $pre1 != length $pre2 ) {
#       $self->throw( "Programming error - cannot called write_line_GenBank_regex with different length pre1 and pre2 tags!");
#   }

   my $subl = $length - (length $pre1) - 2;
   my @lines = ();

   CHUNK: while($line) {
       foreach my $pat ($regex, '[,;\.\/-]\s|'.$regex, '[,;\.\/-]|'.$regex) {
	   if($line =~ m/^(.{1,$subl})($pat)(.*)/) {
	       $line = $3;
	       # be strict about not padding spaces according to 
	       # genbank format
	       my $l = $1.$2;
	       $l =~ s/\s+$//;
	       push(@lines, $l);
	       next CHUNK;
	   }
       }
       # if we get here none of the patterns matched $subl or less chars
       $self->warn("trouble dissecting \"$line\" into chunks ".
		   "of $subl chars or less - this tag won't print right");
       # insert a space char to prevent infinite loops
       $line = substr($line,0,$subl) . " " . substr($line,$subl);
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
   my ($obj,$value) = @_;
   if( defined $value) {
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
   my ($obj,$value) = @_;
   if( defined $value) {
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
   my ($obj,$value) = @_;
   if( defined $value ) {
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
   my ($obj,$value) = @_;
   if( defined $value ) {
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
   my ($obj,$value) = @_;
   if( defined $value ) {      
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
   my ($obj,$value) = @_;
   if( defined $value ) {
      $obj->{'_kw_generation_func'} = $value;
    }
    return $obj->{'_kw_generation_func'};
}

1;
