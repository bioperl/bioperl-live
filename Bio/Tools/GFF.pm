package Bio::Tools::GFF;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::Root::RootI Bio::SeqAnalysisParserI);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($file, $fh, $gff_version) = $self->_rearrange([qw(FILE FH GFF_VERSION)],@args);
  
  if( defined $file && defined $fh ) {
      $self->throw("Cannot define both a file and fh for input");
  }
  if( defined $file ) {
      $fh = Symbol::gensym();
      open ($fh,$file) || $self->throw("Could not open file [$file] $!");
  } elsif( defined $fh ) {
      if (ref $fh !~ /GLOB/)
      { $self->throw("Expecting a GLOB reference, not $fh!"); }
  }
  $self->fh($fh);
  
  $gff_version ||= 2;
  if(($gff_version != 1) && ($gff_version != 2)) {
    $self->throw("Can't build a GFF object with the unknown version".
		 $gff_version);
  }
  $self->gff_version($gff_version);
  return $self;
}

sub next_feature {
  my ($self) = @_;
  
  my $gff_string = $self->_readLine();
  my $feat = Bio::SeqFeature::Generic->new();
  
  if($self->gff_version() == 1)  {
    $self->_from_gff_string($feat, $gff_string);
  } else {
    $self->_from_gff2_string($feat, $gff_string);
  }
  
  return $feat;
}

=head2 _from_gff_string

 Title   : _from_gff_string
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub _from_gff_string {
   my ($gff, $feat, $string) = @_;

   my ($seqname, $source, $primary, $start, $end, $score, $strand, $frame, @group) = split(/\s+/, $string);

   if ( !defined $frame ) {
       $feat->throw("[$string] does not look like GFF to me");
   }
   $frame = 0 unless( $frame =~ /^\d+$/);
   $feat->seqname($seqname);
   $feat->source_tag($source);
   $feat->primary_tag($primary);
   $feat->start($start);
   $feat->end($end);
   $feat->frame($frame);
   if ( $score eq '.' ) {
       #$feat->score(undef);
   } else {
       $feat->score($score);
   }
   if ( $strand eq '-' ) { $feat->strand(-1); }
   if ( $strand eq '+' ) { $feat->strand(1); }
   if ( $strand eq '.' ) { $feat->strand(0); }
   foreach my $g ( @group ) {
       if ( $g =~ /(\S+)=(\S+)/ ) {
	   my $tag = $1;
	   my $value = $2;
	   $feat->add_tag_value($1, $2);
       } else {
	   $feat->add_tag_value('group', $g);
       }
   }
}

=head2 _from_gff2_string

 Title   : _from_gff2_string
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub _from_gff2_string {
   my ($gff, $feat, $string) = @_;
   # according to the Sanger website, GFF2 should be single-tab separated elements, and the
   # free-text at the end should contain text-translated tab symbols but no "real" tabs,
   # so splitting on \t is safe, and $attribs gets the entire attributes field to be parsed later
   my ($seqname, $source, $primary, $start, $end, $score, $strand, $frame, @attribs) = split(/\t+/, $string);
   my $attribs = join '', @attribs;  # just in case the rule against tab characters has been broken
   if ( !defined $frame ) {
       $feat->throw("[$string] does not look like GFF2 to me");
   }
   $feat->seqname($seqname);
   $feat->source_tag($source);
   $feat->primary_tag($primary);
   $feat->start($start);
   $feat->end($end);
   $feat->frame($frame);
   if ( $score eq '.' ) {
       #$feat->score(undef);
   } else {
       $feat->score($score);
   }
   if ( $strand eq '-' ) { $feat->strand(-1); }
   if ( $strand eq '+' ) { $feat->strand(1); }
   if ( $strand eq '.' ) { $feat->strand(0); }

   $attribs =~ s/\#(.*)$//;				 # remove comments field (format:  #blah blah blah...  at the end of the GFF line)
   my @key_vals = split /;/, $attribs;   # attributes are semicolon-delimited

   foreach my $pair ( @key_vals ) {
       my ($blank, $key, $values) = split  /^\s*([\w\d]+)\s/, $pair;	# separate the key from the value based on the = sign

       my @values;								

       while ($values =~ s/"(.*?)"//){          # free text is quoted, so match each free-text block
       		if ($1){push @values, $1};           # and push it on to the list of values (tags may have more than one value...)
       }

       my @othervals = split /\s+/, $values;  # and what is left over should be space-separated non-free-text values
       foreach my $othervalue(@othervals){
       		if (CORE::length($othervalue) > 0){push @values, $othervalue}  # get rid of any empty strings which might result from the split
       }

       foreach my $value(@values){
       	   	$feat->add_tag_value($key, $value);
      }
   }
}

sub write_feature {
  my ($self, $feature) = @_;
  
  if($self->gff_version() == 1) {
    $self->_print($self->_gff_string($feature)."\n");
  } else {
    $self->_print($self->_gff2_string($feature)."\n");
  }
}

sub _gff_string{
   my ($gff, $feat) = @_;
   my ($str,$score,$frame,$name,$strand);

   if( $feat->can('score') ) {
       $score = $feat->score();
   }
   $score = '.' unless defined $score;

   if( $feat->can('frame') ) {
       $frame = $feat->frame();
   }
   $frame = '.' unless defined $frame;

   $strand = $feat->strand();
   if(! $strand) {
       $strand = ".";
   } elsif( $strand == 1 ) {
       $strand = '+';
   } elsif ( $feat->strand == -1 ) {
       $strand = '-';
   }
   
   if( $feat->can('seqname') ) {
       $name = $feat->seqname();
       $name ||= 'SEQ';
   } else {
       $name = 'SEQ';
   }


   $str = join("\t",
                 $name,
		 $feat->source_tag(),
		 $feat->primary_tag(),
		 $feat->start(),
		 $feat->end(),
		 $score,
		 $strand,
		 $frame);

   foreach my $tag ( $feat->all_tags ) {
       foreach my $value ( $feat->each_tag_value($tag) ) {
	   $str .= " $tag=$value";
       }
   }


   return $str;
}

sub _gff2_string{
   my ($gff, $feat) = @_;
   my ($str,$score,$frame,$name,$strand);

   if( $feat->can('score') ) {
       $score = $feat->score();
   }
   $score = '.' unless defined $score;

   if( $feat->can('frame') ) {
       $frame = $feat->frame();
   }
   $frame = '.' unless defined $frame;

   $strand = $feat->strand();
   if(! $strand) {
       $strand = ".";
   } elsif( $strand == 1 ) {
       $strand = '+';
   } elsif ( $feat->strand == -1 ) {
       $strand = '-';
   }

   if( $feat->can('seqname') ) {
       $name = $feat->seqname();
       $name ||= 'SEQ';
   } else {
       $name = 'SEQ';
   }


   $str = join("\t",
                 $name,
		 $feat->source_tag(),
		 $feat->primary_tag(),
		 $feat->start(),
		 $feat->end(),
		 $score,
		 $strand,
		 $frame);

   # the routine below is the only modification I made to the original
   # ->gff_string routine (above) as on November 17th, 2000, the
   # Sanger webpage describing GFF2 format reads: "From version 2
   # onwards, the attribute field must have a tag value structure
   # following the syntax used within objects in a .ace file,
   # flattened onto one line by semicolon separators. Tags must be
   # standard identifiers ([A-Za-z][A-Za-z0-9_]*).  Free text values
   # must be quoted with double quotes".

   # MW

   my $valuestr;
   if ($feat->all_tags){  # only play this game if it is worth playing...
        $str .= "\t";     # my interpretation of the GFF2 specification suggests the need for this additional TAB character...??
        foreach my $tag ( $feat->all_tags ) {
            my $valuestr; # a string which will hold one or more values for this tag, with quoted free text and space-separated individual values.
            foreach my $value ( $feat->each_tag_value($tag) ) {
         		if ($value =~ /[^A-Za-z0-9_]/){
         			$value =~ s/\t/\\t/g;         # substitute tab and newline characters
         			$value =~ s/\n/\\n/g;          # to their UNIX equivalents
         			$value = '"' . $value . '" '}  # if the value contains anything other than valid tag/value characters, then quote it
         		$valuestr .= $value;								# with a trailing space in case there are multiple values
         															# for this tag (allowed in GFF2 and .ace format)		
            }
            $str .= "$tag $valuestr ; ";                              # semicolon delimited with no '=' sign
        }
   		chop $str; chop $str  # remove the trailing semicolon and space
    }
   return $str;
}

sub fh {
  my ($self, $value) = @_;
  if(defined $value && ref($value) =~ /GLOB/i ) {
    $self->{'FH'} = $value;
  } 
  return $self->{'FH'};
}

sub gff_version {
  my ($self, $value) = @_;
  if(defined $value && (($value == 1) || ($value == 2))) {
    $self->{'GFF_VERSION'} = $value;
  }
  return $self->{'GFF_VERSION'};
}

sub _print {
  my ($self, $str) = @_;
  
  my $fh = $self->fh() || \*STDOUT;
  print $fh $str;
}

sub _readLine {
  my ($self) = @_;
  
  my $fh = $self->fh() || \*STDIN;
  my $line = <$fh>;
  
  return $line;
}
