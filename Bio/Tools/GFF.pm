# $Id$
#
# BioPerl module for Bio::Tools::GFF
#
# Cared for by the Bioperl core team
#
# Copyright Matthew Pocock
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::GFF - A Bio::SeqAnalysisParserI compliant GFF format parser

=head1 SYNOPSIS

    use Bio::Tool::GFF;

    # specify input via -fh or -file
    my $gffio = Bio::Tools::GFF->new(-fh => \*STDIN, -gff_version => 2);
    my $feature;
    # loop over the input stream
    while($feature = $gffio->next_feature()) {
        # do something with feature
    }
    $gffio->close();

    # you can also obtain a GFF parser as a SeqAnalasisParserI in
    # HT analysis pipelines (see Bio::SeqAnalysisParserI and
    # Bio::Factory::SeqAnalysisParserFactory)
    my $factory = Bio::Factory::SeqAnalysisParserFactory->new();
    my $parser = $factory->get_parser(-input => \*STDIN, -method => "gff");
    while($feature = $parser->next_feature()) {
        # do something with feature
    }

=head1 DESCRIPTION

This class provides a simple GFF parser and writer. In the sense of a
SeqAnalysisParser, it parses an input file or stream into SeqFeatureI
objects, but is not in any way specific to a particular analysis
program and the output that program produces.

That is, if you can get your analysis program spit out GFF, here is
your result parser.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Matthew Pocock

Email mrp@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::GFF;

use vars qw(@ISA);
use strict;

use Bio::Root::IO;
use Bio::SeqAnalysisParserI;
use Bio::SeqFeature::Generic;

@ISA = qw(Bio::Root::Root Bio::SeqAnalysisParserI Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : 
 Function: Creates a new instance. Recognized named parameters are -file, -fh,
           and -gff_version.

 Returns : a new object
 Args    : names parameters


=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($gff_version) = $self->_rearrange([qw(GFF_VERSION)],@args);

  # initialize IO
  $self->_initialize_io(@args);
    
  $gff_version ||= 2;
  if(($gff_version != 1) && ($gff_version != 2)) {
    $self->throw("Can't build a GFF object with the unknown version ".
		 $gff_version);
  }
  $self->gff_version($gff_version);
  return $self;
}

=head2 next_feature

 Title   : next_feature
 Usage   : $seqfeature = $gffio->next_feature();
 Function: Returns the next feature available in the input file or stream, or
           undef if there are no more features.
 Example :
 Returns : A Bio::SeqFeatureI implementing object, or undef if there are no
           more features.
 Args    : none    

=cut

sub next_feature {
    my ($self) = @_;
    
    my $gff_string;

    # be graceful about empty lines or comments, and make sure we return undef
    # if the input's consumed
    while(($gff_string = $self->_readline()) && defined($gff_string)) {
	next if($gff_string =~ /^\#/);
	next if($gff_string =~ /^\s*$/);
	last;
    }
    return undef unless $gff_string;

    my $feat = Bio::SeqFeature::Generic->new();
    $self->from_gff_string($feat, $gff_string);

    return $feat;
}

=head2 from_gff_string

 Title   : from_gff_string
 Usage   : $gff->from_gff_string($feature, $gff_string);
 Function: Sets properties of a SeqFeatureI object from a GFF-formatted
           string. Interpretation of the string depends on the version
           that has been specified at initialization.

           This method is used by next_feature(). It actually dispatches to
           one of the version-specific (private) methods.
 Example :
 Returns : void
 Args    : A Bio::SeqFeatureI implementing object to be initialized
           The GFF-formatted string to initialize it from

=cut

sub from_gff_string {
    my ($self, $feat, $gff_string) = @_;

    if($self->gff_version() == 1)  {
	$self->_from_gff1_string($feat, $gff_string);
    } else {
	$self->_from_gff2_string($feat, $gff_string);
    }
}

=head2 _from_gff1_string

 Title   : _from_gff1_string
 Usage   :
 Function:
 Example :
 Returns : void
 Args    : A Bio::SeqFeatureI implementing object to be initialized
           The GFF-formatted string to initialize it from

=cut

sub _from_gff1_string {
   my ($gff, $feat, $string) = @_;
   chomp $string;
   my ($seqname, $source, $primary, $start, $end, $score, $strand, $frame, @group) = split(/\t/, $string);

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
 Returns : void
 Args    : A Bio::SeqFeatureI implementing object to be initialized
           The GFF2-formatted string to initialize it from


=cut

sub _from_gff2_string {
   my ($gff, $feat, $string) = @_;
   chomp($string);
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

   #  <Begin Inefficient Code from Mark Wilkinson>
   # this routine is necessay to allow the presence of semicolons in quoted text
   # Semicolons are the delimiting character for new tag/value attributes.
   # it is more or less a "state" machine, with the "quoted" flag going up and down
   # as we pass thorugh quotes to distinguish free-text semicolon and hash symbols from GFF control characters

   #$attribs =~ s/\#(.*)$//;				 # remove comments field (format:  #blah blah blah...  at the end of the GFF line)
   my @att = split //, $attribs;         # split into individual characters
   my $num = $#att;                 # count them
   my $flag = 0;
	my @parsed;		# this is needed to hold the characters that have been parsed	
	for (my $a = 0; $a <= $num ; $a +=1){   # run through each character one at a time and check it
		if ($att[$a] eq "\""){$flag=($flag==0)?1:0}  # flag up on entering quoted text, down on leaving it
		if (($att[$a] eq ";") && $flag){$att[$a] = "INSERT_SEMICOLON_HERE"}  # replace semicolon with an unusual message if the quoted text flag is up
		if (($att[$a] eq "#") && !$flag){last}  # an unquoted hash symbol means the beginning of the comments field - discard
		push @parsed, $att[$a]                  # take the parsed character and push it onto the parsed list
	}

	$attribs = join "", @parsed; # rejoin into a single string

   # <End Inefficient Code>   Please feel free to fix this and make it more "perlish"

   my @key_vals = split /;/, $attribs;   # attributes are semicolon-delimited

   foreach my $pair ( @key_vals ) {
       $pair =~ s/INSERT_SEMICOLON_HERE/;/g;  # replace semicolons that were removed from free-text above.
       my ($blank, $key, $values) = split  /^\s*([\w\d]+)\s/, $pair;	# separate the key from the value

       my @values;								

       while ($values =~ s/"(.*?)"//){          # free text is quoted, so match each free-text block and remove it from the $values string
       		push @values, $1;          # and push it on to the list of values (tags may have more than one value... and the value may be undef)
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

=head2 write_feature

 Title   : write_feature
 Usage   : $gffio->write_feature($feature);
 Function: Writes the specified SeqFeatureI object in GFF format to the stream
           associated with this instance.
 Example :
 Returns : 
 Args    : A Bio::SeqFeatureI implementing object to be serialized

=cut

sub write_feature {
    my ($self, $feature) = @_;
    
    $self->_print($self->gff_string($feature)."\n");
}

=head2 gff_string

 Title   : gff_string
 Usage   : $gffstr = $gffio->gff_string($feature);
 Function: Obtain the GFF-formatted representation of a SeqFeatureI object.
           The formatting depends on the version specified at initialization.

           This method is used by write_feature(). It actually dispatches to
           one of the version-specific (private) methods.
 Example :
 Returns : A GFF-formatted string representation of the SeqFeature
 Args    : A Bio::SeqFeatureI implementing object to be GFF-stringified

=cut

sub gff_string{
    my ($self, $feature) = @_;

    if($self->gff_version() == 1) {
	return $self->_gff1_string($feature);
    } else {
	return $self->_gff2_string($feature);
    }
}

=head2 _gff1_string

 Title   : _gff1_string
 Usage   : $gffstr = $gffio->_gff1_string
 Function: 
 Example :
 Returns : A GFF1-formatted string representation of the SeqFeature
 Args    : A Bio::SeqFeatureI implementing object to be GFF-stringified

=cut

sub _gff1_string{
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

=head2 _gff2_string

 Title   : _gff2_string
 Usage   : $gffstr = $gffio->_gff2_string
 Function: 
 Example :
 Returns : A GFF2-formatted string representation of the SeqFeature
 Args    : A Bio::SeqFeatureI implementing object to be GFF2-stringified

=cut

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
         			$value =~ s/\t/\\t/g;          # substitute tab and newline characters
         			$value =~ s/\n/\\n/g;          # to their UNIX equivalents
         			$value = '"' . $value . '" '}  # if the value contains anything other than valid tag/value characters, then quote it
				$value = "\"\"" unless $value;  # if it is completely empty, then just make empty double quotes
         		$valuestr .=  $value . " ";								# with a trailing space in case there are multiple values
         															# for this tag (allowed in GFF2 and .ace format)		
            }
            $str .= "$tag $valuestr ; ";                              # semicolon delimited with no '=' sign
        }
   		chop $str; chop $str  # remove the trailing semicolon and space
    }
   return $str;
}

=head2 gff_version

 Title   : _gff_version
 Usage   : $gffversion = $gffio->gff_version
 Function: 
 Example :
 Returns : The GFF version this parser will accept and emit.
 Args    : none

=cut

sub gff_version {
  my ($self, $value) = @_;
  if(defined $value && (($value == 1) || ($value == 2))) {
    $self->{'GFF_VERSION'} = $value;
  }
  return $self->{'GFF_VERSION'};
}

1;

