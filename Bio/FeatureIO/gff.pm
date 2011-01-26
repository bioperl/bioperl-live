=pod

=head1 NAME

Bio::FeatureIO::gff - read/write GFF feature files

=head1 SYNOPSIS

  my $feature; #get a Bio::SeqFeature::Annotated somehow
  my $featureOut = Bio::FeatureIO->new(
    -format => 'gff',
    -version => 3,
    -fh => \*STDOUT,
    -validate_terms => 1, #boolean. validate ontology terms online?  default 0 (false).
  );
  $featureOut->write_feature($feature);

=head1 DESCRIPTION

 Currently implemented:

 version         read?   write?
 ------------------------------
 GFF 1             N       N
 GFF 2             N       N
 GFF 2.5 (GTF)     N       Y
 GFF 3             Y       Y

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bioperl.org/wiki/Mailing_list  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

 Allen Day, <allenday@ucla.edu>

=head1 CONTRIBUTORS

 Steffen Grossmann, <grossman@molgen.mpg.de>
 Scott Cain, <scain@cpan.org>
 Rob Edwards <rob@salmonella.org>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::FeatureIO::gff;
use strict;

#these are alphabetical, keep them that way.
use Bio::Annotation::DBLink;
use Bio::Annotation::OntologyTerm;
use Bio::Annotation::SimpleValue;
use Bio::Annotation::Target;
use Bio::FeatureIO;
use Bio::Ontology::OntologyStore;
use Bio::OntologyIO;
use Bio::SeqFeature::Annotated;
use Bio::SeqIO;
use URI::Escape;

use base qw(Bio::FeatureIO);

use constant DEFAULT_VERSION => 3;
my $RESERVED_TAGS   = "ID|Name|Alias|Parent|Target|Gap|Derives_from|Note|Dbxref|dbxref|Ontology_term|Index|CRUD";

sub _initialize {
  my($self,%arg) = @_;

  $self->SUPER::_initialize(%arg);

  $self->version( $arg{-version}        || DEFAULT_VERSION);
  $self->validate($arg{-validate_terms} || 0);

  $self->ignore_seq_region($arg{-ignore_seq_region} || 0);

  if ($arg{-file} =~ /^>.*/ ) {
    $self->_print("##gff-version " . $self->version() . "\n");
  }
  else {
    my $directive;
    while(($directive = $self->_readline()) && ($directive =~ /^##/) ){
      $self->_handle_directive($directive);
    }
    $self->_pushback($directive);
  }
  
  #need to validate against SOFA, no SO
  if ($self->validate) {
    $self->so(
              Bio::Ontology::OntologyStore->get_ontology('Sequence Ontology Feature Annotation')
              );
  }
}

=head2 next_feature()

 Usage   : my $feature = $featureio->next_feature();
 Function: reads a feature record from a GFF stream and returns it as an object.
 Returns : a Bio::SeqFeature::Annotated object
 Args    : N/A

=cut

sub next_feature {
  my $self = shift;
  my $gff_string;

  my($f) = $self->_buffer_feature();
  if($f){
    return $f;
  }

  return if $self->fasta_mode();

  # be graceful about empty lines or comments, and make sure we return undef
  # if the input is consumed
  while(($gff_string = $self->_readline()) && defined($gff_string)) {
    next if $gff_string =~ /^\s*$/;   #skip blank lines
    next if $gff_string =~ /^\#[^#]/; #skip comments, but not directives
    last;
  }

  return unless $gff_string;

  # looks like we went into FASTA mode without a directive.
  if($gff_string =~ /^>/){
    $self->_pushback($gff_string);
    $self->fasta_mode(1);
    return;
  }

  # got a directive
  elsif($gff_string =~ /^##/){
    $self->_handle_directive($gff_string);
    # recurse down to  the next line.  this will bottom out on finding a real feature or EOF
    return $self->next_feature();
  }

  # got a feature
  else {
    return $self->_handle_feature($gff_string);
  }
}

=head2 next_feature_group

 Title   : next_feature_group
 Usage   : @feature_group = $stream->next_feature_group
 Function: Reads the next feature_group from $stream and returns it.

           Feature groups in GFF3 files are separated by '###' directives. The
           features in a group might form a hierarchical structure. The
           complete hierarchy of features is returned, i.e. the returned array
           represents only the top-level features.  Lower-level features can
           be accessed using the 'get_SeqFeatures' method recursively.

 Example : # getting the complete hierarchy of features in a GFF3 file
           my @toplevel_features;
           while (my @fg = $stream->next_feature_group) {
               push(@toplevel_features, @fg);
           }
 Returns : an array of Bio::SeqFeature::Annotated objects
 Args    : none

=cut

sub next_feature_group {
  my $self = shift;

  my $feat;
  my %seen_ids;
  my @all_feats;
  my @toplevel_feats;

  $self->{group_not_done} = 1;

  while ($self->{group_not_done} && ($feat = $self->next_feature()) && defined($feat)) {
	# we start by collecting all features in the group and
	# memorizing those which have an ID attribute
    my $anno_ID = $feat->get_Annotations('ID');
	if(ref($anno_ID)) {
      my $attr_ID = $anno_ID->value;
      $self->throw("Oops! ID $attr_ID exists more than once in your file!")
		if (exists($seen_ids{$attr_ID}));
      $seen_ids{$attr_ID} = $feat;
	}
	push(@all_feats, $feat);
  }

  # assemble the top-level features
  foreach $feat (@all_feats) {
	my @parents = $feat->get_Annotations('Parent');
	if (@parents) {
      foreach my $parent (@parents) {
		my $parent_id = $parent->value;
		$self->throw("Parent with ID $parent_id not found!") unless (exists($seen_ids{$parent_id}));
		$seen_ids{$parent_id}->add_SeqFeature($feat);
      }
	} else {
	    push(@toplevel_feats, $feat);
      }
  }

  return @toplevel_feats;
}

=head2 next_seq()

access the FASTA section (if any) at the end of the GFF stream.  note that this method
will return undef if not all features in the stream have been handled

=cut

sub next_seq() {
  my $self = shift;
  return unless $self->fasta_mode();

  #first time next_seq has been called.  initialize Bio::SeqIO instance
  if(!$self->seqio){
    $self->seqio( Bio::SeqIO->new(-format => 'fasta', -fh => $self->_fh()) );
  }
  return $self->seqio->next_seq();
}

=head2 write_feature()

 Usage   : $featureio->write_feature( Bio::SeqFeature::Annotated->new(...) );
 Function: writes a feature in GFF format.  the GFF version used is governed by the
           '-version' argument passed to Bio::FeatureIO->new(), and defaults to GFF
           version 3.
 Returns : ###FIXME
 Args    : a Bio::SeqFeature::Annotated object.

=cut

sub write_feature {
  my($self,$feature) = @_;
  if (!$feature) {
    $self->throw("gff.pm cannot write_feature unless you give a feature to write.\n");
  }
  $self->throw("only Bio::SeqFeature::Annotated objects are writeable") unless $feature->isa('Bio::SeqFeature::Annotated');

  if($self->version == 1){
    return $self->_write_feature_1($feature);
  } elsif($self->version == 2){
    return $self->_write_feature_2($feature);
  } elsif($self->version == 2.5){
    return $self->_write_feature_25($feature);
  } elsif($self->version == 3){
    return $self->_write_feature_3($feature);
  } else {
    $self->throw(sprintf("don't know how to write GFF version %s",$self->version));
  }
}

################################################################################

=head1 ACCESSORS

=cut

=head2 fasta_mode()

 Usage   : $obj->fasta_mode($newval)
 Function: 
 Example : 
 Returns : value of fasta_mode (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

Side effect when setting: rewind the file handle a little bit to get the last
carriage return that was swallowed when the previous line was processed.

=cut

sub fasta_mode {
  my($self,$val) = @_;

  $self->{'fasta_mode'} = $val if defined($val);

  if ($val && $val == 1) {
  #  seek $self->_fh(), -1, 1; #rewind 1 byte to get the previous line's \n
    $self->_pushback("\n");
  }

  return $self->{'fasta_mode'};
}

=head2 seqio()

 Usage   : $obj->seqio($newval)
 Function: holds a Bio::SeqIO instance for handling the GFF3 ##FASTA section.
 Returns : value of seqio (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub seqio {
  my($self,$val) = @_;
  $self->{'seqio'} = $val if defined($val);
  return $self->{'seqio'};
}

=head2 sequence_region()

 Usage   :
 Function: ###FIXME
 Returns : 
 Args    :


=cut

sub sequence_region {
  my ($self,$k,$v) = @_;
  if(defined($k) && defined($v)){
    $self->{'sequence_region'}{$k} = $v;
    return $v;
  }
  elsif(defined($k)){
    return $self->{'sequence_region'}{$k};
  }
  else {
    return;
  }
}


=head2 so()

 Usage   : $obj->so($newval)
 Function: holds a Sequence Ontology instance
 Returns : value of so (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub so {
  my $self = shift;
  my $val = shift;
  ###FIXME validate $val object's type
  $self->{so} = $val if defined($val);
  return $self->{so};
}

=head2 validate()

 Usage   : $obj->validate($newval)
 Function: true if encountered ontology terms in next_feature()
           mode should be validated.
 Returns : value of validate (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub validate {
  my($self,$val) = @_;
  $self->{'validate'} = $val if defined($val);
  return $self->{'validate'};
}

=head2 version()

 Usage   : $obj->version($newval)
 Function: version of GFF to read/write.  valid values are 1, 2, 2.5, and 3.
 Returns : value of version (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub version {
  my $self = shift;
  my $val = shift;
  my %valid = map {$_=>1} (1, 2, 2.5, 3);
  if(defined $val && $valid{$val}){
    return $self->{'version'} = $val;
  }
  elsif(defined($val)){
    $self->throw('invalid version.  valid versions: '.join(' ', sort keys %valid));
  }
  return $self->{'version'};
}

################################################################################

=head1 INTERNAL METHODS

=cut

=head2 _buffer_feature()

 Usage   :
 Function: ###FIXME
 Returns : 
 Args    :

=cut

sub _buffer_feature {
  my ($self,$f) = @_;

  if ( $f ) {
    push @{ $self->{'buffer'} }, $f;
    return $f;
  }
  elsif ( $self->{'buffer'} ) {
    return shift @{ $self->{'buffer'} };
  }
  else {
    return;
  }
}


=head1 ignore_seq_region

Set this flag to keep FeatureIO from returning 
a feature for a ##sequence-region directive

=cut

sub ignore_seq_region {
  my($self,$val) = @_;
  $self->{'ignore_seq_region'} = $val if defined($val);
  return $self->{'ignore_seq_region'};
}

=head1 _handle_directive()

this method is called for lines beginning with '##'.

=cut

sub _handle_directive {
  my($self,$directive_string) = @_;

  $directive_string =~ s/^##//; #remove escape
  my($directive,@arg) = split /\s+/, $directive_string;

  if($directive eq 'gff-version'){
    my $version = $arg[0];
    if($version != 3){
      $self->throw("this is not a gff version 3 document, it is version '$version'");
    }
  }

  elsif($directive eq 'sequence-region'){
    # RAE: Sequence regions are in the format sequence-region seqid start end
    # for these we want to store the seqid, start, and end. Then when we validate
    # we want to make sure that the features are within the seqid/start/end

    return if $self->ignore_seq_region();

    $self->throw('Both start and end for sequence region should be defined')
      unless $arg[1] && $arg[2];
    my $fta = Bio::Annotation::OntologyTerm->new();
    $fta->name( 'region');

    my $f = Bio::SeqFeature::Annotated->new();
    $f->seq_id( $arg[0] );
    $f->start(  $arg[1] );
    $f->end(    $arg[2] );

    $f->type(   $fta    );

    #cache this in sequence_region(), we may need it for validation later.
    $self->sequence_region($f->seq_id => $f);

    #NOTE: is this the right thing to do -- treat this as a feature? -allenday
    #buffer it to be returned by next_feature()
    $self->_buffer_feature($f);
  }

  elsif($directive eq 'feature-ontology'){
    $self->warn("'##$directive' directive handling not yet implemented");
  }

  elsif($directive eq 'attribute-ontology'){
    $self->warn("'##$directive' directive handling not yet implemented");
  }

  elsif($directive eq 'source-ontology'){
    $self->warn("'##$directive' directive handling not yet implemented");
  }

  elsif($directive eq 'FASTA' or $directive =~ /^>/){
    #next_seq() will take care of this.
    $self->fasta_mode(1);
    return;
  }

  elsif($directive eq '#'){
    #all forward references resolved
    $self->{group_not_done} = 0;
  }

  elsif($directive eq 'organism') {
    my $organism = $arg[0];
    $self->organism($organism);
  }

  else {
    $self->throw("don't know what do do with directive: '##".$directive."'");
  }
}

=head1 _handle_feature()

this method is called for each line not beginning with '#'.  it parses the line and returns a
Bio::SeqFeature::Annotated object.

=cut

sub _handle_feature {
  my($self,$feature_string) = @_;

  my $feat = Bio::SeqFeature::Annotated->new();

  my($seq,$source,$type,$start,$end,$score,$strand,$phase,$attribute_string) = split /\t/, $feature_string;

  $feat->seq_id($seq);
  $feat->source_tag($source);
  $feat->start($start) unless $start eq '.';
  $feat->end($end) unless $end eq '.';
  $feat->strand($strand eq '+' ? 1 : $strand eq '-' ? -1 : 0);
  $feat->score($score);
  $feat->phase($phase);

  my $fta = Bio::Annotation::OntologyTerm->new();

  if($self->validate()){
    # RAE Added a couple of validations based on the GFF3 spec at http://song.sourceforge.net/gff3.shtml
    # 1. Validate the id
    if ($seq =~ /[^a-zA-Z0-9\.\-\:\^\*\$\@\!\+\_\?]/) { # I just escaped everything.
      $self->throw("Validation Error: seqid ($seq) contains characters that are not [a-zA-Z0-9.:^*\$\@!+_?\-] and not escaped");
    }

    if ($seq =~ /\s/) {
      $self->throw("Validation Error: seqid ($seq) contains unescaped whitespace")
    }

    # NOTE i think we're handling this in as a directive, and this test may be removed -allenday
    if ($seq =~ /^>/) {
      $self->throw("Validation Error: seqid ($seq) begins with a >")
    }

    # 2. Validate the starts and stops.
    # these need to be within the region's bounds, and
    # also start <= end.  bail out if either is not true.
    if ($start > $end) {
      $self->throw("Validation Error: start ($start) must be less than or equal to end in $seq");
    }

    my $region = $self->sequence_region($seq);
    # NOTE: we can only validate against sequence-region that are declared in the file.
    # if i reference some region from elsewhere, can't validate.  if we want to be really strict
    # we should bail out here. -allenday
    if ( defined($region) && $start < $region->start() || $end > $region->end() ) {
      $self->throw("Validation Error: sequence location ($seq from $start to $end) does not appear to lie within a defined sequence-region")
    }

    # 3. Validate the strand.
    # In the unvalidated version +=1 and -=-1. Everything else is 0. We just need to warn when it is not [+-.?]
    $self->throw("Validation Error: strand is not one of [+-.?] at $seq") if ($strand =~ /^[^\+\-\.\?]$/);

    # 4. Validate the phase to be one of [.012]
    $self->throw("Validation Error: phase is not one of [.012] at $seq") if ($phase =~ /^[^\.012]$/);

    my $feature_type;
    if($type =~ /^\D+:\d+$/){
      #looks like an identifier
      ($feature_type) = $self->so->find_terms(-identifier => $type);
    } else {
      #looks like a name
      ($feature_type) = $self->so->find_terms(-name => $type);
    }

    if(!$feature_type){
      $self->throw("Validation Error: couldn't find ontology term for '$type'.");
    }
    $fta->term($feature_type);
  } else {

    if($type =~ /^\D+:\d+$/){
      #looks like an identifier
      $fta->identifier($type)
    } else {
      $fta->name($type);
    }
  }

  $feat->type($fta);

  my %attr = ();
  chomp $attribute_string;

  unless ( $attribute_string eq '.' ) {
    my @attributes = split ';', $attribute_string;
    foreach my $attribute (@attributes){
      my($key,$values) = split '=', $attribute;

      # remove leading and trailing quotes from values
      $values =~ s/^["']//;
      $values =~ s/["']$//; #' terminate the quote for emacs

      my @values = map{uri_unescape($_)} split ',', $values;

     #minor hack to allow for multiple instances of the same tag
      if ($attr{$key}) {
        my @tmparray = @{$attr{$key}};
        push @tmparray, @values;
        $attr{$key} = [@tmparray];
      } else {
        $attr{$key} = [@values];
      }
    }
  }

  #Handle Dbxref attributes
  if($attr{Dbxref} or $attr{dbxref}){
    foreach my $value (@{ $attr{Dbxref} }, @{ $attr{dbxref} }){
      my $a = Bio::Annotation::DBLink->new();
      my($db,$accession) = $value =~ /^(.+?):(.+)$/;

      if(!$db or !$accession){ #dbxref malformed
        $self->throw("Error in line:\n$feature_string\nDbxref value '$value' did not conform to GFF3 specification");
        next;
      }

      $a->database($db);
      $a->primary_id($accession);
      $feat->add_Annotation('Dbxref',$a);
    }
  }

  #Handle Ontology_term attributes
  if($attr{Ontology_term}){
    foreach my $id (@{ $attr{Ontology_term} }){
      my $a = Bio::Annotation::OntologyTerm->new();

      if($self->validate()){
        my $ont_name = Bio::Ontology::OntologyStore->guess_ontology($id);
        my $ont = Bio::Ontology::OntologyStore->get_ontology($ont_name);
        my($term) = $ont->find_terms(-identifier => $id);
        $a->term($term);
      } else {
        $a->identifier($id);
      }

      $feat->add_Annotation('Ontology_term',$a);
    }
  }

  #Handle Gap attributes
  if($attr{Gap}){
     for my $value (@{ $attr{Gap} }) {
       my $a = Bio::Annotation::SimpleValue->new();
       $a->value($value);
       $feat->add_Annotation('Gap',$a);
     }
  }

  #Handle Target attributes
  if($attr{Target}){
    my $target_collection = Bio::Annotation::Collection->new();

    foreach my $target_string (@{ $attr{Target} } ) {

      #only replace + for space if + has been used in place of it
      #that is, + could also mean plus strand, and we don't want
      #to accidentally remove it
 
      #presumably you can't use + for space and + for strand in the same string.      
      $target_string =~ s/\+/ /g unless $target_string =~ / /; 

      my ($t_id,$tstart,$tend,$strand,$extra) = split /\s+/, $target_string; 
      if (!$tend || $extra) { # too much or too little stuff in the string
        $self->throw("The value in the Target string, $target_string, does not conform to the GFF3 specification");
      }

      my $a = Bio::Annotation::Target->new(
           -target_id => $t_id,
           -start     => $tstart,
           -end       => $tend,
      );

      if ($strand && $strand eq '+') {
        $strand = 1;
      } elsif ($strand && $strand eq '-') {
        $strand = -1;
      } else {
        $strand = '';
      }

      $a->strand($strand) if $strand;
      $feat->add_Annotation('Target',$a); 
    }
  }

  #Handle ID attribute.  May only have one ID, throw error otherwise

  if($attr{ID}){
    if(scalar( @{ $attr{ID} } ) > 1){
      $self->throw("Error in line:\n$feature_string\nA feature may have at most one ID value");
    }

    #ID's must be unique in the file
    if ($self->{'allIDs'}->{${$attr{ID}}[0]} && $self->validate()) {
      $self->throw("Validation Error: The ID ${$attr{ID}}[0] occurs more than once in the file, but should be unique");
    }
    $self->{'allIDs'}->{${$attr{ID}}[0]} = 1;


    my $a = Bio::Annotation::SimpleValue->new();
    $a->value( @{ $attr{ID} }[0] );
    $feat->add_Annotation('ID',$a);
  }

  #Handle Name attribute.  May only have one Name, throw error otherwise
  if($attr{Name}){
    if(scalar( @{ $attr{Name} } ) > 1){
      $self->throw("Error in line:\n$feature_string\nA feature may have at most one Name value");
    }

    my $a = Bio::Annotation::SimpleValue->new();
    $a->value( @{ $attr{Name} }[0] );
    $feat->add_Annotation('Name',$a);
  }

  foreach my $other_canonical (qw(Alias Parent Note Derives_from Index CRUD)){
    if($attr{$other_canonical}){
      foreach my $value (@{ $attr{$other_canonical} }){
        my $a = Bio::Annotation::SimpleValue->new();
        $a->value($value);
        $feat->add_Annotation($other_canonical,$a);
      }
    }
  }

  my @non_reserved_tags = grep {/^[a-z]/} keys %attr;
  foreach my $non_reserved_tag (@non_reserved_tags) {
    next if ($non_reserved_tag eq 'dbxref');
    foreach my $value (@{ $attr{$non_reserved_tag} }){
      $feat = $self->_handle_non_reserved_tag($feat,$non_reserved_tag,$value);
    }
  }

  my @illegal_tags = grep 
 {!/($RESERVED_TAGS)/} 
 grep {/^[A-Z]/} keys %attr;

  if (@illegal_tags > 0) {
      my $tags = join(", ", @illegal_tags);
      $self->throw("The following tag(s) are illegal and are causing this parser to die: $tags");
  }

  return $feat;
}

=head2 _handle_non_reserved_tag()

 Usage   : $self->_handle_non_reserved_tag($feature,$tag,$value)
 Function: Deal with non-reserved word tags in the ninth column
 Returns : An updated Bio::SeqFeature::Annotated object
 Args    : A Bio::SeqFeature::Annotated and a tag/value pair

Note that this method can be overridden in a subclass to provide
special handling of non-reserved word tags.

=cut

sub _handle_non_reserved_tag {
  my $self = shift;
  my ($feat,$tag,$value) = @_;

  # to customize through subclassing and overriding:
  #if ($tag eq 'someTagOfInterest') {
  #  do something different
  # else { do what is below

  my $a;
  if ($tag eq 'comment') {
    $a = Bio::Annotation::Comment->new();
  }
  else {
    $a = Bio::Annotation::SimpleValue->new();
  }
  $a->value($value); 
  $feat->add_Annotation($tag,$a);
  
  return $feat;
}

=head1 organims

Gets/sets the organims from the organism directive

=cut

sub organism {
    my $self = shift;
    my $organism = shift if defined(@_);
    return $self->{'organism'} = $organism if defined($organism);
    return $self->{'organism'};
}


=head1 _write_feature_1()

write a feature in GFF v1 format.  currently not implemented.

=cut

sub _write_feature_1 {
  my($self,$feature) = @_;
  $self->throw(sprintf("write_feature unimplemented for GFF version %s",$self->version));
}

=head1 _write_feature_2()

write a feature in GFF v2 format.  currently not implemented.

=cut

sub _write_feature_2 {
  my($self,$feature) = @_;
  $self->throw(sprintf("write_feature unimplemented for GFF version %s",$self->version));
}

=head1 _write_feature_25()

write a feature in GFF v2.5 (aka GTF) format.

=cut

sub _write_feature_25 {
  my($self,$feature,$group) = @_;

  #the top-level feature is an aggregate of all subfeatures
  my ($transcript_id, $gene_id) = (($feature->get_Annotations('transcript_id'))[0], ($feature->get_Annotations('gene_id'))[0]);
  if(!defined($group)){
    $group = ($feature->get_Annotations('ID'))[0];
    $transcript_id ||= $group;
    $gene_id ||= $group;
  }
  

  my $seq    = ref($feature->seq_id) ? $feature->seq_id->value : $feature->seq_id;
  my $source = $feature->source->value;
  my $type   = $feature->type->name;
  $type = 'EXON' if $type eq 'exon'; #a GTF peculiarity, incosistent with the sequence ontology.
  my $min    = $feature->start   || '.';
  my $max    = $feature->end     || '.';
  my $strand = $feature->strand == 1 ? '+' : $feature->strand == -1 ? '-' : '.';
  my $score  = defined($feature->score) ? (ref($feature->score) ? $feature->score->value : $feature->score) : '.'; # score is optional
  my $frame  = defined($feature->frame) ? (ref($feature->frame) ? $feature->frame->value : $feature->frame) : (ref($feature->phase) ? $feature->phase->value : $feature->phase);

  #these are the only valid types in a GTF document
  if($type eq 'EXON' or $type eq 'CDS' or $type eq 'start_codon' or $type eq 'stop_codon'){
    my $attr = sprintf('gene_id "%s"; transcript_id "%s";',$gene_id ? $gene_id->value : '',$transcript_id ? $transcript_id->value : '');
    my $outstring = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                            $seq,$source,$type,$min,$max,$score,$strand,$frame eq '.' ? 0 : $frame,$attr);

    $self->_print($outstring);
  }

  foreach my $subfeat ($feature->get_SeqFeatures){
    $self->_write_feature_25($subfeat,$group);
  }
}

=head1 _write_feature_3()

write a feature in GFF v3 format.

=cut

sub _write_feature_3 {
  my($self,$feature) = @_;
  my $seq    = ref($feature->seq_id) ? $feature->seq_id->value : $feature->seq_id;
  my $source;
  if ($feature->source()) {
    $source = $feature->source->value;
  }
  else {
    $source = $feature->source() || "unknownsource";
  }
  my $type;
  if ($feature->type()) { $type = $feature->type->name; }
  else { $type = "region"; }
  my $min    = $feature->start   || '.';
  my $max    = $feature->end     || '.';
  my $strand = $feature->strand == 1 ? '+' : $feature->strand == -1 ? '-' : '.';
  my $score  = defined($feature->score) ? (ref($feature->score) ? $feature->score->value : $feature->score) : undef;
  my $phase  = defined($feature->phase) ? (ref($feature->phase) ? $feature->phase->value : $feature->phase) : undef;

  my @attr;
  if(my @v = ($feature->get_Annotations('Name'))){
    my $vstring = join ',', map {uri_escape($_->value)} @v;
    push @attr, "Name=$vstring";
  }
  if(my @v = ($feature->get_Annotations('ID'))){
    my $vstring = join ',', map {uri_escape($_->value)} @v;
    push @attr, "ID=$vstring";
    $self->throw('GFF3 features may have at most one ID, feature with these IDs is invalid:\n'.$vstring) if scalar(@v) > 1;
  }
  if(my @v = ($feature->get_Annotations('Parent'))){
    my $vstring = join ',', map {uri_escape($_->value)} @v;
    push @attr, "Parent=$vstring";
  }
  if(my @v = ($feature->get_Annotations('dblink'))){
    my $vstring = join ',', map {uri_escape($_->database .':'. $_->primary_id)} @v;
    push @attr, "Dbxref=$vstring";
  }
  if(my @v = ($feature->get_Annotations('ontology_term'))){
    my $vstring = join ',', map {uri_escape($_->identifier)} @v;
    push @attr, "Ontology_term=$vstring";
  }
  if(my @v = ($feature->get_Annotations('comment'))){
    my $vstring = join ',', map {uri_escape($_->text)} @v;
    push @attr, "Note=$vstring";
  }
  if(my @v = ($feature->get_Annotations('Target'))){
    my %strand_map = ( 1=>'+', 0=>'', -1=>'-', '+' => '+', '-' => '-' );
    my $vstring = join ',', map {
      uri_escape($_->target_id).' '.$_->start.' '.$_->end.(defined $_->strand ? ' '.$strand_map{$_->strand} : '')
    } @v;
    push @attr, "Target=$vstring";
  }

  my $attr = join ';', @attr;

  my $outstring = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                          $seq,$source,$type,$min,$max,$score,$strand,$phase,$attr);

  $self->_print($outstring);

  foreach my $subfeat ($feature->get_SeqFeatures){
    $self->_write_feature_3($subfeat);
  }
}




1;
