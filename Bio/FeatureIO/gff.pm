=pod

=head1 NAME

Bio::FeatureIO::gff - DESCRIPTION of Object

=head1 SYNOPSIS

  my $feature; #get a Bio::SeqFeature::Annotated somehow
  my $featureOut = Bio::FeatureIO->new(-format => 'gff',
                                       -version => 3,
                                       -fh => \*STDOUT,
                                       -validate_terms => 1, #boolean. validate ontology
                                                             #terms online?  default false.
                                      );
  $featureOut->write_feature($feature);

=head1 DESCRIPTION

 currently implemented:

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

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Allen Day

Email allenday@ucla.edu

Describe contact details here

=head1 CONTRIBUTORS

 Steffen Grossmann, E<lt>grossman@molgen.mpg.deE<gt>
 Scott Cain, E<lt>cain@cshl.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::FeatureIO::gff;

use strict;
use base qw(Bio::FeatureIO);

use Bio::FeatureIO;

use Bio::SeqFeature::Annotated;
use Bio::OntologyIO;

use Bio::Annotation::DBLink;
use Bio::Annotation::OntologyTerm;
use Bio::Annotation::SimpleValue;
use Bio::Annotation::Target;

use Bio::Ontology::OntologyStore;

use URI::Escape;

use constant DEFAULT_VERSION => 3;

sub _initialize {
  my($self,%arg) = @_;

  $self->SUPER::_initialize(%arg);

  $self->version($arg{-version} || DEFAULT_VERSION);
  $self->validate($arg{-validate_terms} || 0);

  #read headers
  my $directive;
  while(($directive = $self->_readline()) && $directive =~ /^##/){
    $self->_handle_directive($directive)
  }
  $self->_pushback($directive);
  
#  if ($self->mode eq 'w') {   # [SG] I would love to use 'mode' but it does strange things...
  if ($arg{-file} =~ /^>.*/ ) {
      $self->_print("##gff-version " . $self->version() . "\n");
  }

  #need to validate against SOFA, no SO
  $self->so(
            Bio::Ontology::OntologyStore->get_ontology('Sequence Ontology Feature Annotation')
           );
}

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

sub _write_feature_1 {
  my($self,$feature) = @_;
  $self->throw(sprintf("write_feature unimplemented for GFF version %s",$self->version));
}

sub _write_feature_2 {
  my($self,$feature) = @_;
  $self->throw(sprintf("write_feature unimplemented for GFF version %s",$self->version));
}

sub _write_feature_25 {
  my($self,$feature,$group) = @_;

  #the top-level feature is an aggregate of all subfeatures
  if(!defined($group)){
    $group = ($feature->get_Annotations('ID'))[0]->value;
  }

  my $seq    = $feature->seq_id->value;
  my $source = $feature->source->value;
  my $type   = $feature->type->name;
  $type = 'EXON' if $type eq 'exon'; #a GTF peculiarity, incosistent with the sequence ontology.
  my $min    = $feature->start   || '.';
  my $max    = $feature->end     || '.';
  my $strand = $feature->strand == 1 ? '+' : $feature->strand == -1 ? '-' : '.';
  my $score  = $feature->score->value;
  my $phase  = $feature->phase->value;

  #these are the only valid types in a GTF document
  if($type eq 'EXON' or $type eq 'CDS' or $type eq 'start_codon' or $type eq 'stop_codon'){
    my $attr = sprintf('gene_id "%s"; transcript_id "%s";',$group,$group);
    my $outstring = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                            $seq,$source,$type,$min,$max,$score,$strand,$phase,$attr);

    $self->_print($outstring);
  }

  foreach my $subfeat ($feature->get_SeqFeatures){
    $self->_write_feature_25($subfeat,$group);
  }
}

sub _write_feature_3 {
  my($self,$feature) = @_;
  my $seq    = $feature->seq_id->value;
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
  my $score  = $feature->score->value;
  my $phase  = $feature->phase->value;

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
    my $vstring = join ',', map {uri_escape($_->primary_id)} @v;
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

  my $attr = join ';', @attr;

  my $outstring = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                          $seq,$source,$type,$min,$max,$score,$strand,$phase,$attr);

  $self->_print($outstring);

  foreach my $subfeat ($feature->get_SeqFeatures){
    $self->_write_feature_3($subfeat);
  }
}

sub next_feature {
  my $self = shift;
  
  my $feat = $self->_next_feature_or_directive();
  $feat = $self->next_feature() if ($feat eq 'directive'); # we repeat as long as we processed directives

  return $feat;
}

sub _next_feature_or_directive {
  my $self = shift;
  my $gff_string;

  # be graceful about empty lines or comments, and make sure we return undef
  # if the input's consumed
  while(($gff_string = $self->_readline()) && defined($gff_string)) {
    next if $gff_string =~ /^\s*$/;   #skip blank lines
    next if $gff_string =~ /^\#[^#]/; #skip comments, but not directives
    last;
  }

  return undef unless $gff_string;

  if($gff_string =~ /^##/ or $gff_string =~ /^>/){
    $self->_handle_directive($gff_string);
    return 'directive';   # we have to be able to detect when we processed a directive
  } else {
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

    while ($self->{group_not_done} && ($feat = $self->_next_feature_or_directive()) && defined($feat)) {
	next if ($feat eq 'directive');
	# we start by collecting all features in the group and 
	# memorizing those which have an ID attribute
	if(my $anno_ID = $feat->get_Annotations('ID')) {
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
    $self->warn("'##$directive' directive handling not yet implemented");
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

  elsif($directive eq 'FASTA' or $directive =~ /^>(.+)/){
    my $fasta_directive_id = $1 if $1;
    $self->warn("'##$directive' directive handling not yet implemented");
    local $/ = '>';
    while(my $read = $self->_readline()){
       chomp $read;
       my $fasta_id;
       my @seqarray = split /\n/, $read;
       if ($fasta_directive_id) {
         $fasta_id = $fasta_directive_id;
         $fasta_directive_id = '';
       } else {
         $fasta_id = shift @seqarray;
       }
       my $seq = join '', @seqarray;
    }
  }

  elsif($directive eq '#'){
    #all forward references resolved
    $self->{group_not_done} = 0;  
    $self->warn("'##$directive' directive handling not yet implemented");
  }

  else {
    $self->throw("don't know what do do with directive: '##".$directive."'");
  }
}

sub _handle_feature {
  my($self,$feature_string) = @_;

  my $feat = Bio::SeqFeature::Annotated->new();

  my($seq,$source,$type,$start,$end,$score,$strand,$phase,$attribute_string) = split /\t/, $feature_string;

  $feat->seq_id($seq);
  $feat->source($source);
  $feat->start($start) unless $start eq '.';
  $feat->end($end) unless $end eq '.';
  $feat->strand($strand eq '+' ? 1 : $strand eq '-' ? -1 : 0);
  $feat->score($score);
  $feat->phase($phase);

  my $fta = Bio::Annotation::OntologyTerm->new();

  if($self->validate()){

    my $feature_type;
    if($type =~ /^\D+:\d+$/){
      #looks like an identifier
      ($feature_type) = $self->so->find_terms(-identifier => $type);
    } else {
      #looks like a name
      ($feature_type) = $self->so->find_terms(-name => $type);
    }

    if(!$feature_type){
      $self->throw("couldn't find ontology term for '$type'.");
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
  my @attributes = split ';', $attribute_string;
  foreach my $attribute (@attributes){
    my($key,$values) = split '=', $attribute;
    my @values = map{uri_unescape($_)} split ',', $values;
    $attr{$key} = [@values];
  }

  #Handle Dbxref attributes
  if($attr{Dbxref}){
    foreach my $value (@{ $attr{Dbxref} }){
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
    $self->warn("Warning for line:\n$feature_string\nGap attribute handling not yet implemented, skipping it");
  }

  #Handle Target attributes
  if($attr{Target}){
    my $target_collection = Bio::Annotation::Collection->new();

    foreach my $target_string (@{ $attr{Target} } ) {
      $target_string =~ s/\+/ /g; 
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

  foreach my $other_canonical (qw(Alias Parent Note)){
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
    foreach my $value (@{ $attr{$non_reserved_tag} }){
      $feat = $self->_handle_non_reserved_tag($feat,$non_reserved_tag,$value);
    }
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


  my $a = Bio::Annotation::SimpleValue->new();
  $a->value($value);
  $feat->add_Annotation($tag,$a);

  return $feat;
}

=head2 version()

 Usage   : $obj->version($newval)
 Function: version of GFF to read/write.  valid values are 1, 2, 2.5, and 3.
 Returns : value of version (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub version{
    my $self = shift;

    return $self->{'version'} = shift if @_;
    return $self->{'version'};
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


1;
