package Bio::FeatureIO::gff;

use strict;
use base qw(Bio::FeatureIO);
use Bio::SeqFeature::Annotated;
use Bio::OntologyIO;

use Bio::Annotation::OntologyTerm;
use Bio::Annotation::SimpleValue;
use Bio::Annotation::DBLink;

use Bio::Ontology::OntologyStore;

=pod

currently only handles gff version 3.  spec at http://song.sf.net/

=cut

sub _initialize {
  my($self,%arg) = @_;

  $self->SUPER::_initialize(%arg);

  #read headers
  my $directive;
  while(($directive = $self->_readline()) && $directive =~ /^##/){
    $self->_handle_directive($directive)
  }
  $self->_pushback($directive);

  $self->so(
            Bio::Ontology::OntologyStore->get_ontology('Sequence Ontology')
           );
}

sub write_feature {
  my($self,$feature) = @_;
  $self->throw("only Bio::SeqFeature::Annotated objects are writeable") unless $feature->isa('Bio::SeqFeature::Annotated');

  my $seq    = $feature->seq_id || '.';
  my $source = ($feature->annotation->get_Annotations('source'))[0]->value;
  my $type   = ($feature->annotation->get_Annotations('feature_type'))[0]->name;
  my $min    = $feature->start   || '.';
  my $max    = $feature->end     || '.';
  my $strand = $feature->strand == 1 ? '+' : $feature->strand == -1 ? '-' : '.';
  my $score  = $feature->score   || '.';
  my $phase  = $feature->frame   || '.';

  my @attr;
  if(my @v = ($feature->annotation->get_Annotations('Name'))){
    my $vstring = join ',', map {$_->value} @v;
    push @attr, "Name=$vstring";
  }
  if(my @v = ($feature->annotation->get_Annotations('ID'))){
    my $vstring = join ',', map {$_->value} @v;
    push @attr, "ID=$vstring";
  }
  if(my @v = ($feature->annotation->get_Annotations('Parent'))){
    my $vstring = join ',', map {$_->value} @v;
    push @attr, "Parent=$vstring";
  }

  my $attr = join ';', @attr;

  my $outstring = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                          $seq,$source,$type,$min,$max,$score,$strand,$phase,$attr);

  $self->_print($outstring);

  foreach my $subfeat ($feature->get_SeqFeatures){
    $self->write_feature($subfeat);
  }
}

sub next_feature {
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

  if($gff_string =~ /^##/){
    $self->_handle_directive($gff_string);
    return $self->next_feature();
  } else {
    return $self->_handle_feature($gff_string);
  }
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

  elsif($directive eq 'FASTA'){
    $self->warn("'##$directive' directive handling not yet implemented");
  }

  elsif($directive eq '#'){
    #all forward references resolved
    $self->warn("'##$directive' directive handling not yet implemented");
  }

  else {
    $self->throw("don't know what do do with directive: '##".$directive."'");
  }
}

sub _handle_feature {
  my($self,$feature_string) = @_;

  my $feat = Bio::SeqFeature::Annotated->new();
  my $ac = $feat->annotation(); #initialize Bio::Annotation::Collection

  my($seq,$source,$type,$start,$end,$score,$strand,$phase,$attribute_string) = split /\s+/, $feature_string;

  $feat->seq_id($seq);
  $feat->start($start) unless $start eq '.';
  $feat->end($end) unless $end eq '.';
  $feat->strand($strand eq '+' ? 1 : $strand eq '-' ? -1 : 0);
  $feat->score($score) unless $score eq '.';
  $feat->frame($phase);

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
  my $fta = Bio::Annotation::OntologyTerm->new();
  $fta->term($feature_type);

  my %attr = ();
  my @attributes = split ';', $attribute_string;
  foreach my $attribute (@attributes){
    my($key,$values) = split '=', $attribute;
    my @values = split ',', $values;
    $attr{$key} = [@values];
  }

  #Handle Dbxref attributes
  if($attr{Dbxref}){
    foreach my $value (@{ $attr{Dbxref} }){
      my $a = Bio::Annotation::DBLink->new();
      my($db,$accession) = $value =~ /^(.+?):(.+)$/;

      if(!$db or !$accession){ #dbxref malformed
        $self->throw("Error in line:\n$feature_string\nDbxref value '$value' did not conform to GFF3 specification");
      }

      $a->database($db);
      $a->primary_id($accession);
    }
  }

  #Handle Ontology_term attributes
  if($attr{Ontology_term}){
    $self->warn("Warning for line:\n$feature_string\nOntology_term attribute handling not yet implemented, skipping it");
  }

  #Handle Gap attributes
  if($attr{Gap}){
    $self->warn("Warning for line:\n$feature_string\nGap attribute handling not yet implemented, skipping it");
  }

  #Handle Target attributes
  if($attr{Target}){
    $self->warn("Warning for line:\n$feature_string\nTarget attribute handling not yet implemented, skipping it");
  }

  #Handle ID attribute.  May only have one ID, throw error otherwise
  if($attr{ID}){
    if(scalar( @{ $attr{ID} } ) > 1){
      $self->throw("Error in line:\n$feature_string\nA feature may have at most one ID value");
    }

    my $a = Bio::Annotation::SimpleValue->new();
    $a->value( @{ $attr{ID} }[0] );
    $ac->add_Annotation('ID',$a);
  }

  #Handle Name attribute.  May only have one Name, throw error otherwise
  if($attr{Name}){
    if(scalar( @{ $attr{Name} } ) > 1){
      $self->throw("Error in line:\n$feature_string\nA feature may have at most one Name value");
    }

    my $a = Bio::Annotation::SimpleValue->new();
    $a->value( @{ $attr{Name} }[0] );
    $ac->add_Annotation('Name',$a);
  }


  foreach my $other_canonical (qw(Alias Parent Note)){
    if($attr{$other_canonical}){
      foreach my $value (@{ $attr{$other_canonical} }){
        my $a = Bio::Annotation::SimpleValue->new();
        $a->value($value);
        $ac->add_Annotation($other_canonical,$a);
      }
    }
  }

  $ac->add_Annotation('feature_type',$fta);

  return $feat;
}

sub so {
  my $self = shift;
  my $val = shift;
  $self->{so} = $val if defined($val);
  return $self->{so};
}

1;
