package Bio::FeatureIO::gff;

use strict;
use base qw(Bio::FeatureIO);
use Bio::SeqFeature::Generic;
use Bio::OntologyIO;

use Bio::Annotation::OntologyTerm;

=pod

currently only handles gff version 3.  spec at http://song.sf.net/

=cut

sub _initialize {
  my $self = shift;
  my %arg = @_;

  #read headers
  my $directive;
  while(($directive = $self->_readline()) && $directive =~ /^##/){
    $self->_handle_directive($directive)
  }
  $self->_pushback($directive);

  $self->_setup_ontology('Sequence Ontology',
                         "http://cvs.sourceforge.net/viewcvs.py/*checkout*/song/ontology/so.ontology?rev=HEAD",
                         "http://cvs.sourceforge.net/viewcvs.py/*checkout*/song/ontology/so.definition?rev=HEAD",
                        );
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

  my $feat = Bio::SeqFeature::Generic->new();
  my $ac = $feat->annotation(); #initialize Bio::Annotation::Collection

  my($seq,$source,$type,$start,$end,$score,$strand,$phase,$attribute_string) = split /\s+/, $feature_string;

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

  my $fta = Bio::Annotation::OntologyTerm->new();
  $fta->term($feature_type);

  $ac->add_Annotation('feature_type',$fta);

  return $feat;  
}

sub _setup_ontology {
  my($self,$name,$ont_url,$def_url) = @_;

  my $ont_parser = Bio::OntologyIO->new(
                     -format       => 'so',
                     -defs_file    => $def_url,
                     -file         => $ont_url,
                   );
  my $so = $ont_parser->next_ontology();
  $self->so($so);
}

sub so {
  my $self = shift;
  my $val = shift;
  $self->{so} = $val if defined($val);
  return $self->{so};
}

1;
