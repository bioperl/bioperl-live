package Bio::Biblio::IO::ncbi_eutils;
use Bio::Biblio::IO::pubmedxml;
use Bio::Biblio::PubmedLinkout;
use Text::Abbrev;

use vars qw(@ISA @AUTO_ATTRIBUTES @OTHER_ATTRIBUTES %SYNONYMS);
use vars qw($WHAT %WHAT);
use strict;
use Bio::Root::AutoClass;
use Bio::Biblio::IO;
@ISA = qw(Bio::Root::AutoClass Bio::Biblio::IO);
$WHAT='abstracts';
%WHAT=abbrev qw(pubmed_ids pmids ids abstracts links);

BEGIN {
  @AUTO_ATTRIBUTES=qw(data_source format _parser);
  @OTHER_ATTRIBUTES=qw(what);
  %SYNONYMS=();
  Bio::Root::AutoClass::declare(__PACKAGE__,\@AUTO_ATTRIBUTES,\%SYNONYMS);
}

sub _init_self {
  my($self,$class,$args)=@_;
  return unless $class eq __PACKAGE__; # to prevent subclasses from re-running this
  $self->_parser(new Bio::Biblio::IO(-format=>'pubmedxml'));
}

sub parse {
  my($self,@xml)=@_;
  my $what= $self->what;
  return $self->_parse_abstracts(@xml) if $what eq 'abstracts';
  return $self->_parse_ids(@xml) if $what eq 'ids';
  return $self->_parse_links(@xml) if $what eq 'links';
}
sub _parse_abstracts {
  my($self,@xml)=@_;
  my @results;
  my $parser=$self->_parser;
  for my $xml (@xml) {
    $parser->{'_data'}=$xml;
    $parser->_parse;		                   # from Bio::Biblio::IO::medlinexml
    my $result=shift (@{$parser->{'_citations'}}); # ditto
    push(@results,$result);
  }
  wantarray? @results: $results[0];
}
sub _parse_ids {
  my($self,@xml)=@_;
  my @results;
  for my $xml (@xml) {
    my($result)=$xml=~/<Id>(\d+)<\/Id>/s;
    push(@results,$result);
  }
  wantarray? @results: $results[0];
}
sub _parse_links {
  my($self,@xml)=@_;
  my @results;
  my $tags={Id=>'text',ObjUrl=>{Url=>'text',SubjectType=>'text',
				Attribute=>'text',
				Provider=>{Name=>'text',NameAbbr=>'text',Id=>'text',Url=>'text',IconUrl=>'text'}}};
  my $parser=new XML::Parser(Style=>'Tree');
  for my $xml (@xml) {
    my $tree=$parser->parse($xml);
    my $result=new Bio::Biblio::PubmedLinkout(_prune($tree->[1],$tags));
    push(@results,$result);
  }
  wantarray? @results: $results[0];
}

sub _prune {
  my($tree,$tags)=@_;
  my @tags=keys %$tags;
  my(%content,%counts);
  my($attr,@content)=@$tree;
  for (my $i=0;$i<@content;$i+=2) {
    my($tag,$value)=($content[$i],$content[$i+1]);
    if (grep {lc($_) eq lc($tag)} @tags) {
      $value=_prune($value,$tags->{$tag}) if 'HASH' eq ref $tags->{$tag};
      $value=_get_text($value) if 'text' eq $tags->{$tag};
      $content{$tag}=$value if $counts{$tag}==0;
      $content{$tag}=[$content{$tag},$value] if $counts{$tag}==1;
      push(@{$content{$tag}},$value) if $counts{$tag}>1;
      $counts{$tag}++;
    }
  }
  \%content;
}

sub _get_text {
  my($value)=@_;
  my($attr,$tag,$text)=@$value;
  die "Tag should be 0 in text field, not $tag" unless $tag eq '0';
  $text;
}

sub what {
  my $self=shift @_;
  my $what;
  if (@_) {
    $what=_what($_[0]) ||  $self->throw("Invalid what $_[0]");
    $self->{'WHAT'}=$what;
  } else {
    $what=$self->{'WHAT'} || $WHAT;
  }
  $what;
}
sub _what {
  my($what_param)=@_;
  my $what=$WHAT{lc($what_param)};
  $what='ids' if $what=~/pubmed_ids|pmids|ids/;
  $what;
}
1;
