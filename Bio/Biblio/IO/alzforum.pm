package Bio::Biblio::IO::alzforum;
use XML::Parser;
use Text::Abbrev;

use vars qw(@ISA @AUTO_ATTRIBUTES @OTHER_ATTRIBUTES %SYNONYMS);
use vars qw($WHAT %WHAT);
use strict;
use Bio::Root::AutoClass;
use Bio::Biblio::IO;
use Bio::Biblio::AlzforumCitation;
@ISA = qw(Bio::Root::AutoClass Bio::Biblio::IO);
$WHAT='cits';
%WHAT=abbrev qw(cits);

BEGIN {
  @AUTO_ATTRIBUTES=qw(data_source format);
  @OTHER_ATTRIBUTES=qw(what);
  %SYNONYMS=();
  Bio::Root::AutoClass::declare(__PACKAGE__,\@AUTO_ATTRIBUTES,\%SYNONYMS);
}

sub _init_self {
  my($self,$class,$args)=@_;
  return unless $class eq __PACKAGE__; # to prevent subclasses from re-running this
}
sub parse {
  my($self,@xml)=@_;
  my $what= $self->what;
  return $self->_parse_citations(@xml) if $what eq 'cits';
}
sub _parse_citations {
  my($self,@xml)=@_;
  my @results;
  my $parser=new XML::Parser(Style=>'Tree');
  for my $xml (@xml) {
    my $tree=$parser->parse($xml);
    my $result=new Bio::Biblio::AlzforumCitation($tree->[1]->[0]);
    push(@results,$result);
  }
  wantarray? @results: $results[0];
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
  my $what;
  $what='cits' if $what_param=~/^c/i;
  $what or $what=$WHAT{lc($what_param)};
  $what;
}
1;

