package Bio::DB::Biblio::ncbi_eutils;
use LWP::Simple;
use URI::Escape;
use Text::Abbrev;
use Bio::Biblio::IO;
use overload 
  q{""} => sub {'#'.($_[0]->query_key);},
  'cmp' => sub {my $res=(overload::StrVal($_[0]) cmp overload::StrVal($_[1]));
		$_[2]? -$res: $res},
  'fallback'=>1,
  'nomethod' => sub {$_[0]};

use vars qw(@ISA @AUTO_ATTRIBUTES @OTHER_ATTRIBUTES %SYNONYMS);
use vars qw($EUTILS $WHAT %WHAT $CHUNK_SIZE);
use strict;
use Bio::Root::AutoClass;
use Bio::Biblio;
@ISA = qw(Bio::Root::AutoClass Bio::Biblio);
$EUTILS='http://www.ncbi.nlm.nih.gov/entrez/eutils/';
$WHAT='refs';
%WHAT=abbrev qw(pubmed_ids pmids ids abstracts links references refs ref);
$CHUNK_SIZE=100;

BEGIN {
  @AUTO_ATTRIBUTES=qw(access format count query auto_parse chunk_size parser session query_key webenv 
		      _current _next_buffer _next_fetch);
  @OTHER_ATTRIBUTES=qw(what file);
  %SYNONYMS=(get_collection_id=>'collection_id');
  Bio::Root::AutoClass::declare(__PACKAGE__,\@AUTO_ATTRIBUTES,\%SYNONYMS);
}

sub _init_self {
  my($self,$class,$args)=@_;
  return unless $class eq __PACKAGE__; # to prevent subclasses from re-running this
  if ($self->file) {		# file access?
				# nothing else to do here
  } elsif (!$self->session) {	# new session?
    if ($self->query) {		# new session with query?
				# yes: make session w/o query
      my $session=$self->new(-access=>$self->access);
      $self->session($session);	# connect session to current object
      $self->_esearch;		# and run query
    }
    else {			# new session w/o query
      $self->session($self);	# session points to itself
      $self->query_key(0);	# so stringification will look okay
    }
  }
  $self->auto_parse(1) unless defined $self->auto_parse;
  $self->chunk_size($CHUNK_SIZE) unless defined $self->chunk_size;
  $self->parser(new Bio::Biblio::IO(-format=>'ncbi_eutils')) if $self->auto_parse && !defined $self->parser;
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
  $what='refs' if $what_param=~/^r/i;
  $what or $what=$WHAT{lc($what_param)};
  $what='ids' if $what=~/pubmed_ids|pmids|ids/;
#  $what='refs' if $what=~/ref/;
  $what;
}
sub file {
  $_[0]->throw("Illegal to set file parameter when using online access") if @_>1 && $_[1];
  $_[0]->{'FILE'};
}
sub collection_id {
  my $self=shift @_;
  $self->{'COLLECTION_ID'}=$_[0] if @_;
  defined $self->{'COLLECTION_ID'}? $self->{'COLLECTION_ID'}: $self.' '.$self->session->webenv;
}
sub find {
  my($self,$query)=@_;
#  $self->throw("find only legal when using online access") unless $self->online;
  my $term=$query;
  $query="($query) AND $self" if $self->is_query; # AND previous query to new one
  my $result=$self->new(-access=>$self->access,
			-session=>$self->session,
			-query=>$query);
  $result->_esearch;
  $result;
}
sub get_all {
  my $self=shift @_;
  my $what=@_? _what($_[0]) || $self->throw("Invalid what $_[0]"): $self->what;
  my @results;
  if ($what eq 'refs') {
    @results=$self->_get_all_refs;
  } else {
    @results=$self->_get_all_abstracts if $what eq 'abstracts';
    @results=$self->_get_all_ids if $what eq 'ids';
    @results=$self->_get_all_links if $what eq 'links';
    @results=$self->_parse($what,@results) if $self->auto_parse;
  }
  @results;
};
sub _get_all_refs {
  my($self)=@_;
  $self->throw("Can only get complete references if auto_parse is set") unless $self->auto_parse;
  my @abstracts=$self->_abstracts_query->get_all;
  my @links=$self->_links_query->get_all;
  $self->_merge(\@abstracts,\@links); # returns list of results
}
sub _get_all_abstracts {
  my($self)=@_;
  my @results;
  for(my $start=0;$start<$self->count;$start+=$self->chunk_size) {
    my $xml=$self->_efetch($start,$self->chunk_size,'abstract');
    my($header,@xml)=$self->_split_xml($xml,'PubmedArticleSet','PubmedArticle');
    push(@results,@xml);
  }
  @results;
}
sub _get_all_ids {
  my($self)=@_;
  my @results;
  for(my $start=0;$start<$self->count;$start+=$self->chunk_size) {
    my $xml=$self->_efetch($start,$self->chunk_size,'uilist');
    my($header,@xml)=$self->_split_xml($xml,'IdList','Id');
    push(@results,@xml);
  }
  @results;
}
sub _get_all_links {
  my($self)=@_;
  my $xml=$self->_elink('prlinks');
  my($header,@xml)=$self->_split_xml($xml,'IdUrlList','IdUrlSet');
  @xml
}
sub _merge {
  my($self,$abstracts,$links)=@_;
  for (my $i=0;$i<@$abstracts;$i++) {
    my $abstract=$abstracts->[$i];
    my $link=$links->[$i];
    $abstract->linkout($link);
  }
  @$abstracts;
}
sub get_by_id {
  my($self,@ids)=@_;
  my $session=$self->session;
  my $query=$session->find(join(' ',@ids));
  my @results=$query->get_all;
  wantarray? @results: $results[0];
}
sub get_by_ids {my $self=shift @_; $self->get_by_id(@_);}

sub has_next {
  my($self)=@_;
  $self->count-$self->_current;
}
sub get_next {
  my $self=shift @_;
  my @results;
  if ($self->what eq 'refs') {
    my @abstracts=$self->_abstracts_query->get_next(@_);
    my @links=$self->_links_query->get_next(@_);
    @results=$self->_merge(\@abstracts,\@links);
  } else {
    my $n=@_? $_[0]: 1;
    $n=min($self->count-$self->_current,$n);
    if ($n>0) {
      $self->_shift_buffer($self->_next_buffer); 
      $self->_fill_buffer($n) unless $n<=$self->_buffer_count;
      @results=$self->_use_buffer($n);
      $self->_next_buffer($n);
      $self->_current($self->_current+$n);
    }
    @results=$self->_parse($self->what,@results) if $self->auto_parse;
  }
  wantarray? @results: $results[0];
}
sub get_this {
  my $self=shift @_;
  my @results;
  if ($self->what eq 'refs') {
    my @abstracts=$self->_abstracts_query->get_this(@_);
    my @links=$self->_links_query->get_this(@_);
    @results=$self->_merge(\@abstracts,\@links);
  } else {
    my $n=@_? $_[0]: 1;
    $n=min($self->count-$self->_current+$self->_next_buffer,$n);
    if ($n>0) {
      $self->_fill_buffer($n) unless $n<=$self->_buffer_count;
      @results=$self->_use_buffer($n);
    }
    @results=$self->_parse($self->what,@results) if $self->auto_parse;
  }
  wantarray? @results: $results[0];
}
sub reset {
  my $self=shift @_;
  @_ and $self->what($_[0]);
  if ($self->what eq 'refs') {
    $self->_abstracts_query->reset;
    $self->_links_query->reset;
  } else {
    $self->_current(0);
    $self->_next_fetch(0);
    $self->_next_buffer(0);
    $self->_buffer([]);
  }
  $self;
}
sub has_more {my $self=shift @_; $self->has_next(@_);}
sub get_more {my $self=shift @_; $self->get_next(@_);}
sub reset_retrieval {my $self=shift @_; $self->reset(@_);}

sub _buffer {
  @_>1 and $_[0]->{'BUFFER'}=$_[1];
  $_[0]->{'BUFFER'} || ($_[0]->{'BUFFER'}=[]);
}
sub _buffer_count {{@{$_[0]->_buffer}+0};}

sub _fill_buffer {
  my $self=shift @_;
  my $n=@_? $_[0]: $self->chunk_size;
  # fill as much as user wants if > chunk_size, else chunk_size, but
  # not more than exists 
  $n=min($self->count-$self->_next_fetch,
	 max($n-$self->_buffer_count,$self->chunk_size));
  return unless $n>0;
  my $what=$self->what;
  $self->_fill_abstracts($n) if $what eq 'abstracts';
  $self->_fill_ids($n) if $what eq 'ids';
  $self->_fill_links($n) if $what eq 'links';
}
sub _shift_buffer {
  my $self=shift @_;
  my $n=@_? $_[0]: 1;
  my $buffer=$self->_buffer;
  return unless $n>0;
  splice(@$buffer,0,$n);
}
sub _use_buffer {
  my $self=shift @_;
  my $n=@_? $_[0]: 1;
  my $buffer=$self->_buffer;
  return () unless $n>0;
  @$buffer[0..$n-1];
}
sub _fill_abstracts {
  my($self,$n)=@_;
  my $buffer=$self->_buffer;
  my $start=$self->_next_fetch;
  my $end=$start+$n;
  for(;$start<$end;$start+=$self->chunk_size) {
    my $xml=$self->_efetch($start,$self->chunk_size,'abstract');
    my($header,@xml)=$self->_split_xml($xml,'PubmedArticleSet','PubmedArticle');
    push(@$buffer,@xml);
  }
  $self->_next_fetch($self->_current+$self->_buffer_count);
}
sub _fill_ids {
  my($self,$n)=@_;
  my $buffer=$self->_buffer;
  my $start=$self->_next_fetch;
  my $end=$start+$n;
  for(;$start<$end;$start+=$self->chunk_size) {
    my $xml=$self->_efetch($start,$self->chunk_size,'uilist');
    my($header,@xml)=$self->_split_xml($xml,'IdList','Id');
    push(@$buffer,@xml);
  }
  $self->_next_fetch($self->_current+$self->_buffer_count);
}
sub _fill_links {
  my($self,$n)=@_;
  my $buffer=$self->_buffer;
  my $xml=$self->_elink('prlinks');
  my($header,@xml)=$self->_split_xml($xml,'IdUrlList','IdUrlSet');
  push(@$buffer,@xml);
  $self->_next_fetch($self->_current+$self->_buffer_count);
}
sub _abstracts_query {
  my($self)=@_;
  my $query=$self->{'_abstracts_query'} ||
    ($self->{'_abstracts_query'}=
     $self->new(-query_key=>$self->query_key,-session=>$self->session,-count=>$self->count,
		-auto_parse=>1,-chunk_size=>$self->chunk_size,
		-what=>'abstracts'));
  $query;
}
sub _links_query {
  my($self)=@_;
  my $query=$self->{'_links_query'} || 
    ($self->{'_links_query'}=
     $self->new(-query_key=>$self->query_key,-session=>$self->session,-count=>$self->count,
		-auto_parse=>1,-chunk_size=>$self->chunk_size,
		-what=>'links'));
  $query;
}

# Based on NCBI's eutils_example.pl by Oleg Khovayko 
sub _esearch {
  my($self)=@_;
  my $webenv=$self->session->webenv;
  my $query=$self->query;
  my $url=_eutil_url($EUTILS,qw(esearch db=Pubmed retmax=0 usehistory=y),
		     "term=$query",
		     $webenv? "WebEnv=$webenv": undef);			 
  my $result=get($url);
  my($count,$query_key,$new_webenv)=$result=~m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;
  $self->count($count);
  $self->query_key($query_key);
  $self->session->webenv($new_webenv) unless $webenv;
}
sub _efetch {
  my($self,$start,$chunk,$rettype)=@_;
  my $url=_eutil_url($EUTILS,qw(efetch db=Pubmed usehistory=y retmode=xml),
		     "retstart=$start","retmax=$chunk","rettype=$rettype",
		     "query_key=".$self->query_key,"WebEnv=".$self->session->webenv);
#  print "$url\n";
  my $result=get($url);
#  print "$result\n";
  $result;
}
sub _elink {
  my($self,$cmd)=@_;
  my $url=_eutil_url($EUTILS,qw(elink dbfrom=Pubmed usehistory=y retmode=xml),
		     "cmd=$cmd",
		     "query_key=".$self->query_key,"WebEnv=".$self->session->webenv);
#  print "$url\n";
  my $result=get($url);
#  print "$result\n";
  $result;
}
sub _eutil_url {
  my($location,$command,@params)=@_;
  @params=grep {defined $_} @params;
  $location='http://'.$location unless $location=~/^http:\/\//;
  $location.='/' unless $location=~/\/$/;
  $command.='.fcgi' unless $command=~/\.fcgi/;
  my $params=join('&',@params);
  uri_escape($location.$command.'?'.$params);
}
sub _split_xml {
  my($self,$xml,$set_type,$element_type)=@_;
  my($xml_version)=$xml=~/(<\?xml\s+version.*?\?>)/s;
  my($doctype)=$xml=~/(<\!DOCTYPE.*?>)/s;
  $xml=~s/^.*<$set_type>\s*//s;
  $xml=~s/\s*<\/$set_type>.*$//s;
  my $delim="<\/$element_type>\\s*";
  my(@elements)=split(/$delim/,$xml);
  @elements=map {$_."<\/$element_type>\n"} @elements;
  ($xml_version.$doctype,@elements);
}
sub _parse {
  my($self,$what,@results)=@_;
  $self->parser->what($what);
  $self->parser->parse(@results);
}
sub is_session {$_[0] eq $_[0]->session;}
sub is_query {$_[0]->session && $_[0] ne $_[0]->session;}
sub is_file {$_[0]->file;}
sub is_online {$_[0]->session? 1: undef;}
sub is_offline {$_[0]->session? undef: 1;}

sub is_connected {
  my($self)=@_;
  $_[0]->session->webenv;
}

sub min {
  if ($#_==0) {@_=@{$_[0]} if 'ARRAY' eq ref $_[0];}
  @_=grep {defined $_} @_; 
  return undef unless @_;
  if ($#_==1) {my($x,$y)=@_; return ($x<=$y?$x:$y);}
  my $min=shift @_;
  map {$min=$_ if $_<$min} @_;
  $min;
}

sub max {
  if ($#_==0) {@_=@{$_[0]} if 'ARRAY' eq ref $_[0];}
  return undef unless @_;
  if ($#_==1) {my($x,$y)=@_; return ($x>=$y?$x:$y);}
  my $max=shift @_;
  map {$max=$_ if $_>$max} @_;
  $max;
}

# Code below is from Bio::Biblio::IO::medlinexml by Martin Senger
sub fh {
  my $self = shift;
  my $class = ref($self) || $self;
  my $s = Symbol::gensym;
  tie $$s,$class,$self;
  return $s;
}
sub TIEHANDLE {
  my ($class,$val) = @_;
  return bless {'biblio' => $val}, $class;
}
sub READLINE {
  my $self = shift;
  return $self->{'biblio'}->get_next() unless wantarray;
  my (@list, $obj);
  push @list, $obj while $obj = $self->{'biblio'}->get_next();
  return @list;
}
sub PRINT {
  my $self = shift;
  $self->throw("Sorry, can't print to ncbi_eutils. It's read-only.");
}
sub CLOSE {}
sub UNTIE {}
1;
