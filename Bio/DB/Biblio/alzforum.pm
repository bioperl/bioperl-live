package Bio::DB::Biblio::alzforum;
use LWP::Simple;
use URI::Escape;
use Text::Abbrev;
use Bio::Biblio::IO;

use vars qw(@ISA @AUTO_ATTRIBUTES @OTHER_ATTRIBUTES %SYNONYMS);
use vars qw($ALZ_LOCATION $ALZ_SEARCH $WHAT %WHAT);
use strict;
use Bio::Root::AutoClass;
use Bio::Biblio;
@ISA = qw(Bio::Root::AutoClass Bio::Biblio);
$ALZ_LOCATION='http://www.alzforum.org/pap/';
$ALZ_SEARCH='hdpowsearch5.asp';
$WHAT='cits';
%WHAT=abbrev qw(cits);

BEGIN {
  @AUTO_ATTRIBUTES=qw(access format query auto_parse parser session
		      _next _this _is_filled);
  @OTHER_ATTRIBUTES=qw(what file);
  %SYNONYMS=(get_collection_id=>'collection_id');
  Bio::Root::AutoClass::declare(__PACKAGE__,\@AUTO_ATTRIBUTES,\%SYNONYMS);
}

sub _init_self {
  my($self,$class,$args)=@_;
  return unless $class eq __PACKAGE__; # to prevent subclasses from re-running this
  if ($self->file) {		# file access?
    $self->_fill_buffer;	# fill buffer from file
  } elsif ($self->query) {	# it's a query object
    if (!$self->session) {	# but it needs a session
      my $session=$self->new(-access=>$self->access);
      $self->session($session);	# connect session to current object
    }
    $self->_fill_buffer;	# run query and fill our buffer
  } else {			# new session w/o query
    $self->session($self);	# session points to itself
  }
  $self->auto_parse(1) unless defined $self->auto_parse;
  $self->parser(new Bio::Biblio::IO(-format=>'alzforum')) if $self->auto_parse && !defined $self->parser;
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
sub file {
  $_[0]->throw("Illegal to set file parameter when using online access") if @_>1 && $_[1];
  $_[0]->{'FILE'};
}
sub collection_id {
  my $self=shift @_;
  $self->{'COLLECTION_ID'}=$_[0] if @_;
  defined $self->{'COLLECTION_ID'}? $self->{'COLLECTION_ID'}: "$self";
}
sub count {
  my $self=shift @_;
  return $self->{'COUNT'} if defined $self->{'COUNT'};
  # else have to calculate it...
  $self->_fill_buffer unless $self->_is_filled;	# now redundant -- new always fills
  return $self->{'COUNT'}=$self->_buffer_count;
}
sub find {
  my($self,$query)=@_;
#  $self->throw("find only legal when using online access") unless $self->online;
  my $result=$self->new(-access=>$self->access,
			-session=>$self->session,
			-query=>$query);
  $result;
}
sub get_all {
  my $self=shift @_;
  my $what=@_? _what($_[0]) || $self->throw("Invalid what $_[0]"): $self->what;
  my @results=$self->_use_buffer;
  @results=$self->_parse($self->what,@results) if $self->auto_parse;
  @results;
};
#sub get_by_id {
#  my($self,@ids)=@_;
#  my $session=$self->session;
#  my $query=$session->find(join(' ',@ids));
#  my @results=$query->get_all;
#  wantarray? @results: $results[0];
#}
#sub get_by_ids {my $self=shift @_; $self->get_by_id(@_);}

sub has_next {
  my($self)=@_;
  $self->count-$self->_next;
}
sub get_next {
  my $self=shift @_;
  my @results;
  my $n=@_? $_[0]: 1;
  $n=min($self->count-$self->_next,$n);
  if ($n>0) {
    @results=$self->_use_buffer($self->_next,$n);
    $self->_this($self->_next);	# set _this before bumping _next!
    $self->_next($self->_next+$n);
  }
  @results=$self->_parse($self->what,@results) if $self->auto_parse;
  wantarray? @results: $results[0];
}
sub get_this {
  my $self=shift @_;
  my @results;
  my $n=@_? $_[0]: 1;
  $n=min($self->count-$self->_this,$n);
  if ($n>0) {
    @results=$self->_use_buffer($self->_this,$n);
  }
  @results=$self->_parse($self->what,@results) if $self->auto_parse;
  wantarray? @results: $results[0];
}
sub reset {
  my $self=shift @_;
  @_ and $self->what($_[0]);
  $self->_next(0);
  $self->_this(0);
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
  # with Alzforum, have to get all results in one shot
  my $what=$self->what;
  $self->_fill_cits if $what eq 'cits';
  $self->_is_filled(1);
}
sub _use_buffer {
  my $self=shift @_;
  my $buffer=$self->_buffer;
  my($start,$n)=@_;
  return @$buffer[$start..$start+$n-1] if $start || $n;
  return @$buffer;
}
sub _fill_cits {
  my($self)=@_;
  my $buffer=[];
  my $xml=$self->_search;
  my($header,@$buffer)=$self->_split_xml($xml,'Rows');
  $self->_buffer($buffer);
}
sub _search {
  my($self)=@_;
  my $query=$self->query;
  my $url=_eutil_url($ALZ_LOCATION,$ALZ_SEARCH,"$query");
  my $xml=get($url);
  $xml;
}
sub _eutil_url {
  my($location,$command,@params)=@_;
  @params=grep {defined $_} @params;
  $location='http://'.$location unless $location=~/^http:\/\//;
  $location.='/' unless $location=~/\/$/;
  my $params=join('&',@params);
  uri_escape($location.$command.'?'.$params);
}
sub _split_xml {
  my($self,$xml,$set_type)=@_;
  my($xml_version)=$xml=~/(<\?xml\s+version.*?\?>)/s;
  my($doctype)=$xml=~/(<\!DOCTYPE.*?>)/s;
  $xml=~s/^.*<$set_type>\s*//s;
  $xml=~s/\s*<\/$set_type>.*$//s;
  my $delim="\/>\s*";
  my(@elements)=split(/$delim/,$xml);
  @elements=map {$_."\/>\n"} @elements;
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
  $self->throw("Sorry, can't print to alzforum. It's read-only.");
}
sub CLOSE {}
sub UNTIE {}
1;
