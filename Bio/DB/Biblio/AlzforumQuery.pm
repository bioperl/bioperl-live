package Bio::DB::Biblio::AlzforumQuery;
use overload 
  q{""} => sub {$_[0]->as_string;},
  'cmp' => sub {my $res=(overload::StrVal($_[0]) cmp overload::StrVal($_[1]));
		$_[2]? -$res: $res},
  'fallback'=>1,
  'nomethod' => sub {$_[0]};

use vars qw(@ISA @AUTO_ATTRIBUTES @OTHER_ATTRIBUTES %SYNONYMS);
use strict;
use Bio::Root::AutoClass;
@ISA = qw(Bio::Root::AutoClass);

BEGIN {
  @AUTO_ATTRIBUTES=qw(
		      phrase1
		      phrase2
		      andOr
		      author
		      commentator
		      journal
		      PubDateLowMonth
		      PubDateLowYear
		      PubDateHighMonth
		      PubDateHighYear
		      commented
		      ARFRec
		      mileStone
		      sort
		     );
  @OTHER_ATTRIBUTES=qw();
  %SYNONYMS=();
  Bio::Root::AutoClass::declare(__PACKAGE__,\@AUTO_ATTRIBUTES,\%SYNONYMS,'lower');
}

sub _init_self {
  my($self,$class,$args)=@_;
  return unless $class eq __PACKAGE__; # to prevent subclasses from re-running this
}
sub as_string {
  my($self)=@_;
  my @query;
  my @fields=@AUTO_ATTRIBUTES;
  for my $field (@fields) {
    my $value=$self->$field;
    next unless defined $value;
    push(@query,"$field=$value");
  }
  my $query=join('&',@query);
  $query;
}
1;
