package Bio::Biblio::PubmedLinkout;

use vars qw(@ISA @AUTO_ATTRIBUTES @OTHER_ATTRIBUTES %SYNONYMS);
use strict;
use Bio::Root::AutoClass;
@ISA = qw(Bio::Root::AutoClass);

BEGIN {
  @AUTO_ATTRIBUTES=qw();
  @OTHER_ATTRIBUTES=qw();
  %SYNONYMS=(id=>'Id',pmid=>'Id',links=>'ObjUrl');
  Bio::Root::AutoClass::declare(__PACKAGE__,\@AUTO_ATTRIBUTES,\%SYNONYMS);
}

sub _init_self {
  my($self,$class,$args)=@_;
  return unless $class eq __PACKAGE__; # to prevent subclasses from re-running this
  $self->{'Id'}=$args->{'ID'};
  my $links=$args->{'OBJURL'} ||[];
  $links=[$links] unless 'ARRAY' eq ref $links;
  $self->{'ObjUrl'}=$links;
}

sub Id {$_[0]->{'Id'};}
sub id {$_[0]->Id;}
sub pmid {$_[0]->Id;}

sub ObjUrl {$_[0]->{'ObjUrl'};}
sub links {$_[0]->ObjUrl;}

sub full_text {
  my($self)=@_;
  my $full_text=$self->_full_text;
  $full_text and $full_text->{'Url'};
}
sub has_full_text {$_[0]->full_text? 1: undef;}

sub subject_type {
  my($self)=@_;
  $self->_full_text and 'publishers/providers';
}
sub provider {
  my($self)=@_;
  my $full_text=$self->_full_text;
  $full_text and $full_text->{'Provider'};
}
sub provider_name {
  my($self)=@_;
  my $provider=$self->provider;
  $provider and $provider->{'Name'};
}
sub provider_abbr {
  my($self)=@_;
  my $provider=$self->provider;
  $provider and $provider->{'NameAbbr'};
}
sub provider_id {
  my($self)=@_;
  my $provider=$self->provider;
  $provider and $provider->{'Id'};
}
sub provider_url {
  my($self)=@_;
  my $provider=$self->provider;
  $provider and $provider->{'Url'};
}
sub provider_icon {
  my($self)=@_;
  my $provider=$self->provider;
  $provider and $provider->{'IconUrl'};
}
sub Url {$_[0]->full_text;}
sub HasFullText {$_[0]->has_full_text;}
sub SubjectType {$_[0]->subject_type;}
sub Provider {$_[0]->provider;}
sub ProviderName {$_[0]->provider_name;}
sub ProviderNameAbbr {$_[0]->provider_abbr;}
sub ProviderId {$_[0]->provider_id;}
sub ProviderUrl {$_[0]->provider_url;}
sub ProviderIconUrl {$_[0]->provider_icon;}

sub _full_text {
  my($self)=@_;
  return $self->{'_FULL_TEXT'} if $self->{'_FULL_TEXT'};
  my $links=$self->links;
  for my $link (@$links) {
    next unless $link->{'SubjectType'} eq 'publishers/providers';
    return $self->{'_FULL_TEXT'}=$link;
  }
  undef;
}

1;
