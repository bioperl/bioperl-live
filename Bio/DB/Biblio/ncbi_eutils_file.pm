package Bio::DB::Biblio::ncbi_eutils_file;

use vars qw(@ISA @AUTO_ATTRIBUTES @OTHER_ATTRIBUTES %SYNONYMS);
use strict;
use Bio::Root::AutoClass;
use FileHandle;
use Bio::DB::Biblio::ncbi_eutils;
@ISA = qw(Bio::DB::Biblio::ncbi_eutils);

BEGIN {
  @AUTO_ATTRIBUTES=qw(filehandle);
  @OTHER_ATTRIBUTES=qw();
  %SYNONYMS=();
  Bio::Root::AutoClass::declare(__PACKAGE__,\@AUTO_ATTRIBUTES,\%SYNONYMS);
}
sub _init_self {
  my($self,$class,$args)=@_;
  return unless $class eq __PACKAGE__; # to prevent subclasses from re-running this
  $self->what('abstracts');	# change default, since 'refs' doesn't work here
}
sub what {
  my $self=shift @_;
  my $what=$self->SUPER::what(@_);
  $self->throw("Cannot set what to 'refs' unless using online access") if $what=~/ref/;
  $what;
}
sub file {
  my $self=shift @_;
  if (@_) {
    my $file=shift @_;
    $self->throw("Cannot read file parameter $file") if $file && !-r $file;
    $self->{'FILE'}=$file;
  }
  $self->{'FILE'};
}
sub count {
  my $self=shift @_;
  return $self->{'COUNT'}=$_[0] if @_;
  return $self->{'COUNT'} if defined $self->{'COUNT'};
  # else have to get the data before we can count it
  $self->_fill;			# this sets count
  $self->{'COUNT'};
}
sub find {
  my($self,$query)=@_;
  $self->throw("Cannot run find unless using online access");
}
sub reset {
  my $self=shift @_;
  @_ and $self->what($_[0]);
  $self->_fill;
  $self;
}

sub _get_all_abstracts {
  my($self)=@_;
  $self->_fill_abstracts;
  @{$self->_buffer};
}
sub _get_all_ids {
  my($self)=@_;
  $self->_fill_ids;
  @{$self->_buffer};
}
sub _get_all_links {
  my($self)=@_;
  $self->_fill_links;
  @{$self->_buffer};
}

# for file access, load entire file into buffer
sub _fill_abstracts {$_[0]->_fill('PubmedArticleSet','PubmedArticle');}
sub _fill_ids {$_[0]->_fill('IdList','Id');}
sub _fill_links {$_[0]->_fill('IdUrlList','IdUrlSet');}
sub _fill {
  my($self,$set_type,$element_type)=@_;
  if (@_==1) {
    ($set_type,$element_type)=
      $self->what eq 'abstracts'? ('PubmedArticleSet','PubmedArticle'):
	$self->what eq 'ids'? ('IdList','Id'):
	  $self->what eq 'links'? ('IdUrlList','IdUrlSet'): 
	    $self->throw("Unrecognized 'what' ".$self->what.": should have been caught earlier!");
  }
  my $buffer=[];
  my $fh;
  if ($self->file) {
    $fh=new FileHandle $self->file,"r";
    $self->filehandle($fh);
  } else {
    $fh=$self->filehandle;
    my $fh_okay=1;
    if ('GLOB' eq ref $fh) {
      $fh_okay=seek($fh,0,0) if eof($fh);
    } elsif (UNIVERSAL::isa($fh,'FileHandle')) { 
      $fh_okay=$fh->seek(0,0) if $fh->eof()
    }
    $self->throw("Unable to reset filehandle") unless $fh_okay;
  }
  my $xml=join('',<$fh>);
  my($header,@$buffer)=$self->_split_xml($xml,$set_type,$element_type);
  $self->_buffer($buffer);
  $self->_next_fetch($self->_buffer_count);
  $self->_current(0);
  $self->_next_buffer(0);
  $self->count($self->_buffer_count);
}
1;
