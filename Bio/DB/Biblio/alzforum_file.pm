package Bio::DB::Biblio::alzforum_file;

use vars qw(@ISA @AUTO_ATTRIBUTES @OTHER_ATTRIBUTES %SYNONYMS);
use strict;
use Bio::Root::AutoClass;
use FileHandle;
use Bio::DB::Biblio::alzforum;
@ISA = qw(Bio::DB::Biblio::alzforum);

BEGIN {
  @AUTO_ATTRIBUTES=qw(filehandle);
  @OTHER_ATTRIBUTES=qw();
  %SYNONYMS=();
  Bio::Root::AutoClass::declare(__PACKAGE__,\@AUTO_ATTRIBUTES,\%SYNONYMS);
}
sub _init_self {
  my($self,$class,$args)=@_;
  return unless $class eq __PACKAGE__; # to prevent subclasses from re-running this
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
sub find {
  my($self,$query)=@_;
  $self->throw("Cannot run find unless using online access");
}
sub _search {
  my($self)=@_;
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
  $xml;
}
1;
