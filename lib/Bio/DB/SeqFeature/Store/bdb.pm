package Bio::DB::SeqFeature::Store::bdb;


=head1 NAME

Bio::DB::SeqFeature::Store::bdb - fetch and store objects from a BerkeleyDB

=head1 DESCRIPTION

This is a partial implementation -- just enough has been implemented so that we can
fetch and store objects. It is used as a temporary failsafe store by the GFF3Loader module

=cut

use strict;
use base 'Bio::DB::SeqFeature::Store';
use Bio::DB::GFF::Util::Rearrange 'rearrange';
use DB_File;
use Fcntl qw(O_RDWR O_CREAT);
use File::Temp 'tempdir';
use File::Path 'rmtree';

###
# object initialization
#
sub init {
  my $self          = shift;
  my ($directory,
      $is_temporary) = rearrange([['DSN','DIR','DIRECTORY'],
				 ['TMP','TEMP','TEMPORARY']
				 ],@_);
  $directory ||= $is_temporary ? File::Spec->tmpdir : '.';
  $directory = tempdir(__PACKAGE__.'_XXXXXX',TMPDIR=>1,CLEANUP=>1,DIR=>$directory) if $is_temporary;
  -d $directory && -w _ or $self->throw("Can't write into the directory $directory");
  $self->default_settings;
  $self->directory($directory);
  $self->temporary($is_temporary);

  my %h;
  tie (%h,'DB_File',$self->path,O_RDWR|O_CREAT,0666,$DB_HASH) or $self->throw("Couldn't tie: $!");
  $self->db(\%h);
  $h{'.next_id'} ||= 1;
}

sub _store {
  my $self = shift;
  my $indexed = shift;
  my $db   = $self->db;
  my $count = 0;
  for my $obj (@_) {
    my $primary_id = $obj->primary_id;
    $primary_id    = $db->{'.next_id'}++ unless defined $primary_id;
    $db->{$primary_id} = $self->freeze($obj);
    $obj->primary_id($primary_id);
    $count++;
  }
  $count;
}

sub _update {
  my $self = shift;
  my ($object,$primary_id) = @_;
  my $db = $self->db;
  $self->throw("$object is not in database") unless exists $db->{$primary_id};
  $db->{$primary_id} = $self->freeze($object);
}

sub _fetch {
  my $self = shift;
  my $id   = shift;
  my $db = $self->db;
  my $obj = $self->thaw($db->{$id},$id);
  $obj;
}

sub db {
  my $self = shift;
  my $d = $self->setting('db');
  $self->setting(db=>shift) if @_;
  $d;
}

sub directory {
  my $self = shift;
  my $d = $self->setting('directory');
  $self->setting(directory=>shift) if @_;
  $d;
}

sub temporary {
  my $self = shift;
  my $d = $self->setting('temporary');
  $self->setting(temporary=>shift) if @_;
  $d;
}

sub path {
  my $self = shift;
  return $self->directory .'/' . 'feature.bdb';
}

sub DESTROY {
  my $self = shift;
  my $db   = $self->db;
  untie %$db;
  rmtree($self->directory,0,1) if $self->temporary;
}

1;
