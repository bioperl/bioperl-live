=head1 NAME

Bio::DB::GFF::Typename -- The name of a feature type

=head1 SYNOPSIS

  use Bio::DB::GFF;

  my $type = Bio::DB::GFF::Typename->new(similarity => 'BLAT_EST_GENOME');
  my $segment = $segment->features($type);

=head1 DESCRIPTION

Bio::DB::GFF::Typename objects encapsulate the combination of feature
method and source used by the GFF flat file format.  They can be used
in the Bio::DB::GFF modules wherever a feature type is called for.

Since there are relatively few types and many features, this module
maintains a memory cache of unique types so that two features of the
same type will share the same Bio::DB::GFF::Typename object.

=head1 METHODS

=cut

package Bio::DB::GFF::Typename;

use strict;
use overload 
  '""'     => 'asString',
  fallback => 1;

use vars '$VERSION';
$VERSION=1.1;

# cut down on the number of equivalent objects we have to create
my %OBJECT_CACHE;

=head2 new

 Title   : new
 Usage   : $type = Bio::DB::GFF::Typename->new($method,$source)
 Function: create a new Bio::DB::GFF::Typename object
 Returns : a new Bio::DB::GFF::Typename object
 Args    : method and source
 Status  : Public

=cut

sub new    {
  my $package = shift;
  my ($method,$source) = @_;
  $method ||= '';
  $source ||= '';
  return $OBJECT_CACHE{"$method:$source"} ||= bless [$method,$source],$package;
}

=head2 method

 Title   : method
 Usage   : $method = $type->method([$newmethod])
 Function: get or set the method
 Returns : a method name
 Args    : new method name (optional)
 Status  : Public

=cut

sub method {
  my $self = shift;
  my $d = $self->[0];
  $self->[0] = shift if @_;
  $d;
}


=head2 source

 Title   : source
 Usage   : $source = $type->source([$newsource])
 Function: get or set the source
 Returns : a source name
 Args    : new source name (optional)
 Status  : Public

=cut

sub source {
  my $self = shift;
  my $d = $self->[1];
  $self->[1] = shift if @_;
  $d;
}

=head2 asString

 Title   : asString
 Usage   : $string = $type->asString
 Function: get the method and source as a string
 Returns : a string in "method:source" format
 Args    : none
 Status  : Public

This method is used by operator overloading to overload the '""'
operator.

=cut

sub asString { join ':',@{shift()}};

=head2 clone

 Title   : clone
 Usage   : $new_clone = $type->clone;
 Function: clone this object
 Returns : a new Bio::DB::GFF::Typename object
 Args    : none
 Status  : Public

This method creates an exact copy of the object.

=cut

sub clone {
  my $self = shift;
  return bless [@$self],ref $self;
}

1;

=head1 BUGS

This module is still under development.

=head1 SEE ALSO

L<bioperl>, L<Bio::DB::GFF>, L<Bio::DB::RelSegment>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

1;
