=head1 NAME

Bio::DB::GFF::Homol -- A segment of DNA that is homologous to another

=head1 SYNOPSIS

See L<Bio::DB::GFF>.

=head1 DESCRIPTION

Bio::DB::GFF::Homol is a named subtype of Bio::DB::GFF::Segment.  It
inherits all the methods of its parent, and was created primarily to
allow for isa() queries and for compatibility with
Ace::Sequence::Homol.  

A Homol object is typically returned as the method result of the
Bio::DB::GFF::Feature-E<gt>target() method.

=head1 METHODS

=cut

package Bio::DB::GFF::Homol;
use strict;

use base qw(Bio::DB::GFF::Segment);

=head2 name

 Title   : name
 Usage   : $name = $homol->name
 Function: get the ID of the homology object
 Returns : a string
 Args    : none
 Status  : Public

=cut

sub name     { shift->refseq }

=head2 asString

 Title   : asString
 Usage   : $name = $homol->asString
 Function: same as name(), for operator overloading
 Returns : a string
 Args    : none
 Status  : Public

=cut

sub asString { shift->name }


=head2 id

 Title   : id
 Usage   : $id = $homol->id
 Function: get database ID in class:id format
 Returns : a string
 Args    : none
 Status  : Public

=cut

sub id       {
  my $self = shift;
  return "$self->{class}:$self->{name}";
}

sub new_from_segment {
  my $package   = shift;
  $package      = ref $package if ref $package;
  my $segment   = shift;
  my $new = {};
  @{$new}{qw(factory sourceseq start stop strand class ref refstart refstrand)}
    = @{$segment}{qw(factory sourceseq start stop strand class ref refstart refstrand)};
  return bless $new,__PACKAGE__;
}

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
