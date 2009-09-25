package Bio::DB::GFF::Adaptor::ace;

=head1 NAME

Bio::DB::GFF::Adaptor::ace -- ace interface (for multiple inheritance)

=head1 SYNOPSIS

Pending

See L<Bio::DB::GFF> and L<Bio::DB::GFF::Adaptor::dbi::mysql>

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bioperl>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2002 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

use strict;
use Ace;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()

sub dna_db      { 
  my $self = shift;
  my $d = $self->{dna_db};
  $self->{dna_db} = shift if @_;
  $d;
}
sub acedb      { 
  my $self = shift;
  my $d = $self->{acedb};
  $self->{acedb} = shift if @_;
  $d;
}

=head2 freshen_ace

 Title   : freshen
 Usage   : $flag = Bio::DB::GFF->freshen_ace;
 Function: Refresh internal acedb handle
 Returns : flag if correctly freshened
 Args    : none
 Status  : Public

ACeDB has an annoying way of timing out, leaving dangling database
handles.  This method will invoke the ACeDB reopen() method, which
causes dangling handles to be refreshed.  It has no effect if you are
not using ACeDB to create ACeDB objects.

=cut

sub freshen_ace {
  my $acedb = shift->acedb or return;
  $acedb->reopen();
}

1;
