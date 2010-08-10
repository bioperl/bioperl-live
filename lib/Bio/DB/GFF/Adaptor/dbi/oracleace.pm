package Bio::DB::GFF::Adaptor::dbi::oracleace;

=head1 NAME

Bio::DB::GFF::Adaptor::dbi::oracleace -- Unholy union between oracle GFF database and acedb database

=head1 SYNOPSIS

Pending

See L<Bio::DB::GFF> and L<Bio::DB::GFF::Adaptor::dbi::oracle>

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bioperl>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2002 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

use strict;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()

use base qw(Bio::DB::GFF::Adaptor::ace Bio::DB::GFF::Adaptor::dbi::oracle);

# Create a new Bio::DB::GFF::Adaptor::dbi object
sub new {
  my $class = shift;
  my $self  = $class->SUPER::new(@_);
  my ($dna_db,$acedb) = rearrange([[qw(DNADB DNA FASTA FASTA_DIR)],'ACEDB'],@_);
  if ($dna_db) {
    if (!ref($dna_db)) {
      require Bio::DB::Fasta;
      my $fasta_dir = $dna_db;
      $dna_db = Bio::DB::Fasta->new($fasta_dir);
      $dna_db or $class->throw("new(): Failed to create new Bio::DB::Fasta from files in $fasta_dir");
    } else {
      $dna_db->isa('Bio::DB::Fasta') or $class->throw("new(): $dna_db is not a Bio::DB::Fasta object");
    }
    $self->dna_db($dna_db);
  }

  if ($acedb) {
    $acedb->isa('Ace') or $class->throw("$acedb is not an acedb accessor object");
    $self->acedb($acedb);
  }
  $self;
}

sub make_object {
  my $self = shift;
  my ($class,$name,$start,$stop) = @_;

  if (my $db = $self->acedb) {

    # for Notes we just return a text, no database associated
    return $class->new(Text=>$name) if $class eq 'Note';

    # for homols, we create the indicated Protein or Sequence object
    # then generate a bogus Homology object (for future compatability??)
    if ($start ne '') {
      require Ace::Sequence::Homol;
      return Ace::Sequence::Homol->new_homol($class,$name,$db,$start,$stop);
    }

    # General case:
    my $obj = $db->class->new($class=>$name,$self->acedb);

    return $obj if defined $obj;

    # Last resort, return a Text
    return $class->new(Text=>$name);
  }

  return $self->SUPER::make_object($class,$name,$start,$stop);
}

sub get_dna {
  my $self = shift;
  my ($ref,$start,$stop,$class) = @_;
  my $dna_db = $self->dna_db or return $self->SUPER::get_dna(@_);
  return $dna_db->seq($ref,$start,$stop,$class);
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

1;
