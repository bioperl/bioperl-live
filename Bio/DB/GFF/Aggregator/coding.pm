=head1 NAME

Bio::DB::GFF::Aggregator::coding -- The Coding Region Aggregator

=head1 SYNOPSIS

  use Bio::DB::GFF;

  # Open the sequence database
  my $db      = Bio::DB::GFF->new( -adaptor => 'dbi:mysql',
                                   -dsn     => 'dbi:mysql:elegans42',
				   -aggregator => ['coding'],
				 );

 ------------------------------------------------------------------------
 Aggregator method: coding
 Main method:       mRNA
 Sub methods:       CDS
 ------------------------------------------------------------------------

=head1 DESCRIPTION

Bio::DB::GFF::Aggregator::coding aggregates "CDS" features into a
feature called "coding" and was written to be compatible with the
Sequence Ontology canonical gene.  The CDS features are expected to
belong to a parent of type "mRNA," but the aggregator will work even
if this isn't the case.

=cut

package Bio::DB::GFF::Aggregator::coding;

use strict;

use base qw(Bio::DB::GFF::Aggregator);

=head2 method

 Title   : method
 Usage   : $aggregator->method
 Function: return the method for the composite object
 Returns : the string "coding"
 Args    : none
 Status  : Public

=cut

sub method { 'coding' }

=head2 part_names

 Title   : part_names
 Usage   : $aggregator->part_names
 Function: return the methods for the sub-parts
 Returns : the list (CDS cds)
 Args    : none
 Status  : Public

=cut

sub part_names {
  return qw(CDS cds);
}

=head2 main_name

 Title   : main_name
 Usage   : $aggregator->main_name
 Function: return the method for the main component
 Returns : the string "mRNA"
 Args    : none
 Status  : Public

=cut

sub main_name {
  return 'mRNA';
}

1;
__END__

=head1 BUGS

None reported.


=head1 SEE ALSO

L<Bio::DB::GFF>, L<Bio::DB::GFF::Aggregator>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

