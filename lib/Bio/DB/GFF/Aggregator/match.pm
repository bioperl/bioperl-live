=head1 NAME

Bio::DB::GFF::Aggregator::match -- Match aggregator

=head1 SYNOPSIS

  use Bio::DB::GFF;

  # Open the sequence database
  my $db      = Bio::DB::GFF->new( -adaptor => 'dbi:mysql',
                                   -dsn     => 'dbi:mysql:elegans42',
				   -aggregator => ['match'],
				 );

 -------------------------------------------------
 Aggregator method: match
 Main method:       match
 Sub methods:       similarity HSP
 -------------------------------------------------

=head1 DESCRIPTION

This aggregator is used for Sequence Ontology-compatible gapped
alignments, in which there is a single top-level alignment called
"match" and a series of subalignments called either "similarity" or
"HSP".

Also see the "alignment" aggregator.

=cut

package Bio::DB::GFF::Aggregator::match;

use strict;

use base qw(Bio::DB::GFF::Aggregator);

=head2 method

 Title   : method
 Usage   : $aggregator->method
 Function: return the method for the composite object
 Returns : the string "match"
 Args    : none
 Status  : Public

=cut

sub method { 'match' }

=head2 part_names

 Title   : part_names
 Usage   : $aggregator->part_names
 Function: return the methods for the sub-parts
 Returns : the list "similarity", "HSP"
 Args    : none
 Status  : Public

=cut

sub part_names {
  return qw(similarity HSP);
}

=head2 main_name

 Title   : main_name
 Usage   : $aggregator->main_name
 Function: return the method for the main component
 Returns : the string "match"
 Args    : none
 Status  : Public

=cut

sub main_name {
  return 'match';
}

sub require_whole_object {1}

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

