=head1 NAME

Bio::DB::GFF::Aggregator::orf -- An aggregator for orf regions

=head1 SYNOPSIS

  use Bio::DB::GFF;

  # Open the sequence database
  my $db      = Bio::DB::GFF->new( -adaptor => 'dbi:mysql',
                                   -dsn     => 'dbi:mysql:elegans42',
				   -aggregator => ['orf','clone'],
				 );

 ---------------------------
 Aggregator method: orf
 Main method:       -none-
 Sub methods:       ORF
 ---------------------------

=head1 DESCRIPTION

Bio::DB::GFF::Aggregator::orf was written to work with the "cds"
glyph.  GFF files.  It aggregates raw "ORF" features into "coding"
features. This is basically identical to the "coding" aggregator,
except that it looks for features of type "ORF" rather than "cds".

=cut

package Bio::DB::GFF::Aggregator::orf;

use strict;
use Bio::DB::GFF::Aggregator;

use base qw(Bio::DB::GFF::Aggregator);

=head2 method

 Title   : method
 Usage   : $aggregator->method
 Function: return the method for the composite object
 Returns : the string "orf"
 Args    : none
 Status  : Public

=cut

sub method { 'orf' }

# sub require_whole_object { 1; }

=head2 part_names

 Title   : part_names
 Usage   : $aggregator->part_names
 Function: return the methods for the sub-parts
 Returns : the list "CDS"
 Args    : none
 Status  : Public

=cut

sub part_names {
  return qw(ORF);
}

1;
__END__

=head1 BUGS

None reported.


=head1 SEE ALSO

L<Bio::DB::GFF>, L<Bio::DB::GFF::coding>, 
L<Bio::DB::GFF::Aggregator>, L<Bio::Graphics::Glyph::cds>


=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

