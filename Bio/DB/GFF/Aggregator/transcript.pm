=head1 NAME

Bio::DB::GFF::Aggregator::transcript -- Transcript aggregator

=head1 SYNOPSIS

  use Bio::DB::GFF;

  # Open the sequence database
  my $db      = Bio::DB::GFF->new( -adaptor => 'dbi:mysql',
                                   -dsn     => 'dbi:mysql:elegans42',
				   -aggregator => ['transcript','clone'],
				 );


=head1 DESCRIPTION

Bio::DB::GFF::Aggregator::transcript is one of the default
aggregators, and was written to be compatible with the C elegans GFF
files.  It aggregates raw "Sequence", "exon", "CDS", "5'UTR", "3'UTR",
"polyA" and "TSS" features into "transcript" features.

=cut

package Bio::DB::GFF::Aggregator::transcript;

use strict;
use Bio::DB::GFF::Aggregator;

use vars qw($VERSION @ISA);
@ISA = qw(Bio::DB::GFF::Aggregator);

$VERSION = '0.10';

=head2 method

 Title   : method
 Usage   : $aggregator->method
 Function: return the method for the composite object
 Returns : the string "transcript"
 Args    : none
 Status  : Public

=cut

sub method { 'transcript' }

=head2 part_names

 Title   : part_names
 Usage   : $aggregator->part_names
 Function: return the methods for the sub-parts
 Returns : the list "intron", "exon" and "CDS"
 Args    : none
 Status  : Public

=cut

sub part_names {
  return qw(exon CDS 5'UTR 3'UTR TSS PolyA);
}

=head2 main_name

 Title   : main_name
 Usage   : $aggregator->main_name
 Function: return the method for the main component
 Returns : the string "Sequence"
 Args    : none
 Status  : Public

=cut

sub main_name {
#  return 'Sequence:curated';
  return 'Sequence';
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

