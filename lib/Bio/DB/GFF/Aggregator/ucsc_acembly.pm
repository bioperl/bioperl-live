=head1 NAME

Bio::DB::GFF::Aggregator::ucsc_acembly -- UCSC acembly aggregator

=head1 SYNOPSIS

  use Bio::DB::GFF;

  # Open the sequence database
  my $db      = Bio::DB::GFF->new( -adaptor => 'dbi:mysql',
                                   -dsn     => 'dbi:mysql:elegans42',
				   -aggregator => ['transcript','clone'],
				 );

 -------------------------------------------------
 Aggregator method: transcript
 Main method:       transcript
 Sub methods:       exon CDS 5'UTR 3'UTR TSS PolyA
 -------------------------------------------------

=head1 DESCRIPTION

L<Bio::DB::GFF::Aggregator::transcript>

=cut

package Bio::DB::GFF::Aggregator::ucsc_acembly;

use strict;

use base qw(Bio::DB::GFF::Aggregator);


=head2 method

 Title   : method
 Usage   : $aggregator->method
 Function: return the method for the composite object
 Returns : the string "acembly"
 Args    : none
 Status  : Public

=cut

sub method { 'acembly' }

=head2 part_names

 Title   : part_names
 Usage   : $aggregator->part_names
 Function: return the methods for the sub-parts
 Returns : empty list
 Args    : none
 Status  : Public

=cut

sub part_names {
  return ();
}

=head2 main_name

 Title   : main_name
 Usage   : $aggregator->main_name
 Function: return the method for the main component
 Returns : the string "transcript:acembly"
 Args    : none
 Status  : Public

=cut

sub main_name {
  return 'transcript:acembly';
}

1;
__END__

=head1 BUGS

None reported.


=head1 SEE ALSO

L<Bio::DB::GFF>, L<Bio::DB::GFF::Aggregator>

=head1 AUTHOR

Allen Day E<lt>allenday@ucla.eduE<gt>.

Copyright (c) 2002 Allen Day, University of California, Los Angeles.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

