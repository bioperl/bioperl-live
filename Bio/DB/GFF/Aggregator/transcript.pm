=head1 NAME

Bio::DB::GFF::Aggregator::transcript -- Transcript aggregator

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

Bio::DB::GFF::Aggregator::transcript is one of the default
aggregators, and was written to be compatible with the C elegans GFF
files.  It aggregates raw ""exon", "CDS", "5'UTR", "3'UTR", "polyA"
and "TSS" features into "transcript" features.  For compatibility with
the idiosyncrasies of the Sanger GFF format, it expects that the full
range of the transcript is contained in a main feature of type
"Transcript" (notice the capital "T").

Internally this module is very simple.  To override it with one that
recognizes a main feature named "gene", simply follow this
template:

 my $db = Bio::DB::GFF->new(...etc...)
 my $aggregator = Bio::DB::GFF::Aggregator->new(-method => 'transcript',
 					        -main_method => 'gene',
					        -sub_parts => ['exon','CDS']);
 $db->add_aggregator($aggregator);

=cut

package Bio::DB::GFF::Aggregator::transcript;

use strict;

use base qw(Bio::DB::GFF::Aggregator);

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
 Returns : the string "transcript"
 Args    : none
 Status  : Public

=cut

sub main_name {
  return 'transcript';
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

