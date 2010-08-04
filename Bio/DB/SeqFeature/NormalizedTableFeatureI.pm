package Bio::DB::SeqFeature::NormalizedTableFeatureI;


=head1 NAME

Bio::DB::SeqFeature::NormalizedTableFeatureI -- Interface for normalized features whose hierarchy is stored in a table

=head1 SYNOPSIS

none

=head1 DESCRIPTION

This is an extremely simple interface that contains a single method,
subfeatures_are_stored_in_a_table(). This method returns a true value.

Bio::DB::SeqFeature::Store feature classes will inherit this interface
to flag that in addition to being able to store features in a
normalized way, they will use the Bio::DB::SeqFeature::Store database
to record their parent/child relationships. A class that inherits from
NormalizedTableFeatureI will also inherit from NormalizedFeatureI, as
the first is a subclass of the second.

=head1 BUGS

None, but the whole class design might be flawed.

=head1 SEE ALSO

L<bioperl>,
L<Bio::DB::SeqFeature>,
L<Bio::DB::SeqFeature::Store>,
L<Bio::DB::SeqFeature::Segment>,
L<Bio::DB::SeqFeature::GFF3Loader>,
L<Bio::DB::SeqFeature::Store::DBI::mysql>,
L<Bio::DB::SeqFeature::Store::bdb>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

use base 'Bio::DB::SeqFeature::NormalizedFeatureI';

sub subfeatures_are_stored_in_a_table { 1 }

1;
