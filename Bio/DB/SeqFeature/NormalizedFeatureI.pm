package Bio::DB::SeqFeature::NormalizedFeatureI;


=head1 NAME

Bio::DB::SeqFeature::NormalizedFeatureI -- Interface for normalized features

=head1 SYNOPSIS

none

=head1 DESCRIPTION

This is an extremely simple interface that contains a single method,
subfeatures_are_normalized(). This method returns a true value.

Bio::DB::SeqFeature::Store feature classes will inherit this interface
to flag that they are able to store subfeatures in a normalized way
such that the subfeature is actually contained in the
Bio::DB::SeqFeature::Store database and the parent feature contains
only the subfeatures primary ID.

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



sub subfeatures_are_normalized { 1 }

1;
