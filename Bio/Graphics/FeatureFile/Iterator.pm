package Bio::Graphics::FeatureFile::Iterator;

# $Id$

=head1 NAME

Bio::Graphics::FeatureFile::Iterator -- Iterator across a Bio::Graphics::FeatureFile

=head1 SYNOPSIS

 use Bio::Graphics::FeatureFile;
 my $data  = Bio::Graphics::FeatureFile->new(-file => 'features.txt');
 my $iterator = $data->get_seq_stream;
 while (my $feature = $iterator->next_seq) {
   print $feature->display_id,"\t",$feature->start,"\t",$feature->end,"\n";
 }

=head1 DESCRIPTION

This is a Bio::SeqIO-like object that recognizes the next_seq() and
next_feature() methods.  The two methods are synonymous.

There is also a rewind() method which will start iterating from the
beginning again.

=head1 SEE ALSO

L<Bio::Graphics::Panel>,
L<Bio::Graphics::Glyph>,
L<Bio::Graphics::Feature>,
L<Bio::Graphics::FeatureFile>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

sub new {
  my $package = shift;
  return bless {features => shift,
		index    => 0},$package;
}

sub next_seq {
  my $self = shift;
  return unless $self->{features};
  return $self->{features}[$self->{index}++];
}

*next_features = \&next_seq;

sub rewind {
  my $self = shift;
  $self->{index} = 0;
}

1;
