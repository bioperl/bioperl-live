package Bio::Graphics;

use Bio::Graphics::Panel;
use strict;

use vars '$VERSION';
$VERSION = '1.04';

1;

=head1 NAME

Bio::Graphics - Generate GD images of Bio::Seq objects

=head1 SYNOPSIS

  use Bio::Graphics;
  use Bio::DB::BioFetch;  # or some other Bio::SeqI generator
  # get a Bio::SeqI object somehow
  my $bf     = Bio::DB::BioFetch->new;
  my $cosmid = $bf->getSeq_by_id('CEF58D5');

  my @features = $seq->all_SeqFeatures;
  my @CDS      = grep {$_->primary_tag eq 'CDS'}  @features;
  my @gene     = grep {$_->primary_tag eq 'gene'} @features;
  my @tRNAs    = grep {$_->primary_tag eq 'tRNA'} @features;
  # let the drawing begin...
  my $panel = Bio::Graphics::Panel->new(
				      -segment => $cosmid,
				      -width  => 800
				     );

  $panel->add_track(arrow => $cosmid,
	  	   -bump => 0,
		   -double=>1,
		   -tick => 2);

  $panel->add_track(transcript  => \@gene,
		   -bgcolor    =>  'blue',
		   -fgcolor    =>  'black',
		   -key        => 'Genes',
		   -bump       =>  +1,
		   -height     =>  10,
		   -label      => 1,
		   -description=> 1
		 ) ;

  $panel->add_track(transcript2  => \@CDS,
		    -bgcolor    =>  'cyan',
		    -fgcolor    =>  'black',
		    -key        => 'CDS',
		    -bump       =>  +1,
		    -height     =>  10,
		    -label      => \&cds_label,
		    -description=> \&cds_description,
		 );

  $panel->add_track(generic    => \@tRNAs,
		    -bgcolor   =>  'red',
		    -fgcolor   =>  'black',
		    -key       => 'tRNAs',
		    -bump      =>  +1,
		    -height    =>  8,
		    -label      => 1,
		   );

  my $gd = $panel->gd;
  print $gd->can('png') ? $gd->png : $gd->gif;

  # these are callbacks used to generate nice labels and descriptions for
  # the features...
  sub cds_label {
    my $feature = shift;
    my @notes;
    foreach (qw(product gene)) {
      next unless $feature->has_tag($_);
      @notes = $feature->each_tag_value($_);
      last;
    }
    $notes[0];
  }

  sub cds_description {
    my $feature = shift;
    my @notes = $feature->each_tag_value('notes')
                if $feature->has_tag('notes');
    return unless @notes;
    substr($notes[0],30) = '...' if length $notes[0] > 30;
    $notes[0];
  }

=head1 DESCRIPTION

Please see L<Bio::Graphics::Panel> for the full API.

=head1 SEE ALSO

L<Bio::Graphics::Panel>,
L<Bio::Graphics::Glyph>,
L<Bio::SeqI>,
L<Bio::SeqFeatureI>,
L<Bio::Das>,
L<Bio::DB::GFF::Feature>,
L<Ace::Sequence>,
L<GD>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

