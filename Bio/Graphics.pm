package Bio::Graphics;

use Bio::Graphics::Panel;
use strict;

use vars '$VERSION';
$VERSION = '0.98';


1;

=head1 NAME

Bio::Graphics - Generate GD images of Bio::Seq objects

=head1 SYNOPSIS

  use Bio::Graphics;

  # get a set of Bio::SeqFeature objects, in this case from AcePerl
  use Ace::Sequence;
  my $db     = Ace->connect(-host=>'brie2.cshl.org',-port=>2005) or die;
  my $cosmid = Ace::Sequence->new(-seq=>'Y16B4A',
				  -db=>$db,-start=>-15000,-end=>15000) or die;
  my @transcripts = $cosmid->transcripts;

  # let the drawing begin...
  my $panel = Bio::Graphics::Panel->new(
				      -segment => $cosmid,
				      -width  => 800
				     );


  $panel->add_track(arrow => $cosmid,
 		  -bump => 0,
 		  -tick=>2);

  $panel->add_track(transcript => \@transcripts,
 		    -bgcolor   =>  'wheat',
 		    -fgcolor   =>  'black',
                    -key       => 'Curated Genes',
 		    -bump      =>  +1,
 		    -height    =>  10,
 		    -label     =>  1);

  my $boxes = $panel->boxes;
  print $panel->png;

=head1 DESCRIPTION

Please see Bio::Graphics::Panel for the full API.

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

Lincoln Stein <lstein@cshl.org>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

