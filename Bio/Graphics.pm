package Bio::Graphics;

use Bio::Graphics::Panel;
use strict;

use vars '$VERSION';
$VERSION = 1.2003;

1;

=head1 NAME

Bio::Graphics - Generate GD images of Bio::Seq objects

=head1 SYNOPSIS

 # This script generates a PNG picture of a 10K region containing a
 # set of red features and a set of blue features. Call it like this:
 #         red_and_blue.pl > redblue.png
 # you can now view the picture with your favorite image application


 # This script parses a GenBank or EMBL file named on the command
 # line and produces a PNG rendering of it.  Call it like this:
 # render.pl my_file.embl | display -

 use strict;
 use Bio::Graphics;
 use Bio::SeqIO;

 my $file = shift                       or die "provide a sequence file as the argument";
 my $io = Bio::SeqIO->new(-file=>$file) or die "could not create Bio::SeqIO";
 my $seq = $io->next_seq                or die "could not find a sequence in the file";

 my @features = $seq->all_SeqFeatures;

 # sort features by their primary tags
 my %sorted_features;
 for my $f (@features) {
   my $tag = $f->primary_tag;
   push @{$sorted_features{$tag}},$f;
 }

 my $wholeseq = Bio::SeqFeature::Generic->new(-start=>1,-end=>$seq->length);

 my $panel = Bio::Graphics::Panel->new(
				      -length    => $seq->length,
 				      -key_style => 'between',
 				      -width     => 800,
 				      -pad_left  => 10,
 				      -pad_right => 10,
 				      );
 $panel->add_track($wholeseq,
 		  -glyph => 'arrow',
 		  -bump => 0,
 		  -double=>1,
 		  -tick => 2);

 $panel->add_track($seq,
 		  -glyph  => 'generic',
 		  -bgcolor => 'blue',
 		  -label  => 1,
 		 );

 # general case
 my @colors = qw(cyan orange blue purple green chartreuse magenta yellow aqua);
 my $idx    = 0;
 for my $tag (sort keys %sorted_features) {
   my $features = $sorted_features{$tag};
   $panel->add_track($features,
 		    -glyph    =>  'generic',
 		    -bgcolor  =>  $colors[$idx++ % @colors],
 		    -fgcolor  => 'black',
 		    -font2color => 'red',
 		    -key      => "${tag}s",
 		    -bump     => +1,
 		    -height   => 8,
 		    -label    => 1,
 		    -description => 1,
 		   );
 }

 print $panel->png;
 exit 0;

=head1 DESCRIPTION

Please see L<Bio::Graphics::Panel> for the full interface.

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

