#!/usr/bin/perl

use strict;
use lib '.','./blib/lib';
use Bio::DB::BioFetch;
use Bio::Graphics;

my $accession = shift or die <<END;
Usage: $0 <accession> [start] [stop]
   Render a GenBank/EMBL accession into drawable form.
   Return as a GIF or PNG image on standard output.

   If start and stop are specified, then that segment
   will be displayed.

   To view, pipe to a viewer program as shown below.

Example to try:
   render_entry CEF58D5 | display -
END
;

my ($start,$stop) = @ARGV;

my $bf = eval {require Bio::DB::FileCache}
  ? Bio::DB::FileCache->new(-seqdb=>Bio::DB::BioFetch->new,
			    -file =>'/usr/tmp/my_seq_cache',
			    -keep =>1)
  : Bio::DB::Biofetch->new;

my $seq = $bf->get_Seq_by_id($accession);

my @features = $seq->all_SeqFeatures;
my @CDS      = grep {$_->primary_tag eq 'CDS'}  @features;
my @gene     = grep {$_->primary_tag eq 'gene'} @features;
my @tRNAs    = grep {$_->primary_tag eq 'tRNA'} @features;

warn "rendering...\n";
$start = $seq->start unless defined $start;
$stop  = $seq->end   unless defined $stop;

my $panel = Bio::Graphics::Panel->new(
				      -offset  => $start,
				      -length  => $stop - $start + 1,
				      -width   => 1000,
				      );
$panel->add_track(arrow => $seq,
		  -bump => 0,
		  -double=>1,
		  -tick => 2);

$panel->add_track(transcript2  => \@gene,
		  -bgcolor    =>  'blue',
		  -fgcolor    =>  'black',
		  -key        => 'Genes',
		  -bump       =>  +1,
		  -height     =>  10,
		  -label      => \&gene_label,
		  -description=> \&gene_description,
		 );

$panel->add_track(transcript2  => \@CDS,
		  -bgcolor    =>  'cyan',
		  -fgcolor    =>  'black',
		  -key        => 'CDS',
		  -bump       =>  +1,
		  -height     =>  10,
		  -label      => \&gene_label,
		  -description=> \&gene_description,
		 );

$panel->add_track(generic    => \@tRNAs,
		  -bgcolor   =>  'red',
		  -fgcolor   =>  'black',
		  -key       => 'tRNAs',
		  -bump      =>  +1,
		  -height    =>  8,
		  -label      => \&gene_label,
		 );

my $gd = $panel->gd;

print $gd->can('png') ? $gd->png : $gd->gif;

sub gene_label {
  my $feature = shift;
  my @notes;
  foreach (qw(product gene)) {
    next unless $feature->has_tag($_);
    @notes = $feature->each_tag_value($_);
    last;
  }
  $notes[0];
}
sub gene_description {
  my $feature = shift;
  my @notes;
  foreach (qw(note)) {
    next unless $feature->has_tag($_);
    @notes = $feature->each_tag_value($_);
    last;
  }
  return unless @notes;
  substr($notes[0],30) = '...' if length $notes[0] > 30;
  $notes[0];
}
