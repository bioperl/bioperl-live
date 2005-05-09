#!/usr/bin/perl

use strict;
use lib '.','../blib/lib';
use Bio::DB::BioFetch;
use Bio::Graphics;

my $accession = shift;
if (!defined $accession || $accession =~ /^-/) { die <<END; }
Usage: $0 <accession> [start] [stop]
   Render a GenBank/EMBL accession into drawable form.
   Return as a GIF or PNG image on standard output.

   If start and stop are specified, then that segment
   will be displayed.

   To view, pipe to a viewer program as shown below.

Example to try:
   render_sequence.pl CEF58D5 | display -

By default, will look for accession in the "embl" namespace.  To
choose other namespaces, use these formats:

  swall:CEF58D5
  refseq:NC_001320

END
;

my ($start,$stop) = @ARGV;

my $db = 'embl';
if ($accession =~ /^(embl|swall|refseq):(.+)/) {
  $db        = $1;
  $accession = $2;
}

my $bf = eval {require Bio::DB::FileCache}
  ? Bio::DB::FileCache->new(-seqdb=>Bio::DB::BioFetch->new(-db=>$db),
			    -file =>'/usr/tmp/my_seq_cache',
			    -keep =>1)
  : Bio::DB::BioFetch->new(-db=>$db);

warn "fetching...\n";
my $seq = $bf->get_Seq_by_id($accession);

my @features = $seq->all_SeqFeatures;
my @CDS      = grep {$_->primary_tag eq 'CDS'}  @features;
my @gene     = grep {$_->primary_tag eq 'gene'} @features;
my @tRNAs    = grep {$_->primary_tag eq 'tRNA'} @features;

warn "rendering...\n";
$start = 1 unless defined $start;
$stop  = $seq->length   unless defined $stop;

my $panel = Bio::Graphics::Panel->new(
				      -offset  => $start,
				      -length  => $stop - $start + 1,
				      -width   => 1000,
				      );
$panel->add_track(arrow => 
		  Bio::Graphics::Feature->new(-start => $start,
					      -stop   => $stop,
					      -name   => $seq->display_id),
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
      next unless $feature->can('has_tag') && $feature->has_tag($_);
      @notes = $feature->each_tag_value($_);
      last;
  }
  $notes[0];
}
sub gene_description {
  my $feature = shift;
  my @notes;
  
  foreach (qw(note)) {
      next unless $feature->can('has_tag') && $feature->has_tag($_);
      @notes = $feature->each_tag_value($_);
      last;
  }
  return unless @notes;
  substr($notes[0],30) = '...' if length $notes[0] > 30;
  $notes[0];
}
