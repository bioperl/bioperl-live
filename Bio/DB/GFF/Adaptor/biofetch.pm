package Bio::DB::GFF::Adaptor::biofetch;

=head1 NAME

Bio::DB::GFF::Adaptor::biofetch -- Cache BioFetch objects in a Bio::DB::GFF database

=head1 SYNOPSIS

Proof of principle.  Not for production use.

=head1 DESCRIPTION

This adaptor is a proof-of-principle.  It is used to fetch BioFetch
sequences into a Bio::DB::GFF database (currently uses a hard-coded
mysqlopt database) as needed.  This allows the Generic Genome Browser
to be used as a Genbank/EMBL browser.

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright 2002 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

use strict;
use Bio::DB::GFF::Adaptor::dbi::mysqlopt;
use Bio::DB::BioFetch;

use vars qw($VERSION @ISA);
@ISA = qw(Bio::DB::GFF::Adaptor::dbi::mysqlopt);
$VERSION = 0.10;

sub segment {
  my $self = shift;
  my @segments = $self->SUPER::segment(@_);

  if (!@segments) {
    my $refclass = $self->refclass;

    my %args = $self->setup_segment_args(@_);
    if ($args{-class} && $args{-class} =~ /$refclass/oi) {
      return unless $self->load_from_embl('embl'=>$args{-name});
      @segments = $self->SUPER::segment(@_);
    } elsif ($args{-class} && $args{-class} =~ /refseq|swall|embl/i) { #hack to get refseq names
      return unless $self->load_from_embl(lc($args{-class})=>$args{-name});
      $args{-class} = $self->refclass;
      @segments = $self->SUPER::segment(%args);
    }
  }

  $self->_multiple_return_args(@segments);
}

# default is to return 'Sequence' as the class of all references
sub refclass {
  my $self = shift;
  my $refname = shift;
  'Accession';
}

sub load_from_embl {
  my $self = shift;
  my $db   = shift;
  my $acc  = shift or $self->throw('Must provide an accession ID');
  my $biofetch = $self->{_biofetch}{$db} ||= Bio::DB::BioFetch->new(-db=>$db);
  my $seq  = eval {$biofetch->get_Seq_by_id($acc)} or return;
  my $refclass = $self->refclass;

  # begin loading
  $self->setup_load();

  # first synthesize the entry for the top-level feature
  my @aliases;
  foreach ($seq->accession,$seq->get_secondary_accessions) {
    next if lc($_) eq lc($acc);
    push @aliases,[Alias => $_];
  }
  $self->load_gff_line(
		       {
			ref    => $acc,
			class  => $refclass,
			source => 'EMBL',
			method => 'origin',
			start  => $seq->start,
			stop   => $seq->end,
			score  => undef,
			strand => '.',
			phase  => '.',
			gclass => $self->refclass,
			gname  => $acc,
			tstart => undef,
			tstop  => undef,
			attributes  => [[Note => $seq->desc],@aliases],
		       }
		      );
  # now load each feature in turn
  for my $feat ($seq->all_SeqFeatures) {
    my $attributes = $self->get_attributes($feat);
    my $first = (shift @$attributes);

    my $location = $feat->location;
    my @segments = map {[$_->start,$_->end]} 
      $location->can('sub_Location') ? $location->sub_Location : $location;

    for my $segment (@segments) {

      $self->load_gff_line( {
			     ref    => $acc,
			     class  => $refclass,
			     source => 'EMBL',
			     method => $feat->primary_tag,
			     start  => $segment->[0],
			     stop   => $segment->[1],
			     score  => $feat->score || undef,
			     strand => $feat->strand > 0 ? '+' : ($feat->strand < 0 ? '-' : '.'),
			     phase  => $feat->frame || '.',
			     gclass => $first->[0],
			     gname  => $first->[1],
			     tstart => undef,
			     tstop  => undef,
			     attributes  => $attributes,
			    }
			  );
    }
  }

  # finish loading
  $self->finish_load();

  # now load the DNA
  $self->load_sequence_string($acc,$seq->seq);

  1;
}

sub get_attributes {
  my $self = shift;
  my $seq  = shift;

  my @tags = $seq->all_tags or return;
  my @result;
  foreach my $tag (@tags) {
    foreach my $value ($seq->each_tag_value($tag)) {
      push @result,[$tag=>$value];
    }
  }
  \@result;
}



1;
