package Bio::DB::GFF::Adaptor::biofetch;
#$Id$
=head1 NAME

Bio::DB::GFF::Adaptor::biofetch -- Cache BioFetch objects in a Bio::DB::GFF database

=head1 SYNOPSIS

Proof of principle.  Not for production use.

=head1 DESCRIPTION

This adaptor is a proof-of-principle.  It is used to fetch BioFetch
sequences into a Bio::DB::GFF database (currently uses a hard-coded
EMBL database) as needed.  This allows the Generic Genome Browser to
be used as a Genbank/EMBL browser.

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright 2002 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

use strict;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use Bio::DB::BioFetch;
use Bio::SeqIO;

use vars qw(%preferred_tags);

# THIS IS WRONG: biofetch should delegate to an underlying
# database adaptor, and not inherit from one.
use base qw(Bio::DB::GFF::Adaptor::dbi::mysql);

# priority for choosing names of CDS tags, higher is higher priority
%preferred_tags = (
		      strain        => 10,
		      organism      => 20,
		      protein_id    => 40,
		      locus_tag     => 50,
		      locus         => 60,
		      gene          => 70,
		      standard_name => 80,
		      );

=head2 new

 Title   : new
 Usage   : $db = Bio::DB::GFF->new(-adaptor=>'biofetch',@args)
 Function: create a new adaptor
 Returns : a Bio::DB::GFF object
 Args    :   -adaptor : required.  Which adaptor to use; biofetch for mysql, biofetch_oracle for Oracle
             -preferred_tags : optional.  A hash of {classname => weight,...}
                               used to determine the class and name of the feature
                               when a choice of possible feature classes is available
                               (e.g. a feature has both a 'gene' and a 'locus' tag).
                               Common defaults are provided that work well for eukaryotic
                               features (but not well for viral/prokaryotic)
              see below for additional arguments.
 Status  : Public

This is the constructor for the adaptor.  It is called automatically
by Bio::DB::GFF-E<gt>new.  In addition to arguments that are common among
all adaptors, the following class-specific arguments are recgonized:

  Argument       Description
  --------       -----------

  -dsn           the DBI data source, e.g. 'dbi:mysql:ens0040'

  -user          username for authentication

  -pass          the password for authentication

  -proxy         [['http','ftp'],'http://proxy:8080']

  -source        source to use for loaded features ('EMBL')

-dsn,-user and -pass indicate the local database to cache results in,
and as are per Bio::DB::GFF::Adaptor::dbi.  The -proxy argument allows
you to set the biofetch web proxy, and uses the same syntax described
for the proxy() method of L<Bio::DB::WebDBSeqI>, except that the
argument must be passed as an array reference.

=cut

sub new {
  my $class = shift;
  my $self  = $class->SUPER::new(@_);
  my ($preferred,$proxy,$source) = rearrange(['PREFERRED_TAGS','PROXY','SOURCE'],@_);

  # if the caller sent their own preferences, then use these, otherwise use defaults.
  $self->_preferred_tags($preferred ? $preferred : \%preferred_tags);
  $self->_source($source || 'EMBL');

  if ($proxy) {
    my @args = ref($proxy) ? @$proxy : eval $proxy;
    $self->{_proxy} = \@args if @args;
  }
  $self;
}

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
  'Sequence';
}

sub load_from_embl {
  my $self = shift;
  my $db   = shift;
  my $acc  = shift or $self->throw('Must provide an accession ID');

  my $biofetch;
  if ($self->{_biofetch}{$db}) {
    $biofetch = $self->{_biofetch}{$db};
  } else {
    $biofetch = $self->{_biofetch}{$db} = Bio::DB::BioFetch->new(-db=>$db);
    $biofetch->retrieval_type('tempfile');
    $biofetch->proxy(@{$self->{_proxy}}) if $self->{_proxy};
  }

  my $seq  = eval {$biofetch->get_Seq_by_id($acc)} or return;
  $self->_load_embl($acc,$seq);
  1;
}

sub load_from_file {
  my $self = shift;
  my $file = shift;

  my $format = $file =~ /\.(gb|genbank|gbk)$/i ? 'genbank' : 'embl';

  my $seqio = Bio::SeqIO->new( '-format' => $format, -file => $file);
  my $seq   = $seqio->next_seq;

  $self->_load_embl($seq->accession,$seq);
  1;
}

sub _load_embl {
  my $self = shift;
  my $acc  = shift;
  my $seq  = shift;
  my $refclass = $self->refclass;
  my $locus    = $seq->id;
  my $source   = $self->_source;

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
			source => $source,
#			method => 'origin',
			method => 'region',
			start  => 1,
			stop   => $seq->length,
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
  my ($transcript_version,$mRNA_version) = (0,0);
  for my $feat ($seq->all_SeqFeatures) {
    my $attributes = $self->get_attributes($feat);
    my $name       = $self->guess_name($attributes);

    my $location = $feat->location;
    my @segments = map {[$_->start,$_->end,$_->seq_id]}
      $location->can('sub_Location') ? $location->sub_Location : $location;

# this changed CDS to coding, but that is the wrong thing to do, since
# CDS is in SOFA and coding is not
#    my $type     =   $feat->primary_tag eq 'CDS'   ? 'coding'
#                   : $feat->primary_tag;
    my $type=  $feat->primary_tag;
    next if (lc($type) eq 'contig');
#    next if (lc($type) eq 'variation');

    if (lc($type) eq 'variation' and $feat->length == 1) {
      $type = 'SNP';
    } elsif (lc($type) eq 'variation' ) {
      $type = 'chromosome_variation';
    }

    if ($type  eq 'source') {
      $type = 'region';
    }

    if ($type =~ /misc.*RNA/i) {
      $type = 'RNA';
    }

    if ($type eq 'misc_feature' and $name->[1] =~ /similar/i) {
      $type = 'computed_feature_by_similarity';
    } elsif ($type eq 'misc_feature') {
      warn "skipping a misc_feature\n";
      next;
    }

    my $parttype =  $feat->primary_tag eq 'mRNA'   ? 'exon' : $feat->primary_tag;

    if ($type eq 'gene') {
      $transcript_version = 0;
      $mRNA_version       = 0;
    } elsif ($type eq 'mRNA') {
      $name->[1] = sprintf("%s.t%02d",$name->[1],++$transcript_version);
    } elsif ($type eq 'CDS') {
      $name->[0] = 'mRNA';
      $name->[1] = sprintf("%s.t%02d",$name->[1],$transcript_version);
    }

    my $strand = $feat->strand;
    my $str    = defined $strand ?
                                     ($strand > 0 ? '+' : '-')
				   : '.';
    $self->load_gff_line( {
			   ref    => $acc,
			   class  => $refclass,
			   source => $source,
			   method => $type,
			   start  => $location->start,
			   stop   => $location->end,
			   score  => $feat->score || undef,
			   strand => $str,
			   phase  => $feat->frame || '.',
			   gclass => $name->[0],
			   gname  => $name->[1],
			   tstart => undef,
			   tstop  => undef,
			   attributes  => $attributes,
			  }
			) if ($type &&
                           ($type ne 'CDS'||($type eq 'CDS'&&@segments==1) ) );

    @$attributes = ();

    next if @segments == 1;
    for my $segment (@segments) {

      my $strand = $feat->strand;
      my $str    = defined $strand ?
                                     ($strand > 0 ? '+' : '-')
				   : '.';
      $self->load_gff_line( {
			     ref    => $segment->[2] eq $locus ? $acc : $segment->[2],
			     class  => $refclass,
			     source => $source,
			     method => $parttype,
			     start  => $segment->[0],
			     stop   => $segment->[1],
			     score  => $feat->score || undef,
			     strand => $str,
			     phase  => $feat->frame || '.',
			     gclass => $name->[0],
			     gname  => $name->[1],
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

sub guess_name {
  my $self = shift;
  my $attributes = shift;
# remove this fix when Lincoln fixes it properly
  return ["Misc" => "Misc"] unless ($attributes);  # these are arbitrary, and possibly destructive defaults
  my @ordered_attributes = sort {($self->_preferred_tags->{$a->[0]} || 0) <=> ($self->_preferred_tags->{$b->[0]} || 0)} @$attributes;
  my $best = pop @ordered_attributes;
  @$attributes = @ordered_attributes;
  return $best;
}


sub _preferred_tags {
  my $self = shift;
  $self->{preferred_tags} = shift if @_;
  return $self->{preferred_tags};
}

sub _source {
  my $self = shift;
  $self->{source} = shift if @_;
  $self->{source};
}

1;
