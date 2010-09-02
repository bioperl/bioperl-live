#
# bioperl module for Bio::LiveSeq::Gene
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Joseph Insana <insana@ebi.ac.uk> <jinsana@gmx.net>
#
# Copyright Joseph Insana
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::LiveSeq::Gene - Range abstract class for LiveSeq

=head1 SYNOPSIS

  # documentation needed

=head1 DESCRIPTION

This is used as storage for all object references concerning a particular gene.

=head1 AUTHOR - Joseph A.L. Insana

Email:  Insana@ebi.ac.uk, jinsana@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LiveSeq::Gene;
use strict;
use Carp;
use Bio::LiveSeq::Prim_Transcript; # needed to create maxtranscript obj

=head2 new

  Title   : new
  Usage   : $gene = Bio::LiveSeq::Gene->new(-name => "name",
                                            -features => $hashref
                                            -upbound => $min
                                            -downbound => $max);

  Function: generates a new Bio::LiveSeq::Gene
  Returns : reference to a new object of class Gene
  Errorcode -1
  Args    : one string and one hashreference containing all features defined
            for the Gene and the references to the LiveSeq objects for those
            features.
            Two labels for defining boundaries of the gene. Usually the
            boundaries will reflect max span of transcript, exon... features,
            while the DNA sequence will be created with some flanking regions
            (e.g. with the EMBL_SRS::gene2liveseq routine).
            If these two labels are not given, they will default to the start
            and end of the DNA object.
  Note    : the format of the hash has to be like
               DNA => reference to LiveSeq::DNA object
               Transcripts => reference to array of transcripts objrefs
               Transclations => reference to array of transcripts objrefs
               Exons => ....
               Introns => ....
               Prim_Transcripts => ....
               Repeat_Units => ....
               Repeat_Regions => ....
            Only DNA and Transcripts are mandatory

=cut

sub new {
  my ($thing, %args) = @_;
  my $class = ref($thing) || $thing;
  my ($i,$self,%gene);

  my ($name,$inputfeatures,$upbound,$downbound)=($args{-name},$args{-features},$args{-upbound},$args{-downbound});

  unless (ref($inputfeatures) eq "HASH") {
    carp "$class not initialised because features hash not given";
    return (-1);
  }

  my %features=%{$inputfeatures}; # this is done to make our own hash&ref, not
  my $features=\%features; # the ones input'ed, that could get destroyed
  
  my $DNA=$features->{'DNA'};
  unless (ref($DNA) eq "Bio::LiveSeq::DNA") {
    carp "$class not initialised because DNA feature not found";
    return (-1);
  }

  my ($minstart,$maxend);# used to calculate Gene->maxtranscript from Exon, Transcript (CDS) and Prim_Transcript features

  my ($start,$end);

  my @Transcripts=@{$features->{'Transcripts'}};

  my $strand;
  unless (ref($Transcripts[0]) eq "Bio::LiveSeq::Transcript") {
    $self->warn("$class not initialised: first Transcript not a LiveSeq object");
    return (-1);
  } else {
    $strand=$Transcripts[0]->strand; # for maxtranscript consistency check
  }

  for $i (@Transcripts) {
    ($start,$end)=($i->start,$i->end);
    unless ((ref($i) eq "Bio::LiveSeq::Transcript")&&($DNA->valid($start))&&($DNA->valid($end))) {
      $self->warn("$class not initialised because of problems in Transcripts feature");
      return (-1);
    } else {
    }
    unless($minstart) { $minstart=$start; } # initialize
    unless($maxend) { $maxend=$end; } # initialize
    if ($i->strand != $strand) {
      $self->warn("$class not initialised because exon-CDS-prim_transcript features do not share the same strand!");
      return (-1);
    }
    if (($strand == 1)&&($start < $minstart)||($strand == -1)&&($start > $minstart)) { $minstart=$start; }
    if (($strand == 1)&&($end > $maxend)||($strand == -1)&&($end < $maxend)) { $maxend=$end; }
  }  
  my @Translations; my @Introns; my @Repeat_Units; my @Repeat_Regions;
  my @Prim_Transcripts; my @Exons;
  if (defined($features->{'Translations'})) {
    @Translations=@{$features->{'Translations'}}; }
  if (defined($features->{'Exons'})) {
    @Exons=@{$features->{'Exons'}}; }
  if (defined($features->{'Introns'})) {
    @Introns=@{$features->{'Introns'}}; }
  if (defined($features->{'Repeat_Units'})) {
    @Repeat_Units=@{$features->{'Repeat_Units'}}; }
  if (defined($features->{'Repeat_Regions'})) {
    @Repeat_Regions=@{$features->{'Repeat_Regions'}}; }
  if (defined($features->{'Prim_Transcripts'})) {
    @Prim_Transcripts=@{$features->{'Prim_Transcripts'}}; }

  
  if (@Translations) {
    for $i (@Translations) {
      ($start,$end)=($i->start,$i->end);
      unless ((ref($i) eq "Bio::LiveSeq::Translation")&&($DNA->valid($start))&&($DNA->valid($end))) {
	$self->warn("$class not initialised because of problems in Translations feature");
	return (-1);
      }
    }
  }
  if (@Exons) {
    for $i (@Exons) {
      ($start,$end)=($i->start,$i->end);
      unless ((ref($i) eq "Bio::LiveSeq::Exon")&&($DNA->valid($start))&&($DNA->valid($end))) {
	$self->warn("$class not initialised because of problems in Exons feature");
	return (-1);
      }
      if ($i->strand != $strand) {
	$self->warn("$class not initialised because exon-CDS-prim_transcript features do not share the same strand!");
	return (-1);
      }
      if (($strand == 1)&&($start < $minstart)||($strand == -1)&&($start > $minstart)) { $minstart=$start; }
      if (($strand == 1)&&($end > $maxend)||($strand == -1)&&($end < $maxend)) { $maxend=$end; }
    }
  }
  if (@Introns) {
    for $i (@Introns) {
      ($start,$end)=($i->start,$i->end);
      unless ((ref($i) eq "Bio::LiveSeq::Intron")&&($DNA->valid($start))&&($DNA->valid($end))) {
	$self->warn("$class not initialised because of problems in Introns feature");
	return (-1);
      }
    }
  }
  if (@Repeat_Units) {
    for $i (@Repeat_Units) {
      ($start,$end)=($i->start,$i->end);
      unless ((ref($i) eq "Bio::LiveSeq::Repeat_Unit")&&($DNA->valid($start))&&($DNA->valid($end))) {
	$self->warn("$class not initialised because of problems in Repeat_Units feature");
	return (-1);
      }
    }
  }
  if (@Repeat_Regions) {
    for $i (@Repeat_Regions) {
      ($start,$end)=($i->start,$i->end);
      unless ((ref($i) eq "Bio::LiveSeq::Repeat_Region")&&($DNA->valid($start))&&($DNA->valid($end))) {
	$self->warn("$class not initialised because of problems in Repeat_Regions feature");
	return (-1);
      }
    }
  }
  if (@Prim_Transcripts) {
    for $i (@Prim_Transcripts) {
      ($start,$end)=($i->start,$i->end);
      unless ((ref($i) eq "Bio::LiveSeq::Prim_Transcript")&&($DNA->valid($start))&&($DNA->valid($end))) {
	$self->warn("$class not initialised because of problems in Prim_Transcripts feature");
	return (-1);
      }
      if ($i->strand != $strand) {
	$self->warn("$class not initialised because exon-CDS-prim_transcript features do not share the same strand!");
	return (-1);
      }
      if (($strand == 1)&&($start < $minstart)||($strand == -1)&&($start > $minstart)) { $minstart=$start; }
      if (($strand == 1)&&($end > $maxend)||($strand == -1)&&($end < $maxend)) { $maxend=$end; }
    }
  }

  # create an array containing all obj references for all Gene Features
  # useful for _set_Gene_in_all
  my @allfeatures;
  push (@allfeatures,$DNA,@Transcripts,@Translations,@Exons,@Introns,@Repeat_Units,@Repeat_Regions,@Prim_Transcripts);

  # create hash holding numbers for Gene Features
  my %multiplicity; 
  my $key; my @array;
  foreach $key (keys(%features)) {
    unless ($key eq "DNA") {
      @array=@{$features{$key}};
      $multiplicity{$key}=scalar(@array);
    }
  }
  $multiplicity{DNA}=1;

  # create maxtranscript object. It's a Prim_Transcript with start as the
  # minimum start and end as the maximum end.
  # usually these start and end will be the same as the gene->upbound and
  # gene->downbound, but maybe there could be cases when this will be false
  # (e.g. with repeat_units just before the prim_transcript or first exon,
  # but still labelled with the same /gene qualifier)

  my $maxtranscript=Bio::LiveSeq::Prim_Transcript->new(-start => $minstart, -end => $maxend, -strand => $strand, -seq => $DNA);


  # check the upbound downbound parameters
  if (defined($upbound)) {
    unless ($DNA->valid($upbound)) {
      $self->warn("$class not initialised because upbound label not valid");
      return (-1);
    }
  } else {
    $upbound=$DNA->start;
  }
  if (defined($downbound)) {
    unless ($DNA->valid($downbound)) {
      $self->warn("$class not initialised because downbound label not valid");
      return (-1);
    }
  } else {
    $downbound=$DNA->end;
  }

  %gene = (name => $name, features => $features,multiplicity => \%multiplicity,
          upbound => $upbound, downbound => $downbound, allfeatures => \@allfeatures, maxtranscript => $maxtranscript);
  $self = \%gene;
  $self = bless $self, $class;
  _set_Gene_in_all($self,@allfeatures);
  return $self;
}

# this sets the "gene" objref in all the objects "belonging" to the Gene,
# i.e. in all its Features.
sub _set_Gene_in_all {
  my $Gene=shift;
  my $self;
  foreach $self (@_) {
    $self->gene($Gene);
  }
}

# you can get or set the name of the gene
sub name {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'name'} = $value;
  }
  unless (exists $self->{'name'}) {
    return "unknown";
  } else {
    return $self->{'name'};
  }
}

# gets the features hash
sub features {
  my $self=shift;
  return ($self->{'features'});
}
sub get_DNA {
  my $self=shift;
  return ($self->{'features'}->{'DNA'});
}
sub get_Transcripts {
  my $self=shift;
  return ($self->{'features'}->{'Transcripts'});
}
sub get_Translations {
  my $self=shift;
  return ($self->{'features'}->{'Translations'});
}
sub get_Prim_Transcripts {
  my $self=shift;
  return ($self->{'features'}->{'Prim_Transcripts'});
}
sub get_Repeat_Units {
  my $self=shift;
  return ($self->{'features'}->{'Repeat_Units'});
}
sub get_Repeat_Regions {
  my $self=shift;
  return ($self->{'features'}->{'Repeat_Regions'});
}
sub get_Introns {
  my $self=shift;
  return ($self->{'features'}->{'Introns'});
}
sub get_Exons {
  my $self=shift;
  return ($self->{'features'}->{'Exons'});
}
sub featuresnum {
  my $self=shift;
  return ($self->{'multiplicity'});
}
sub upbound {
  my $self=shift;
  return ($self->{'upbound'});
}
sub downbound {
  my $self=shift;
  return ($self->{'downbound'});
}
sub printfeaturesnum {
  my $self=shift;
  my ($key,$value);
  my %hash=%{$self->featuresnum};
  foreach $key (keys(%hash)) {
    $value=$hash{$key};
    print "\t$key => $value\n";
  }
}
sub maxtranscript {
  my $self=shift;
  return ($self->{'maxtranscript'});
}

sub delete_Obj {
  my $self = shift;
  my @values= values %{$self};
  my @keys= keys %{$self};

  foreach my $key ( @keys ) {
    delete $self->{$key};
  }
  foreach my $value ( @values ) {
    if (index(ref($value),"LiveSeq") != -1) { # object case
      eval {
	# delete $self->{$value};
	$value->delete_Obj;
      };
    } elsif (index(ref($value),"ARRAY") != -1) { # array case
      my @array=@{$value};
      my $element;
      foreach $element (@array) {
	eval {
	  $element->delete_Obj;
	};
      }
    } elsif (index(ref($value),"HASH") != -1) { # object case
      my %hash=%{$value};
      my $element;
      foreach $element (%hash) {
	eval {
	  $element->delete_Obj;
	};
      }
    }
  }
  return(1);
}


=head2 verbose

 Title   : verbose
 Usage   : $self->verbose(0)
 Function: Sets verbose level for how ->warn behaves
           -1 = silent: no warning
            0 = reduced: minimal warnings
            1 = default: all warnings
            2 = extended: all warnings + stack trace dump
            3 = paranoid: a warning becomes a throw and the program dies

           Note: a quick way to set all LiveSeq objects at the same verbosity
           level is to change the DNA level object, since they all look to
           that one if their verbosity_level attribute is not set.
           But the method offers fine tuning possibility by changing the
           verbose level of each object in a different way.

           So for example, after $loader= and $gene= have been retrieved
           by a program, the command $gene->verbose(0); would
           set the default verbosity level to 0 for all objects.

 Returns : the current verbosity level
 Args    : -1,0,1,2 or 3

=cut


sub verbose {
  my $self=shift;
  my $value = shift;
  return $self->{'features'}->{'DNA'}->verbose($value);
}

sub warn {
  my $self=shift;
  my $value = shift;
  return $self->{'features'}->{'DNA'}->warn($value);
}



1;
