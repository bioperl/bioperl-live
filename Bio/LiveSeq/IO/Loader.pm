#
# bioperl module for Bio::LiveSeq::IO::Loader
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

Bio::LiveSeq::IO::Loader - Parent Loader for LiveSeq

=head1 SYNOPSIS

  #documentation needed

=head1 DESCRIPTION

This package holds common methods used by BioPerl and file loaders.
It contains methods to create LiveSeq objects out of entire entries or from a
localized sequence region surrounding a particular gene.

=head1 AUTHOR - Joseph A.L. Insana

Email:  Insana@ebi.ac.uk, jinsana@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LiveSeq::IO::Loader;

use strict;
use Carp qw(cluck croak carp);
use Bio::LiveSeq::DNA;
use Bio::LiveSeq::Exon;
use Bio::LiveSeq::Transcript ;
use Bio::LiveSeq::Translation;
use Bio::LiveSeq::Gene;
use Bio::LiveSeq::Intron;
use Bio::LiveSeq::Prim_Transcript;
use Bio::LiveSeq::Repeat_Region;
use Bio::LiveSeq::Repeat_Unit;
use Bio::LiveSeq::AARange;
use Bio::Tools::CodonTable;

=head2 entry2liveseq

  Title   : entry2liveseq
  Usage   : @translationobjects=$loader->entry2liveseq();
          : @translationobjects=$loader->entry2liveseq(-getswissprotinfo => 0);
  Function: creates LiveSeq objects from an entry previously loaded
  Returns : array of references to objects of class Translation
  Errorcode 0
  Args    : optional boolean flag to avoid the retrieval of SwissProt
            information for all Transcripts containing SwissProt x-reference
            default is 1 (to retrieve those information and create AARange
            LiveSeq objects)
  Note    : this method can get really slow for big entries. The lightweight
            gene2liveseq method is recommended

=cut

sub entry2liveseq {
  my ($self, %args) = @_;
  my ($getswissprotinfo)=($args{-getswissprotinfo});
  if (defined($getswissprotinfo)) {
    if (($getswissprotinfo ne 0)&&($getswissprotinfo ne 1)) {
      carp "-getswissprotinfo argument can take only boolean (1 or 0) values. Setting it to 0, i.e. not trying to retrieve SwissProt information....";
      $getswissprotinfo=0;
    }
  } else {
    $getswissprotinfo=1;
  }
  my $hashref=$self->{'hash'};
  unless ($hashref) { return (0); }
  my @translationobjects=$self->hash2liveseq($hashref,$getswissprotinfo);
  my $test_transl=0;
  if ($test_transl) { $self->test_transl($hashref,\@translationobjects);}
  return @translationobjects;
}

=head2 novelaasequence2gene

  Title   : novelaasequence2gene
  Usage   : $gene=$loader->novelaasequence2gene(-aasequence => "MGLAAPTRS*");
          : $gene=$loader->novelaasequence2gene(-aasequence => "MGLAAPTRS*");
                                             -taxon => 9606,
                                             -gene_name => "tyr-kinase");

  Function: creates LiveSeq objects from a novel amino acid sequence,
            using codon usage database to choose codons according to
            relative frequencies.
            If a taxon ID is not specified, the default is to use the human
            one (taxonomy ID 9606).
  Returns : reference to a Gene object containing references to LiveSeq objects
  Errorcode 0
  Args    : string containing an amino acid sequence
            integer (optional) with a taxonomy ID
            string specifying a gene name

=cut

=head2 gene2liveseq

  Title   : gene2liveseq
  Usage   : $gene=$loader->gene2liveseq(-gene_name => "gene name");
          : $gene=$loader->gene2liveseq(-gene_name => "gene name",
                                        -flanking => 64);
          : $gene=$loader->gene2liveseq(-gene_name => "gene name",
                                        -getswissprotinfo => 0);
          : $gene=$loader->gene2liveseq(-position => 4);

  Function: creates LiveSeq objects from an entry previously loaded
            It is a "light weight" creation because it creates
            a LiveSequence just for the interesting region in an entry
            (instead than for the total entry, like in entry2liveseq) and for
            the flanking regions up to 500 nucleotides (default) or up to
            the specified amount of nucleotides (given as argument) around the
            Gene.
  Returns : reference to a Gene object containing possibly alternative
            Transcripts.
  Errorcode 0
  Args    : string containing the gene name as in the EMBL feature qualifier
            integer (optional) "flanking": amount of flanking bases to be kept
            boolean (optional) "getswissprotinfo": if set to "0" it will avoid
             trying to fetch information from a crossreference to a SwissProt
             entry, avoding the process of creation of AARange objects
             It is "1" (on) by default

            Alternative to a gene_name, a position can be given: an
            integer (1-) containing the position of the desired CDS in the
            loaded entry

=cut

sub gene2liveseq {
  my ($self, %args) = @_;
  my ($gene_name,$flanking,$getswissprotinfo,$cds_position)=($args{-gene_name},$args{-flanking},$args{-getswissprotinfo},$args{-position});
  my $input;
  unless (($gene_name)||($cds_position)) {
    carp "Gene_Name or Position not specified for gene2liveseq loading function";
    return (0);
  }
  if (($gene_name)&&($cds_position)) {
    carp "Gene_Name and Position cannot be given together";
    return (0);
  } elsif ($gene_name) {
    $input=$gene_name;
  } else {
    $input="cds-position:".$cds_position;
  }

  if (defined($getswissprotinfo)) {
    if (($getswissprotinfo ne 0)&&($getswissprotinfo ne 1)) {
      carp "-getswissprotinfo argument can take only boolean (1 or 0) values. Setting it to 0, i.e. not trying to retrieve SwissProt information....";
      $getswissprotinfo=0;
    }
  } else {
    $getswissprotinfo=1;
  }

  if (defined($flanking)) {
    unless ($flanking >= 0) {
      carp "No sense in specifying a number < 0 for flanking regions to be created for gene2liveseq loading function";
      return (0);
    }
  } else {
    $flanking=500; # the default flanking length
  }
  my $hashref=$self->{'hash'};
  unless ($hashref) { return (0); }
  my $gene=$self->hash2gene($hashref,$input,$flanking,$getswissprotinfo);
  unless ($gene) { # if $gene == 0 it means problems in hash2gene
    carp "gene2liveseq produced error";
    return (0);
  }
  return $gene;
}

# TODO: update so that it will work even if CDS is not only accepted FEATURE!!
# this method is for now deprecated and not supported
sub test_transl {
  my ($self,$entry)=@_;
  my @features=@{$entry->{'Features'}};
  my @translationobjects=@{$_[1]};
  my ($i,$translation);
  my ($obj_transl,$hash_transl);
  my @cds=@{$entry->{'CDS'}};
  foreach $translation (@translationobjects) {
    $obj_transl=$translation->seq;
    $hash_transl=$cds[$i]->{'qualifiers'}->{'translation'};
    #before seq was changed in Translation 1.4# chop $obj_transl; # to remove trailing "*"
    unless ($obj_transl eq $hash_transl) {
      cluck "Serious error: Translation from the Entry does not match Translation from object's seq for CDS at position $i";
      carp "\nEntry's transl: ",$hash_transl,"\n";
      carp "\nObject's transl: ",$obj_transl,"\n";
      exit;
    }
    $i++;
  }
}

# argument: hashref containing the EMBL entry datas,
#           getswissprotinfo boolean flag
# creates the liveseq objects
# returns: an array of Translation object references
sub hash2liveseq {
  my ($self,$entry,$getswissprotinfo)=@_;
  my $i;
  my @transcripts;
  my $dna=Bio::LiveSeq::DNA->new(-seq => $entry->{'Sequence'} );
  $dna->alphabet(lc($entry->{'Molecule'}));
  $dna->display_id($entry->{'ID'});
  $dna->accession_number($entry->{'AccNumber'});
  $dna->desc($entry->{'Description'});
  my @cds=@{$entry->{'CDS'}};
  my ($swissacc,$swisshash); my @swisshashes;
  for $i (0..$#cds) {
    #my @transcript=@{$cds[$i]->{'range'}};
    #$transcript=\@transcript;
    #push (@transcripts,$transcript);
    push (@transcripts,$cds[$i]->{'range'});
    if ($getswissprotinfo) {
      $swissacc=$cds[$i]->{'qualifiers'}->{'db_xref'};
      $swisshash=$self->get_swisshash($swissacc);
      #$self->printswissprot($swisshash); # DEBUG
      push (@swisshashes,$swisshash);
    }
  }
  my @translations=($self->transexonscreation($dna,\@transcripts));
  my $translation; my $j=0;
  foreach $translation (@translations) {
    if ($swisshashes[$j]) { # if not 0
      $self->swisshash2liveseq($swisshashes[$j],$translation);
    }
    $j++;
  }
  return (@translations);
}

# only features pertaining to a specified gene are created
# only the sequence of the gene and appropriate context flanking regions
# are created as chain
# arguments: hashref, gene_name (OR: cds_position), length_of_flanking_sequences, getswissprotinfo boolean flag
# returns: reference to Gene object
#
# Note: if entry contains just one CDS, all the features get added
#       this is useful because often the features in these entries do not
#       carry the /gene qualifier
#
# errorcode: 0
sub hash2gene {
  my ($self,$entry,$input,$flanking,$getswissprotinfo)=@_;
  my $entryfeature;
  my $genefeatureshash;

  my @cds=@{$entry->{'CDS'}};

  # checking if a position has been given instead than a gene_name
  if (index($input,"cds-position:") == 0 ) {
    my $cds_position=substr($input,13); # extracting the cds position
    if (($cds_position >= 1)&&($cds_position <= scalar(@cds))) {
      $genefeatureshash=$self->_findgenefeatures($entry,undef,$cds_position,$getswissprotinfo);
    }
  } else {
    $genefeatureshash=$self->_findgenefeatures($entry,$input,undef,$getswissprotinfo);
  }

  unless (($genefeatureshash)&&(scalar(@{$genefeatureshash->{'genefeatures'}}))) { # array empty, no gene features found
    my @genes=$self->genes($entry);
    my $cds_number=scalar(@cds);
    warn "Warning! Not even one genefeature found for /$input/....
    The genes present in this entry are:\n\t@genes\n
    The number of CDS in this entry is:\n\t$cds_number\n";
    return(0);
  }

  # get max and min, check flankings
  my ($min,$max)=$self->rangeofarray(@{$genefeatureshash->{'labels'}}); # gene "boundaries"
  my $seqlength=$entry->{'SeqLength'};
  my ($mindna,$maxdna); # some flanking region next to the gene "boundaries"
  if ($min-$flanking < 1) {
    $mindna=1;
  } else {
    $mindna=$min-$flanking;
  }
  if ($max+$flanking > $seqlength) {
    $maxdna=$seqlength;
  } else {
    $maxdna=$max+$flanking;
  }
  my $subseq=substr($entry->{'Sequence'},$mindna-1,$maxdna-$mindna+1);

  # create LiveSeq objects

  # create DNA
  my $dna=Bio::LiveSeq::DNA->new(-seq => $subseq, -offset => $mindna);
  $dna->alphabet(lc($entry->{'Molecule'}));
  $dna->source($entry->{'Organism'});
  $dna->display_id($entry->{'ID'});
  $dna->accession_number($entry->{'AccNumber'});
  $dna->desc($entry->{'Description'});

  my @transcripts=@{$genefeatureshash->{'transcripts'}};
  # create Translations, Transcripts, Exons out of the CDS
  unless (@transcripts) {
    cluck "no CDS feature found for /$input/....";
    return(0);
  }
  my @translationobjs=$self->transexonscreation($dna,\@transcripts);
  my @transcriptobjs;

  # get the Transcript obj_refs
  my $translation;
  my $j=0;
  my @ttables=@{$genefeatureshash->{'ttables'}};
  my @swisshashes=@{$genefeatureshash->{'swisshashes'}};
  foreach $translation (@translationobjs) {
    push(@transcriptobjs,$translation->get_Transcript);
    if ($ttables[$j]) { # if not undef
      $translation->get_Transcript->translation_table($ttables[$j]);
    #} else { # DEBUG
    #  print "\n\t\tno translation table information....\n";
    }
    if ($swisshashes[$j]) { # if not 0
      $self->swisshash2liveseq($swisshashes[$j],$translation);
    }
    $j++;
  }

  my %gene; # this is the hash to store created object references
  $gene{DNA}=$dna;
  $gene{Transcripts}=\@transcriptobjs;
  $gene{Translations}=\@translationobjs;

  my @exonobjs; my @intronobjs;
  my @repeatunitobjs; my @repeatregionobjs;
  my @primtranscriptobjs;

  my ($object,$range,$start,$end,$strand);

  my @exons=@{$genefeatureshash->{'exons'}};
  my @exondescs=@{$genefeatureshash->{'exondescs'}};
  if (@exons) {
    my $exoncount = 0;
    foreach $range (@exons) {
      ($start,$end,$strand)=@{$range};
      $object = Bio::LiveSeq::Exon->new(-seq=>$dna,-start=>$start,-end=>$end,-strand=>$strand);
      if ($object != -1) {
	$object->desc($exondescs[$exoncount]) if defined $exondescs[$exoncount];
	$exoncount++;
	push (@exonobjs,$object);
      } else {
	$exoncount++;
      }
    }
    $gene{Exons}=\@exonobjs;
  }
  my @introns=@{$genefeatureshash->{'introns'}};
  my @introndescs=@{$genefeatureshash->{'introndescs'}};
  if (@introns) {
    my $introncount = 0;
    foreach $range (@introns) {
      ($start,$end,$strand)=@{$range};
      $object=Bio::LiveSeq::Intron->new(-seq=>$dna,-start=>$start,-end=>$end,-strand=>$strand);
      if ($object != -1) {
	$object->desc($introndescs[$introncount]);
	$introncount++;
	push (@intronobjs,$object);
      } else {
	$introncount++;
      }
    }
    $gene{Introns}=\@intronobjs;
  }
  my @prim_transcripts=@{$genefeatureshash->{'prim_transcripts'}};
  if (@prim_transcripts) {
    foreach $range (@prim_transcripts) {
      ($start,$end,$strand)=@{$range};
      $object=Bio::LiveSeq::Prim_Transcript->new(-seq=>$dna,-start=>$start,-end=>$end,-strand=>$strand);
      if ($object != -1) { push (@primtranscriptobjs,$object); }
    }
    $gene{Prim_Transcripts}=\@primtranscriptobjs;
  }
  my @repeat_regions=@{$genefeatureshash->{'repeat_regions'}};
  my @repeat_regions_family=@{$genefeatureshash->{'repeat_regions_family'}};
  if (@repeat_regions) {
    my $k=0;
    foreach $range (@repeat_regions) {
      ($start,$end,$strand)=@{$range};
      $object=Bio::LiveSeq::Repeat_Region->new(-seq=>$dna,-start=>$start,-end=>$end,-strand=>$strand);
      if ($object != -1) {
	$object->desc($repeat_regions_family[$k]);
	$k++;
	push (@repeatregionobjs,$object);
      } else {
	$k++;
      }
    }
    $gene{Repeat_Regions}=\@repeatregionobjs;
  }
  my @repeat_units=@{$genefeatureshash->{'repeat_units'}};
  my @repeat_units_family=@{$genefeatureshash->{'repeat_units_family'}};
  if (@repeat_units) {
    my $k=0;
    foreach $range (@repeat_units) {
      ($start,$end,$strand)=@{$range};
      $object=Bio::LiveSeq::Repeat_Unit->new(-seq=>$dna,-start=>$start,-end=>$end,-strand=>$strand);
      if ($object != -1) {
	$object->desc($repeat_units_family[$k]);
	$k++;
	push (@repeatunitobjs,$object);
      } else {
	$k++;
      }
    }
    $gene{Repeat_Units}=\@repeatunitobjs;
  }

  # create the Gene
  my $gene_name=$genefeatureshash->{'gene_name'}; # either a name or a cdspos
  return (Bio::LiveSeq::Gene->new(-name=>$gene_name,-features=>\%gene,
                                  -upbound=>$min,-downbound=>$max));
}

# maybe this function will be moved to general utility package
# argument: array of numbers
# returns: (min,max) numbers in the array
sub rangeofarray {
  my $self=shift;
  my @array=@_;
  #print "\n-=-=-=-=-=-=-=-=-=-=array: @array\n";
  my ($max,$min,$element);
  $min=$max=shift(@array);
  foreach $element (@array) {
      $element = 0 unless defined $element;
    if ($element < $min) {
      $min=$element;
    }
    if ($element > $max) {
      $max=$element;
    }
  }
  #print "\n-=-=-=-=-=-=-=-=-=-=min: $min\tmax: $max\n";
  return ($min,$max);
}


# argument: reference to DNA object, reference to array of transcripts
# returns: an array of Translation object references
sub transexonscreation {
  my $self=shift;
  my $dna=$_[0];
  my @transcripts=@{$_[1]};

  my (@transexons,$start,$end,$strand,$exonref,$exonobj,$transcript,$transcriptobj);
  my $translationobj;
  my @translationobjects;
  foreach $transcript (@transcripts) {
    foreach $exonref (@{$transcript}) {
      ($start,$end,$strand)=@{$exonref};
      #print "Creating Exon: start $start end $end strand $strand\n";
      $exonobj=Bio::LiveSeq::Exon->new(-seq=>$dna,-start=>$start,-end=>$end,-strand=>$strand);
      #push (@exonobjects,$exonobj);
      push (@transexons,$exonobj);
    }
    $transcriptobj=Bio::LiveSeq::Transcript->new(-exons => \@transexons );
    if ($transcriptobj != -1) {
      $translationobj=Bio::LiveSeq::Translation->new(-transcript=>$transcriptobj);
      @transexons=(); # cleans it
      #push (@transcriptobjects,$transcriptobj);
      push (@translationobjects,$translationobj);
    }
  }
  return (@translationobjects);
}

#sub printgene {
# deleted. Some functionality placed in Gene->printfeaturesnum

=head2 printswissprot

  Title   : printswissprot
  Usage   : $loader->printswissprot($hashref);
  Function: prints out all information loaded from a database entry into the
            loader. Mainly used for testing purposes.
  Args    : a hashref containing the SWISSPROT entry datas
  Note    : the hashref can be obtained with a call to the method
               $loader->get_swisshash()      (BioPerl via Bio::DB::EMBL.pm)
	    that takes as argument a string like "SWISS-PROT:P10275"

=cut

# argument: hashref containing the SWISSPROT entry datas
# prints out that hash, showing the information loaded
sub printswissprot {
  my ($self,$entry)=@_;
  unless ($entry) {
    return;
  }
  printf "ID: %s\n",
      $entry->{'ID'};
  printf "ACC: %s\n",
      $entry->{'AccNumber'};
  printf "GENE: %s\n",
      $entry->{'Gene'};
  printf "DES: %s\n",
      $entry->{'Description'};
  printf "ORG: %s\n",
      $entry->{'Organism'};
  printf "SEQLN: %s\n",
      $entry->{'SeqLength'};
  printf "SEQ: %s\n",
      substr($entry->{'Sequence'},0,64);
  if ($entry->{'Features'}) {
    my @features=@{$entry->{'Features'}};
    my $i;
    for $i (0..$#features) {
      print "|",$features[$i]->{'name'},"|";
      print " at ",$features[$i]->{'location'},": ";
      print "",$features[$i]->{'desc'},"\n";
    }
  }
}

=head2 printembl

  Title   : printembl
  Usage   : $loader->printembl();
  Function: prints out all information loaded from a database entry into the
            loader. Mainly used for testing purposes.
  Args    : none

=cut

# argument: hashref containing the EMBL entry datas
# prints out that hash, showing the information loaded
sub printembl {
  my ($self,$entry)=@_;
  unless ($entry) {
    $entry=$self->{'hash'};
  }
  my ($i,$featurename);
  printf "ID: %s\n",
      $entry->{'ID'};
  printf "ACC: %s\n",
      $entry->{'AccNumber'};
  printf "MOL: %s\n",
      $entry->{'Molecule'};
  printf "DIV: %s\n",
      $entry->{'Division'};
  printf "DES: %s\n",
      $entry->{'Description'};
  printf "ORG: %s\n",
      $entry->{'Organism'};
  printf "SEQLN: %s\n",
      $entry->{'SeqLength'};
  printf "SEQ: %s\n",
      substr($entry->{'Sequence'},0,64);
  my @features=@{$entry->{'Features'}};
  my @cds=@{$entry->{'CDS'}};
  print "\nFEATURES\nNumber of CDS: ",scalar(@cds)," (out of ",$entry->{'FeaturesNumber'}, " total features)\n";
  my ($exonref,$transcript);
  my ($qualifiernumber,$qualifiers,$key);
  my ($start,$end,$strand);
  my $j=0;
  for $i (0..$#features) {
    $featurename=$features[$i]->{'name'};
    if ($featurename eq "CDS") {
      print "|CDS| number $j at feature position: $i\n";
      #print $features[$i]->{'location'},"\n";
      $transcript=$features[$i]->{'range'};
      foreach $exonref (@{$transcript}) {
	($start,$end,$strand)=@{$exonref};
	print "\tExon: start $start end $end strand $strand\n";
      }
      $j++;
    } else {
      print "|$featurename| at feature position: $i\n";
      print "\trange: ";
      print join("\t",@{$features[$i]->{'range'}}),"\n";
      #print $features[$i]->{'location'},"\n";
    }
    $qualifiernumber=$features[$i]->{'qual_number'};
    $qualifiers=$features[$i]->{'qualifiers'}; # hash
    foreach $key (keys (%{$qualifiers})) {
    print "\t\t",$key,": ";
      print $qualifiers->{$key},"\n";
    }
  }
}

=head2 genes

  Title   : genes
  Usage   : $loader->genes();
  Function: Returns an array of gene_names (strings) contained in the loaded
            entry.
  Args    : none

=cut

# argument: entryhashref
# returns: array of genenames found in the entry
sub genes {
  my ($self,$entry)=@_;
  unless ($entry) {
    $entry=$self->{'hash'};
  }
  my @entryfeatures=@{$entry->{'Features'}};
  my ($genename,$genenames,$entryfeature);
  for $entryfeature (@entryfeatures) {
    $genename=$entryfeature->{'qualifiers'}->{'gene'};
    if ($genename) {
      if (index($genenames,$genename) == -1) { # if name is new
	$genenames .= $genename . " "; # add the name
      }
    }
  }
  return (split(/ /,$genenames)); # assumes no space inbetween each genename
}

# arguments: swisshash, translation objref
# adds information to the Translation, creates AARange objects, sets the
# aa_range attribute on the Translation, pointing to those objects
sub swisshash2liveseq {
  my ($self,$entry,$translation)=@_;
  my $translength=$translation->length;
  $translation->desc($translation->desc . $entry->{'Description'});
  $translation->display_id("SWISSPROT:" . $entry->{'ID'});
  $translation->accession_number("SWISSPROT:" . $entry->{'AccNumber'});
  $translation->name($entry->{'Gene'});
  $translation->source($entry->{'Organism'});
  my @aarangeobjects;
  my ($start,$end,$aarangeobj,$feature);
  my @features; my @newfeatures;
  if ($entry->{'Features'}) {
    @features=@{$entry->{'Features'}};
  }
  my $cleavedmet=0;
  # check for cleaved Met
  foreach $feature (@features) {
    if (($feature->{'name'} eq "INIT_MET")&&($feature->{'location'} eq "0 0")) {
      $cleavedmet=1;
      $translation->{'offset'}="1"; # from swissprot to liveseq protein sequence
    } else {
      push(@newfeatures,$feature);
    }
  }

  my $swissseq=$entry->{'Sequence'};
  my $liveseqtransl=$translation->seq;
  chop $liveseqtransl; # to take away the trailing STOP "*"
  my $translated=substr($liveseqtransl,$cleavedmet);

  my ($liveseq_aa,$swiss_aa,$codes_aa)=$self->_get_alignment($translated,$swissseq); # alignment after cleavage of possible initial met

  if ((index($liveseq_aa,"-") != -1)||(index($swiss_aa,"-") != -1)) { # there are gaps, how to proceed?
    print "LIVE-SEQ=\'$liveseq_aa\'\nIDENTITY=\'$codes_aa\'\nSWS-PROT=\'$swiss_aa\'\n";
    carp "Nucleotides translation and SwissProt translation are different in size, cannot attach the SwissSequence to the EMBL one, cannot add any AminoAcidRange object/Domain information!";
    return;
  }

  #my $i=0; # debug
  @features=@newfeatures;
  foreach $feature (@features) {
    #print "Processing SwissProtFeature: $i\n"; # debug
    ($start,$end)=split(/ /,$feature->{'location'});
    # Note: cleavedmet is taken in account for updating numbering
    $aarangeobj=Bio::LiveSeq::AARange->new(-start => $start+$cleavedmet, -end => $end+$cleavedmet, -name => $feature->{'name'}, -description => $feature->{'description'}, -translation => $translation, -translength => $translength);
    if ($aarangeobj != -1) {
      push(@aarangeobjects,$aarangeobj);
    }
    # $i++; # debug
  }
  $translation->{'aa_ranges'}=\@aarangeobjects;
}

# if there is no SRS support, the default will be to return 0
# i.e. this function is overridden in SRS package
sub get_swisshash {
  return (0);
}

# Args: $entry hashref, gene_name OR cds_position (undef is used to
# choose between the two), getswissprotinfo boolean flag
# Returns: an hash holding various arrayref used in the hash2gene method
# Function: examines the nucleotide entry, identifying features belonging
# to the gene (defined either by its name or by the position of its CDS in
# the entry)

sub _findgenefeatures {
  my ($self,$entry,$gene_name,$cds_position,$getswissprotinfo)=@_;

  my @entryfeatures=@{$entry->{'Features'}};
  my @exons; my @introns; my @prim_transcripts; my @transcripts;
  my @repeat_units; my @repeat_regions;
  my @repeat_units_family; my @repeat_regions_family; my $rpt_family;
  my $entryfeature; my @genefeatures;
  my $desc; my @exondescs; my @introndescs;

  # for swissprot xreference
  my ($swissacc,$swisshash); my @swisshashes;

  # for translation_tables
  my @ttables;

  # to create labels
  my ($name,$exon);
  my @range; my @cdsexons; my @labels;

  # maybe here also could be added special case when there is no CDS feature
  # in the entry (e.g. tRNA entry -> TOCHECK).
  # let's deal with the special case in which there is just one gene per entry
  # usually without /gene qualifier
  my @cds=@{$entry->{'CDS'}};

  my $skipgenematch=0;
  if (scalar(@cds) == 1) {
    #carp "Note: only one CDS in this entry. Treating all features found in entry as Gene features.";
    $skipgenematch=1;
  }

  my ($cds_begin,$cds_end,$proximity);
  if ($cds_position) { # if a position has been requested
    my @cds_exons=@{$cds[$cds_position-1]->{'range'}};
    ($cds_begin,$cds_end)=($cds_exons[0]->[0],$cds_exons[-1]->[1]); # begin and end of CDS
    $gene_name=$cds[$cds_position-1]->{'qualifiers'}->{'gene'};
    # DEBUG
    unless ($skipgenematch) {
      carp "--DEBUG-- cdsbegin $cds_begin cdsend $cds_end--------";
    }
    $proximity=100; # proximity CONSTANT to decide whether a feature "belongs" to the CDS
  }

  for $entryfeature (@entryfeatures) { # get only features for the desired gene
    if (($skipgenematch)||(($cds_position)&&($self->_checkfeatureproximity($entryfeature->{'range'},$cds_begin,$cds_end,$proximity)))||(!($cds_position)&&($entryfeature->{'qualifiers'}->{'gene'} eq "$gene_name"))) {
      push(@genefeatures,$entryfeature);

      my @range=@{$entryfeature->{'range'}};
      $name=$entryfeature->{'name'};
      my %qualifierhash=%{$entryfeature->{'qualifiers'}};
      if ($name eq "CDS") { # that has range containing array of exons

	# swissprot crossindexing (if without SRS support it will fill array
	# with zeros and do nothing
	if ($getswissprotinfo) {
	  $swissacc=$entryfeature->{'qualifiers'}->{'db_xref'};
	  $swisshash=$self->get_swisshash($swissacc);
	  #$self->printswissprot($swisshash); # DEBUG
	  push (@swisshashes,$swisshash);
	}

	push (@ttables,$entryfeature->{'qualifiers'}->{'transl_table'}); # undef if not specified
	
	# create labels array
	for $exon (@range) {
	  push(@labels,$exon->[0],$exon->[1]); # start and end of every exon of the CDS
	}
	push (@transcripts,$entryfeature->{'range'});
      } else {
	# "simplifying" the joinedlocation features. I.e. changing them from
	# multijoined ones to simple plain start-end features, taking only
	# the start of the first "exon" and the end of the last "exon" as
	# start and end of the entire feature
	if ($entryfeature->{'locationtype'} && $entryfeature->{'locationtype'} eq "joined") { # joined location
	  @range=($range[0]->[0],$range[-1]->[1]);
	}
	push(@labels,$range[0],$range[1]); # start and end of every feature
	if ($name eq "exon") {
	  $desc=$entryfeature->{'qualifiers'}->{'number'};
	  if ($entryfeature->{'qualifiers'}->{'note'}) {
	    if ($desc) {
	      $desc .= "|" . $entryfeature->{'qualifiers'}->{'note'};
	    } else {
	      $desc = $entryfeature->{'qualifiers'}->{'note'};
	    }
	  }
	  push (@exondescs,$desc || "unknown");
	  push(@exons,\@range);
	}
	if ($name eq "intron") {
 	  $desc=$entryfeature->{'qualifiers'}->{'number'};
	  if ($desc) {
	    $desc .= "|" . $entryfeature->{'qualifiers'}->{'note'};
	  } else {
	    $desc = $entryfeature->{'qualifiers'}->{'note'};
	  }
	  push (@introndescs,$desc || "unknown"); 
	  push(@introns,\@range);
	}
	if (($name eq "prim_transcript")||($name eq "mRNA")) { push(@prim_transcripts,\@range); }
	if ($name eq "repeat_unit") { push(@repeat_units,\@range);
	  $rpt_family=$entryfeature->{'qualifiers'}->{'rpt_family'};
	  push (@repeat_units_family,$rpt_family || "unknown");
	}
	if ($name eq "repeat_region") { push(@repeat_regions,\@range);
	  $rpt_family=$entryfeature->{'qualifiers'}->{'rpt_family'};
	  push (@repeat_regions_family,$rpt_family || "unknown");
	}
      }
    }
  }
  unless ($gene_name) { $gene_name="cds-position:".$cds_position; }
  my %genefeatureshash;
  $genefeatureshash{gene_name}=$gene_name;
  $genefeatureshash{genefeatures}=\@genefeatures;
  $genefeatureshash{labels}=\@labels;
  $genefeatureshash{ttables}=\@ttables;
  $genefeatureshash{swisshashes}=\@swisshashes;
  $genefeatureshash{transcripts}=\@transcripts;
  $genefeatureshash{exons}=\@exons;
  $genefeatureshash{exondescs}=\@exondescs;
  $genefeatureshash{introns}=\@introns;
  $genefeatureshash{introndescs}=\@introndescs;
  $genefeatureshash{prim_transcripts}=\@prim_transcripts;
  $genefeatureshash{repeat_units}=\@repeat_units;
  $genefeatureshash{repeat_regions}=\@repeat_regions;
  $genefeatureshash{repeat_units_family}=\@repeat_units_family;
  $genefeatureshash{repeat_regions_family}=\@repeat_regions_family;
  return (\%genefeatureshash);
}


# used by _findgenefeatures, when a CDS at a certain position is requested,
# to retrieve only features quite close to the wanted CDS.
# Args: range hashref, begin and end positions of the CDS, $proximity
# $proximity holds the maximum distance between the extremes of the CDS
# and of the feature under exam.
# Returns: boolean
sub _checkfeatureproximity {
  my ($self,$range,$cds_begin,$cds_end,$proximity)=@_;
  my @range=@{$range};
  my ($begin,$end,$strand);
  if (ref($range[0]) eq "ARRAY") { # like in CDS, whose range equivals to exons
    ($begin,$end,$strand)=($range[0]->[0],$range[-1]->[1],$range[0]->[2]);
  } else {
    ($begin,$end,$strand)=@range;
  }
  if ($cds_begin > $cds_end) { # i.e. reverse strand CDS
    ($cds_begin,$cds_end)=($cds_end,$cds_begin); # swap boundaries
  }
  if ($strand == -1) { # reverse strand
    ($begin,$end)=($end,$begin); # swap boundaries
  }
  if (($cds_begin-$end)>$proximity) {
    carp "--DEBUG-- feature rejected: begin $begin end $end -------";
    return (0);
  }
  if (($begin-$cds_end)>$proximity) {
    carp "--DEBUG-- feature rejected: begin $begin end $end -------";
    return (0);
  }
  carp "--DEBUG-- feature accepted: begin $begin end $end -------";
  return (1); # otherwise ok, feature considered next to CDS
}


# function that calls the external program "align" (on the fasta2 package)
# to create an alignment between two sequences, returning the aligned
# strings and the codes for the identity (:: ::::)

sub _get_alignment {
  my ($self,$seq1,$seq2)=@_;

  my $null = ($^O =~ m/mswin/i) ? 'NUL' : '/dev/null';
  my $fastafile1="/tmp/tmpfastafile1";
  my $fastafile2="/tmp/tmpfastafile2";
  my $grepcut='egrep -v "[[:digit:]]|^ *$|sequences" | cut -c8-'; # grep/cut
  my $alignprogram="/usr/local/etc/bioinfo/fasta2/align -s /usr/local/etc/bioinfo/fasta2/idnaa.mat $fastafile1 $fastafile2 2>$null | $grepcut"; # ALIGN
  open my $TMPFASTAFILE1, '>', $fastafile1 or croak "Could not write file '$fastafile1' for aa alignment: $!";
  open my $TMPFASTAFILE2, '>', $fastafile2 or croak "Could not write file '$fastafile2' for aa alignment: $!";
  print $TMPFASTAFILE1 ">firstseq\n$seq1\n";
  print $TMPFASTAFILE2 ">secondseq\n$seq2\n";
  close $TMPFASTAFILE1;
  close $TMPFASTAFILE2;
  my $alignment=`$alignprogram`;
  my @alignlines=split(/\n/,$alignment);
  my ($linecount,$seq1_aligned,$seq2_aligned,$codes);
  for ($linecount=0; $linecount < @alignlines; $linecount+=3) {  
    $seq1_aligned .= $alignlines[$linecount];
    $codes .= $alignlines[$linecount+1];
    $seq2_aligned .= $alignlines[$linecount+2];
  }
  return ($seq1_aligned,$seq2_aligned,$codes);
}

# common part of the function to create a novel liveseq gene structure
# from an amino acid sequence, using codon usage frequencies
# args: codon_usage_array transltableid aasequence gene_name
sub _common_novelaasequence2gene {
  my ($species_codon_usage,$ttabid,$aasequence,$gene_name)=@_;
  my @species_codon_usage=@{$species_codon_usage};
  my @codon_usage_label= 
      qw (cga cgc cgg cgt aga agg cta ctc ctg ctt tta ttg tca tcc tcg
      tct agc agt aca acc acg act cca ccc ccg cct gca gcc gcg gct gga
      ggc ggg ggt gta gtc gtg gtt aaa aag aac aat caa cag cac cat gaa
      gag gac gat tac tat tgc tgt ttc ttt ata atc att atg tgg taa tag
      tga);
  my ($i,$j);
  my %codon_usage_value;
  my $aa_codon_total;
  for ($i=0;$i<64;$i++) {
    $codon_usage_value{$codon_usage_label[$i]}=$species_codon_usage[$i];
  }

  my $CodonTable  = Bio::Tools::CodonTable->new ( -id => $ttabid );
  my @aminoacids = split(//,uc($aasequence));
  my @alt_codons; my ($relativeusage,$dnasequence,$chosen_codon,$dice,$partial,$thiscodon);
  for $i (@aminoacids) {
    @alt_codons = $CodonTable->revtranslate($i);
    unless (@alt_codons) {
      carp "No reverse translation possible for aminoacid \'$i\'";
      $dnasequence .= "???";
    } else {
      $aa_codon_total=0;
      for $j (@alt_codons) {
	$aa_codon_total+=$codon_usage_value{$j};
      }
      # print "aminoacid $i, codonchoice: "; # verbose
      #$partial=0;
      #for $j (@alt_codons) {
	#printf "%s %.2f ",$j,$partial+$codon_usage_value{$j}/$aa_codon_total;
	#$partial+=($codon_usage_value{$j}/$aa_codon_total);
      #}
      #print "\n";
      $dice=rand(1);
      #print "roulette: $dice\n"; # verbose
      $partial=0;
      $chosen_codon="";
      CODONCHOICE:
      for $j (0..@alt_codons) { # last one not accounted
	$thiscodon=$alt_codons[$j];
	$relativeusage=($codon_usage_value{$thiscodon}/$aa_codon_total);
	if ($dice < $relativeusage+$partial) {
	  $chosen_codon=$thiscodon;
	  last CODONCHOICE;
	} else {
	  $partial += $relativeusage;
	}
      }
      unless ($chosen_codon) {
	$chosen_codon = $alt_codons[-1]; # the last one
      }
      # print ".....adding $chosen_codon\n"; # verbose
      $dnasequence .= $chosen_codon;
    }
  }

  my $dna = Bio::LiveSeq::DNA->new(-seq => $dnasequence);
  my $min=1;
  my $max=length($dnasequence);
  my $exon = Bio::LiveSeq::Exon->new(-seq => $dna, -start => $min, -end => $max, -strand => 1);
  my @exons=($exon);
  my $transcript = Bio::LiveSeq::Transcript->new(-exons => \@exons);
  $transcript->translation_table($ttabid);
  my @transcripts=($transcript);
  my $translation = Bio::LiveSeq::Translation->new(-transcript => $transcript);
  my @translations=($translation);
  my %features=(DNA => $dna, Transcripts => \@transcripts, Translations => \@translations);
  my $gene = Bio::LiveSeq::Gene->new(-name => $gene_name, -features => \%features, -upbound => $min, -downbound => $max);

  # creation of gene
  unless ($gene) { # if $gene == 0 it means problems in hash2gene
    carp "Error in Gene creation phase";
    return (0);
  }
  return $gene;
}

1;
