# $Id$
#
# bioperl module for Bio::LiveSeq::IO::SRS
#
# Cared for by Joseph Insana <insana@ebi.ac.uk> <jinsana@gmx.net>
#
# Copyright Joseph Insana
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

  Bio::LiveSeq::IO::SRS - Loader for LiveSeq from EMBL entries with SRS

=head1 SYNOPSIS

  my $db="EMBL";
  my $acc_id="M20132";
  my $query="embl-acc:$acc_id";

  my $loader=Bio::LiveSeq::IO::SRS->load(-db=>"EMBL", -query=>"$query");

  my @translationobjects=$loader->entry2liveseq();

  my $gene="AR";
  my $gene=$loader->gene2liveseq("gene");

  NOTE: The only -db now supported is EMBL. Hence it defaults to EMBL.

=head1 DESCRIPTION

This package uses the SRS (Sequence Retrieval System) to fetch a sequence
database entry, analyse it and create LiveSeq objects out of it.

An embl-acc ID has to be passed to this package which will return references
to all translation objects created from the EMBL entry.
Transcription, DNA and Exon objects' references can all be retrieved departing
from these.

Alternatively, a specific "gene" name can be specified, together with the
embl-acc ID. This will create a LiveSeq::Gene object with all relevant gene
features attached/created.

=head1 AUTHOR - Joseph A.L. Insana

Email:  Insana@ebi.ac.uk, jinsana@gmx.net

Address: 

     EMBL Outstation, European Bioinformatics Institute
     Wellcome Trust Genome Campus, Hinxton
     Cambs. CB10 1SD, United Kingdom 

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LiveSeq::IO::SRS;
$VERSION=2.4;

# Version history:
# Wed Apr  5 13:06:43 BST 2000 v 1.0 restarted as a child of Loader.pm
# Wed Apr  5 15:57:22 BST 2000 v 1.1 load() created
# Thu Apr  6 02:52:11 BST 2000 v 1.2 new field "range" in hash
# Thu Apr  6 03:11:29 BST 2000 v 1.22 transition from $location to @range
# Thu Apr  6 03:40:04 BST 2000 v 1.25 added Division
# Tue Apr 18 17:15:26 BST 2000 v 2.0 started coding swissprot2hash
# Thu Apr 20 02:17:28 BST 2000 v 2.1 mRNA added to valid_feature_names
# Wed Jun  7 02:08:57 BST 2000 v 2.2 added support for joinedlocation features
# Thu Jun 29 19:07:59 BST 2000 v 2.22 Gene stripped of possible newlines (horrible characteristic of some entries!!!!)
# Fri Jun 30 14:08:21 BST 2000 v 2.3 SRS querying routines now conveniently wrapped in eval{} blocks to avoid SRS crash the whole LiveSeq
# Tue Jul  4 14:07:52 BST 2000 v 2.31 note&number added in val_qual_names
# Mon Sep  4 17:46:42 BST 2000 v 2.4 novelaasequence2gene() added

use strict;
use Carp qw(cluck croak carp);
use vars qw($VERSION @ISA);
use lib $ENV{SRSEXE};
use srsperl;
use Bio::Tools::CodonTable; # for novelaasequence2gene

use Bio::LiveSeq::IO::Loader 2.2;

@ISA=qw(Bio::LiveSeq::IO::Loader);

# This package can in the future host other databases loading subroutines.
# e.g. ensembl2hash

=head1 load

  Title   : load
  Usage   : my $acc_id="M20132";
            my $query="embl-acc:$acc_id";
            $loader=Bio::LiveSeq::IO::SRS->load(-db=>"EMBL", -query=>"$query");

  Function: loads an entry with SRS from a database into a hash
  Returns : reference to a new object of class IO::SRS holding an entry
  Errorcode 0
  Args    : an SRS query resulting in one fetched EMBL (by default) entry

=cut

sub load {
  my ($thing, %args) = @_;
  my $class = ref($thing) || $thing;
  my ($obj,%loader);

  my ($db,$query)=($args{-db},$args{-query});

  if (defined($db)) {
    unless ($db eq "EMBL") {
      carp "Note: only EMBL now supported!";
      return(0);
    }
  } else {
    $db="EMBL";
  }

  my $hashref;
  if ($db eq "EMBL") {
    my $test_transl=0; # change to 0 to avoid comparison of translation

    # these can be changed for future needs
    my @embl_valid_feature_names=qw(CDS exon prim_transcript intron repeat_unit repeat_region mRNA);
    my @embl_valid_qual_names=qw(gene codon_start db_xref product note number rpt_family transl_table);

    # dunno yet how to implement test_transl again....
    # probably on a one-on-one basis with each translation?
    if ($test_transl) {
      push (@embl_valid_qual_names,"translation"); # needed for test_transl
    }

    # not to have the whole program die because of SRS death
    eval {
      $hashref=&embl2hash("$query",\@embl_valid_feature_names,\@embl_valid_qual_names);
    };
    my $errormsg;
    if ($@) {
      foreach $errormsg (split(/\n/,$@)) {
	if (index($errormsg,"in cleanup") == -1) {
	  carp "SRS EMBL Loader: $@";
	}
      }
    }
  }
  unless ($hashref) { return (0); }

  %loader = (db => $db, query => $query, hash => $hashref);
  $obj = \%loader;
  $obj = bless $obj, $class;
  return $obj;
}

=head1 embl2hash

  Title   : embl2hash
  Function: retrieves with SRS an EMBL entry, parses it and creates
            a hash that contains all the information.
  Returns : a reference to a hash
  Errorcode: 0
  Args    : an SRS query resulting in one fetched EMBL entry
              i.e. a string in a form like "embl-acc:accession_number"
	    two array references to skip features and qualifiers (for
	    performance)
  Example: @valid_features=qw(CDS exon prim_transcript mRNA);
           @valid_qualifiers=qw(gene codon_start db_xref product rpt_family);
           $hashref=&embl2hash("$query",\@valid_features,\@valid_qualifiers);

=cut

# this has to be defined here as it is the only thing really proper to
# the /SRS/ loader. Every loader has to sport his own "entry2hash" function
# the main functions will be in the Loader.pm
# arguments: embl SRS query resulting in one fetched EMBL entry
# to skip features and qualifiers (for performance), two array
# references must be passed (this can change into string arguments to
# be passed....)
# returns: a reference to a hash containing the important features requested
sub embl2hash {
  my ($i,$j);
  my $query=$_[0];
  my %valid_features; my %valid_names;
  if ($_[1]) {
    %valid_features = map {$_, 1} @{$_[1]}; # to skip features
  }
  if ($_[2]) {
    %valid_names = map {$_, 1} @{$_[2]}; # to skip qualifiers
  }
  # SRS API used to fetch all relevant fields
  my $sess = new Session;
  my $set = $sess->query("[$query]", "");
  my $numEntries=$set->size();
  if ($numEntries < 1) {
    carp "No sequence entry found for the query $query";
    return (0);
  } elsif ($numEntries > 1) {
    carp "Too many entries found for the input given";
    return (0);
  } else {
    my $entry = $set->getEntry(0);
    my ($entry_ID,$entry_AccNumber,$entry_Molecule,$entry_Description,$entry_Organism,$entry_SeqLength,$entry_Sequence,$entry_Division);
    # Fetch what we can fetch without the loader
    $entry_ID = $entry->fieldToString("id","","","");
    $entry_AccNumber = $entry->fieldToString("acc","","","");
    $entry_Molecule = $entry->fieldToString("mol","","","");
    $entry_Division = $entry->fieldToString("div","","","");
    $entry_Description = $entry->fieldToString("des","","","");
    $entry_Description =~ s/\n/ /g;
    $entry_Organism = $entry->fieldToString("org","","","");
    $entry_SeqLength = $entry->fieldToString("sl","","","");
    # Now use the loader
    my $loadedentry = $entry->load("EMBL");
    # Fetch the rest via the loader
    $entry_Sequence = $loadedentry->attrStr("sequence");
    $entry_Sequence =~ s/\n//g; # from plain format to raw string

    # put into the hash
    my %entryhash;
    $entryhash{ID}=$entry_ID;
    $entryhash{AccNumber}=$entry_AccNumber;
    $entryhash{Molecule}=$entry_Molecule;
    $entryhash{Division}=$entry_Division;
    $entryhash{Description}=$entry_Description;
    $entryhash{Organism}=$entry_Organism;
    $entryhash{Sequence}=$entry_Sequence;
    $entryhash{SeqLength}=$entry_SeqLength;

    # create features array
    my $features = $loadedentry->attrObjList("features");
    my $featuresnumber= $features->size();
    $entryhash{FeaturesNumber}=$featuresnumber;
    my ($feature,$feature_name,$feature_location);
    my ($feature_qual_names,$feature_qual_values);
    my ($feature_qual_name,$feature_qual_value);
    my ($feature_qual_number,$feature_qual_number2);
    my @features;
    for $i (0..$featuresnumber-1) {
      my %feature;
      $feature = $features->get($i);
      $feature_name = $feature->attrStr("FtKey");

      #unless ($feature_name eq "CDS") {
      unless ($valid_features{$feature_name}) {
	#print "not valid feature: $feature_name\n";
	next;
      }
      #print "now processing feature: $feature_name\n";
      $feature_location = $feature->attrStr("FtLocation");
      $feature_location =~ s/\n//g;
      $feature_qual_names= $feature->attrStrList("FtQualNames");
      $feature_qual_values= $feature->attrStrList("FtQualValues");
      $feature_qual_number= $feature_qual_names->size();
      $feature_qual_number2= $feature_qual_values->size(); # paranoia check
      if ($feature_qual_number > $feature_qual_number2) {
	carp ("Warning with Feature at position $i ($feature_name): There are MORE qualifier names than qualifier values.");
	# this can happen, e.g. "/partial", let's move on, just issue a warning
	#return (0);
      } elsif ($feature_qual_number < $feature_qual_number2) {
	carp ("Error with Feature at position $i ($feature_name): There are LESS qualifier names than qualifier values. Stopped");
	return (0);
      }
      #} else {print "NUMBER OF QUALIFIERS: $feature_qual_number \n";} # DEBUG

      # Put things inside hash
      $feature{position}=$i;
      $feature{name}=$feature_name;

      # new range field in hash
      my @range;
      if ($feature_name eq "CDS") {
	@range=cdslocation2transcript($feature_location);
	$feature{locationtype}="joined";
      } else {
	if (index($feature_location,"join") == -1) {
	  @range=location2range($feature_location);
	} else { # complex location in feature different than CDS
	  @range=joinedlocation2range($feature_location);
	  $feature{locationtype}="joined";
	}
      }
      $feature{range}=\@range;
      $feature{location}="deprecated";
      my %feature_qualifiers;
      for $j (0..$feature_qual_number-1) {
	$feature_qual_name=$feature_qual_names->get($j);
	$feature_qual_name =~ s/^\///; # eliminates heading "/"

	# skip all not interesting (for now) qualifiers
        unless ($valid_names{$feature_qual_name}) {
	  #print "not valid name: $feature_qual_name\n";
	  next;
	}
	#print "now processing: $feature_qual_name\n";
	$feature_qual_value=$feature_qual_values->get($j);
	$feature_qual_value =~ s/^"//; # eliminates heading "
	$feature_qual_value =~ s/"$//; # eliminates trailing "
	$feature_qual_value =~ s/\n//g; # eliminates all new lines
	$feature_qualifiers{$feature_qual_name}=$feature_qual_value;
      }
      $feature{qualifiers}=\%feature_qualifiers;
      push (@features,\%feature); # array of features
    }
    $entryhash{Features}=\@features; # put this also into the hash
    my @cds; # array just of CDSs
    for $i (0..$#features) {
      if ($features[$i]->{name} eq "CDS") {
	push(@cds,$features[$i]);
      }
    }
    $entryhash{CDS}=\@cds; # put this also into the hash
    return (\%entryhash);
  }
}

# argument: location field of an EMBL feature
# returns: array with correct $start $end and $strand to create LiveSeq obj
sub location2range {
  my ($location)=@_;
  my ($start,$end,$strand);
  if (index($location,"complement") == -1) { # forward strand
    $strand=1;
  } else { # reverse strand
    $location =~ s/complement\(//g;
    $location =~ s/\)//g;
    $strand=-1;
  }
  $location =~ s/\<//g;
  $location =~ s/\>//g;
  my @range=split(/\.\./,$location);
  if (scalar(@range) == 1) { # special case of range with just one position (e.g. polyA_site EMBL features
    $start=$end=$range[0];
  } else {
    if ($strand == 1) {
      ($start,$end)=@range;
    } else { # reverse strand
      ($end,$start)=@range;
    }
  }
  return ($start,$end,$strand);
}

# argument: location field of a CDS feature
# returns: array of exons arrayref, useful to create LiveSeq Transcript obj
sub cdslocation2transcript {
  my ($location)=@_;
  my @exonlocs;
  my $exonloc;
  my @exon;
  my @transcript=();
  $location =~ s/^join\(//;
  $location =~ s/\)$//;
  @exonlocs = split (/\,/,$location);
  for $exonloc (@exonlocs) {
    my @exon=location2range($exonloc);
    push (@transcript,\@exon);
  }
  return (@transcript);
}

# argument: location field of a CDS feature
# returns: array of exons arrayref, useful to create LiveSeq Transcript obj
sub joinedlocation2range {
  &cdslocation2transcript;
}


=head1 get_swisshash

  Title   : get_swisshash
  Usage   : $loader->get_swisshash();
  Example : $swisshash=$loader->swissprot2hash("SWISS-PROT:P10275")
  Function: retrieves with SRS a SwissProt entry, parses it and creates
            a hash that contains all the information.
  Returns : a reference to a hash
  Errorcode: 0
  Args    : the db_xref qualifier's value from an EMBL CDS Feature
            i.e. a string in the form "SWISS-PROT:accession_number"
  Note    : this can be modified (adding a second argument) to retrieve
            and parse SWTREMBL, SWALL... entries

=cut

# argument: db_xref qualifier's value from EMBL CDS
# errorcode: 0
# returns hashref
sub get_swisshash {
  my ($self,$query)=@_;
  if (index($query,"SWISS-PROT") == -1) {
    return (0);
  }
  $query =~ s/SWISS-PROT/swissprot-acc/;
  my $hashref;
  eval {
    $hashref=&swissprot2hash("$query");
  };
  my $errormsg;
  if ($@) {
    foreach $errormsg (split(/\n/,$@)) {
      if (index($errormsg,"in cleanup") == -1) {
	carp "SRS Swissprot Loader: $@";
      }
    }
  }
  unless ($hashref) {
    return (0);
  } else {
    return ($hashref);
  }
}

=head1 swissprot2hash

  Title   : swissprot2hash
  Usage   : $loader->swissprot2hash();
  Example : $swisshash=$loader->swissprot2hash("swissprot-acc:P10275")
  Function: retrieves with SRS a SwissProt entry, parses it and creates
            a hash that contains all the information.
  Returns : a reference to a hash
  Errorcode: 0
  Args    : an SRS query resulting in one entry from SwissProt database
            i.e. a string in the form "swissprot-acc:accession_number"
  Note    : this can be modified (adding a second argument) to retrieve
            and parse SWTREMBL, SWALL... entries

=cut

# arguments: swissprot SRS query resulting in one fetched swissprot entry
# returns: a reference to a hash containing the important features requested
sub swissprot2hash {
  my ($i,$j);
  my $query=$_[0];
  # SRS API used to fetch all relevant fields
  my $sess = new Session;
  my $set = $sess->query("[$query]", "");
  my $numEntries = $set->size();
  if ($numEntries < 1) {
    carp "No sequence entry found for the query $query";
    return (0);
  } elsif ($numEntries > 1) {
    carp "Too many entries found for the input given";
    return (0);
  } else {
    my $entry = $set->getEntry(0);
    my ($entry_ID,$entry_AccNumber,$entry_Molecule,$entry_Description,$entry_Organism,$entry_SeqLength,$entry_Sequence,$entry_Gene);
    # Fetch what we can fetch without the loader
    $entry_ID = $entry->fieldToString("id","","","");
    $entry_AccNumber = $entry->fieldToString("acc","","","");
    $entry_Gene = $entry->fieldToString("gen","","","");
    $entry_Gene =~ s/\n/ /g;
    $entry_Description = $entry->fieldToString("des","","","");
    $entry_Description =~ s/\n/ /g;
    $entry_Organism = $entry->fieldToString("org","","","");
    chop $entry_Organism;
    $entry_SeqLength = $entry->fieldToString("sl","","","");
    # Now use the loader
    my $loadedentry = $entry->load("Swissprot");
    # Fetch the rest via the loader
    $entry_Sequence = $loadedentry->attrStr("Sequence");
    $entry_Sequence =~ s/\n//g; # from plain format to raw string

    # put into the hash
    my %entryhash;
    $entryhash{ID}=$entry_ID;
    $entryhash{AccNumber}=$entry_AccNumber;
    $entryhash{Molecule}=$entry_Molecule;
    $entryhash{Gene}=$entry_Gene;
    $entryhash{Description}=$entry_Description;
    $entryhash{Organism}=$entry_Organism;
    $entryhash{Sequence}=$entry_Sequence;
    $entryhash{SeqLength}=$entry_SeqLength;

    # create features array
    my $features = $loadedentry->attrObjList("Features");
    my $featuresnumber= $features->size();
    $entryhash{FeaturesNumber}=$featuresnumber;
    my ($feature,$feature_name,$feature_description,$feature_location);
    my @features;
    for $i (0..$featuresnumber-1) {
      my %feature;
      $feature = $features->get($i);
      $feature_name = $feature->attrStr("FtKey");
      $feature_location = $feature->attrStr("FtLocation");
      $feature_location =~ s/ +/ /g;
      $feature_description = $feature->attrStr("FtDescription");
      chop $feature_description;
      $feature_description =~ s/\nFT                                / /g;

      # Put things inside hash
      $feature{position}=$i;
      $feature{name}=$feature_name;
      $feature{location}=$feature_location;
      $feature{description}=$feature_description;

      push (@features,\%feature); # array of features
    }
    $entryhash{Features}=\@features; # put this also into the hash
    return (\%entryhash);
  }
}

=head1 novelaasequence2gene

  Title   : novelaasequence2gene
  Usage   : $gene=Bio::LiveSeq::IO::SRS->novelaasequence2gene(-aasequence => "MGLAAPTRS*");
          : $gene=Bio::LiveSeq::IO::SRS->novelaasequence2gene(-aasequence => "MGLAAPTRS*",
                                             -genome => "Homo sapiens");
          : $gene=Bio::LiveSeq::IO::SRS->novelaasequence2gene(-aasequence => "MGLAAPTRS*",
                                             -genome => "Mitochondrion Homo sapiens",
                                             -gene_name => "tyr-kinase");

  Function: creates LiveSeq objects from a novel amino acid sequence,
            using codon usage database to choose codons according to
            relative frequencies.
            If a genome latin name is not specified, the default is to use
            the Homo sapiens' one (taxonomy ID 9606).
  Returns : reference to a Gene object containing references to LiveSeq objects
  Errorcode 0
  Args    : string containing an amino acid sequence
            string (optional) with a species/genome latin name
            string specifying a gene name
  Note    : SRS access to TAXON and CODONUSAGE databases is required

=cut

sub novelaasequence2gene {
  my ($self, %args) = @_;
  my ($gene_name,$species_name,$aasequence)=($args{-gene_name},$args{-genome},$args{-aasequence});
  unless ($aasequence) {
    carp "aasequence not given";
    return (0);
  }
  unless ($gene_name) {
    $gene_name="Novel Unknown";
  }
  unless ($species_name) {
    $species_name="Homo sapiens";
  }
  
  my $sess = new Session;
  my ($e,$numEntries,$set);

  # codonusage lookup
  my $codonusage_usage;
  my @species_codon_usage;
  $set = $sess->query("[codonusage:'$species_name']", "");
  $numEntries = $set->size();
  if ($numEntries > 0) {
    $e = $set->getEntry(0);
    $species_name = $e->fieldToString("id","","","");
    $codonusage_usage = $e->fieldToString("usage","","","");
    @species_codon_usage=split(/\s/,$codonusage_usage); # spaces or tabs
    if ($numEntries > 1) {
      carp "Search in Codon Usage DB resulted in $numEntries results. Taking the first one: $species_name";
    }
  } else {
      carp "Genome not found in codon usage DB.";
      return (0);
  }  

  # taxonomy lookup
  my $mito_flag = 0;
  my $species_origin;
  if ($species_name =~ /^Mitochondrion /) {
    $mito_flag = 1;
    $species_origin = $species_name;
    $species_origin =~ s/^Mitochondrion //;
    $set = $sess->query("[taxonomy-species:'$species_origin']", "");
  } elsif ($species_name =~ /^Chloroplast |^Kinetoplast |^Chromoplast /) {
    $species_origin = $species_name;
    $species_origin =~ s/^Chromoplast //;
    $species_origin =~ s/^Kinetoplast //;
    $species_origin =~ s/^Chloroplast //;
    $set = $sess->query("[taxonomy-species:'$species_origin']", "");
  } else {
    $set = $sess->query("[taxonomy-species:'$species_name']", "");
  }
  $numEntries = $set->size();
  my ($taxonomy_id,$taxonomy_gc,$taxonomy_mgc,$taxonomy_name);
  if ($numEntries > 0) {
    $e = $set->getEntry(0);
    $taxonomy_id = $e->fieldToString("id","","","");
    $taxonomy_name = $e->fieldToString("species","","","");
    $taxonomy_gc = $e->fieldToString("gc","","","");
    $taxonomy_mgc = $e->fieldToString("mgc","","","");
    if ($numEntries > 1) {
      carp "Note! More than one entry found in taxonomy DB for the genome query given. Using the first one: $taxonomy_name ($taxonomy_id)";
    }
  } else {
    carp "Genome not found in taxonomy DB.";
    return (0);
  }
  # Lookup appropriate translation table
  my $ttabid;
  if ($mito_flag) {
    $ttabid = $taxonomy_mgc;
  } else {
    $ttabid = $taxonomy_gc;
  }

  my $gene=Bio::LiveSeq::IO::Loader::_common_novelaasequence2gene(\@species_codon_usage,$ttabid,$aasequence,$gene_name);
  return ($gene);
}

 
