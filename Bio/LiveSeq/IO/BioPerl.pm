# $Id$
#
# bioperl module for Bio::LiveSeq::IO::BioPerl
#
# Cared for by Joseph Insana <insana@ebi.ac.uk> <jinsana@gmx.net>
#
# Copyright Joseph Insana
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::LiveSeq::IO::BioPerl - Loader for LiveSeq from EMBL entries with BioPerl

=head1 SYNOPSIS

  my $db="EMBL";
  my $file="../data/M20132";

  my $loader=Bio::LiveSeq::IO::BioPerl->load(-db=>"EMBL", -file=>"$file");

  my @translationobjects=$loader->entry2liveseq();

  my $gene="AR";
  my $gene=$loader->gene2liveseq("gene");

  NOTE: The only -db now supported is EMBL. Hence it defaults to EMBL.

=head1 DESCRIPTION

This package uses BioPerl (SeqIO) to fetch a sequence database entry,
analyse it and create LiveSeq objects out of it.

An filename has to be passed to this package which will return
references to all translation objects created from the EMBL
entry. References to Transcription, DNA and Exon objects can all be
retrieved departing from these.

Alternatively, a specific "gene" name can be specified, together with
the embl-acc ID. This will create a LiveSeq::Gene object with all
relevant gene features attached/created.

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

package Bio::LiveSeq::IO::BioPerl;
$VERSION=2.3;

# Version history:
# Thu Apr  6 00:25:46 BST 2000 v 1.0 begun
# Thu Apr  6 03:40:04 BST 2000 v 1.25 added Division
# Thu Apr  6 03:40:36 BST 2000 v 2.0 working
# Thu Apr 20 02:17:28 BST 2000 v 2.1 mRNA added to valid_feature_names
# Tue Jul  4 14:07:52 BST 2000 v 2.11 note&number added in val_qual_names
# Fri Sep 15 15:41:02 BST 2000 v 2.22 novelaasequence2gene now works without SRS
# Mon Jan 29 17:40:06 EST 2001 v 2.3 made it work with the new split_location of BioPerl 0.7

# TODO->TOCHECK
# each_secondary_access not working
# why array from each_tag_value($qual) ? When will there be more than one
#                                        element in such array?
# what is the annotation object? ($seqobj->annotation)
# unsatisfied by both BioPerl binomial and SRS "org" to retrieve Organism info

use strict;
use Carp qw(cluck croak carp);
use vars qw($VERSION @ISA);
use Bio::SeqIO;

use Bio::LiveSeq::IO::Loader 2.0;

@ISA=qw(Bio::LiveSeq::IO::Loader);

# This package can in the future host other databases loading subroutines.
# e.g. ensembl2hash

=head2 load

  Title   : load
  Usage   : my $filename="../data/M20132";
            $loader=Bio::LiveSeq::IO::BioPerl->load(-db=>"EMBL", -file=>"$filename");

  Function: loads an entry with BioPerl from a database into a hash
  Returns : reference to a new object of class IO::BioPerl holding an entry
  Errorcode 0
  Args    : an BioPerl query resulting in one fetched EMBL (by default) entry

=cut

sub load {
  my ($thing, %args) = @_;
  my $class = ref($thing) || $thing;
  my ($obj,%loader);

  my ($db,$filename)=($args{-db},$args{-file});

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
    my @embl_valid_feature_names=qw(CDS CDS_span exon prim_transcript intron repeat_unit repeat_region mRNA);
    my @embl_valid_qual_names=qw(gene codon_start db_xref product note number rpt_family transl_table);

    # dunno yet how to implement test_transl again....
    # probably on a one-on-one basis with each translation?
    if ($test_transl) {
      push (@embl_valid_qual_names,"translation"); # needed for test_transl
    }
    $hashref=&embl2hash("$filename",\@embl_valid_feature_names,\@embl_valid_qual_names);
  }
  unless ($hashref) { return (0); }

  %loader = (db => $db, filename => $filename, hash => $hashref);
  $obj = \%loader;
  $obj = bless $obj, $class;
  return $obj;
}

=head2 embl2hash

  Title   : embl2hash
  Function: retrieves with BioPerl an EMBL entry, parses it and creates
            a hash that contains all the information.
  Returns : a reference to a hash
  Errorcode: 0
  Args    : a filename pointing to a file containing one EMBL entry
	    two array references to skip features and qualifiers (for
	    performance)
  Example: @valid_features=qw(CDS exon prim_transcript mRNA);
           @valid_qualifiers=qw(gene codon_start db_xref product rpt_family);
           $hashref=&embl2hash("$file",\@valid_features,\@valid_qualifiers);

=cut

# arguments: embl filename containing one EMBL entry
# to skip features and qualifiers (for performance), two array
# references must be passed (this can change into string arguments to
# be passed....)
# returns: a reference to a hash containing the important features requested
sub embl2hash {
  my $filename=$_[0];
  my %valid_features; my %valid_names;
  if ($_[1]) {
    %valid_features = map {$_, 1} @{$_[1]}; # to skip features
  }
  if ($_[2]) {
    %valid_names = map {$_, 1} @{$_[2]}; # to skip qualifiers
  }

  my $stream = Bio::SeqIO->new('-file' => $filename, '-format' => 'EMBL');
 
  my $seqobj = $stream->next_seq();
  my $annobj = $seqobj->annotation(); # what's this?

  my $entry_Sequence = lc($seqobj->seq()); # SRS returns lowercase

  my $entry_ID = $seqobj->display_id;
  my $entry_AccNumber = $seqobj->accession; # or maybe accession_number ?
  my $secondary_acc; # to fetch the other acc numbers
  foreach $secondary_acc ($seqobj->get_secondary_accessions) { # not working!
    $entry_AccNumber .= " $secondary_acc";
  }
  my $entry_Molecule = $seqobj->molecule; # this alone returns molec+division
  my $entry_Division = $seqobj->division;
  # fixed: now Molecule works in BioPerl, no need for next lines
  #my @Molecule=split(" ",$entry_Molecule);
  #my $entry_Division = pop(@Molecule); # only division
  #$entry_Molecule = join(" ",@Molecule); # only molecule
  my $entry_Description = $seqobj->desc;

  my $speciesobj = $seqobj->species;
  my $entry_Organism = $speciesobj->binomial;

  my $entry_SeqLength = $seqobj->length;
  
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

  my @topfeatures=$seqobj->top_SeqFeatures();
  # create features array
  my $featuresnumber= scalar(@topfeatures);
  $entryhash{FeaturesNumber}=$featuresnumber;
  my $feature_name;
  my @feature_qual_names; my @feature_qual_value;
  my ($feature_qual_name,$feature_qual_number);
  my @features;

  my ($feat,$qual,$subfeat);
  my @subfeat;
  my $i=0;
  foreach $feat (@topfeatures) {
      my %feature;
      $feature_name = $feat->primary_tag;
      unless ($valid_features{$feature_name}) {
	  #print "skipping $feature_name\n";
	  next;
      }
# works ok with 0.6.2
#    if ($feature_name eq "CDS_span") { # case of CDS with various exons 0.6.2
#      $feature_name="CDS"; # 0.6.2
      my $featlocation=$feat->location; # 0.7
      if (($feature_name eq "CDS")&&($featlocation->isa('Bio::Location::SplitLocationI'))) { # case of CDS with various exons BioPerl 0.7
#      @subfeat=$feat->sub_SeqFeature; # 0.6.2
	  @subfeat=$featlocation->sub_Location(); # 0.7
	  my @transcript;
	  foreach $subfeat (@subfeat) {
	      my @range;
	      if ($subfeat->strand == -1) {
		  @range=($subfeat->end,$subfeat->start,$subfeat->strand);
	      } else {
		  @range=($subfeat->start,$subfeat->end,$subfeat->strand);
	      }
	      push (@transcript,\@range);
	  }
	  $feature{range}=\@transcript;
      } else {
	  my @range;
	  ($feat->strand == -1) ? (@range = ($feat->end, $feat->start, $feat->strand) ) :
	      (@range = ( $feat->start,$feat->end,$feat->strand) );
# works ok with 0.6.2
	  if ($feature_name eq "CDS") { # case of single exon CDS (CDS name but not split location)
	      my @transcript=(\@range);
	      $feature{range}=\@transcript;
	  } else { # all other range features
	      $feature{range}=\@range;
	  }
      }
      $feature{location}="deprecated";
      
      $feature{position}=$i;
      $feature{name}=$feature_name;
      
      @feature_qual_names= $feat->all_tags();
      $feature_qual_number= scalar(@feature_qual_names);
      
      $feature{qual_number}=$feature_qual_number;
      
      my %feature_qualifiers;
      for $qual (@feature_qual_names) {
	  $feature_qual_name=$qual;
	  unless ($valid_names{$feature_qual_name}) {
	      next;
	  }
      @feature_qual_value=$feat->each_tag_value($qual);
	  #print "$qual => @feature_qual_value \n";
	  $feature_qualifiers{$feature_qual_name}=$feature_qual_value[0]; # ?
      # maybe the whole array should be entered, not just the 1st element?
	  # what could be the other elements? TOCHECK!
      }
      $feature{qualifiers}=\%feature_qualifiers;
      push (@features,\%feature); # array of features
      $i++;
  }
  $entryhash{Features}=\@features; # put this also into the hash
  
  my @cds; # array just of CDSs
  for $i (0..$#features) {
      if ($features[$i]->{'name'} eq "CDS") {
	  push(@cds,$features[$i]);
      }
  }
  $entryhash{CDS}=\@cds; # put this also into the hash
  return (\%entryhash);
}

=head2 novelaasequence2gene

  Title   : novelaasequence2gene
  Usage   : $gene=Bio::LiveSeq::IO::BioPerl->novelaasequence2gene(-aasequence => "MGLAAPTRS*");
          : $gene=Bio::LiveSeq::IO::BioPerl->novelaasequence2gene(-aasequence => "MGLAAPTRS*",
                                             -cusg_data => "58 44 7 29 3 3 480 267 105 143 122 39 144 162 14 59 53 25 233 292 19 113 88 246 28 68 161 231 27 102 128 151 67 60 138 131 48 61 153 19 233 73 150 31 129 38 147 71 138 43 181 81 44 15 255 118 312 392 236 82 20 10 14 141");
          : $gene=Bio::LiveSeq::IO::BioPerl->novelaasequence2gene(-aasequence => "MGLAAPTRS*",
                                             -cusg_data => "58 44 7 29 3 3 480 267 105 143 122 39 144 162 14 59 53 25 233 292 19 113 88 246 28 68 161 231 27 102 128 151 67 60 138 131 48 61 153 19 233 73 150 31 129 38 147 71 138 43 181 81 44 15 255 118 312 392 236 82 20 10 14 141",
                                             -translation_table => "2",
                                             -gene_name => "tyr-kinase");

  Function: creates LiveSeq objects from a novel amino acid sequence,
            using codon usage information (loaded from a file) to choose
            codons according to relative frequencies.
            If a codon_usage information is not specified,
            the default is to use Homo sapiens data (taxonomy ID 9606).
            If a translation_table ID is not specified, it will default to 1
            (standard code).
  Returns : reference to a Gene object containing references to LiveSeq objects
  Errorcode 0
  Args    : string containing an amino acid sequence
	    string (optional) with codon usage data (64 integer numbers)
            string (optional) specifying a gene_name
            integer (optional) specifying a translation_table ID

=cut

sub novelaasequence2gene {
  my ($self, %args) = @_;
  my ($gene_name,$cusg_data,$aasequence,$ttabid)=($args{-gene_name},$args{-cusg_data},$args{-aasequence},$args{-translation_table});

  my @species_codon_usage;
  unless ($aasequence) {
    carp "aasequence not given";
    return (0);
  }
  unless ($gene_name) {
    $gene_name="Novel Unknown";
  }
  unless ($ttabid) {
    $ttabid=1;
  }
  unless ($cusg_data) {
    @species_codon_usage=
	qw(68664 118404 126679 51100 125600 123646 75667 210903 435317
	139009 79303 135218 128429 192616 49456 161556 211962 131222
	162837 213626 69346 140780 182506 219428 76684 189374 173010
	310626 82647 202329 180955 250410 180001 118798 76398 160764
	317359 119013 262630 359627 218376 186915 130857 377006 162826
	113684 317703 441298 287040 245435 174805 133427 134523 108740
	225633 185619 78463 240138 174021 244236 142435 8187 5913
	14381); # updated 21Jul2000
  } else {
    @species_codon_usage=split(/ /,$cusg_data);
  }
  
  my $gene=Bio::LiveSeq::IO::Loader::_common_novelaasequence2gene(\@species_codon_usage,$ttabid,$aasequence,$gene_name);
  return ($gene);
}

1;
