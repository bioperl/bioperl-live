# Working with development sources

BioPerl uses [Dist::Zilla](http://dzil.org/) to author releases.  You
will also need the `Dist::Zilla::PluginBundle::BioPerl` installed as
well as its dependencies.  Then, you can run the following commands:

    dzil test
    dzil install

# The Directory Structure

The bioperl-live repository structure is organized as follows:

* `lib/` - BioPerl modules

* `examples/` - Scripts demonstrating the many uses of BioPerl

* `scripts/` - Useful production-quality scripts with POD documentation

* `t/` - Perl built-in tests, tests are divided into subdirectories
  based on the specific classes being tested

* `t/data/` - Data files used for the tests, provides good example data

* `travis_scripts/` - script to customize Travis

## `Bio::` namespace summary

The BioPerl project is split over multiple Perl module distributions.
The BioPerl distribution is the BioPerl core distribution, including a
selection of modules and namespaces but not all.  For example, the
entire Bio::Biblio is not included in the BioPerl distribution.
Similarly, while many Bio::SearchIO modules in the BioPerl
distribution, there also Bio::SearchIO modules in other distributions
such as Bio-SearchIO-blastxml.

This section describes most of the Bio:: namespaces developed by the
BioPerl project, including those which are not part of the BioPerl
distribution.  For example, the Bio::Biblio and Bio::Assembly are
documented here but are not part of the BioPerl distribution.

* `Bio::Seq` is for *Sequences* (protein and DNA).
    * `Bio::PrimarySeq` is a plain sequence (sequence data +
      identifiers)
    * `Bio::Seq` is a fancier `PrimarySeq`, in that it has annotation
      (via `Bio::Annotation::Collection`) and sequence features (via
      `Bio::SeqFeatureI` objects, attached via `Bio::FeatureHolderI`).
    * `Bio::Seq::RichSeq` is all of the above, plus it has slots for
      extra information specific to GenBank/EMBL/SwissProt files.
    * `Bio::Seq::LargeSeq` is for sequences which are too big for
      fitting into memory.

* `Bio::SeqIO` is for *reading and writing Sequences*.  It is a front
  end module for separate driver modules supporting the different
  sequence formats.

* `Bio::SeqFeature` represent start/stop/strand-based localised
  annotations (features) of sequences
    * `Bio::SeqFeature::Generic` is basic catchall
    * `Bio::SeqFeature::Similarity` a similarity sequence feature
    * `Bio::SeqFeature::FeaturePair` a sequence feature which is
      pairwise such as query/hit pairs

* `Bio::SearchIO` is for reading and writing pairwise alignment
  reports, like BLAST or FASTA.

* `Bio::Search` is where the alignment objects for `SearchIO` are
  defined
    * `Bio::Search::Result::GenericResult` is the result object (a
      blast query is a `Result` object)
    * `Bio::Search::Hit::GenericHit` is the `Hit` object (a query will
      have 0 to many hits in a database)
    * `Bio::Search::HSP::GenericHSP` is the High-scoring Segment Pair
      object defining the alignment(s) of the query and hit.

* `Bio::SimpleAlign` is for multiple sequence alignments

* `Bio::AlignIO` is for reading and writing multiple sequence
  alignment formats

* `Bio::Assembly` provides the start of an infrastructure for assemblies and
  `Bio::Assembly::IO` *IO converters* for them

* `Bio::DB` is the namespace for database query classes
    * `Bio::DB::GenBank/GenPept` are two modules which query NCBI
      entrez for sequences.
    * `Bio::DB::SwissProt/EMBL` query various EMBL and SwissProt
      repositories for a sequences.
    * `Bio::DB::GFF` is Lincoln Stein's fast, lightweight feature and
      sequence database which is the backend to his
      [GBrowse](www.gmod.org) system.
    * `Bio::DB::Flat` is a fast implementation of the OBDA flat-file
      indexing system (cross-language and cross-platform supported by
      O|B|F projects see http://obda.open-bio.org).
    * `Bio::DB::BioFetch/DBFetch` for OBDA, Web (HTTP) access to
      remote databases.
    * `Bio::DB::InMemoryCache/FileCache` (fast local caching of
      sequences from remote dbs to speed up your access).
    * `Bio::DB::Registry` interface to the OBDA specification for
      remote data sources.
    * `Bio::DB::Biblio` for access to remote bibliographic databases.
    * `Bio::DB::EUtilities` is the initial set of modules used for
      generic queried using NCBI's eUtils.

* `Bio::Annotation` collection of annotation objects (comments,
  DBlinks, References, and misc key/value pairs)

* `Bio::Coordinate`** is a system for mapping between different
  coordinate systems such as DNA to protein or between assemblies

* `Bio::Index` is for locally indexed flatfiles with BerkeleyDB

* `Bio::Tools` contains many *miscellaneous parsers and functions* for different
  bioinformatics needs such as:
    * Gene prediction parser (Genscan, MZEF, Grail, Genemark)
    * Annotation format (GFF)
    * Enumerate codon tables and valid sequences symbols (CodonTable,
      IUPAC)
    * Phylogenetic program parsing (PAML, Molphy, Phylip)

* `Bio::Map` represents genetic and physical map representations

* `Bio::Structure` parse and represent protein structure data

* `Bio::TreeIO` is for reading and writing Tree formats

* `Bio::Tree` is the namespace for all associated Tree classes
    * `Bio::Tree::Tree` is the basic tree object
    * `Bio::Tree::Node` are the nodes which make up the tree
    * `Bio::Tree::Statistics` is for computing statistics for a tree
    * `Bio::Tree::TreeFunctionsI` is where specific tree functions are
      implemented (like `is_monophyletic` and `lca`)

* `Bio::Biblio` is where bibliographic data and database access
  objects are kept

* `Bio::Variation` represent sequences with mutations and variations
  applied so one can compare and represent wild-type and mutation
  versions of a sequence.

* `Bio::Root` are basic objects for the internals of BioPerl


# Releases

BioPerl currently uses a [semantic versioning](https://semver.org/)
scheme for version numbers.  Basically, a version has three numbers in
the form `MAJOR.MINOR.PATH`, each of which changes when:

1. `MAJOR` --- incompatible API changes,
2. `MINOR` --- new functionality in a backwards-compatible manner,
3. `PATCH` --- backwards-compatible bug fixes.

## 1.7 releases

Before 1.7 release, the BioPerl project had a single distribution with
all of BioPerl modules.  During the 1.7 release series, subsets of the
modules were extracted into separate distribution.

## Pre 1.7 releases

From version 1.0 until 1.6, even numbers (e.g. version 1.4) indicated
stable releases.  Stable releases were well tested and recommended for
most uses.  Odd numbers (e.g. version 1.3) were development releases
which one would only use if interested in the latest features.  The
final number (e.g. in `1.2.1`) is the point or patch release. The
higher the number the more bug fixes has been incorporated. In theory
you can upgrade from one point or patch release to the next with no
changes to your own code (for production cases, obviously check things
out carefully before you switch over).
