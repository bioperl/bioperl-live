# Getting Started

Please see the the `INSTALL` or `INSTALL.WIN` documents for installation
instructions.

# About BioPerl

BioPerl is a package of public domain Perl tools for computational molecular
biology.

Our website (http://bioperl.org/) provides an online resource of modules,
scripts, and web links for developers of Perl-based software for life science
research.

# Contact info

BioPerl mailing list: bioperl-l@bioperl.org

There's quite a variety of tools available in BioPerl, and more are added all
the time. If the tool you're looking for isn't described in the documentation
please write us, it could be undocumented or in process.

* Project website : http://bioperl.org/

* Bug reports : https://github.com/bioperl/bioperl-live/issues

Please send us bugs, in particular about documentation which you think is
unclear or problems in installation. We are also very interested in functions
which don't work the way you think they do!

# The directory structure

The BioPerl directory structure is organized as follows:

* **`Bio/`** - BioPerl modules

* **`doc/`** - Documentation utilities

* **`examples/`** - Scripts demonstrating the many uses of BioPerl

* **`ide/`** - files for developing BioPerl using an IDE

* **`maintenance/`** - BioPerl housekeeping scripts

* **`models/`** - DIA drawing program generated OO UML for BioPerl classes
  (these are quite out-of-date)

* **`scripts/`** - Useful production-quality scripts with POD documentation

* **`t/`** - Perl built-in tests, tests are divided into subdirectories
  based on the specific classes being tested

* **`t/data/`** - Data files used for the tests, provides good example data

# Documentation

For documentation on BioPerl see the **HOWTO** documents and tutorials online at
http://bioperl.org.

Useful documentation in the form of example code can also be found in the
**`examples/`** and **`scripts/`** directories. The current collection includes
scripts that run BLAST, index flat files, parse PDB structure files, make
primers, retrieve ESTs based on tissue, align protein to nucleotide sequence,
run GENSCAN on multiple sequences, and much more! See `bioscripts.pod` for a
complete listing.

Individual `*.pm` modules have their own embedded POD documentation as well. A
complete set of hyperlinked POD, or module, documentation is available at
http://www.bioperl.org/.

Remember that '`perldoc`' is your friend. You can use it to read any file
containing POD formatted documentation without needing any type of translator
(e.g. '`perldoc Bio::SeqIO`').

If you used the Build.PL installation, and depending on your platform, you may
have documentation installed as man pages, which can be accessed in the usual
way.

# Releases

BioPerl releases are always available from the website at
http://www.bioperl.org/DIST or in CPAN. The latest code can be found at
https://github.com/bioperl.

* BioPerl currently uses a sematic numbering scheme to indicate stable release
  series vs. development release series. A release number is a three digit
  number like `1.2.0`.
  * The *first digit indicates the major release*, the idea being that all the
    API calls in a major release are reasonably consistent.
  * The *second number is the release series*. This is probably the most
    important number, and represents added functionality that is
    backwards-compatible.
  * The *third number is the point or patch release* and represents mainly bug
    fixes or additional code that doesn't add significant functionality to the
    code base.

* From the **1.0 release until the 1.6 release**, even numbers (`1.0`, `1.2`, etc)
  indicated stable releases. Stable releases were well tested and recommended
  for most uses. Odd numbers (`1.1`, `1.3`, etc) were development releases which one
  would only use if one were interested in the latest and greatest features. The
  final number (e.g. `1.2.0`, `1.2.1`) is the bug fix release. The higher the number
  the more bug fixes has been incorporated. In theory you can upgrade from one
  bug fix release to the next with no changes to your own code (for production
  cases, obviously check things out carefully before you switch over).

* The upcoming **1.7 release** will be the last release series to utilize the
  alternating 'stable'/'developer' convention. Starting immediately after the
  final 1.6 branch, we will start splitting BioPerl into several smaller
  easier-to-manage distributions. These will have independent versions, all
  likely starting with v1.7.0. **We do not anticipate major API changes in the
  1.7.x release series*, merely that the code will be restructured in a way to
  make maintenance more feasible. We anticipate retaining semantic versioning
  until the **v2.x** release.

# Caveats and warnings

When you run the tests ("`./Build test`") some tests may issue warnings messages
or even fail. Sometimes this is because we didn't have anyone to test the test
system on the combination of your operating system, version of perl, and
associated libraries and other modules. Because BioPerl depends on several
outside libraries we may not be able to test every single combination so if
there are warnings you may find that the package is still perfectly useful.

If you install the bioperl-run system and run tests when you don't have the
program installed you'll get messages like '`program XXX not found, skipping
tests`'. That's okay, BioPerl is doing what it is supposed to do. If you wanted
to run the program you'd need to install it first.

Not all scripts in the `examples/` directory are correct and up-to-date. We need
volunteers to help maintain these so if you find they do not submit a bug report
to https://github.com/bioperl/bioperl-live/issues and consider helping out in
their maintenance.

If you are confused about what modules are appropriate when you try and solve a
particular issue in bioinformatics we urge you to look at HOWTO documents first.

# A simple module summary

Here is a quick summary of many of the useful modules and how the toolkit is
laid out:

All modules are in the **`Bio/`** namespace,

* **`Perl`** is for *new users*, and gives a functional interface to the main
  parts of the package.

* **`Seq`** is for *Sequences* (protein and DNA).
    * `Bio::PrimarySeq` is a plain sequence (sequence data + identifiers)
    * `Bio::Seq` is a fancier `PrimarySeq`, in that it has annotation (via
    `Bio::Annotation::Collection`) and sequence features (via `Bio::SeqFeatureI` objects, attached via
    `Bio::FeatureHolderI`).
    * `Bio::Seq::RichSeq` is all of the above, plus it has slots for extra information specific to GenBank/EMBL/SwissProt files.
    * `Bio::Seq::LargeSeq` is for sequences which are too big for
    fitting into memory.

* **`SeqIO`** is for *reading and writing Sequences*. It is a front end module
  for separate driver modules supporting the different sequence formats

* **`SeqFeature`** represent *start/stop/strand-based localized annotations (features) of sequences*
    * **`Bio::SeqFeature::Generic`** is basic catchall
    * **`Bio::SeqFeature::Similarity`** a similarity sequence feature
    * **`Bio::SeqFeature::FeaturePair`** a sequence feature which is pairwise
    such as query/hit pairs

* **`SearchIO`** is for *reading and writing pairwise alignment reports*, like
  BLAST or FASTA

* **`Search`** is where the *alignment objects for `SearchIO` are defined*
    * **`Bio::Search::Result::GenericResult`** is the result object (a blast
    query is a `Result` object)
    * **`Bio::Search::Hit::GenericHit`** is the `Hit` object (a query will have
    0 to many hits in a database)
    * **`Bio::Search::HSP::GenericHSP`** is the High-scoring Segment Pair
    object defining the alignment(s) of the query and hit.

* **`SimpleAlign`** is for *multiple sequence alignments*

* **`AlignIO`** is for *reading and writing multiple sequence alignment
  formats*

* **`Assembly`** provides the start of an *infrastructure for assemblies* and
  **`Assembly::IO`** *IO converters* for them

* **`DB`** is the namespace for *all the database query classes*
    * **`Bio::DB::GenBank/GenPept`** are two modules which query NCBI entrez for
      sequences
    * **`Bio::DB::SwissProt/EMBL`** query various EMBL and SwissProt
      repositories for a sequences
    * **`Bio::DB::GFF`** is Lincoln Stein's fast, lightweight feature and
      sequence database which is the backend to his GBrowse system (see
      www.gmod.org)
    * **`Bio::DB::Flat`** is a fast implementation of the OBDA flat-file
      indexing system (cross-language and cross-platform supported by O|B|F
      projects see http://obda.open-bio.org).
    * **`Bio::DB::BioFetch/DBFetch`** for OBDA, Web (HTTP) access to remote
      databases.
    * **`Bio::DB::InMemoryCache/FileCache`** (fast local caching of sequences
      from remote dbs to speed up your access).
    * **`Bio::DB::Registry`** interface to the OBDA specification for remote
      data sources
    * **`Bio::DB::Biblio`** for access to remote bibliographic databases.
    * **`Bio::DB::EUtilities`** is the initial set of modules used for generic
      queried using NCBI's eUtils.

* **`Annotation`** collection of *annotation objects* (comments, DBlinks,
  References, and misc key/value pairs)

* **`Coordinate`** is a system for *mapping between different coordinate systems*
  such as DNA to protein or between assemblies

* **`Index`** is for *locally indexed flatfiles* with BerkeleyDB

* **`Tools`** contains many *miscellaneous parsers and functions* for different
  bioinformatics needs
    * Gene prediction parser (Genscan, MZEF, Grail, Genemark)
    * Annotation format (GFF)
    * Enumerate codon tables and valid sequences symbols (CodonTable,
    IUPAC)
    * Phylogenetic program parsing (PAML, Molphy, Phylip)

* **`Map`** represents *genetic and physical map representations*

* **`Structure`** - parse and represent *protein structure data*

* **`TreeIO`** is for reading and writing *Tree formats*

* **`Tree`** is the namespace for **all associated Tree classes**
    * **`Bio::Tree::Tree`** is the basic tree object
    * **`Bio::Tree::Node`** are the nodes which make up the tree
    * **`Bio::Tree::Statistics`** is for computing statistics for a tree
    * **`Bio::Tree::TreeFunctionsI`** is where specific tree functions are
      implemented (like `is_monophyletic` and `lca`)

* **`Bio::Biblio`** is where *bibliographic data and database access objects*
  are kept

* **`Variation`** represent *sequences with mutations and variations* applied so
  one can compare and represent wild-type and mutation versions of a sequence.

* **`Root`**, basic objects for the *internals of BioPerl*

# Upgrading from an older version

If you have a previously installed version of BioPerl on your system some of
these notes may help you.

* Some modules have been removed because they have been superceded by new
  development efforts. They are documented in the **`DEPRECATED`** file that is
  included in the release.

* Some methods, or the Application Programming Interface (API), have changed or
  been removed. You may find that scripts which worked with BioPerl 1.4 may give
  you warnings or may not work at all (although we have tried very hard to
  minimize this!). Send an email to the list and we'll be happy to give you
  pointers.
