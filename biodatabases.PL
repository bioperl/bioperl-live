#!/usr/bin/perl
use strict;
use LocalConfig;
use File::Basename qw(basename dirname);
my $file;
($file = basename($0)) =~ s/\.PL$/\.pod/;
open OUT,">$file" or die "Can't create $file: $!";

print OUT<<'!NO!GROK!';
#
# Documentation about the backend of Bioperl
# 
# This is a quick overview about how to use databases with Bioperl
#

=head1 NAME

Bioperl databases - how to use databases with Bioperl

=head1 SYNOPSIS

Not really appropiate for a synopsis. Read on...

=head1 DESCRIPTION

This document is designed to let you use Bioperl with databases.
Bioperl can work with a number of sequence formats and allows users
to interconvert easily amongst these formats. Bioperl can also index
flat files in some of these standard formats, allowing very fast retrieval.
Other layers above flat file indexing are provided, allowing users to
retrieve sequences from remote databases via the Web or build and query
in-house relational databases (RDBs).

Different scripts are provided to get you started with the Bioperl
backend, for example:

=over 2

=item scripts/bpfetch.pl 

Fetches sequences from a Database

=item scripts/bpindex.pl

Builds indexes for flat files databases which are easily accessible
by bpfetch.pl. These scripts are discussed in more detail below.

=back

The core of the backend system is found in following modules

=over 2

=item Bio::DB::*

Generic access to databases, whether flat file, web or RDB. At the
moment, this provides random access retrieval, on the basis of ids or
accession numbers, but does not provide the ability to loop over the
entire database, nor does it provide any complex querying ability. For
these advanced features see L<USING BIOPERL WITH A RELATIONAL DATABASE (bioperl-db)>.

Bio::DB::BioSeqI is the abstract interface for the databases (hence
the I).  Bio::DB::GenBank and Bio::DB::GenPept are concrete
implementations for network access to the GenBank and GenPept
databases at NCBI, and others, using HTTP as a protocol. See
L<Bio::DB::GenBank>, L<Bio::DB::GenPept>, and L<REMOTE SEQUENCE RETRIEVAL (Bio::DB::*)> for more information.

=item Bio::Index::* and Bio::DB::Fasta

Flat file indexing system, for read-only, flat file distributions. These
provide for format-specific, generic type access, but the underlying
machinery can be customised for any number of different flat file systems.
The indices allow extremely fast retrieval based on unique keys.

The Index modules EMBL and Fasta, as they are designed as sequence databases,
conform to the Bio::DB::BioSeqI interface, meaning they can be used whereever
the Bio::DB::BioSeqI is expected.

Flat file indexing of Fasta files is also provided by Bio::DB::Fasta,
please see L<Bio::DB::Fasta> for more information - this module has some
useful features not contained in Bio::Index::Fasta.

=item Bio::SeqIO::*

Conversion systems for Bio::Seq objects, either to or from sequence
streams. The movement of things into SeqIO prevents the Bio::Seq object
bloating up with format code, and the SeqIO system has the benefit
of being very easy to extend to new formats. Please see L<Bio::SeqIO>
for more information.

=item bioperl-db

The bioperl-db package is installed separately from Bioperl but is
tightly integrated with Bioperl's objects. The bioperl-db RDB contains
a set of tables that map to Bioperl's Seq and SeqFeature
objects, and sequences in familiar formats are easily loaded into
the schema in relational form, and retrieved as Seq objects as well.
The RDB is also designed to contain genetic marker and map information.
Currently only mysql is supported, but other RDBMs are planned. See
L<USING BIOPERL WITH A RELATIONAL DATABASE (bioperl-db)> below.

=back

=head1 REMOTE SEQUENCE RETRIEVAL (Bio::DB::*)

Bioperl offers a number of modules for retrieving sequences over the
network. The available remote databases include GenBank, GenPept, GDB,
EMBL, SwissProt, XEMBL, and remote Ace servers. A typical method is

   get_Seq_by_id($id)

which returns a Seq object, or

   get_Stream_by_id($ref_to_array_of_ids)

which returns a SeqIO object. See L<Bio::DB::GenBank>, L<Bio::DB::GenPept>,
L<Bio::DB::GDB>, L<Bio::DB::EMBL>, L<Bio::DB::SwissProt>, L<Bio::DB::XEMBL>,
and L<Bio::DB::Ace>, and L<Bio::SeqIO> for more information.


=head1 SETTING UP BIOPERL INDICES (Bio::Index::*)

If you want to use Bioperl indicies of Fasta, EMBL/SwissProt .dat files,
SwissPfam, GenBank, or Blast files then the bpfetch.pl and
bpindex.pl scripts are great ways to start off (and also reading the
scripts shows you how to use the Bioperl indexing stuff). bpfetch.pl and
bpindex.pl coordinate using two environment variables

  BIOPERL_INDEX - directory where the indices are kept

  BIOPERL_INDEX_TYPE - type of DBM file to use for the index

The basic way of indexing a database, once BIOPERL_INDEX has been
set up, is to go

  bpindex.pl <index-name> <filenames as full path>

eg, for Fasta files

  bpindex.pl est /nfs/somewhere/fastafiles/est*.fa

Or, for EMBL/Swissprot files

  bpindex.pl -fmt=EMBL swiss /nfs/somewhere/swiss/swissprot.dat

To retrieve sequences from the index go

  bpfetch.pl <index-name>:<id>

eg,

  bpfetch.pl est:AA01234

or

  bpfetch.pl swiss:VAV_HUMAN


bpfetch.pl also has other options to connect to Genbank across the network.


=head1 CHECKLIST

   mkdir /nfs/datadisk/bioperlindex/

or any other directory

   setenv BIOPERL_INDEX /nfs/datadisk/bioperlindex/

in .cshrc or .tcshrc (or B<set> and B<export> in bash and its .bashrc)

go 

   bpindex.pl swissprot /nfs/datadisk/swiss/swissprot.dat

etc. You are now ready to use bpfetch.pl. See L<Bio::Index::Fasta>, 
L<Bio::Index::Genbank>, L<Bio::Index::Blast>, L<Bio::Index::EMBL>,
L<Bio::Index::SwissPfam>, and L<Bio::Index::Swissprot> for more.	

Flat file indexing of Fasta files is also provided by Bio::DB::Fasta,
please see L<Bio::DB::Fasta> for more information - this module provides
some functionality absent from Bio::Index::Fasta.


=head1 USING BIOPERL WITH A RELATIONAL DATABASE (bioperl-db)

First download and install the bioperl-db package
(ftp://bioperl.org/pub/DIST/bioperl-db-0.1.tar.gz). The current version
relies on the prior installation of the free mysql server, from
http://www.mysql.com. The bioperl-db package integrates neatly with
Bioperl's objects and contains tables for sequences, entries, sequence
features, taxonomic information, references, keywords, a wide variety of
genetic markers, and more. With a Seq object and a built-in database
"adaptor" one can load the object's data in one step, as in 

  $adaptor->store($newid,$seqobj)

Similarly, one can query the database with id's or using query objects
and retrieve arrays of Seq objects. For example

  @seqobjs = $adaptor->fetch_by_query( -constraints =>
            ["or", 
                  ["and", # find all human topoisomerases
                         "species=Human",
                         "keywords like *topoisomerase*"],
                  ["and", # and Escherichia gyrases
                         "species like Escherichia*",
                         "keywords like *gyrase*"],
            ]);

See L<Bio::DB::SQL::SeqAdaptor>, L<Bio::DB::SQL::QueryConstraint>, 
L<Bio::DB::SQL::CommentAdaptor>, and L<Bio::DB::SQL::BioQuery> for
more examples.

=head1 BIO::DB::GFF: INDEXING AND RETRIEVAL OF SEQUENCE FEATURE DATA

An alternative view of sequence feature data is provided by the
Bio::DB::GFF module, part of the core package.  This module takes a
set of sequence features stored in GFF (gene-finding format) format
and loads them into a relational database optimized for positional
queries.  This allows a variety of data mining operations, such as
finding sequence features that are within a certain distance of each
other or which overlap. The module can also store large contiguous
sequences and extract subsequences rapidly.

See L<Bio::DB::GFF>, L<Bio::DB::GFF::RelSegment>,
L<Bio::DB::GFF::Feature>, and L<Bio::DB::GFF::Adaptor::dbi>.

!NO!GROK!
