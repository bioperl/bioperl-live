#!/usr/bin/perl
use strict;
use LocalConfig;
use File::Basename qw(basename dirname);
my $file;
($file = basename($0)) =~ s/\.PL$/\.pod/;
open OUT,">$file" or die "Can't create $file: $!";

print OUT<<'!NO!GROK!';

=head1 NAME

Bioscripts

=head1 SYNOPSIS

A list of the scripts in the Bioperl package

=head1 DESCRIPTION

These scripts have been contributed by the developers and users of
Bioperl. They are organized into directories roughly mirroring those
in the Bioperl Bio/ directory. There are 2 directories for these 
scripts, scripts/ and examples/. The scripts in scripts/ are 
production quality scripts, the scripts in examples/ are useful
examples of Bioperl code.

Please contact bioperl-l@bioperl.org if you are interested in 
contributing your own script.

=head1 PRODUCTION SCRIPTS

=head2 scripts/Bio-DB-GFF/generate_histogram.PLS

Create a GFF-formatted histogram of the density of the indicated set
of feature types.

=head2 scripts/Bio-DB-GFF/process_gadfly.PLS

Transforms Gadfly GFF files into correct format.

=head2 scripts/Bio-DB-GFF/process_ncbi_human.PLS

Trnasforms NCBI's chromosome annotations into correct format.

=head2 scripts/Bio-DB-GFF/process_sgd.PLS

Transform SGD format annotations into GFF format.

=head2 scripts/Bio-DB-GFF/process_wormbase.PLS

Transforms Wormbase's GFF files into correct format. Requires Ace.pm.

=head2 scripts/DB/flanks.PLS

Fetch a sequence, find the sequences flanking a variant or SNP in the
sequence given its position.

=head2 scripts/DB/biofetch_genbank_proxy.PLS

A CGI scripts that queries NCBI's eutils to provide database access
according to the BioFetch protocol. Requires Cache::FileCache.

=head2 scripts/DB/biogetseq.PLS

Sequence retrieval the OBDA registry.

=head2 scripts/graphics/feature_draw.PLS

Script that accepts files in GFF or tab-delimited format and creates
corresponding PNG image files. See L<Bio::Graphics> and
L<Bio::Graphics::FeatureFile> for more information.

=head2 scripts/index/bpfetch.PLS

Fetch sequences from local indexed database or over the network and
reformat using Bio::Index* and Bio::DB*.

=head2 scripts/index/bpindex.PLS

Indexes local databases, partners with bpfetch.pl.

=head2 scripts/seq/extract_feature_seq.PLS

Extract the sequence for a specified feature type.

=head2 scripts/taxa/taxid4species.PLS

Retrieve the NCBI Tada ID for a given species.

=head2 scripts/tree/blast2tree.PLS

Builds a phylogenetic tree based on a sequence search (Fasta,
BLAST, HMMER).

=head2 scripts/utilities/remote_blast.PLS

This script executes a remote Blast search using RemoteBlast. See 
L<Bio::Tools::Run::RemoteBlast> for more information.

=head2 scripts/utilities/bp_mrtrans.PLS

Perl implementation of Bill Pearson's mrtrans to project protein 
alignment back into cDNA coordinates.

=head2 scripts/utilities/bp_nrdb.PLS

Make a non-redundant database based on sequence, not id.
Requires Digest::MD5.

=head2 scripts/utilities/bp_sreformat.PLS

Perl implementation of Sean Eddy's sreformat, a sequence and 
alignment converter.

=head2 scripts/utilities/pairwise_kaks.PLS

Takes cDNA sequences as input, aligns them as proteins, 
projects the alignment back into cDNA and estimates 
the Ka (non-synonymous) and Ks (synonymous) substitutions.

=head2 scripts/utilities/seq_length.PLS

Reports the total number of residues and total number of individual 
sequences in a specified sequence database file.

=head1 EXAMPLE SCRIPTS

=head2 examples/align/align_on_codons.pl

Aligns nucleotide sequences based on codons in a specified reading frame.

=head2 examples/align/aligntutorial.pl

Examples using EMBOSS, pSW, Clustalw, TCoffee, and Blast to align sequences.

=head2 examples/align/clustalw.pl

A demonstration of the various uses of Alignment::Clustalw. See 
L<Bio::Tools::Run::Alignment::Clustalw> for more.

=head2 examples/align/simplealign.pl

A script that demonstrates some uses of AlignIO. Please see
L<Bio::AlignIO> for more information.

=head2 examples/biblio/biblio.pl

A script that shows how to query bibliographic databases, such as
Medline, using ids, keywords, and other fields. See L<Bio::Biblio> for 
details.

=head2 examples/biblio/biblio_soap.pl

Connect to and test a SOAP server using a Bio::Biblio object.

=head2 examples/bioperl.pl

A Bioperl shell!

=head2 examples/cluster/dbsnp.pl

How to parse a dbsnp XML file. See L<Bio::ClusterIO> for details.

=head2 examples/contributed/nmrpdb_parse.pl

Extracts individual conformers from an NMR-derived PDB file.

=head2 examples/contributed/prosite2perl.pl

Convert Prosite motifs to Perl regular expressions.

=head2 examples/contributed/rebase2list.pl

Script to convert rebase file to format compatible with 
Bio::Tools::RestrictionEnzyme.

=head2 examples/das/das_server.pl

Sets up a minimal DAS annotation server, requires Apache::DBI and
Bio::DB::GFF. See L<Bio::DB::GFF> for details.

=head2 examples/db/dbfetch

Creates a Web page to query a local SRS server and fetch sequences.

=head2 examples/db/est_tissue_query.pl

Fetch EST sequences from local files or Genbank filtered by tissue
using Bio::DB* or Bio::Index*.

=head2 examples/db/gb2features.pl

Shows how to extract all the features from a Genbank file. See 
L<Bio::Seq> for more information on features.

=head2 examples/db/getGenBank.pl

Retrieving Genbank entries over the Web using DB::GenBank. See
L<Bio::DB::GenBank> for more information.

=head2 examples/db/get_seqs.pl

Fetches and formats sequences from GenBank, EMBL, or SwissProt over
the network using Bio::DB*.

=head2 examples/db/gff/*

Scripts that reformat sequence to GFF and load GFF format files into
an indexed database - see L<Bio::DB::GFF> for more information.

=head2 examples/db/rfetch.pl

A script that uses Bio::DB::Registry to retrieve sequences from EMBL,
reformat them, and print them. See L<Bio::DB::Registry>.

=head2 examples/db/use_registry.pl

Script that shows how to use Bio::DB::Registry, part of Bioperl's
integration with OBDA, the Open Bio Database Access registry
scheme. See L<Bio::DB::Registry> for more information.

=head2 examples/exceptions/*

Scripts that demonstrate how to throw and catch Error.pm objects.

=head2 examples/generate_random_seq.pl

Writes random RNA, DNA, or protein sequence of given length.

=head2 examples/graphics/render_sequence.pl

This scripts fetchs a sequence from a remote database, extracts its
features (CDS, gene, tRNA), and creates a graphic representation of
the sequence in PNG or GIF format. See L<Bio::DB::BioFetch>
and L<Bio::Graphics>.

=head2 examples/liveseq/change_gene.pl

A script showing how to use LiveSeq::Mutator and
LiveSeq::Mutation. Please see L<Bio::LiveSeq::Mutator> and 
L<Bio::LiveSeq::Mutation> for more information.

=head2 examples/longorf.pl

A script that finds the longest ORF in one or more nucleotide sequences.

=head2 examples/make_mrna_protein.pl

Translate a cDNA or ORF to protein using Bio::Seq's translate() method.

=head2 examples/make_primers.pl

Design PCR primers given a sequence and the positions of the start and
stop codons in the sequence's ORF.

=head2 examples/rev_and_trans.pl

Examples using Bio::Seq.pm for reversing and translating
sequences. See L<Bio::Seq> for more information.

=head2 examples/revcom_dir.pl

Eeturn reverse complement sequences of all sequences in the current
directory and save them in the same directory.

=head2 examples/root_object/*

Scripts that demonstrate uses of the Bio::Root modules.

=head2 examples/root_object/vector/vector.pl

Script to test Bio::Root::Vector.

=head2 examples/searchio/*

A number of scripts illustrating the use of Bio::SearchIO for parsing
Blast and PSI-Blast results. See L<Bio::SearchIO> for more information.

=head2 examples/searchio/waba2gff.pl

Convert raw WABA output to one type of GFF.

=head2 examples/seq/*

Example code for working with multiple sequence files, including
formatting, printing, and filtering based on length or description or ID.

=head2 examples/seqstats/aacomp.pl

Calculate amino acid composition of a protein using Tools::CodonTable
and Tools::IUPAC. See L<Bio::Tools::IUPAC> and
L<Bio::Tools::CodonTable> for more information.

=head2 examples/seqstats/chaos_plot.pl

Produce a PNG or JPEG chaos plot given a DNA sequence using GD.pm.

=head2 examples/seqstats/gccalc.pl

Calculate %GC given a DNA sequence using Tools::SeqStats. See
L<Bio::Tools::SeqStats> for more information.

=head2 examples/seqstats/oligo_count.pl

Calculates oligomer frequencies given an oligomer length and a sequence.

=head2 examples/structure/struct_example*

Scripts that show how to examine details of the 3D structure of a
protein by parsing a PDB file. See L<Bio::Structure::IO> for more information.

=head2 examples/subsequence.cgi

CGI script to fetch a sequence from Genbank and extract a subsequence
using DB::GenBank. See L<Bio::DB::GenBank>.

=head2 examples/tk/gsequence

Create a Protein Sequence Control Panel GUI with Gtk.

=head2 examples/tk/hitdisplay.pl

Create a GUI for displaying Blast results using Tk::HitDisplay.
Please see L<Bio::Tk::HitDisplay> for more information.

=head2 examples/tools/gb_to_gff.pl

Extracts top-level sequence features from Genbank-formatted sequence
files using Tools::GFF. See L<Bio::Tools::GFF>.

=head2 examples/tools/gff2ps.pl

Takes an input file in GFF format and draws its genes and features as
Postscript using Tools::GFF. See L<Bio::Tools::GFF>.

=head2 examples/tools/parse_codeml.pl

Script that parses output from codeml, one of the PAML programs. See
L<Bio::Tools::Phylo::PAML>.

=head2 examples/tools/psw.pl

Example code for using the XS extensions for comparing proteins 
using Smith-Waterman.

=head2 examples/tools/restriction.pl

Example code for using the RestrictionEnzyme module. See
L<Bio::Tools::RestrictionEnzyme> for more information.

=head2 examples/tools/run_genscan.pl

Run GENSCAN on multiple sequences and create output sequence files
using Tools::Genscan. Please see L<Bio::Tools::Genscan> for more information.

=head2 examples/tools/seq_pattern.pl

A script that shows how to use sequences as regular expressions using
Tools::SeqPattern. Please see L<Bio::Tools::SeqPattern> for more information.

=head2 examples/tools/state-machine.pl

A demonstration of how to create a state machine using
StateMachine::AbstractStateMachine. Please see
L<Bio::Tools::StateMachine::AbstractStateMachine> for more information.

=head2 examples/tools/test-genscan.pl

Script for testing and demonstrating Bio::Tools::Genscan.

=head2 examples/tree/paup2phylip.pl

Convert a PAUP tree block to Phylip format.

=cut
!NO!GROK!
