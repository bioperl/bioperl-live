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

A list of the scripts in Bioperl

=head1 DESCRIPTION

These scripts have been contributed by the developers and users of
Bioperl. They are organized into directories roughly mirroring those
in the Bioperl Bio directory. Please contact bioperl-l@bioperl.org if
you are interested in contributing your own script.

=head2 scripts/align/align_on_codons.pl

Aligns nucleotide sequences based on codons in a specified reading frame.

=head2 scripts/align/aligntutorial.pl

Examples using EMBOSS, pSW, Clustalw, TCoffee, and Blast to align sequences.

=head2 scripts/align/clustalw.pl

A demonstration of the various uses of Alignment::Clustalw. See 
L<Bio::Tools::Run::Alignment::Clustalw> for more.

=head2 scripts/align/simplealign.pl

A script that demonstrates some uses of AlignIO. Please see
L<Bio::AlignIO> for more information.

=head2 scripts/biblio/biblio.pl

A script that shows how to query bibliographic databases, such as
Medline, using ids, keywords, and other fields. See L<Bio::Biblio> for 
details.

=head2 scripts/biblio/biblio_soap.pl

Connect to and test a SOAP server using a Bio::Biblio object.

=head2 scripts/bioperl.pl

A Bioperl shell!

=head2 scripts/cluster/dbsnp.pl

How to parse a dbsnp XML file. See L<Bio::ClusterIO> for details.

=head2 scripts/contributed/nmrpdb_parse.pl

Extracts individual conformers from an NMR-derived PDB file.

=head2 scripts/contributed/prosite2perl.pl

Convert Prosite motifs to Perl regular expressions.

=head2 scripts/contributed/rebase2list.pl

Script to convert rebase file to format compatible with 
Bio::Tools::RestrictionEnzyme.

=head2 scripts/das/das_server.pl

Sets up a minimal DAS annotation server, requires Apache::DBI and
Bio::DB::GFF. See L<Bio::DB::GFF> for details.

=head2 scripts/db/dbfetch

Creates a Web page to query a local SRS server and fetch sequences.

=head2 scripts/db/est_tissue_query.pl

Fetch EST sequences from local files or Genbank filtered by tissue
using Bio::DB* or Bio::Index*.

=head2 scripts/db/flanks.pl

Fetch a sequence, find the sequences flanking a variant or SNP in the
sequence given its position.

=head2 scripts/db/gb2features.pl

Shows how to extract all the features from a Genbank file. See 
L<Bio::Seq> for more information on features.

=head2 scripts/db/getGenBank.pl

Retrieving Genbank entries over the Web using DB::GenBank. See
L<Bio::DB::GenBank> for more information.

=head2 scripts/db/get_seqs.pl

Fetches and formats sequences from GenBank, EMBL, or SwissProt over
the network using Bio::DB*.

=head2 scripts/db/gff/*

Scripts that reformat sequence to GFF and load GFF format files into
an indexed database - see L<Bio::DB::GFF> for more information.

=head2 scripts/db/rfetch.pl

A script that uses Bio::DB::Registry to retrieve sequences from EMBL,
reformat them, and print them. See L<Bio::DB::Registry>.

=head2 scripts/db/use_registry.pl

Script that shows how to use Bio::DB::Registry, part of Bioperl's
integration with OBDA, the Open Bio Database Access registry
scheme. See L<Bio::DB::Registry> for more information.

=head2 scripts/exceptions/*

Scripts that demonstrate how to throw and catch Error.pm objects.

=head2 scripts/generate_random_seq.pl

Writes random RNA, DNA, or protein sequence of given length.

=head2 scripts/graphics/feature_draw.pl

Script that accepts files in GFF or tab-delimited format and creates
corresponding PNG image files. See L<Bio::Graphics> and
L<Bio::Graphics::FeatureFile> for more information.

=head2 scripts/graphics/render_sequence.pl

This scripts fetchs a sequence from a remote database, extracts its
features (CDS, gene, tRNA), and creates a graphic representation of
the sequence in PNG or GIF format. See L<Bio::DB::BioFetch>
and L<Bio::Graphics>.

=head2 scripts/index/bpfetch.pl

Fetch sequences from local indexed database or over the network and
reformat using Bio::Index* and Bio::DB*.

=head2 scripts/index/bpindex.pl

Indexes local databases, partners with bpfetch.pl.

=head2 scripts/liveseq/change_gene.pl

A script showing how to use LiveSeq::Mutator and
LiveSeq::Mutation. Please see L<Bio::LiveSeq::Mutator> and 
L<Bio::LiveSeq::Mutation> for more information.

=head2 scripts/longorf.pl

A script that finds the longest ORF in one or more nucleotide sequences.

=head2 scripts/make_mrna_protein.pl

Translate a cDNA or ORF to protein using Bio::Seq's translate() method.

=head2 scripts/make_primers.pl

Design PCR primers given a sequence and the positions of the start and
stop codons in the sequence's ORF.

=head2 scripts/rev_and_trans.pl

Examples using Bio::Seq.pm for reversing and translating
sequences. See L<Bio::Seq> for more information.

=head2 scripts/revcom_dir.pl

Eeturn reverse complement sequences of all sequences in the current
directory and save them in the same directory.

=head2 scripts/root_object/*

Scripts that demonstrate uses of the Bio::Root modules.

=head2 scripts/root_object/vector/vector.pl

Script to test Bio::Root::Vector.

=head2 scripts/searchio/*

A number of scripts illustrating the use of Bio::SearchIO for parsing
Blast and PSI-Blast results. See L<Bio::SearchIO> for more information.

=head2 scripts/searchio/waba2gff.pl

Convert raw WABA output to one type of GFF.

=head2 scripts/seq/*

Example code for working with multiple sequence files, including
formatting, printing, and filtering based on length or description or ID.

=head2 scripts/seqstats/aacomp.pl

Calculate amino acid composition of a protein using Tools::CodonTable
and Tools::IUPAC. See L<Bio::Tools::IUPAC> and
L<Bio::Tools::CodonTable> for more information.

=head2 scripts/seqstats/chaos_plot.pl

Produce a PNG or JPEG chaos plot given a DNA sequence using GD.pm.

=head2 scripts/seqstats/gccalc.pl

Calculate %GC given a DNA sequence using Tools::SeqStats. See
L<Bio::Tools::SeqStats> for more information.

=head2 scripts/seqstats/oligo_count.pl

Calculates oligomer frequencies given an oligomer length and a sequence.

=head2 scripts/structure/struct_example*

Scripts that show how to examine details of the 3D structure of a
protein by parsing a PDB file. See L<Bio::Structure::IO> for more information.

=head2 scripts/subsequence.cgi

CGI script to fetch a sequence from Genbank and extract a subsequence
using DB::GenBank. See L<Bio::DB::GenBank>.

=head2 scripts/tk/gsequence

Create a Protein Sequence Control Panel GUI with Gtk.

=head2 scripts/tk/hitdisplay.pl

Create a GUI for displaying Blast results using Tk::HitDisplay.
Please see L<Bio::Tk::HitDisplay> for more information.

=head2 scripts/tools/gb_to_gff.pl

Extracts top-level sequence features from Genbank-formatted sequence
files using Tools::GFF. See L<Bio::Tools::GFF>.

=head2 scripts/tools/gff2ps.pl

Takes an input file in GFF format and draws its genes and features as
Postscript using Tools::GFF. See L<Bio::Tools::GFF>.

=head2 scripts/tools/parse_codeml.pl

Script that parses output from codeml, one of the PAML programs. See
L<Bio::Tools::Phylo::PAML>.

=head2 scripts/tools/psw.pl

Example code for using the XS extensions for comparing proteins 
using Smith-Waterman.

=head2 scripts/tools/remote_blast.pl

This script executes remote Blast using RemoteBlast. See 
L<Bio::Tools::Run::RemoteBlast> for more information.

=head2 scripts/tools/restriction.pl

Example code for using the RestrictionEnzyme module. See
L<Bio::Tools::RestrictionEnzyme> for more information.

=head2 scripts/tools/run_genscan.pl

Run GENSCAN on multiple sequences and create output sequence files
using Tools::Genscan. Please see L<Bio::Tools::Genscan> for more information.

=head2 scripts/tools/seq_pattern.pl

A script that shows how to use sequences as regular expressions using
Tools::SeqPattern. Please see L<Bio::Tools::SeqPattern> for more information.

=head2 scripts/tools/state-machine.pl

A demonstration of how to create a state machine using
StateMachine::AbstractStateMachine. Please see
L<Bio::Tools::StateMachine::AbstractStateMachine> for more information.

=head2 scripts/tools/test-genscan.pl

Script for testing and demonstrating Bio::Tools::Genscan.

=head2 scripts/tree/paup2phylip.pl

Convert a PAUP tree block to Phylip format.

=cut
!NO!GROK!
