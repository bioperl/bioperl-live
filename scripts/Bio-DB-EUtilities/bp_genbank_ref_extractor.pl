#!/usr/bin/perl
## Copyright (C) 2011 Carnë Draug <carandraug+dev@gmail.com>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>.

use 5.010;                      # Use Perl 5.10
use warnings;                   # Replacement for the -w flag, but lexically scoped
use strict;                     # Enforce some good programming rules
use Getopt::Long;               # Parse program arguments
use Cwd;                        # Determines current working directory
use File::Spec;                 # Perform operation on file names
use Bio::SeqIO;                 # Handler for SeqIO Formats
use Bio::DB::EUtilities;        # Retrieve entries from Entrez

=head1 NAME

bp_genbank_ref_extractor - retrieves all related sequences for a list of searches on Entrez gene

=head1 SYNOPSIS

B<bp_genbank_ref_extractor> [options] [Entrez Gene Queries]

=head1 DESCRIPTION

This script searches on I<Entrez Gene> database and retrieves not only the gene sequence but
also the related transcript and protein sequences.

The gene UIDs of multiple searches are collected before attempting to retrieve them so each gene
will only be analyzed once even if appearing as result on more than one search.

Note that I<by default no sequences are saved> (see options and examples).

=head1 OPTIONS

Several options can be used to fine tune the script behaviour. It is possible to obtain extra
base pairs upstream and downstream of the gene, control the naming of files and genome assembly to use.

See the section bugs for problems when using default values of options.

=over

=cut

=item B<--assembly>

When retrieving the sequence, a specific assemly can be defined. The value expected
is a regex that will be case-insensitive. If it matches more than one assembly, it will
use the first match. It defauls to C<(primary|reference) assembly>.

=cut
my $assembly_regex = '(primary|reference) assembly';

=item B<--debug>

If set, even more output will be printed that may help on debugging. Unlike the messages
from B<--verbose> and B<--very-verbose>, these will not appear on the log file
unless this option is selected. This option also sets B<--very-verbose>.

=cut
my $debug         = 0;

=item B<--downstream>, B<--down>

Specifies the number of extra base pairs to be retrieved downstream of the gene.
This extra base pairs will only affect the gene sequence, not the transcript or proteins.

=cut
my $downstream    = 0;

=item B<--format>

Specifies the format that the sequences will be saved. Defaults to I<genbank> format.
Valid formats are 'genbank' or 'fasta'.

=cut
my $format        = 'genbank';

=item B<--genes>

Specifies the name for gene file. By default, they are not saved. If no value is given
defaults to its UID. Possible values are 'uid', 'name', 'symbol' (the official symbol or
nomenclature).

=cut
my $genes         = '';
sub genes_option_parsing {
  given ($_[1]) {
    when (/^(u)?id$/i)    { $genes = 'uid'; }
    when (/^sym(bol)?$/i) { $genes = 'symbol'; }
    when (/^name$/i)      { $genes = 'name'; }
    ## default is set here, when value is empty
    when ('')             { $genes = 'uid' }
    default               { die "Invalid identifier '$_[1]' for gene files."; }
  }
}

=item B<--limit>

When making a query, limit the result to these first specific results. This is to
prevent the use of specially unspecific queries and a warning will be given if a
query returns more results than the limit. The default value is 200. Note that
this limit is for I<each> search.

=cut
my $limit         = 200;

=item B<--non-coding>, B<--nonon-coding>

Some protein coding genes have transcripts that are non-coding. By default, these sequences are
saved as well. B<--nonon-coding> can be used to ignore those transcripts.

=cut
my $get_noncoding = 1;

=item B<--proteins>

Specifies the name for proteins file. By default, they are not saved. If no value is given
defaults to its accession. Possible values are 'accession', 'description', 'gene' (the corresponding
gene ID) and 'transcript' (the corresponding transcript accesion).

Note that if not using 'accession' is possible for files to be overwritten. It is possible for the same gene
to encode more than one protein or different proteins to have the same description.

=cut
my $proteins      = '';
sub proteins_option_parsing {
  given ($_[1]) {
    when (/^acc(ession)?$/i)      { $proteins = 'accession'; }
    when (/^desc(ription)?$/i)    { $proteins = 'description'; }
    when (/^gene$/i)              { $proteins = 'gene'; }
    when (/^(transcript|mrna)$/i) { $proteins = 'transcript'; }
    ## default is set here, when value is empty
    when ('')                     { $proteins = 'accession' }
    default                       { die "Invalid identifier '$_[1]' for protein files."; }
  }
}

=item B<--pseudo>, B<--nopseudo>

By default, sequences of pseudo genes will be saved. B<--nopseudo> can be used to ignore those genes.

=cut
my $get_pseudo    = 1;

=item B<--save>

Specifies the path for the directory where the sequence and log files will be saved. If the
directory does not exist it will be created altough the path to it must exist. Files on the
directory may be rewritten if necessary. If unspecified, a directory named F<extracted sequences>
on the current directory will be used.

=cut
my $save          = File::Spec->catfile (getcwd, 'extracted sequences');

=item B<--save-data>

This options saves the data (gene UIDs, description, product accessions, etc) to
a file. As an optional value, the file format can be specified. Defaults to CSV.

Currently only CSV is supported.

Saving the data structure as a CSV file, requires the installation of the Text::CSV module.

=cut
my $save_data     = '';
sub save_data_option_parsing {
  given ($_[1]) {
    when (/^csv$/i) { $save_data = 'csv'; require Text::CSV; }
    when ('')       { $save_data = 'csv'; require Text::CSV; } ## Do nothing. If not set, use default
    default         { die "Specified format to save data '$save_data' is not valid."; }
  }
}

=item B<--transcripts>, B<--mrna>

Specifies the name for transcripts file. By default, they are not saved. If no value is given
defaults to its accession. Possible values are 'accession', 'description', 'gene' (the corresponding
gene ID) and 'protein' (the protein the transcript encodes).

Note that if not using 'accession' is possible for files to be overwritten. It is possible for the same gene
to have more than one transcript or different transcripts to have the same description. Also, non-coding
transcripts will create problems if using 'protein'.

=cut
my $transcripts   = '';
sub transcripts_option_parsing {
  given ($_[1]) {
    when (/^acc(ession)?$/i)    { $transcripts = 'accession'; }
    when (/^desc(ription)?$/i)  { $transcripts = 'description'; }
    when (/^gene$/i)            { $transcripts = 'gene'; }
    when (/^protein$/i)         { $transcripts = 'protein'; }
    ## default is set here, when value is empty
    when ('')                   { $transcripts = 'accession' }
    default                     { die "Invalid identifier '$_[1]' for transcript files."; }
  }
}

=item B<--upstream>, B<--up>

Specifies the number of extra base pairs to be extracted upstream of the gene.
This extra base pairs will only affect the gene sequence, not the transcript or proteins.

=cut
my $upstream      = 0;

=item B<--verbose>, B<--v>

If set, program becomes verbose. For an extremely verbose program, use B<--very-verbose> instead.

=cut
my $verbose       = '';

=item B<--very-verbose>, B<--vv>

If set, program becomes extremely verbose. Setting this option, automatically sets B<--verbose> as well.
For help in debugging, consider using B<--debug>

=cut
my $very_verbose  = '';

=back

=head1 EXAMPLES

=over

=item C<bp_genbank_ref_extractor --transcripts=accession '"homo sapiens"[organism] AND H2B'>

Search Entrez gene with the query C<'"homo sapiens"[organism] AND H2B'>, and
save their transcripts sequences. Note that default value of B<--limit> may only extract
some of the hits.

=item C<bp_genbank_ref_extractor --transcripts=accession --proteins=accession --format=fasta '"homo sapiens"[organism] AND H2B' '"homo sapiens"[organism] AND MCPH1'>

Same as first example but also searches for C<'"homo sapiens"[organism] AND MCPH1'>,
proteins sequences, and saves them in the fasta format.

=item C<bp_genbank_ref_extractor --genes --up=100 --down=500 '"homo sapiens"[organism] AND H2B'>

Same search as first example but saves the genomic sequences instead including
100 and 500 bp upstream and downstream.

=item C<bp_genbank_ref_extractor --genes --asembly='Alternate HuRef' '"homo sapiens"[organism] AND H2B'>

Same search as first example but saves genomic sequences and from the Alternate HuRef genome assembly instead.

=item C<bp_genbank_ref_extractor --save-data=CSV '"homo sapiens"[organism] AND H2B'>

Same search as first example but does not save any sequence but saves all the results in a CSV file.

=item C<bp_genbank_ref_extractor --save='search results' --genes=name --upstream=200 downstream=500 --nopseudo --nonnon-coding  --transcripts --proteins  --format=fasta --save-data=CSV '"homo sapiens"[organism] AND H2B' '"homo sapiens"[organism] AND MCPH1'>

Searches on Entrez gene for both C<'"homo sapiens"[organism] AND H2B'> and C<'"homo sapiens"[organism] AND MCPH1'>
and saves the gene sequences of all hits (not passing the default limit and ignoring pseudogenes) plus 200 and 500bp
upstream and downstream of them. It will also save the sequences of all transcripts and proteins of each gene (but
ignoring non-coding transcripts). It will save the sequences in the fasta format, inside a directory C<search results>,
and save the results in a CSV file

=back

=head1 BUGS

If you find any bug, or have a feature request, please report these at L<https://redmine.open-bio.org/projects/bioperl> or
e-mail L<mailto:bioperl-l@lists.open-bio.org>

=over

=item *

When supplying options, it's possible to not supply a value and use their default. However,
when the expected value is a string, the next argument may be confused as value for the
option. For example, when using the following command:

C<bp_genbank_ref_extractor --transcripts 'H2A AND homo sapiens'>

we mean to search for 'H2A AND homo sapiens' saving only the transcripts and using the default
as base for the filename. However, the search terms will be interpreted as the base for the
filenames (but since it's not a valid identifier, it will return an error). To prevent
this, you can either specify the values:

C<bp_genbank_ref_extractor --transcripts 'accession' 'H2A AND homo sapiens'>

C<bp_genbank_ref_extractor --transcripts='accession' 'H2A AND homo sapiens'>

or you can use the double hash to stop processing options. Note that this should only be used
after the last option. All arguments supplied after the double dash will be interpreted as search terms

C<bp_genbank_ref_extractor --transcripts -- 'H2A AND homo sapiens'>

=back

=head1 NOTES ON USAGE

=over

=cut

################################################################################
## Parse options, check and create files and directories needed
################################################################################

GetOptions(
            'assembly:s'          => \$assembly_regex,
            'debug'               => \$debug,
            'down|downstream=i'   => \$downstream,
            'format=s'            => \$format,
            'genes:s'             => \&genes_option_parsing,
            'limit=i'             => \$limit,
            'non-coding!'         => \$get_noncoding,
            'proteins:s'          => \&proteins_option_parsing,
            'pseudo!'             => \$get_pseudo,
            'save=s'              => \$save,
            'save-data:s'         => \&save_data_option_parsing,
            'transcripts|mrna:s'  => \&transcripts_option_parsing,
            'up|upstream=i'       => \$upstream,
            'verbose|v'           => \$verbose,
            'very-verbose|vv'     => \$very_verbose,
          ) or die "Error processing options";
## It is necessary to check success of GetOptions since:
## ''GetOptions returns true to indicate success. It returns false when the function
## detected one or more errors during option parsing. These errors are signalled
## using warn() and can be trapped with $SIG{__WARN__}''

## set verbosity level
my $verbosity;
if ($debug) {
  $verbosity = 9;
} elsif ($very_verbose) {
  $verbosity = 3;
} elsif ($verbose) {
  $verbosity = 2;
} else {
  $verbosity = 1;
}

my $gene_dir  = File::Spec->catfile ($save, 'genes');
my $mrna_dir  = File::Spec->catfile ($save, 'transcripts');
my $prot_dir  = File::Spec->catfile ($save, 'proteins');
check_dir($_) foreach ($save, $gene_dir, $mrna_dir, $prot_dir);
my $log_file  = File::Spec->catfile ($save, 'extractor.log');
open (LOG, ">", $log_file) or die "Couldn't open file $log_file for writing: $!";

## TODO should also say the version
log_it (1, "This is bp_genbank_ref_extractor on ". &get_time);

################################################################################
## Everything is ready. Start accessing the database
################################################################################

my $data = Structure->new;
say "Searching on Entrez gene...";
my @uids;
push (@uids, gb_search ($_)) foreach (@ARGV);
{
  my $start = scalar(@uids);
  clean_array(\@uids);
  my $diff  = $start - scalar(@uids);
  log_it (2, "Entrez gene: removed $diff UIDs from the search results for being repeated.") if $diff > 0;
  log_it (3, "Entrez gene: list of retrieved IDs is: @uids");
}
say "Fetching gene info...";
analyze_entrez_genes ($data, \@uids);

if ($genes)       { say "Fetching gene sequences...";       get_genes($data); }
if ($transcripts) { say "Fetching transcript sequences..."; get_products('transcript', $data); }
if ($proteins)    { say "Fetching protein sequences...";    get_products('protein', $data); }

if ($save_data) { save_structure($data); }
if ($debug) { use Data::Dumper; print Dumper $data; }
exit;

################################################################################
## genbank search subs
################################################################################

sub gb_search {
  log_it (2, "Entrez gene: searching with '$_[0]'");
  my $searcher  = Bio::DB::EUtilities->new(
                                            -eutil   => 'esearch',
                                            -db      => 'gene',
                                            -term    => $_[0],
                                            -retmax  => $limit,
                                            );
  log_it (3, "Entrez gene: query $_[0] translated into '" . $searcher->get_query_translation . "'");
  log_it (2, "Entrez gene: found " . $searcher->get_count . " UIDS");
  if ($searcher->get_count > $limit) {
    my $w_message = "Entrez gene: search returned more ids than the set limit of $limit. Retrieving only the first $limit genes.";
    log_it (2, $w_message);
    warn $w_message;
  }
  return $searcher->get_ids;
}

## we are not using esummary because it doesn't tell us the products of the gene thus forcing us
## to download the sequence and analyze it from that. We are also not using elink because it fails
## too frequently. Also, entrezgene gives currently the most information so it'll be easier to implement
## new stuff at a later time
sub analyze_entrez_genes {
  my $struct  = shift;
  my $uid_ref = shift;  # a reference for the array containing the list of gene UID

  ## TODO may be a good idea to limit this and download only a few sequences rather
  ## than what can possibly be thousands of them.
  my $fetcher = Bio::DB::EUtilities->new(
                                          -eutil   => 'efetch',
                                          -db      => 'gene',
                                          -id      => $uid_ref,
                                          -retmode => 'text',
                                          -rettype => 'asn1',
                                         );
  my $response = $fetcher->get_Response->content;
  open(my $seq_fh, "<", \$response) or die "Could not open sequences string for reading: $!";

  ## TODO when this bug is fixed https://redmine.open-bio.org/issues/3261
  ## this should be fixed to use Bio::SeqIO with format=> 'entrezgene'
  ## then we could use methods to access the data
  use Bio::ASN1::EntrezGene;

  my $parser = Bio::ASN1::EntrezGene->new(
                                          -fh => $seq_fh,
                                          );

  SEQ: while(my $result = $parser->next_seq){
    ## it's possible that when analyzing genes, if a gene has a RefSeq status of secondary, it will point
    ## to another gene, also in the list. To prevent analyzing the same gene twice here, this hash keeps
    ## track of the analyzed gene UIDs (even the ones who are pseudo and may not be on $struct)
    state %analyzed_genes;
    $result = $result->[0] if(ref($result) eq 'ARRAY');
    ## Data::Dumper can be used to look into the structure and find where things are
#      use Data::Dumper;
#      print Dumper ($result);
#      exit;

    my $uid = $result->{'track-info'}->[0]->{'geneid'};
    if ($analyzed_genes{$uid}) {
      log_it (9, "DEBUG: skipping analysis of gene with UID='$uid' since it's already been done.");
      next SEQ;
    }
    $analyzed_genes{$uid} = 1;

    my ($symbol, $name);
    foreach my $p (@{$result->{'properties'}}){
      $p = $p->[0] if(ref($p) eq 'ARRAY');
      next unless ($p->{'label'} && $p->{'label'} eq 'Nomenclature');
      foreach my $pp (@{$p->{'properties'}}){
        $pp     = $pp->[0] if(ref($pp) eq 'ARRAY');
        $name   = $pp->{'text'} if ($pp->{'label'} && $pp->{'label'} eq 'Official Full Name');
        $symbol = $pp->{'text'} if ($pp->{'label'} && $pp->{'label'} eq 'Official Symbol');
      }
    }
    ## if couldn't find the name and symbol on 'properties', try 'gene'
    if (!$symbol || !$name) {
      foreach my $g (@{$result->{'gene'}}){
        $g      = $g->[0] if(ref($g) eq 'ARRAY');
        $name   = $g->{'desc'}  if (!$name   && $g->{'desc'});
        $symbol = $g->{'locus'} if (!$symbol && $g->{'locus'});
      }
    }
    ## if still couldn't find the name and symbol try 'rna' and then 'prot'
    if (!$name) {
      foreach my $r (@{$result->{'rna'}}){
        $r    = $r->[0] if(ref($r) eq 'ARRAY');
        $name = $r->{'ext'}->[0]->{'name'} if (!$name && $r->{'ext'}->[0]->{'name'});
      }
    }
    if (!$name) {
      foreach my $p (@{$result->{'prot'}}){
        $p    = $p->[0] if(ref($p) eq 'ARRAY');
        $name = $p->{'name'}->[0] if (!$name && $p->{'name'}->[0]);
      }
    }

    my ($ensembl);
    foreach my $gdb (@{$result->{'gene'}->[0]->{'db'}}){
      $gdb = $gdb->[0] if(ref($gdb) eq 'ARRAY');
      next unless ($gdb->{'db'} && $gdb->{'db'} eq 'Ensembl');
      $ensembl = $gdb->{'tag'}->[0]->{'str'};
      last;
    }
    ## values for the gene-status (different from RefSeq status)
    ## live           good??
    ## secondary      synonym with merged
    ## discontinued   'deleted', still index and display to public
    ## newentry-      for GeneRif submission
    my $status = $result->{'track-info'}->[0]->{'status'};

    given ($status) {
      when ('discontinued') {
        log_it (3, "Discontinued gene: UID='$uid', symbol='$symbol', name='$name'. Forgetting about it...");
        next SEQ;
      }
      when ('secondary') {
        ## recursivity! UUUUUUUUUU!
        log_it (3, "Secondary gene: UID='$uid', symbol='$symbol', name='$name'. Attempting to find its current UID...");
        my $current_id;
        foreach my $c (@{$result->{'track-info'}->[0]->{'current-id'}}) {
          next unless $c->{'db'} eq 'GeneID';
          $current_id = $c->{'tag'}->[0]->{'id'};
          next unless $current_id;
          log_it (3, "Update: found current UID '$current_id' of secondary gene with UID='$uid'");
          analyze_entrez_genes ($struct, [$current_id]);
        }
        log_it (3, "Update: could not find current UID of secondary gene with UID='$uid'") unless $current_id;
        next SEQ;
      }
      default {
        if (!$status) {
          log_it (1, "WARNING: couldn't find gene status for gene with UID='$uid'. Assuming value of 'live'");
          $status = 'live-assumed';
        }
        my @extra_arguments;
        my $ng_message = "New gene: UID='$uid', gene status='$status', symbol='$symbol', name='$name'";
        if ($ensembl) {
          $ng_message = $ng_message . ", EnsEMBL ID='$ensembl'";
          push (@extra_arguments, ensembl => $ensembl);
        }
        log_it (3, "$ng_message");
        $struct->add_gene(
                          uid     => $uid,
                          status  => $status,
                          name    => $name,
                          symbol  => $symbol,
                          @extra_arguments,
                          );
      }
    }

    ## value is 'pseudo' for pseudo genes, should be 'protein-coding' otherwise
    if ($result->{'type'} eq 'pseudo') {
      log_it (3, "Update: gene with UID='$uid' is '". $result->{'type'} ."' gene. Marking as pseudo...");
      $struct->add_gene(uid => $uid, pseudo  => 1);
      unless ($get_pseudo) {
        log_it (3, "Update: removing gene with UID='$uid' for being pseudo...");
        $struct->remove_gene($uid);
        next SEQ;
      }
    } else {
      log_it (3, "Update: gene with UID='$uid' is '". $result->{'type'} ."' gene. Marking as protein-coding...");
      $struct->add_gene(uid => $uid, pseudo  => 0);
    }

    foreach my $l (@{$result->{'locus'}}){
      $l = $l->[0] if(ref($l) eq 'ARRAY');
      next unless ($l->{'heading'} && $l->{'heading'} =~ m/$assembly_regex/i);
      my $assembly  = $l->{'heading'};
      my $ChrAccVer = $l->{'accession'};
      my $ChrStart  = $l->{'seqs'}->[0]->{'int'}->[0]->{'from'};
      my $ChrStop   = $l->{'seqs'}->[0]->{'int'}->[0]->{'to'};
      my $ChrStrand = $l->{'seqs'}->[0]->{'int'}->[0]->{'strand'};
      if ($ChrStrand eq 'plus') {
        $ChrStrand   = 1;
        $ChrStart   += 1 - $upstream;
        $ChrStop    += 1 + $downstream;
      } else {
        $ChrStrand   = 2;
        $ChrStart   += 1 - $downstream;
        $ChrStop    += 1 + $upstream;
      }
      log_it (3, "Update: gene with UID='$uid' has Accesion number '$ChrAccVer' between coordinates $ChrStart ... $ChrStop on strand $ChrStrand.");
      $struct->add_gene(
                        uid       => $uid,
                        assembly  => $assembly,
                        ChrAccVer => $ChrAccVer,
                        ChrStart  => $ChrStart,
                        ChrStop   => $ChrStop,
                        ChrStrand => $ChrStrand,
                        );
      last; # if we got here once, no point in looking on the others
    }

    ## we will look for products accessions on the comments section instead of the
    ## locus section because locus doesn't say the RefSeq status of the products
    foreach my $c (@{$result->{'comments'}}){
      $c = $c->[0] if(ref($c) eq 'ARRAY');
      
      ## get RefSeq status
      if ($c->{'heading'} && $c->{'heading'} eq 'RefSeq Status') {
        my $refseq_status = $c->{'label'};
        log_it (3, "Update: gene with UID='$uid' has RefSeq status='$refseq_status'");
        $struct->add_gene(uid => $uid, 'RefSeq status' => $refseq_status);
      }

      ## the rest of the block only makes sense if it's not a pseudo gene
      if ($struct->get_info('gene', $uid, 'pseudo') ) {
        log_it (9, "DEBUG: finished analyzing gene with UID='$uid' earlier since it's pseudo gene " . by_caller_and_location('here') );
        next SEQ;
      }
      next unless ($c->{'heading'} && $c->{heading} eq 'NCBI Reference Sequences (RefSeq)');
      foreach my $cc (@{$c->{'comment'}}){
        $cc = $cc->[0] if(ref($cc) eq 'ARRAY');
        next unless $cc->{'heading'} eq 'RefSeqs maintained independently of Annotated Genomes';
        foreach my $ccp (@{$cc->{'products'}}){
          $ccp = $ccp->[0] if(ref($ccp) eq 'ARRAY');
          my $mRNA_acc  = $ccp->{'accession'};
          my $prot_acc  = $ccp->{'products'}->[0]->{'accession'};
          ## for the RefSeq status needs to be on a loop since there's a list of comments
          ## About RefSeq status of products:
          ## http://www.ncbi.nlm.nih.gov/entrez/query/static/help/genefaq.html#faq_g7.2
          ## XXX what can we do with this product RefSeq status?
          my $ref_stat;
          foreach my $ccpc (@{$ccp->{'comment'}}){
            next unless ($ccpc->{'label'} && $ccpc->{'label'} eq 'RefSeq Status');
            $ref_stat = $ccpc->{'text'};
          }
          ## some transcripts are non-coding. In those cases $prot_acc hence the
          ## need for this. For those cases, we want to give no value for 'protein'
          ## and not an empty string. The opposite (protein without corresponding
          ## transcript) should not be possible
          if (!$prot_acc) {
            log_it (3, "Update: gene with UID='$uid' found to encode transcrip='$mRNA_acc' which is non-coding transcript.");
            ## TODO some transcripts, do not encode a protein. This will create errors
            ## must fix for genes such as UID 100507436
            $struct->add_product(
                                  type       => 'transcript',
                                  accession  => $mRNA_acc,
                                  gene       => $uid,
                                  coding     => 0,
                                  );
          } else {
            log_it (3, "Update: gene with UID='$uid' found to encode transcrip='$mRNA_acc' and protein='$prot_acc' with product RefStatus of '$ref_stat'.");
            ## TODO some transcripts, do not encode a protein. This will create errors
            ## must fix for genes such as UID 100507436
            $struct->add_product(
                                  type       => 'transcript',
                                  accession  => $mRNA_acc,
                                  gene       => $uid,
                                  protein    => $prot_acc,
                                  coding     => 1,
                                  );
            $struct->add_product(
                                  type       => 'protein',
                                  accession  => $prot_acc,
                                  gene       => $uid,
                                  transcript => $mRNA_acc,
                                  );
          }
        }
      }
    }
    unless (scalar($struct->get_product_list('transcript', $uid)) >= 1 && scalar($struct->get_product_list('protein', $uid)) >= 1) {
      log_it (1, "WARNING: non-pseudo gene with UID='$uid' returned no protein and transcript.");
    }
=item *

Genes that are marked as 'live' and 'protein-coding' should have at least one
transcript. However, This is not always true due to mistakes on annotation. Such
cases will throw a warning. When faced with this, be nice and write to the entrez
RefSeq maintainers L<http://www.ncbi.nlm.nih.gov/RefSeq/update.cgi>.

=cut

  }
}

sub create_fetcher {
  my $searched  = shift;
  my $rettype;
  given ($format) {
    when ('genbank') {$rettype = 'gb';}
    when ('fasta')   {$rettype = 'fasta';}
    default          {$rettype = 'native'; log_it (1, "WARNING: couldn't convert format '$format' to rettype. Using native as default."); }
  }
  my $fetcher = Bio::DB::EUtilities->new(
                                          -eutil    => 'efetch',
                                          -db       => $searched,
                                          -retmode  => 'text',
                                          -rettype  => $rettype,
                                          );
  return $fetcher;
}

sub get_filehandle {
  my $type        = shift;
  my $product_key = shift;
  my $name_key    = shift;
  my $struct      = shift;
  my $base_dir;
  given ($type) {
    when ('gene')       { $base_dir = $gene_dir; }
    when ('transcript') { $base_dir = $mrna_dir; }
    when ('protein')    { $base_dir = $prot_dir; }
    default { die "Found a bug. Unknow type provided '$type' to generate filename ". by_caller_and_location('before') ." Please report."; }
  }
  my $filename = fix_filename ( $struct->get_info($type, $product_key, $name_key) );
  my $filepath = File::Spec->catfile ($base_dir, $filename . file_extension_for($format));
  open (my $filehandle, ">", $filepath) or die "Couldn't open '$filepath' for writing: $!";
  return $filehandle;
}

sub get_genes {
  my $struct    = shift;
  my @gene_uids = $struct->get_list('gene');
  my $fetcher   = create_fetcher('nucleotide');
  foreach my $gene_uid (@gene_uids) {
    log_it (2, "Fetching gene: trying to fetch gene with UID='$gene_uid'...");
    unless ($struct->get_info('gene', $gene_uid, 'ChrAccVer')) {
      log_it (2, "Update: found no genomic info for gene with UID='$gene_uid'...");
      next;
    }
    my $filehandle = get_filehandle('gene', $gene_uid, $genes, $struct);
    $fetcher->set_parameters (
                              -id         => $struct->get_info('gene', $gene_uid, 'ChrAccVer'),
                              -seq_start  => $struct->get_info('gene', $gene_uid, 'ChrStart'),
                              -seq_stop   => $struct->get_info('gene', $gene_uid, 'ChrStop'),
                              -strand     => $struct->get_info('gene', $gene_uid, 'ChrStrand'),
                             );
    log_it (3, "Update: fetching gene sequence for gene with UID='$gene_uid'...");
    print $filehandle $fetcher->get_Response->content;
    close $filehandle or warn $! ? "WARNING: error closing filehandle: $!" : "WARNING: exit status $? from filehandle";
  }
}

sub get_products {
  my $product       = shift;
  my $struct        = shift;
  my @product_acc   = $struct->get_list($product);
  
  my ($fetcher, $base_name);
  if ($product eq 'transcript') {
    $fetcher    = create_fetcher ('nuccore');
    $base_name  = $transcripts;
  } elsif ($product eq 'protein') {
    $fetcher    = create_fetcher ('protein');
    $base_name  = $proteins;
  } else {
    die "Bug found. Invalid product $product argument given ". by_caller_and_location('before') .". Please report it.";
  }

  ## ideally, there would be no loop, and we'd use $fetcher to get all the sequences in
  ## one go. However, that would force us to get the sequences in Bio::Seq files which
  ## can be different from the actually downloaded file. Better to not take the chance
  foreach my $product_acc (@product_acc) {
    if ( !$get_noncoding && $product eq 'transcript' && !$struct->get_info($product, $product_acc, 'coding') ) {
      log_it (2, "Fetching $product: skipping $product with accession='$product_acc' since it's non-coding...");
    } else {
      log_it (2, "Fetching $product: trying to fetch $product with accession='$product_acc'...");
    }
    $fetcher->set_parameters (
                              -id => $product_acc,
                               );
    my $response = $fetcher->get_Response->content;
    open(my $seq_fh, "<", \$response) or die "Could not open sequence string for reading: $!";
    my $parser = Bio::SeqIO->new(
                                 -fh      => $seq_fh,
                                 -format  => $format,
                                 );
    while (my $seq = $parser->next_seq) {
      my $product_desc = $seq->description;
      $struct->add_product(
                            type        => $product,
                            accession   => $product_acc,
                            description => $product_desc,
                            );
    }
    my $filehandle = get_filehandle($product, $product_acc, $base_name, $struct);
    print $filehandle $response;
    close $filehandle;
  }
}

################################################################################
## other small useful subs
################################################################################

=item *
When creating the directories to save the files, if the directory already exists it will be used and no error
or warning will be issued unless B<--debug> as been set. If a non-directory file already exists with that name
bp_genbank_ref_extractor exits with an error.

=cut

## checks if directories exist and creates them as needed
sub check_dir {
  # TODO must check for permissions as well
  if (!-e $_[0]) {
    mkdir $_[0] or die "Could not mkdir '$_[0]': $!";
    log_it (9, "DEBUG: directory '$_[0]' created.");
  } elsif (!-d $_[0]) {
    die "Directory '$_[0]' to save output already exists as non-directory.";
  } else {
    log_it (9, "DEBUG: directory '$_[0]' NOT created since it already exists.");
  }
}

=item *
On the subject of verbosity, all messages are saved on the log file. The options
B<--verbose> and B<--very-verbose> only affect their printing to standard
output. Debug messages are different as they will only show up (and be logged)
if requested with B<--debug>.

=cut

sub log_it {
  my $level   = shift;
  my $message = shift;
  say STDOUT $message unless $level > $verbosity;
  ## debug messages (level 9) will only be logged if on debug mode
  if ($level == 9) {
    say LOG $message if $debug;
  } else {
    say LOG $message;
  }
}

## Removes repeated elements from an array. Does not respect original order
sub clean_array {
  my %hash;
  foreach (@{$_[0]}) {
    if ($hash{$_}) {
      log_it (9, "DEBUG: value '$_' removed from array " . by_caller_and_location('here') . " called " . by_caller_and_location('before') );
    } else {
      $hash{$_} = 1;
    }
  }
  @{$_[0]} = keys %hash;
}

## Returns a pretty string about current time
sub get_time {
  my ($second, $minute, $hour, $day, $month, $year) = (localtime)[0,1,2,3,4,5];
  return sprintf ("[%04d-%02d-%02d %02d:%02d:%02d]", $year+1900, $month+1, $day, $hour, $minute, $second);
}

=item *
When saving a file, to avoid problems with limited filesystems such as NTFS or FAT, only some
characters are allowed. All other characters will be replaced by an underscore. Allowed characters
are:

B<a-z 0-9 - +  . , () {} []'>

=cut

## Tries to sanitize a filename
sub fix_filename {
  my $file = $_[0];
  $file =~ s/[^a-z0-9\-\+ \.,\(\){}\[\]']/_/ig;
  log_it (9, "DEBUG: filepath '$_[0]' was converted to '$file' " . by_caller_and_location('here') . " called " . by_caller_and_location('before') );
  return $file;
}

sub by_caller_and_location {
  my $level;
  if (!@_ || $_[0] eq 'here') {
    $level = 1;
  } elsif ($_[0] eq 'before'){
    $level = 2;
  } elsif ($_[0] =~ /^[0-9]+$/){
    $level = 1 + $_[0];
  } else {
    die "Bug found when calculating level for caller function. Please report.";
  }
  my $deeper = shift;
  return "by " . (caller($level))[3] . " at line " . (caller($level))[2];
}

=item *

B<bp_genbank_ref_extractor> tries to use the same file extensions that bioperl
would expect when saving the file. If unable it will use the '.seq' extension.

=cut

sub file_extension_for {
  ## TODO in some cases, extension changes whether it's protein or DNA or whatever
  ## and this should be supported
  ## XXX there must be a more elegant to handle the formats on this scripts

  ## to update this list, look in the _guess_format method, inside SeqIO.pm of bioperl
  given ($_[0]) {
    when (/embl/i)       {return '.embl';}
    when (/entrezgene/i) {return '.asn';}
    when (/fasta/i)      {return '.fasta';} # fasta|fast|fas|seq|fa|fsa|nt|aa|fna|faa
    when (/fastq/i)      {return '.fastq';}
    when (/gcg/i)        {return '.gcg';}
    when (/genbank/i)    {return '.gb';} # gb|gbank|genbank|gbk|gbs
    when (/swiss/i)      {return '.swiss';} # swiss|sp
    default {
      log_it (9, "DEBUG: couldn't find the right extension for the requested format. Using '.seq' as default.");
      return ".seq";
    }
  }
}

sub save_structure {
  given ($save_data) {
    when ('csv') { create_csv($_[0]); }
  }
}

sub create_csv {
  my $struct = shift;
  my $csv = Text::CSV->new ({
                              binary => 1,
                              eol => $/,
                              }) or die "Cannot use Text::CSV: ". Text::CSV->error_diag ();

  my $csv_file  = File::Spec->catfile ($save, 'data.csv');
  open (my $fh, ">", $csv_file) or die "Couldn't open file $csv_file for writing: $!";

  $csv->print ($fh, ['gene symbol', 'gene UID', 'EnsEMBL ID', 'gene name', 'pseudo', 'transcript accession','protein accession', 'chromosome accession', 'chromosome start coordinates', 'chromosome stop coordinates', 'assembly'] );

  my @uids = $struct->get_list('gene');
  foreach my $uid(@uids) {
    my @lines;
    my @mRNA_acc = $struct->get_product_list('transcript', $uid);
    if (!@mRNA_acc) { @mRNA_acc = (''); }   # this allows the next loop to run once for pseudo genes
    foreach my $mRNA_acc (@mRNA_acc) {
      push(@lines, [
                    $struct->get_info('gene', $uid, 'symbol'),
                    $uid,
                    $struct->get_info('gene', $uid, 'ensembl'),
                    $struct->get_info('gene', $uid, 'name'),
                    $struct->get_info('gene', $uid, 'pseudo'),
                    $mRNA_acc,
                    $struct->get_info('transcript', $mRNA_acc, 'protein'),
                    $struct->get_info('gene', $uid, 'ChrAccVer'),
                    $struct->get_info('gene', $uid, 'ChrStart'),
                    $struct->get_info('gene', $uid, 'ChrStop'),
                    $struct->get_info('gene', $uid, 'assembly'),
                    ]);
    }
    $csv->print ($fh, $_) for @lines;
  }
  close $fh;
}


################################################################################
## Structure methods
################################################################################
package Structure;

## creates a new instance of the object
sub new {
  my $class = shift;
  my $self  = {};
  bless ($self, $class);
  return $self;
}

## adds information to a specific gene and adds the gene to the structure if it doesn't exist
## $object->add_gene (
##                    uid        => gene_uid_of_the_gene,
##                    gene_name  => gene_name,
##                    other_info => corresponding value
##                    );
sub add_gene {
  my $self  = shift;
  my %data  = @_;
  ## remove the value from the hash so it doesn't show up in the loop
  my $gene_uid = delete $data{'uid'};
  unless ($gene_uid) {
    log_it (1, "WARNING: no gene UID supplied when adding new gene". main::by_caller_and_location('before') );
    return;
  }
  ## when adding a gene (instead of updating, this goes first. Can't just hope
  ## to happen by itself on the loop later because it's possible to add a gene
  ## with no info on it
  if ( !exists($self->{'gene'}->{$gene_uid}) ) {
    main::log_it (9, "DEBUG: creating new gene with UID='$gene_uid'.");
    $self->{'gene'}->{$gene_uid} = {};
    $self->{'gene'}->{$gene_uid}->{'uid'} = $gene_uid;  # this is not stupid. Makes it easier to have a general function to create the filename
    $self->{'gene'}->{$gene_uid}->{'transcript'} = [];
    $self->{'gene'}->{$gene_uid}->{'protein'}    = [];
  }
  ## fill it with all the data
  foreach my $key (keys %data) {
    $self->{'gene'}->{$gene_uid}->{$key} = $data{$key};
    main::log_it (9, "DEBUG: added $key='$data{$key}' to gene with UID='$gene_uid'.");
  }
}

## remove genes from the structure given their uid
## $object->remove_gene($uid)
## $object->remove_gene(@uids)
sub remove_gene {
  my $self  = shift;
  foreach my $uid (@_) {
    delete $self->{'gene'}->{$uid};
    main::log_it (9, "DEBUG: removed gene with UID='$uid'");
  }
}

sub add_product {
  my $self  = shift;
  my %data  = @_;
  ## remove these values from the hash so they don't show up in the loop later
  my $product     = delete $data{'type'};
  die "Bug found. Please report this. Product requested was $product ". main::by_caller_and_location('before') unless ($product eq 'protein' || $product eq 'transcript');
  my $product_acc = delete $data{'accession'};
  unless ($product_acc) {
    main::log_it(1, "WARNING: no $product accession supplied when adding new product ". main::by_caller_and_location('before') . " Adding nothing." );
    return;
  }
  ## Since it's possible that a record for a product be generated automatically when
  ## creating it's related product, the only way to know if it's the first time, it's
  ## to check it's accession (it's the same as the key. It looks kinda of redundant
  ## but allows to have a simpler function to generate filename that uses get_info)
  if ( !exists($self->{$product}->{$product_acc}->{'accession'}) ) {
    main::log_it (9, "DEBUG: creating new $product with accession='$product_acc'.");
    $self->{$product}->{$product_acc}->{'accession'} = $product_acc;
  }
  ## fill it with all the data
  foreach my $key (keys %data) {
    $self->{$product}->{$product_acc}->{$key} = $data{$key};
    ## if we're adding info about gene and related products, the array on their
    ## part of structure  needs to be updated
    if ($key eq 'gene') {
      my $products_on_gene = \@{$self->{'gene'}->{$data{$key}}->{$product}};
      push (@$products_on_gene, $product_acc);
      main::clean_array($products_on_gene);
    } elsif ($key eq 'transcript' && $product eq 'protein') {
      my $transcript_acc  = $data{$key};
      my $current         = $self->{'transcript'}->{$transcript_acc}->{'protein'};
      if ($current && $current ne $product_acc) {
        warn "WARNING: replacing accession $current with $product_acc as product of $transcript_acc. Please report this bug.";
        $self->{'transcript'}->{'protein'} = $product_acc;
      }
    } elsif ($key eq 'protein' && $product eq 'transcript') {
      my $protein_acc  = $data{$key};
      my $current      = $self->{'protein'}->{$protein_acc}->{'transcript'};
      if ($current && $current ne $product_acc) {
        warn "WARNING: replacing accession $current with $product_acc as ''template'' of $product_acc. Please report this bug.";
        $self->{'protein'}->{'transcript'} = $product_acc;
      }
    }
  }
}
## get information from the structure
## $value = $object->get_info('gene', $gene_id, 'wanted info');
## $value = $object->get_info('protein', $protein_acc, 'wanted info');
sub get_info {
  my $self  = shift;
  my $type  = shift;
  my $key   = shift;
  my $req   = shift;
  ## can't check the request here and return the key if the request is accession
  ## or id because even though we're using product accessions and gene UID as keys,
  ## gene also have accessions and products have id and if we had support for them
  ## later, it would create a lot of confusion and bugs
  return $self->{$type}->{$key}->{$req};
}

## returns a list of all gene UIDs or product accessions
##   @gene_uids = $structure->get_list('gene')
##   @mRNA_acc  = $structure->get_list('transcript')
sub get_list {
  my $self  = shift;
  my $type  = shift;
  return keys %{$self->{$type}};
}

## for the specified genes UIDs returns list of the requested products accessions
## If no gene_id is specified, returns a list of all accessions of the product requested
## @mRNA_acc  = structure->get_product_list('transcript', $gene_uid)
## @prot_acc  = structure->get_product_list('protein', @gene_uid)
## @prot_acc  = structure->get_product_list('protein')
sub get_product_list {
  my $self    = shift;
  my $product = shift;
  given ( scalar(@_) ) {
    when ([1]) { return @{$self->{'gene'}->{$_[0]}->{$product}}; }
    when ([0]) { return $self->get_list($product); }
    default    { return map { $self->get_product_list($product, $_) } @_; }
  }
}

## this =back closes the last point on the NOTES on usage section

=back

=cut
