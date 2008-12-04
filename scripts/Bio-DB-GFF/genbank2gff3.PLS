#!/usr/bin/perl -w

#$Id$;

=pod

=head1 NAME 

genbank2gff3.pl -- Genbank-E<gt>gbrowse-friendly GFF3

=head1 SYNOPSIS

  genbank2gff3.pl [options] filename(s)

  # process a directory containing GenBank flatfiles
  perl genbank2gff3.pl --dir path_to_files --zip

  # process a single file, ignore explicit exons and introns
  perl genbank2gff3.pl --filter exon --filter intron file.gbk.gz

  # process a list of files 
  perl genbank2gff3.pl *gbk.gz

  # process data from URL, with Chado GFF model (-noCDS), and pipe to database loader
  curl ftp://ftp.ncbi.nih.gov/genomes/Saccharomyces_cerevisiae/CHR_X/NC_001142.gbk \
  | perl genbank2gff3.pl -noCDS -in stdin -out stdout \
  | perl gmod_bulk_load_gff3.pl -dbname mychado -organism fromdata

    Options:
        --dir     -d  path to a list of genbank flatfiles
        --outdir  -o  location to write GFF files (can be 'stdout' or '-' for pipe)
        --zip     -z  compress GFF3 output files with gzip
        --summary -s  print a summary of the features in each contig
        --filter  -x  genbank feature type(s) to ignore
        --split   -y  split output to seperate GFF and fasta files for
                      each genbank record
        --nolump  -n  seperate file for each reference sequence
                      (default is to lump all records together into one 
                       output file for each input file)
        --ethresh -e  error threshold for unflattener
                      set this high (>2) to ignore all unflattener errors
        --[no]CDS -c  Keep CDS-exons, or convert to alternate gene-RNA-protein-exon 
                      model. --CDS is default. Use --CDS to keep default GFF gene model, 
                      use --noCDS to convert to g-r-p-e.
        --format  -f  Input format (SeqIO types): GenBank, Swiss or Uniprot, EMBL work
                      (GenBank is default)
        --GFF_VERSION 3 is default, 2 and 2.5 and other Bio::Tools::GFF versions available
        --quiet       dont talk about what is being processed 
        --typesource  SO sequence type for source (e.g. chromosome; region; contig)
        --help    -h  display this message


=head1 DESCRIPTION

This script uses Bio::SeqFeature::Tools::Unflattener and
Bio::Tools::GFF to convert GenBank flatfiles to GFF3 with gene
containment hierarchies mapped for optimal display in gbrowse.

The input files are assumed to be gzipped GenBank flatfiles for refseq
contigs.  The files may contain multiple GenBank records.  Either a
single file or an entire directory can be processed.  By default, the
DNA sequence is embedded in the GFF but it can be saved into seperate
fasta file with the --split(-y) option.

If an input file contains multiple records, the default behaviour is
to dump all GFF and sequence to a file of the same name (with .gff
appended).  Using the 'nolump' option will create a seperate file for
each genbank record.  Using the 'split' option will create seperate
GFF and Fasta files for each genbank record.


=head2 Notes

=head3 'split' and 'nolump' produce many files

In cases where the input files contain many GenBank records (for
example, the chromosome files for the mouse genome build), a very
large number of output files will be produced if the 'split' or
'nolump' options are selected.  If you do have lists of files E<gt> 6000,
use the --long_list option in bp_bulk_load_gff.pl or
bp_fast_load_gff.pl to load the gff and/ or fasta files.

=head3 Designed for RefSeq

This script is designed for RefSeq genomic sequence entries.  It may
work for third party annotations but this has not been tested.
But see below, Uniprot/Swissprot works, EMBL and possibly EMBL/Ensembl
if you don't mind some gene model unflattener errors (dgg).

=head3 G-R-P-E Gene Model

Don Gilbert worked this over with needs to produce GFF3 suited to
loading to GMOD Chado databases.  Most of the changes I believe are
suited for general use.  One main chado-specific addition is the
  --[no]cds2protein  flag

My favorite GFF is to set the above as ON by default (disable with --nocds2prot)
For general use it probably should be OFF, enabled with --cds2prot.

This writes GFF with an alternate, but useful Gene model,
instead of the consensus model for GFF3 

  [ gene > mRNA> (exon,CDS,UTR) ]

This alternate is

  gene > mRNA > polypeptide > exon 

means the only feature with dna bases is the exon.  The others
specify only location ranges on a genome.  Exon of course is a child
of mRNA and protein/peptide.   

The protein/polypeptide feature is an important one, having all the
annotations of the GenBank CDS feature, protein ID, translation, GO
terms, Dbxrefs to other proteins.

UTRs, introns, CDS-exons are all inferred from the primary exon bases
inside/outside appropriate higher feature ranges.   Other special gene
model features remain the same.

Several other improvements and bugfixes, minor but useful are included

  * IO pipes now work:
    curl ftp://ncbigenomes/... | genbank2gff3 --in stdin --out stdout | gff2chado ...

  * GenBank main record fields are added to source feature, e.g. organism, date,
    and the sourcetype, commonly chromosome for  genomes, is used.

  * Gene Model handling for ncRNA, pseudogenes are added.

  * GFF header is cleaner, more informative.
    --GFF_VERSION flag allows choice of v2 as well as default v3

  * GFF ##FASTA inclusion is improved, and
    CDS translation sequence is moved to FASTA records.

  * FT -> GFF attribute mapping is improved.

  * --format choice of SeqIO input formats (GenBank default). 
    Uniprot/Swissprot and EMBL work and produce useful GFF.

  * SeqFeature::Tools::TypeMapper has a few FT -> SOFA additions
      and more flexible usage.

=head1 TODO

=head2 Are these additions desired?

 * filter input records by taxon (e.g. keep only organism=xxx or taxa level = classYYY
 * handle Entrezgene, other non-sequence SeqIO structures (really should change
    those parsers to produce consistent annotation tags).

=head2 Related bugfixes/tests

These items from Bioperl mail were tested (sample data generating
errors), and found corrected:

 From: Ed Green <green <at> eva.mpg.de>
 Subject: genbank2gff3.pl on new human RefSeq
 Date: 2006-03-13 21:22:26 GMT 
   -- unspecified errors (sample data works now).

 From: Eric Just <e-just <at> northwestern.edu>
 Subject: bp_genbank2gff3.pl
 Date: 2007-01-26 17:08:49 GMT
   -- bug fixed in genbank2gff3 for multi-record handling

This error is for a /trans_splice gene that is hard to handle, and
unflattner/genbank2 doesn't

 From: Chad Matsalla <chad <at> dieselwurks.com> 
 Subject: genbank2gff3.PLS and the unflatenner - Inconsistent	order?
 Date: 2005-07-15 19:51:48 GMT 

=head1 AUTHOR 

Sheldon McKay (mckays@cshl.edu)

Copyright (c) 2004 Cold Spring Harbor Laboratory.

=head2 AUTHOR of hacks for GFF2Chado loading

Don Gilbert (gilbertd@indiana.edu)


=cut

use strict;

use lib "$ENV{HOME}/bioperl-live";
# chad put this here to enable situations when this script is tested
# against bioperl compiled into blib along with other programs using blib
BEGIN {
        unshift(@INC,'blib/lib');
};
use Pod::Usage;
use Bio::Root::RootI;
use Bio::SeqIO;
use File::Spec;
use Bio::SeqFeature::Tools::Unflattener;
use Bio::SeqFeature::Tools::TypeMapper;
use Bio::SeqFeature::Tools::IDHandler;
use Bio::Location::SplitLocationI;
use Bio::Location::Simple;
use Bio::Tools::GFF;
use Getopt::Long;

use vars qw/$split @filter $zip $outdir $help $ethresh
            $file @files $dir $summary $nolump 
            $source_type %proteinfa %exonpar $didheader $verbose $DEBUG $GFF_VERSION 
            $gene_id $rna_id $tnum $ncrna_id $rnum %method %id %seen/;

use constant GM_NEW_TOPLEVEL => 2;
use constant GM_NEW_PART => 1;
use constant GM_DUP_PART => 0;
use constant GM_NOT_PART => -1;

$GFF_VERSION = 3; # allow v2 ...
$verbose = 1; # right default? -nov to turn off

# dgg: change the gene model to  Gene/mRNA/Polypeptide/exons...
my $CDSkeep= 1; # default should be ON (prior behavior), see gene_features()
my $PROTEIN_TYPE = 'polypeptide'; # for noCDSkeep; 
  # protein = flybase chado usage; GMOD Perls use 'polypeptide' with software support

my $FORMAT="GenBank"; # swiss ; embl; genbank ; ** guess from SOURCEID **
my $SOURCEID= $FORMAT; # "UniProt" "GenBank"  "EMBL" should work  
   # other Bio::SeqIO formats may work.  TEST: EntrezGene < problematic tags; InterPro  KEGG 


my %TAG_MAP = (
  db_xref => 'Dbxref',
  name => 'Name',
  note => 'Note', # also pull GO: ids into Ontology_term
  synonym => 'Alias',
  symbol => 'Alias',  # is symbol still used?
  # protein_id => 'Dbxref', also seen Dbxref tags: EC_number 
  # translation: handled in gene_features
);


$| = 1;
my $quiet= !$verbose;
my $ok= GetOptions( 'd|dir|input:s'   => \$dir,
            'z|zip'     => \$zip, 
            'h|help'    => \$help,
            's|summary' => \$summary,
            'o|outdir|output:s'=> \$outdir,
            'x|filter:s'=> \@filter,
            'y|split'   => \$split,
            "ethresh|e=s"=>\$ethresh,
            'c|CDS!'    => \$CDSkeep,
            'f|format=s' => \$FORMAT,
            'typesource=s' => \$source_type,
            'GFF_VERSION=s' => \$GFF_VERSION,
            'quiet!'    => \$quiet, # swap quiet to verbose
            'DEBUG!'    => \$DEBUG,
            'n|nolump'  => \$nolump);

my $lump = 1 unless $nolump || $split;
$verbose= !$quiet;

# look for help request
pod2usage(2) if $help || !$ok;

# keep SOURCEID as-is and change FORMAT for SeqIO types; 
# note SeqIO uses file.suffix to guess type; not useful here
$SOURCEID= $FORMAT; 
$FORMAT  = "swiss" if $FORMAT =~/UniProt|trembl/;
$verbose =1 if($DEBUG);

# initialize handlers
my $unflattener = Bio::SeqFeature::Tools::Unflattener->new; # for ensembl genomes (-trust_grouptag=>1);
$unflattener->error_threshold($ethresh) if $ethresh;
$unflattener->verbose(1) if($DEBUG);
# $unflattener->group_tag('gene') if($FORMAT =~ /embl/i) ; #? ensembl only? 
# ensembl parsing is still problematic, forget this

my $tm  = Bio::SeqFeature::Tools::TypeMapper->new;
my $idh = Bio::SeqFeature::Tools::IDHandler->new;

# dgg
$source_type ||= "region"; # should really parse from FT.source contents below

my $FTSOmap = $tm->FT_SO_map();
# #convert $FTSOmap undefined to valid SO : moved to TypeMapper->map_types( -undefined => "region")

# stringify filter list if applicable
my $filter = join ' ', @filter  if @filter;

# determine input files
my $stdin=0; # dgg: let dir == stdin == '-' for pipe use
if ($dir && ($dir eq '-' || $dir eq 'stdin')) {
  $stdin=1;  $dir=''; @files=('stdin');
  
} elsif ( $dir ) {
    if ( -d $dir ) {
        opendir DIR, $dir or die "could not open $dir for reading: $!";
        @files = map { "$dir/$_";} grep { /\.gb.*/ } readdir DIR;  
        closedir DIR;
    }
    else {
        die "$dir is not a directory\n";
    }
}
else {
    @files = @ARGV;
    $dir = '';
}

# we should have some files by now
pod2usage(2) unless @files;


my $stdout=0; # dgg: let outdir == stdout == '-' for pipe use
if($outdir && ($outdir eq '-' || $outdir eq 'stdout')) {
  warn("std. output chosen: cannot split\n") if($split);
  warn("std. output chosen: cannot zip\n") if($zip);
  warn("std. output chosen: cannot nolump\n") if($nolump);
  $stdout=1; $lump=1; $split= 0; $zip= 0; # unless we pipe stdout thru gzip
  
} elsif ( $outdir && !-e $outdir ) {
    mkdir($outdir) or die "could not create directory $outdir: $!\n";        
}
elsif ( !$outdir ) {
    $outdir = $dir || '.';
}

$outdir .= '/' unless $outdir =~ m|/$|;

for my $file ( @files ) {
    # dgg ; allow 'stdin' / '-' input ?
    chomp $file;
    die "$! $file" unless($stdin || -e $file);
    print "# Input: $file\n" if($verbose);

    my ($lump_fh, $lumpfa_fh, $outfile, $outfa);
    if ($stdout) {
      $lump_fh= *STDOUT; $lump="stdout$$";
      $outfa= "stdout$$.fa";  # this is a temp file ... see below
      open $lumpfa_fh, ">$outfa" or die "Could not create a lump outfile called ($outfa) because ($!)\n";

    } elsif ( $lump ) {
                # this is better, but still should use catfile
        my ($vol,$dirs,$fileonly) = File::Spec->splitpath($file); 
        $lump   = $outdir . $fileonly . '.gff';
        ($outfa= $lump) =~ s/\.gff/\.fa/;
        open $lump_fh, ">$lump" or die "Could not create a lump outfile called ($lump) because ($!)\n";
        open $lumpfa_fh, ">$outfa" or die "Could not create a lump outfile called ($outfa) because ($!)\n";

    }
    
    # open input file, unzip if req'd
    if ($stdin) {
       *FH= *STDIN;   
    } elsif ( $file =~ /\.gz/ ) {
        open FH, "gunzip -c $file |";
    }
    else {
        open FH, "<$file";
    }

    my $in = Bio::SeqIO->new(-fh => \*FH, -format => $FORMAT, -debug=>$DEBUG);
    my $gffio = Bio::Tools::GFF->new( -noparse => 1, -gff_version => $GFF_VERSION );

    while ( my $seq = $in->next_seq ) {
        my $seq_name = $seq->accession_number;
        my $end = $seq->length;
        my @to_print;

        # arrange disposition of GFF output
        $outfile = $lump || $outdir . $seq_name . ".gff";
        my $out;

        if ( $lump ) {
            $outfile = $lump;
            $out = $lump_fh;
        }
        else {
            $outfile = $outdir . $seq_name . ".gff";
            open $out, ">$outfile";
        }

        # filter out unwanted features
        my $source_feat= undef;
        my @source= filter($seq); $source_feat= $source[0];

        ($source_type,$source_feat)= 
          getSourceInfo( $seq, $source_type, $source_feat ) ;
          # always; here we build main prot $source_feat; # if @source;
          
        # abort if there are no features
        warn "$seq_name has no features, skipping\n" and next
            if !$seq->all_SeqFeatures;


        $FTSOmap->{'source'}= $source_type;
        ## $FTSOmap->{'CDS'}= $PROTEIN_TYPE; # handle this in gene_features
        
        # construct a GFF header
        # add: get source_type from attributes of source feature? chromosome=X tag
        # also combine 1st ft line here with source ft from $seq ..
        my($header,$info)= gff_header($seq_name, $end, $source_type, $source_feat);
        print $out $header;
        print "# working on $info\n" if($verbose);
        
        # unflatten gene graphs, apply SO types, etc; this also does TypeMapper ..
        unflatten_seq($seq);


        # Note that we use our own get_all_SeqFeatures function 
        # to rescue cloned exons
        for my $feature ( get_all_SeqFeatures($seq) ) {
            
            my $method = $feature->primary_tag;
            next if($SOURCEID =~/UniProt|swiss|trembl/i && $method ne $source_type);
            
            $feature->seq_id($seq->id) unless($feature->seq_id);
            $feature->source_tag($SOURCEID);
            
            # dgg; need to convert some Genbank to GFF tags: note->Note; db_xref->Dbxref;
            ## also, pull any GO:000 ids from /note tag and put into Ontology_term
            maptags2gff($feature);
            
            # current gene name.  The unflattened gene features should be in order so any
            # exons, CDSs, etc that follow will belong to this gene
            my $gene_name;
            if ( $method eq 'gene' || $method eq 'pseudogene' ) {
              @to_print= print_held($out, $gffio, \@to_print);
              $gene_id = $gene_name= gene_name($feature); 
            } else {
              $gene_name= gene_name($feature);
            }
        
            #?? should gene_name from /locus_tag,/gene,/product,/transposon=xxx
            # be converted to or added as  Name=xxx (if not ID= or as well)
            ## problematic: convert_to_name ($feature); # drops /locus_tag,/gene, tags
            
            ## dgg: extended to protein|polypeptide
            ## this test ($feature->has_tag('gene') ||) is not good: repeat_regions over genes
            ## in yeast have that genbank tag; why?
            ## these include pseudogene ...
            
            ## Note we also have mapped types to SO, so these RNA's are now transcripts:
            # pseudomRNA => "pseudogenic_transcript", 
            # pseudotranscript" => "pseudogenic_transcript", 
            # misc_RNA=>'processed_transcript',

            warn "#at: $method $gene_id/$gene_name\n" if $DEBUG;

            if ( $method =~ /(gene|RNA|CDS|exon|UTR|protein|polypeptide|transcript)/ 
               || ( $gene_id && $gene_name eq $gene_id ) ) {
               
                my $action = gene_features($feature, $gene_id, $gene_name);  # -1, 0, 1, 2 result
                if ($action == GM_DUP_PART) {
                  # ignore, this is dupl. exon with new parent ...
                  
                } elsif ($action == GM_NOT_PART) {
                  add_generic_id( $feature, $gene_name, "nocount");
                  my $gff= $gffio->gff_string($feature);
                  print $out "$gff\n" if $gff;

                } elsif ($action > 0) {
                 # hold off print because exon etc. may get 2nd, 3rd parents
                  @to_print= print_held($out, $gffio, \@to_print) if ($action == GM_NEW_TOPLEVEL);
                  push(@to_print, $feature);
                }
            }
            
            # otherwise handle as generic feats with IDHandler labels 
            else {
                add_generic_id( $feature, $gene_name, "");
                my $gff= $gffio->gff_string($feature);
                print $out "$gff\n" if $gff;
            }
        }

        # dont like doing this after others; do after each new gene id?
        @to_print= print_held($out, $gffio, \@to_print);
          
        # deal with the corresponding DNA
        my ($fa_out,$fa_outfile);
        my $dna = $seq->seq;
        if($dna || %proteinfa) {
          $method{'RESIDUES'} += length($dna);
          $dna    =~ s/(\S{60})/$1\n/g;
          $dna   .= "\n";
          
          if ($split) {
              $fa_outfile = $outfile;
              $fa_outfile =~ s/gff$/fa/;
              open $fa_out, ">$fa_outfile" or die $!; 
              print $fa_out ">$seq_name\n$dna" if $dna;
              foreach my $aid (sort keys %proteinfa) { 
                my $aa= delete $proteinfa{$aid}; 
                $method{'RESIDUES(tr)'} += length($aa);
                $aa =~ s/(\S{60})/$1\n/g; 
                print $fa_out ">$aid\n$aa\n"; 
                }
               
          }
          else {
          ## problem here when multiple GB Seqs in one file; all FASTA needs to go at end of $out
          ## see e.g. Mouse: mm_ref_chr19.gbk has NT_082868 and NT_039687 parts in one .gbk
          ## maybe write this to temp .fa then cat to end of lumped gff $out
              print $lumpfa_fh ">$seq_name\n$dna" if $dna;
              foreach my $aid (sort keys %proteinfa) { 
                my $aa= delete $proteinfa{$aid}; 
                $method{'RESIDUES(tr)'} += length($aa);
                $aa =~ s/(\S{60})/$1\n/g; 
                print $lumpfa_fh ">$aid\n$aa\n"; 
                }
          }
          
        %proteinfa=();
        }
     
        if ( $zip && !$lump ) {
            system "gzip -f $outfile";
            system "gzip -f $fa_outfile" if($fa_outfile);
            $outfile .= '.gz';
            $fa_outfile .= '.gz' if $split;
        }

        # print "\n>EOF\n" if($stdout); #?? need this if summary goes to stdout after FASTA
        print "# GFF3 saved to $outfile" unless( !$verbose || $stdout || $lump);
        print ($split ? "; DNA saved to $fa_outfile\n" : "\n") unless($stdout|| $lump);
        
        # dgg: moved to after all inputs; here it prints cumulative sum for each record
#         if ( $summary ) {
#             print "# Summary:\n# Feature\tCount\n# -------\t-----\n";
#         
#             for ( keys %method ) {
#                 print "# $_  $method{$_}\n";
#             }
#             print "# \n";
#         }       
    
    }

    print "# GFF3 saved to $outfile\n" if( $verbose && $lump);
    if ( $summary ) {
        print "# Summary:\n# Feature\tCount\n# -------\t-----\n";
    
        for ( keys %method ) {
            print "# $_  $method{$_}\n";
        }
        print "# \n";
    }
    
     ## FIXME for piped output w/ split FA files ...
    close($lumpfa_fh);
    if (!$split && $outfa && $lump_fh) {     
      print $lump_fh "##FASTA\n"; # GFF3 spec
      open $lumpfa_fh, $outfa or warn "reading FA $outfa: $!";
      while( <$lumpfa_fh>) { print $lump_fh $_; } # is $lump_fh still open?
      close($lumpfa_fh); unlink($outfa);
      }
        

    if ( $zip && $lump ) {
        system "gzip -f $lump";
    }
    
    close FH;
}





sub typeorder {
  return 1 if ($_[0] =~ /gene/);
  return 2 if ($_[0] =~ /RNA|transcript/);
  return 3 if ($_[0] =~ /protein|peptide/);
  return 4 if ($_[0] =~ /exon|CDS/);
  return 3; # default before exon (smallest part)
}

sub sort_by_feattype {
  my($at,$bt)= ($a->primary_tag, $b->primary_tag);
  return (typeorder($at) <=> typeorder($bt))  
    or ($at cmp $bt);
    ## or ($a->name() cmp $b->name());
}

sub print_held {  
  my($out,$gffio,$to_print)= @_;
  return unless(@$to_print);
  @$to_print = sort sort_by_feattype  @$to_print; # put exons after mRNA, otherwise chado loader chokes
  while ( my $feature = shift @$to_print) {
    my $gff= $gffio->gff_string($feature); # $gff =~ s/\'/./g; # dang bug in encode
    print $out "$gff\n";
    }
  return (); # @to_print
}

sub maptags2gff {
  my $f = shift;
  ## should copy/move locus_tag to Alias, if not ID/Name/Alias already
  # but see below /gene /locus_tag usage
  foreach my $tag (keys %TAG_MAP) {
    if ($f->has_tag($tag)) {
      my $newtag= $TAG_MAP{$tag};
      my @v= $f->get_tag_values($tag);
      $f->remove_tag($tag);
      $f->add_tag_value($newtag,@v);
      
      ## also, pull any GO:000 ids from /note tag and put into Ontology_term
      ## ncbi syntax in CDS /note is now '[goid GO:0005886]' OR '[goid 0005624]'
      if ($tag eq 'note') {  
        map { s/\[goid (\d+)/\[goid GO:$1/g; } @v; 
        my @go= map { m/(GO:\d+)/g } @v; 
        $f->add_tag_value('Ontology_term',@go) if(@go);
        }
      
      }
    }    
}


sub getSourceInfo {
  my ($seq, $source_type, $sf) = @_;  
  
  my $is_swiss= ($SOURCEID =~/UniProt|swiss|trembl/i);
  my $is_gene = ($SOURCEID =~/entrezgene/i);
  my $is_rich = (ref($seq) =~ /RichSeq/);
  my $seq_name= $seq->accession_number();
  
  unless($sf) { # make one
    $source_type=  $is_swiss ? $PROTEIN_TYPE 
        : $is_gene ? "eneg" # "gene"  # "region" # 
        : $is_rich ? $seq->molecule : $source_type;
    $sf = Bio::SeqFeature::Generic->direct_new();
    
    my $len = $seq->length(); $len=1 if($len<1); my $start = 1; ##$start= $len if ($len<1);
    my $loc= $seq->can('location') ? $seq->location() 
            :  new Bio::Location::Simple( -start => $start, -end => $len);
    $sf->location( $loc ); 
    $sf->primary_tag($source_type);
    $sf->source_tag($SOURCEID);
    $sf->seq_id( $seq_name);
    #? $sf->display_name($seq->id()); ## Name or Alias ?
    $sf->add_tag_value( Alias => $seq->id()); # unless id == accession
    $seq->add_SeqFeature($sf);
    ## $source_feat= $sf;
  }
  
  if ($sf->has_tag("chromosome")) {
    $source_type= "chromosome"; 
    my ($chrname) = $sf->get_tag_values("chromosome");
    ## PROBLEM with Name <> ID, RefName for Gbrowse; use Alias instead
    ## e.g. Mouse chr 19 has two IDs in NCBI genbank now
    $sf->add_tag_value( Alias => $chrname );
    }

  # pull GB Comment, Description for source ft ...
  # add reference - can be long, not plain string...             
  warn "# $SOURCEID:$seq_name fields = ", join(",", $seq->annotation->get_all_annotation_keys()),"\n" if $DEBUG;
  # GenBank   fields: keyword,comment,reference,date_changed
  # Entrezgene fields 850293 =ALIAS_SYMBOL,RefSeq status,chromosome,SGD,dblink,Entrez Gene Status,OntologyTerm,LOCUS_SYNONYM

    # is this just for main $seq object or for all seqfeatures ?
  my %AnnotTagMap= ( 
      'gene_name' => 'Alias',   
      'ALIAS_SYMBOL' => 'Alias',  # Entrezgene
      'LOCUS_SYNONYM' => 'Alias', #?
      'symbol' => 'Alias',   
      'synonym' => 'Alias',   
      'dblink' => 'Dbxref',   
      'product' => 'product',
      'Reference' => 'reference',
      'OntologyTerm' => 'Ontology_term',
      'comment'  => 'Note',
      'comment1' => 'Note',
      # various map-type locations
      # gene accession tag is named per source db !??
      # 'Index terms' => keywords ??
      );
  
  
  my ($desc)= $seq->annotation->get_Annotations("desc") || ( $seq->desc() );
  my ($date)= $seq->annotation->get_Annotations("dates") 
      || $seq->annotation->get_Annotations("update-date")
      || $is_rich ? $seq->get_dates() : ();
  my ($comment)= $seq->annotation->get_Annotations("comment");
  my ($species)= $seq->annotation->get_Annotations("species") 
               || ( $seq->can('species') ? $seq->species()->binomial() : undef );
               
  # update source feature with main GB fields
  $sf->add_tag_value( ID => $seq_name ) unless $sf->has_tag('ID');
  $sf->add_tag_value( Note => $desc ) if($desc && ! $sf->has_tag('Note'));
  $sf->add_tag_value( organism => $species ) if($species && ! $sf->has_tag('organism'));
  $sf->add_tag_value( comment1 => $comment ) if(!$is_swiss && $comment && ! $sf->has_tag('comment1'));
  $sf->add_tag_value( date => $date ) if($date && ! $sf->has_tag('date'));

  $sf->add_tag_value( Dbxref => $SOURCEID.':'.$seq_name ) if $is_swiss || $is_gene;

  foreach my $atag (sort keys %AnnotTagMap) {
    my $gtag= $AnnotTagMap{$atag}; next unless($gtag);
    my @anno = map{ 
          ref $_
         ? split( /[,;] */, $_->value) 
         : split( /[,;] */, "$_") if($_);
         } $seq->annotation->get_Annotations($atag);  
    foreach(@anno) { $sf->add_tag_value( $gtag => $_ ); }
    }
    
#   my @genes = map{ split( /[,;] */, "$_"); } $seq->annotation->get_Annotations('gene_name');  
#   $sf->add_tag_value( Alias => $_ ) foreach(@genes);
# 
#   my @dblink= map { "$_"; } $seq->annotation->get_Annotations("dblink"); # add @all 
#   $sf->add_tag_value( Dbxref => $_ ) foreach(@dblink);

  
  return (wantarray)? ($source_type,$sf) : $source_type; #?
}


sub gene_features {
    my ($f, $gene_id, $genelinkID) = @_;
    local $_ = $f->primary_tag;
    $method{$_}++;
    
    if ( /gene/ ) {
        $f->add_tag_value( ID => $gene_id ) unless($f->has_tag('ID')); # check is same value!? 
        $tnum = $rnum= 0; $ncrna_id= $rna_id  = '';
        return GM_NEW_TOPLEVEL;
        
    } elsif ( /mRNA/ ) {  
        return GM_NOT_PART unless $gene_id;
        return GM_NOT_PART if($genelinkID && $genelinkID ne $gene_id);
        ($rna_id  = $gene_id ) =~ s/gene/mRNA/;
        $rna_id   .= '.t0' . ++$tnum;
        $f->add_tag_value( ID => $rna_id );
        $f->add_tag_value( Parent => $gene_id );
        
    } elsif ( /RNA|transcript/) { 
      ## misc_RNA here; missing exons ... flattener problem?
      #  all of {t,nc,sn}RNA can have gene models now
      ## but problem in Worm chr: mRNA > misc_RNA > CDS with same locus tag
      ## CDS needs to use mRNA, not misc_RNA, rna_id ...
      ## also need to fix cases where tRNA,... lack a 'gene' parent: make this one top-level

        if($gene_id) {
          return GM_NOT_PART if($genelinkID && $genelinkID ne $gene_id);
          ($ncrna_id = $gene_id) =~ s/gene/ncRNA/;
          $ncrna_id .= '.r0' . ++$rnum;
          $f->add_tag_value( Parent => $gene_id );
          $f->add_tag_value( ID => $ncrna_id );
        } else {
          unless ($f->has_tag('ID')) {
            if($genelinkID) {
              $f->add_tag_value( ID => $genelinkID  ) ;
            } else {
              $idh->generate_unique_persistent_id($f);
            }
          }
         ($ncrna_id)= $f->get_tag_values('ID'); 
         return GM_NEW_TOPLEVEL;
          # this feat now acts as gene-top-level; need to print @to_print to flush prior exons?
        }
        
    } elsif ( /exon/ ) { # can belong to any kind of RNA
        return GM_NOT_PART unless ($rna_id||$ncrna_id);
        return GM_NOT_PART if($genelinkID && $genelinkID ne $gene_id);
        ## we are getting duplicate Parents here, which chokes chado loader, with reason...
        ## problem is when mRNA and ncRNA have same exons, both ids are active, called twice
        ## check all Parents
        for my $expar ($rna_id, $ncrna_id) { 
          next unless($expar);
          if ( $exonpar{$expar} and $f->has_tag('Parent') ) {
            my @vals = $f->get_tag_values('Parent');
            next if (grep {$expar eq $_} @vals);
            }
          $exonpar{$expar}++;
          $f->add_tag_value( Parent => $expar); 
          # last; #? could be both
        }
      # now we can skip cloned exons
      # dgg note: multiple parents get added and printed for each unique exon
      return GM_DUP_PART if ++$seen{$f} > 1;
        
    } elsif ( /CDS|protein|polypeptide/ ) {
        return GM_NOT_PART unless $rna_id; ## ignore $ncrna_id ??
        return GM_NOT_PART if($genelinkID && $genelinkID ne $gene_id); #??
        (my $pro_id = $rna_id) =~ s/\.t/\.p/;
        
        if( ! $CDSkeep && /CDS/) {  
          $f->primary_tag($PROTEIN_TYPE); 
        
          ## duplicate problem is Location ..
          if ($f->location->isa("Bio::Location::SplitLocationI")) {
            # my($b,$e)=($f->start, $f->end); # is this all we need?
            my($b,$e)=(-1,0);
            foreach my $l ($f->location->each_Location) {
               $b = $l->start if($b<0 || $b > $l->start);
               $e = $l->end if($e < $l->end);
               }
            $f->location( Bio::Location::Simple->new(
                -start => $b, -end => $e, -strand => $f->strand) );
            }

            $f->add_tag_value( Derives_from => $rna_id ); 
          }
          else {
            $f->add_tag_value( Parent => $rna_id );
          }

        $f->add_tag_value( ID => $pro_id );
        
        move_translation_fasta($f, $pro_id);
#         if( $f->has_tag('translation')) {
#           my ($aa) = $f->get_tag_values("translation");
#           $proteinfa{$pro_id}= $aa;
#           $f->remove_tag("translation");
#           $f->add_tag_value("translation","length.".length($aa)); # hack for odd chado gbl problem
#         }
        
    } else {
        return GM_NOT_PART unless $gene_id;
        $f->add_tag_value( Parent => $gene_id );  
    }
    
    ## return GM_DUP_PART if /exon/ && ++$seen{$f} > 1;
    
    return GM_NEW_PART;
}

## was generic_features >  add_generic_id
sub add_generic_id {
    my ($f, $ft_name, $flags) = @_;
    my $method = $f->primary_tag;
    $method{$method}++ unless($flags =~ /nocount/); ## double counts GM_NOT_PART from above
     
    if ($f->has_tag('ID')) {
    
    }
    elsif ( $f->has_tag($method) ) {
        my ($name) = $f->get_tag_values($method);
        $f->add_tag_value( ID => "$method:$name" );
    }
    elsif($ft_name) { # is this unique ?
        $f->add_tag_value( ID => $ft_name ); 
    }
    else {
        $idh->generate_unique_persistent_id($f);
    }

    move_translation_fasta( $f, ($f->get_tag_values("ID"))[0] )
      if($method =~ /CDS/);

    # return $io->gff_string($f);
}

sub move_translation_fasta {
    my ($f, $ft_id) = @_;
    if( $ft_id && $f->has_tag('translation') ) {
      my ($aa) = $f->get_tag_values("translation");
      if($aa && $aa !~ /^length/) {
        $proteinfa{$ft_id}= $aa;
        $f->remove_tag("translation");
        $f->add_tag_value("translation","length.".length($aa)); # hack for odd chado gbl problem
        }
    }
}

sub gff_header {
    my ($name, $end, $source_type, $source_feat) = @_;
    $source_type ||= "region";

    my $info= "$source_type:$name";
    my $head= "##gff-version $GFF_VERSION\n# sequence-region $name 1 $end\n";
    $head .= "# conversion-by bp_genbank2gff3.pl\n";
    if($source_feat) {
    ## dgg: these header comment fields are not useful when have multi-records, diff organisms
    for my $key (qw(organism date Note)) {
      my($value) = $source_feat->get_tag_values($key);
      $head .= "# $key $value\n" if($value);
      $info .= ", $value" if($value);
      }
      $head="" if($didheader);
    } else {
      $head .= "$name\t$SOURCEID\t$source_type\t1\t$end\t.\t.\t.\tID=$name\n";
    }
    $didheader++;
    return (wantarray) ? ($head,$info) : $head;
}

sub unflatten_seq {
    my $seq = shift;

    ## print "# working on $source_type:", $seq->accession, "\n"; 
    my $uh_oh = "Possible gene unflattening error with" .  $seq->accession_number .
                ": consult STDERR\n";
    
    eval {
        $unflattener->unflatten_seq( -seq => $seq, 
                                     -use_magic => 1 );
    };
    
    # deal with unflattening errors
    if ( $@ ) {
        warn $seq->accession_number . " Unflattening error:\n";
        warn "Details: $@\n";
        print "# ".$uh_oh;
    }

    return 0 if !$seq || !$seq->all_SeqFeatures;

    # map feature types to the sequence ontology
    ## $tm->map_types_to_SO( -seq => $seq );
    $tm->map_types( -seq => $seq, -type_map => $FTSOmap, -undefined => "region" ); #dgg
    
    1;
}

sub filter {
    my $seq = shift;
    ## return unless $filter;
    my @feats;
    my @sources; # dgg; pick source features here; only 1 always?
    if ($filter) {
      for my $f ( $seq->remove_SeqFeatures ) {
        my $m = $f->primary_tag;
        push @sources, $f if ($m eq 'source'); # dgg? but leave in @feats ?
        push @feats, $f unless $filter =~ /$m/i;
        }
      $seq->add_SeqFeature(@feats) if @feats;
    } else {
      for my $f ( $seq->get_SeqFeatures ){
        my $m = $f->primary_tag;
        push @sources, $f if ($m eq 'source'); # dgg? but leave in @feats ?
        }
    }

    return @sources;
}


# The default behaviour of Bio::FeatureHolderI:get_all_SeqFeatures
# changed to filter out cloned features.  We have to implement the old
# method.  These two subroutines were adapted from the v1.4 Bio::FeatureHolderI
sub get_all_SeqFeatures  {
    my $seq = shift;
    my @flatarr;

    foreach my $feat ( $seq->get_SeqFeatures ){
        push(@flatarr,$feat);
        _add_flattened_SeqFeatures(\@flatarr,$feat);
    }
    return @flatarr;
}

sub gene_name {
    my $g = shift;
    my $gene_id = ''; # zero it;
    
    if ($g->has_tag('gene')) {
        ($gene_id) = $g->get_tag_values('gene'); 
    }
    elsif ($g->has_tag('locus_tag')) {
        ($gene_id) = $g->get_tag_values('locus_tag');
    }
    elsif ($g->has_tag('ID')) { # for non-Genbank > Entrezgene
        ($gene_id) = $g->get_tag_values('ID');
    }

## See Unflattener comment:
       # on rare occasions, records will have no /gene or /locus_tag
       # but it WILL have /product tags. These serve the same purpose
       # for grouping. For an example, see AY763288 (also in t/data)
    # eg. product=tRNA-Asp ;  product=similar to crooked neck protein
    elsif ($g->has_tag('product')) {
        my ($name)= $g->get_tag_values('product');
        ($gene_id) = $name unless($name =~ / /); # a description not name
    }

  ## dgg; also handle transposon=xxxx ID/name
  # ID=GenBank:repeat_region:NC_004353:1278337:1281302;transposon=HeT-A{}1685;Dbxref=FLYBASE:FBti0059746
    elsif ($g->has_tag('transposon')) {
        my ($name)= $g->get_tag_values('transposon');
        ($gene_id) = $name unless($name =~ / /); # a description not name
    }
  
    return $gene_id;
}

# same list as gene_name .. change tag to generic Name
sub convert_to_name {
    my $g = shift;
    my $gene_id = ''; # zero it;
    
    if ($g->has_tag('gene')) {
        ($gene_id) = $g->get_tag_values('gene'); 
        $g->remove_tag('gene');
        $g->add_tag_value('Name', $gene_id);
    }
    elsif ($g->has_tag('locus_tag')) {
        ($gene_id) = $g->get_tag_values('locus_tag');
        $g->remove_tag('locus_tag');
        $g->add_tag_value('Name', $gene_id);
    }

    elsif ($g->has_tag('product')) {
        my ($name)= $g->get_tag_values('product');
        ($gene_id) = $name unless($name =~ / /); # a description not name
        ## $g->remove_tag('product');
        $g->add_tag_value('Name', $gene_id);
    }

    elsif ($g->has_tag('transposon')) {
        my ($name)= $g->get_tag_values('transposon');
        ($gene_id) = $name unless($name =~ / /); # a description not name
        ## $g->remove_tag('transposon');
        $g->add_tag_value('Name', $gene_id);
    }
    return $gene_id;
}


sub _add_flattened_SeqFeatures  {
    my ($arrayref,$feat) = @_;
    my @subs = ();

    if ($feat->isa("Bio::FeatureHolderI")) {
        @subs = $feat->get_SeqFeatures;
    } 
    elsif ($feat->isa("Bio::SeqFeatureI")) {
        @subs = $feat->sub_SeqFeature;
    }
    else {
        warn ref($feat)." is neither a FeatureHolderI nor a SeqFeatureI. ".
            "Don't know how to flatten.";
    }

    for my $sub (@subs) {
        push(@$arrayref,$sub);
        _add_flattened_SeqFeatures($arrayref,$sub);
    }

}

