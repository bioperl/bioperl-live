#!/usr/bin/perl



=pod

=head1 NAME 

bp_genbank2gff3.pl -- Genbank-E<gt>gbrowse-friendly GFF3

=head1 SYNOPSIS

  bp_genbank2gff3.pl [options] filename(s)

  # process a directory containing GenBank flatfiles
  perl bp_genbank2gff3.pl --dir path_to_files --zip

  # process a single file, ignore explicit exons and introns
  perl bp_genbank2gff3.pl --filter exon --filter intron file.gbk.gz

  # process a list of files 
  perl bp_genbank2gff3.pl *gbk.gz

  # process data from URL, with Chado GFF model (-noCDS), and pipe to database loader
  curl ftp://ftp.ncbi.nih.gov/genomes/Saccharomyces_cerevisiae/CHR_X/NC_001142.gbk \
  | perl bp_genbank2gff3.pl -noCDS -in stdin -out stdout \
  | perl gmod_bulk_load_gff3.pl -dbname mychado -organism fromdata

    Options:
        --noinfer  -r  don't infer exon/mRNA subfeatures
        --conf     -i  path to the curation configuration file that contains user preferences
                       for Genbank entries (must be YAML format)
                       (if --manual is passed without --ini, user will be prompted to 
                        create the file if any manual input is saved)
        --sofile  -l  path to to the so.obo file to use for feature type mapping
                       (--sofile live will download the latest online revision)
        --manual   -m  when trying to guess the proper SO term, if more than
                       one option matches the primary tag, the converter will 
                       wait for user input to choose the correct one
                       (only works with --sofile)
        --dir      -d  path to a list of genbank flatfiles
        --outdir   -o  location to write GFF files (can be 'stdout' or '-' for pipe)
        --zip      -z  compress GFF3 output files with gzip
        --summary  -s  print a summary of the features in each contig
        --filter   -x  genbank feature type(s) to ignore
        --split    -y  split output to separate GFF and fasta files for
                       each genbank record
        --nolump   -n  separate file for each reference sequence
                       (default is to lump all records together into one 
                       output file for each input file)
        --ethresh  -e  error threshold for unflattener
                       set this high (>2) to ignore all unflattener errors
        --[no]CDS  -c  Keep CDS-exons, or convert to alternate gene-RNA-protein-exon 
                       model. --CDS is default. Use --CDS to keep default GFF gene model, 
                       use --noCDS to convert to g-r-p-e.
        --format   -f  Input format (SeqIO types): GenBank, Swiss or Uniprot, EMBL work
                       (GenBank is default)
        --GFF_VERSION  3 is default, 2 and 2.5 and other Bio::Tools::GFF versions available
        --quiet        don't talk about what is being processed 
        --typesource   SO sequence type for source (e.g. chromosome; region; contig)
        --help     -h  display this message


=head1 DESCRIPTION

This script uses Bio::SeqFeature::Tools::Unflattener and
Bio::Tools::GFF to convert GenBank flatfiles to GFF3 with gene
containment hierarchies mapped for optimal display in gbrowse.

The input files are assumed to be gzipped GenBank flatfiles for refseq
contigs.  The files may contain multiple GenBank records.  Either a
single file or an entire directory can be processed.  By default, the
DNA sequence is embedded in the GFF but it can be saved into separate
fasta file with the --split(-y) option.

If an input file contains multiple records, the default behaviour is
to dump all GFF and sequence to a file of the same name (with .gff
appended).  Using the 'nolump' option will create a separate file for
each genbank record.  Using the 'split' option will create separate
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
    curl ftp://ncbigenomes/... | bp_genbank2gff3 --in stdin --out stdout | gff2chado ...

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
 Subject: genbank2gff3.pl
 Date: 2007-01-26 17:08:49 GMT
   -- bug fixed in genbank2gff3 for multi-record handling

This error is for a /trans_splice gene that is hard to handle, and
unflattner/genbank2 doesn't

 From: Chad Matsalla <chad <at> dieselwurks.com> 
 Subject: genbank2gff3.PLS and the unflatenner - Inconsistent order?
 Date: 2005-07-15 19:51:48 GMT 

=head1 AUTHOR 

Sheldon McKay (mckays@cshl.edu)

Copyright (c) 2004 Cold Spring Harbor Laboratory.

=head2 AUTHOR of hacks for GFF2Chado loading

Don Gilbert (gilbertd@indiana.edu)


=cut

use strict;
use warnings;

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
use List::Util qw(first);
use Bio::OntologyIO;
use YAML qw(Dump LoadFile DumpFile);
use File::Basename; 

use vars qw/$split @filter $zip $outdir $help $ethresh
            $ONTOLOGY %FEATURES %DESCENDANTS @RETURN $MANUAL @GFF_LINE_FEAT
            $CONF $YAML $TYPE_MAP $SYN_MAP $noinfer $SO_FILE 
            $file @files $dir $summary $nolump 
            $source_type %proteinfa %exonpar $didheader $verbose $DEBUG $GFF_VERSION 
            $gene_id $rna_id $tnum $ncrna_id $rnum %method %id %seen/;

use constant SO_URL => 
    'http://song.cvs.sourceforge.net/viewvc/*checkout*/song/ontology/so.obo';
use constant ALPHABET => [qw(a b c d e f g h i j k l m n o p q r s t u v w x y z)];
use constant ALPHABET_TO_NUMBER => {
    a => 0, b => 1, c => 2, d => 3, e => 4, f => 5, g => 6, h => 7, i => 8, 
    j => 9, k => 10, l => 11, m => 12, n => 13, o => 14, p => 15, q => 16,
    r => 17, s => 18, t => 19, u => 20, v => 21, w => 22, x => 23, y => 24,
    z => 25,
    };
use constant ALPHABET_DIVISOR => 26;
use constant GM_NEW_TOPLEVEL => 2;
use constant GM_NEW_PART => 1;
use constant GM_DUP_PART => 0;
use constant GM_NOT_PART => -1;

# Options cycle in multiples of 2 because of left side/right side pairing. 
# You can make this number odd, but displayed matches will still round up
use constant OPTION_CYCLE => 6; 

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
            'r|noinfer' => \$noinfer,
            'i|conf=s' => \$CONF,
            'sofile=s' => \$SO_FILE,
            'm|manual' => \$MANUAL,
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

#my $FTSOmap = $tm->FT_SO_map();
my $FTSOmap;
my $FTSOsynonyms;

if (defined($SO_FILE) && $SO_FILE eq 'live') {
    print "\nDownloading the latest SO file from ".SO_URL."\n\n";
    use LWP::UserAgent;
    my $ua = LWP::UserAgent->new(timeout => 30);
    my $request = HTTP::Request->new(GET => SO_URL);
    my $response = $ua->request($request);


    if ($response->status_line =~ /200/) {
        use File::Temp qw/ tempfile /;
        my ($fh, $fn) = tempfile();
        print $fh $response->content;
        $SO_FILE = $fn;
    } else {
        print "Couldn't download SO file online...skipping validation.\n" 
            . "HTTP Status was " . $response->status_line . "\n" 
            and undef $SO_FILE
    }
}

if ($SO_FILE) {


    my (%terms, %syn);

    my $parser = Bio::OntologyIO->new( -format => "obo", -file => $SO_FILE );
    $ONTOLOGY = $parser->next_ontology();

    for ($ONTOLOGY->get_all_terms) { 
        my $feat = $_;

        $terms{$feat->name} = $feat->name;
        #$terms{$feat->name} = $feat;

        my @syn = $_->each_synonym;

        push @{$syn{$_}}, $feat->name for @syn;
        #push @{$syn{$_}}, $feat for @syn;
    }

    $FTSOmap = \%terms;
    $FTSOsynonyms = \%syn;

    my %hardTerms = %{ $tm->FT_SO_map() };
    map { $FTSOmap->{$_} ||= $hardTerms{$_} } keys %hardTerms;

} else { 
    my %terms = %{ $tm->FT_SO_map() };
    while (my ($k,$v) = each %terms) {
        $FTSOmap->{$k} = ref($v) ? shift @$v : $v;
    }
}

$TYPE_MAP = $FTSOmap;
$SYN_MAP = $FTSOsynonyms;


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
        my ($vol,$dirs,$fileonly) = File::Spec->splitpath($file); 
        $lump = File::Spec->catfile($outdir, $fileonly.'.gff');
        ($outfa = $lump) =~ s/\.gff/\.fa/;
        open $lump_fh, ">$lump" or die "Could not create a lump outfile called ($lump) because ($!)\n";
        open $lumpfa_fh, ">$outfa" or die "Could not create a lump outfile called ($outfa) because ($!)\n";

    }
    
    # open input file, unzip if req'd
    if ($stdin) {
       *FH = *STDIN;
    } elsif ( $file =~ /\.gz/ ) {
        open FH, "gunzip -c $file |";
    }
    else {
        open FH, '<', $file;
    }

    my $in = Bio::SeqIO->new(-fh => \*FH, -format => $FORMAT, -debug=>$DEBUG);
    my $gffio = Bio::Tools::GFF->new( -noparse => 1, -gff_version => $GFF_VERSION );

    while ( my $seq = $in->next_seq() ) {
        my $seq_name = $seq->accession_number;
        my $end = $seq->length;
        my @to_print;

        # arrange disposition of GFF output
        $outfile = $lump || File::Spec->catfile($outdir, $seq_name.'.gff');
        my $out;

        if ( $lump ) {
            $outfile = $lump;
            $out = $lump_fh;
        }
        else {
            $outfile = File::Spec->catfile($outdir, $seq_name.'.gff');
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


        $FTSOmap->{'source'} = $source_type;
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

        @GFF_LINE_FEAT = ();
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
            convert_to_name($feature); 
            
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
                  my $gff = $gffio->gff_string($feature);
                  push @GFF_LINE_FEAT, $feature;
                  #print $out "$gff\n" if $gff;

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
                push @GFF_LINE_FEAT, $feature;
                #print $out "$gff\n" if $gff;
            }
        }

        # don't like doing this after others; do after each new gene id?
        @to_print= print_held($out, $gffio, \@to_print);

        gff_validate(@GFF_LINE_FEAT);

        for my $feature (@GFF_LINE_FEAT) {
            my $gff= $gffio->gff_string($feature);
            print $out "$gff\n" if $gff;
        }

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
        #if ( $summary ) {
        #      print "# Summary:\n# Feature\tCount\n# -------\t-----\n";
        #  
        #      for ( keys %method ) {
        #          print "# $_  $method{$_}\n";
        #      }
        #      print "# \n";
        #  }       
    
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
    close($lumpfa_fh) if $lumpfa_fh;
    if (!$split && $outfa && $lump_fh) {     
        print $lump_fh "##FASTA\n"; # GFF3 spec
        open $lumpfa_fh, $outfa or warn "reading FA $outfa: $!";
        while( <$lumpfa_fh>) {
            print $lump_fh $_;
        } # is $lump_fh still open?
        close($lumpfa_fh);
        unlink($outfa);
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
    push @GFF_LINE_FEAT, $feature;
    #print $out "$gff\n";
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
  my ($species)= $seq->annotation->get_Annotations("species");
  if (!$species 
       && $seq->can('species') 
       && defined $seq->species() 
       && $seq->species()->can('binomial') ) {
    $species= $seq->species()->binomial();
  }

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
          if (ref $_ && $_->can('get_all_values')) { 
              split( /[,;] */, join ";", $_->get_all_values) 
          }
          elsif (ref $_ && $_->can('display_text')) { 
              split( /[,;] */, $_->display_text) 
          }
          elsif (ref $_ && $_->can('value')) { 
              split( /[,;] */, $_->value) 
          } 
          else {
             ();
          }
      } $seq->annotation->get_Annotations($atag);  
      foreach(@anno) { $sf->add_tag_value( $gtag => $_ ); }
  }
    
  #my @genes = map{ split( /[,;] */, "$_"); } $seq->annotation->get_Annotations('gene_name');  
  #$sf->add_tag_value( Alias => $_ ) foreach(@genes);
  # 
  #my @dblink= map { "$_"; } $seq->annotation->get_Annotations("dblink"); # add @all 
  #$sf->add_tag_value( Dbxref => $_ ) foreach(@dblink);
  
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
        #if( $f->has_tag('translation')) {
        #   my ($aa) = $f->get_tag_values("translation");
        #    $proteinfa{$pro_id}= $aa;
        #    $f->remove_tag("translation");
        #    $f->add_tag_value("translation","length.".length($aa)); # hack for odd chado gbl problem
        #}
    } elsif ( /region/ ) {       
        $f->primary_tag('gene_component_region');
        $f->add_tag_value( Parent => $gene_id ); 
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

    my $info = "$source_type:$name";
    my $head = "##gff-version $GFF_VERSION\n".
               "##sequence-region $name 1 $end\n".
               "# conversion-by bp_genbank2gff3.pl\n";
    if ($source_feat) {
        ## dgg: these header comment fields are not useful when have multi-records, diff organisms
        for my $key (qw(organism Note date)) {
            my $value;
            if ($source_feat->has_tag($key)) { 
                ($value) = $source_feat->get_tag_values($key);
            }
            if ($value) {
                $head .= "# $key $value\n";
                $info .= ", $value";
            }
        }
        $head = "" if $didheader;
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
                                     -noinfer => $noinfer,
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
    #$tm->map_types( -seq => $seq, -type_map => $FTSOmap, -undefined => "region" ); #dgg

    map_types( 
        $tm, 
        -seq => $seq, 
        -type_map  => $FTSOmap, 
        -syn_map  => $FTSOsynonyms, 
        -undefined => "region" 
    ); #nml

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
      $seq->add_SeqFeature($_) foreach @feats;
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

    if ($g->has_tag('locus_tag')) {
        ($gene_id) = $g->get_tag_values('locus_tag');
    }

    elsif ($g->has_tag('gene')) {
        ($gene_id) = $g->get_tag_values('gene'); 
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
    elsif ($g->has_tag('ID')) { 
        my ($name)= $g->get_tag_values('ID');
        $g->add_tag_value('Name', $name);
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

sub map_types {

    my ($self, @args) = @_;

    my($sf, $seq, $type_map, $syn_map, $undefmap) =
        $self->_rearrange([qw(FEATURE
                    SEQ
                    TYPE_MAP
                    SYN_MAP
                    UNDEFINED
                    )],
                @args);

    if (!$sf && !$seq) {
        $self->throw("you need to pass in either -feature or -seq");
    }

    my @sfs = ($sf);
    if ($seq) {
        $seq->isa("Bio::SeqI") || $self->throw("$seq NOT A SeqI");
        @sfs = $seq->get_all_SeqFeatures;
    }
    $type_map = $type_map || $self->typemap; # dgg: was type_map;
    foreach my $feat (@sfs) {

        $feat->isa("Bio::SeqFeatureI") || $self->throw("$feat NOT A SeqFeatureI");
        $feat->isa("Bio::FeatureHolderI") || $self->throw("$feat NOT A FeatureHolderI");

        my $primary_tag = $feat->primary_tag;

        #if ($primary_tag =~ /^pseudo(.*)$/) {
        #    $primary_tag = $1;
        #    $feat->primary_tag($primary_tag);
        #}

        my $mtype = $type_map->{$primary_tag};
        if ($mtype) {
            if (ref($mtype)) {
                if (ref($mtype) eq 'ARRAY') {
                    my $soID;
                    ($mtype, $soID) = @$mtype;

                    if ($soID && ref($ONTOLOGY)) {
                        my ($term) = $ONTOLOGY->find_terms(-identifier => $soID);
                        $mtype = $term->name if $term;
                    } 
                    # if SO ID is undefined AND we have an ontology to search, we want to delete 
                    # the feature type hash entry in order to force a fuzzy search
                    elsif (! defined $soID && ref($ONTOLOGY)) {
                        undef $mtype;
                        delete $type_map->{$primary_tag};
                    } 
                    elsif ($undefmap && $mtype eq 'undefined') { # dgg
                        $mtype= $undefmap;
                    }

                    $type_map->{$primary_tag} = $mtype if $mtype;
                }
                elsif (ref($mtype) eq 'CODE') {
                    $mtype = $mtype->($feat);
                }
                else {
                    $self->throw('must be scalar or CODE ref');
                }
            }
            elsif ($undefmap && $mtype eq 'undefined') { # dgg
                $mtype= $undefmap;
            }
            $feat->primary_tag($mtype);
        }

        if ($CONF) {
            conf_read(); 
            my %perfect_matches;
            while (my ($p_tag,$rules) = each %$YAML) {
                RULE:
                for my $rule (@$rules) {
                    for my $tags (@$rule) {
                        while (my ($tag,$values) = each %$tags) {
                            for my $value (@$values) {
                                if ($feat->has_tag($tag)) {
                                    for ($feat->get_tag_values($tag)) {
                                        next RULE unless $_ =~ /\Q$value\E/;
                                    }
                                } elsif ($tag eq 'primary_tag') {
                                    next RULE unless $value eq
                                        $feat->primary_tag; 
                                } elsif ($tag eq 'location') {
                                    next RULE unless $value eq
                                        $feat->start.'..'.$feat->end;
                                } else { next RULE }
                            }
                        }
                    }
                    $perfect_matches{$p_tag}++;
                }
            }
            if (scalar(keys %perfect_matches) == 1) {
                $mtype = $_ for keys %perfect_matches;
            } elsif (scalar(keys %perfect_matches) > 1) {
                warn "There are conflicting rules in the config file for the" .
                     " following types: ";
                warn "\t$_\n" for keys %perfect_matches;
                warn "Until conflict resolution is built into the converter," .
                     " you will have to manually edit the config file to remove the" .
                     " conflict. Sorry :(. Skipping user preference for this entry";
                sleep(2);
            }
        } 

        if ( ! $mtype  && $syn_map) {
            if ($feat->has_tag('note')) {

                my @all_matches;

                my @note = $feat->each_tag_value('note');

                for my $k (keys %$syn_map) {

                    if ($k =~ /"(.+)"/) {

                        my $syn = $1;

                        for my $note (@note) {

                            # look through the notes to see if the description
                            # is an exact match for synonyms
                            if ( $syn eq $note ) { 

                                my @map = @{$syn_map->{$k}};

                                
                                my $best_guess = $map[0];

                                unshift @{$all_matches[-1]}, [$best_guess];

                                $mtype = $MANUAL
                                    ? manual_curation($feat, $best_guess, \@all_matches)
                                    : $best_guess;

                                print '#' x 78 . "\nGuessing the proper SO term for GenBank"
                                    . " entry:\n\n" . GenBank_entry($feat) . "\nis:\t$mtype\n" 
                                    . '#' x 78 . "\n\n";

                            } else {
                                # check both primary tag and and note against
                                # SO synonyms for best matching description

                                SO_fuzzy_match( $k, $primary_tag, $note, $syn, \@all_matches); 
                            }

                        }
                    } 
                }

                #unless ($mtype) {
                for my $note (@note) {
                    for my $name (values %$type_map) {
                    # check primary tag against SO names for best matching
                    # descriptions //NML also need to check against
                    # definition && camel case split terms

                        SO_fuzzy_match($name, $primary_tag, $note, $name, \@all_matches);
                    }
                }
                #}

                if (scalar(@all_matches) && !$mtype) {

                    my $top_matches = first { defined $_ } @{$all_matches[-1]}; 

                    my $best_guess = $top_matches->[0];



                    # if guess has quotes, it is a synonym term. we need to 
                    # look up the corresponding name term
                    # otherwise, guess is a name, so we can use it directly
                    if ($best_guess =~ /"(.+)"/) {

                        $best_guess = $syn_map->{$best_guess}->[0];

                    } 

                    @RETURN = @all_matches;
                    $mtype = $MANUAL
                        ? manual_curation($feat, $best_guess, \@all_matches)
                        : $best_guess;

                    print '#' x 78 . "\nGuessing the proper SO term for GenBank"
                        . " entry:\n\n" . GenBank_entry($feat) . "\nis:\t$mtype\n" 
                        . '#' x 78 . "\n\n";

                }
            }
            $mtype ||= $undefmap;
            $feat->primary_tag($mtype);
        } 
    }


}

sub SO_fuzzy_match {

    my $candidate = shift;
    my $primary_tag = shift;
    my $note = shift;
    my $SO_terms = shift;
    my $best_matches_ref = shift;
    my $modifier = shift; 

    $modifier ||= '';

        my @feat_terms;

    for ( split(" |_", $primary_tag) ) {
        #my @camelCase = /(?:[A-Z]|[a-z])(?:[A-Z]+|[a-z]*)(?=$|[A-Z])/g;
        my @camelCase = /(?:[A-Z]|[a-z])(?:[A-Z]+|[a-z]*)(?=$|[A-Z]|[;:.,])/g;
        push @feat_terms, @camelCase;
    }

    for ( split(" |_", $note) ) {
        #my @camelCase = /(?:[A-Z]|[a-z])(?:[A-Z]+|[a-z]*)(?=$|[A-Z])/g;
        #my @camelCase = /(?:[A-Z]|[a-z])(?:[A-Z]+|[a-z]*)(?=$|[A-Z]|[;:.,])/g;
        (my $word = $_) =~ s/[;:.,]//g;
        push @feat_terms, $word;
    }


    my @SO_terms = split(" |_", $SO_terms);

    # fuzzy match works on a simple point system. When 2 words match,
    # the $plus counter adds one. When they don't, the $minus counter adds
    # one. This is used to sort similar matches together. Better matches
    # are found at the end of the array, near the top.

    # NML: can we improve best match by using synonym tags
    # EXACT,RELATED,NARROW,BROAD?

    my ($plus, $minus) = (0, 0); 
    my %feat_terms;
    my %SO_terms;

    #unique terms
    map {$feat_terms{$_} = 1} @feat_terms;
    map {$SO_terms{$_} = 1} @SO_terms;

    for my $st (keys %SO_terms) {
        for my $ft (keys %feat_terms) {

            ($st =~ m/$modifier\Q$ft\E/) ? $plus++ : $minus++;

        }
    }

    push @{$$best_matches_ref[$plus][$minus]}, $candidate if $plus;

}

sub manual_curation {

    my ($feat, $default_opt,  $all_matches) = @_; 

    my @all_matches = @$all_matches;

    # convert all SO synonyms into names and filter
    # all matches into unique term list because
    # synonyms can map to multiple duplicate names

    my (@unique_SO_terms, %seen);
    for (reverse @all_matches) {
        for (@$_) {
            for (@$_) {
                #my @names;
                if ($_ =~ /"(.+)"/) {
                    for (@{$SYN_MAP->{$_}}) {
                        push @unique_SO_terms, $_ unless $seen{$_};
                        $seen{$_}++;
                    }
                } else { 
                    push @unique_SO_terms, $_ unless $seen{$_};
                    $seen{$_}++;
                }
            }
        }
    }

    my $s = scalar(@unique_SO_terms);

    my $choice = 0;

    my $more = 
        "[a]uto   : automatic input (selects best guess for remaining entries)\r" .
        "[f]ind   : search for other SO terms matching your query (e.g. f gene)\r" . 
        "[i]nput  : add a specific term\r" .
        "[r]eset  : reset to the beginning of matches\r" .
        "[s]kip   : skip this entry (selects best guess for this entry)\r"
    ;

    $more .= 
        "[n]ext   : view the next ".OPTION_CYCLE." terms\r" .
        "[p]rev   : view the previous ".OPTION_CYCLE." terms" if ($s > OPTION_CYCLE);

    my $msg = #"\n\n" . '-' x 156 . "\n"
         "The converter found $s possible matches for the following GenBank entry: ";

    my $directions = 
        "Type a number to select the SO term that best matches"
        . " the genbank entry, or use any of the following options:\r" . '_' x 76 . "\r$more"; 


    # lookup filtered list to pull out definitions
    my @options = map { 
        my $term = $_;
        my %term;
        for (['name', 'name'], ['def', 'definition'], ['synonym',
                'each_synonym']) {
            my ($label, $method) = @$_;
            $term{$label} = \@{[$term->$method]};
        }
        [++$choice, $_->name, ($_->definition || 'none'), \%term, $_->each_synonym ];
    }  map { $ONTOLOGY->find_terms(-name => $_) } @unique_SO_terms;


    my $option = options_cycle(0, OPTION_CYCLE, $msg, $feat, $directions,
            $default_opt, @options);

    if ($option eq 'skip') { return $default_opt 
    } elsif ($option eq 'auto') {
        $MANUAL = 0;
        return $default_opt;
    } else { return $option }

}

sub options_cycle {

    my ($start, $stop, $msg, $feat, $directions, $best_guess, @opt) = @_;

    #NML: really should only call GenBank_entry once. Will need to change
    #method to return array & shift off header
    my $entry = GenBank_entry($feat, "\r");

    my $total = scalar(@opt);

    ($start,$stop) = (0, OPTION_CYCLE) 
        if ( ($start < 0) && ($stop > 0) );

    ($start,$stop) = (0, OPTION_CYCLE) 
        if ( ( ($stop - $start) < OPTION_CYCLE ) && $stop < $total);

    ($start,$stop) = ($total - OPTION_CYCLE, $total) if $start < 0;
    ($start,$stop) = (0, OPTION_CYCLE) if $start >= $total;

    $stop = $total if $stop > $total;

    my $dir_copy = $directions;
    my $msg_copy = $msg;
    my $format = "format STDOUT = \n" .
        '-' x 156 . "\n" . 
        '^' . '<' x 77 .  '| Available Commands:' . "\n" .
        '$msg_copy' . "\n" .
        '-' x 156 . "\n" . 
        ' ' x 78 . "|\n" .
        '^' . '<' x 77 . '| ^' . '<' x 75 . '~' . "\n" .
        '$entry' . ' ' x 74 . '$dir_copy,' . "\n" .
        (' ' x 20 . '^' . '<' x 57 . '| ^' . '<' x 75 . '~' . "\n" .
        ' ' x 20 . '$entry,' . ' ' x 53 . '$dir_copy,' . "\n") x 1000  . ".\n";

    {
        # eval throws redefined warning that breaks formatting. 
        # Turning off warnings just for the eval to fix this.
        no warnings 'redefine';
        eval $format;
    }

    write;

    print '-' x 156 . "\n" .
        'Showing results ' . ( $stop ? ( $start + 1 ) : $start ) . 
        " - $stop of possible SO term matches: (best guess is \"$best_guess\")" .
        "\n" . '-' x 156 . "\n"; 

    for  (my $i = $start; $i < $stop; $i+=2) {

        my ($left, $right) = @opt[$i,$i+1];

        my ($nL, $nmL, $descL, $termL, @synL) = @$left;

        #odd numbered lists can cause fatal undefined errors, so check
        #to make sure we have data
        
        my ($nR, $nmR, $descR, $termR, @synR) = ref($right) ? @$right : (undef, undef, undef);


        my $format = "format STDOUT = \n";

        $format .=
            ' ' x 78 . "|\n" .

            '@>>: name: ^' . '<' x 64 . '~' . ' |' .
                ( ref($right) ? ('@>>: name: ^' . '<' x 64 . '~' ) : '' ) .  "\n" .
            '$nL,' . ' ' x 7 . '$nmL,' .
                ( ref($right) ? (' ' x 63 . '$nR,' .  ' ' x 7 .  "\$nmR,") : '' ) . "\n" .

            ' ' x 11 . '^' . '<' x 61 . '...~' . ' |' . 
                (ref($right) ? ('           ^' . '<' x 61 .  '...~') : '') . "\n" .
            ' ' x 11 . '$nmL,' . 
                (ref($right) ? (' ' x 74 . '$nmR,') : '') . "\n" .
            #' ' x 78 . '|' . "\n" .


            '     def:  ^' . '<' x 65 . ' |' . 
                (ref($right) ? ('     def:  ^' . '<' x 64 . '~') : '') . "\n" .
            ' ' x 11 . '$descL,' . 
                (ref($right) ? (' ' x 72 . '$descR,') : '') . "\n" .


            ('           ^' . '<' x 65 . ' |' . 
                (ref($right) ? ('           ^' . '<' x 64 . '~') : '') . "\n" .
            ' ' x 11 . '$descL,' . 
                (ref($right) ? (' ' x 72 . '$descR,') : '') . "\n") x 5 .


            '           ^' . '<' x 61 . '...~ |' . 
                (ref($right) ? ('           ^' . '<' x 61 . '...~') : '') . "\n" .
            ' ' x 11 . '$descL,' . 
                (ref($right) ? (' ' x 72 . '$descR,') : '') . "\n" .

            ".\n";

        {
            # eval throws redefined warning that breaks formatting. 
            # Turning off warnings just for the eval to fix this.
            no warnings 'redefine';
            eval $format;
        }
        write;

    }   
    print '-' x 156 . "\nenter a command:";

    while (<STDIN>) {

        (my $input = $_) =~ s/\s+$//;

        if ($input =~ /^\d+$/) {
            if ( $input && defined $opt[$input-1] ) {
                return $opt[$input-1]->[1]
            } else {
                print "\nThat number is not an option. Please enter a valid number.\n:";
            }
        } elsif ($input =~ /^n/i | $input =~ /next/i ) {
            return options_cycle($start + OPTION_CYCLE, $stop + OPTION_CYCLE, $msg, 
                    $feat, $directions, $best_guess, @opt)
        } elsif ($input =~ /^p/i | $input =~ /prev/i ) {
            return options_cycle($start - OPTION_CYCLE, $stop - OPTION_CYCLE, $msg,
                    $feat, $directions, $best_guess, @opt)
        } elsif ( $input =~ /^s/i || $input =~ /skip/i ) { return 'skip' 
        } elsif ( $input =~ /^a/i || $input =~ /auto/i ) { return 'auto' 
        } elsif ( $input =~ /^r/i || $input =~ /reset/i ) { 
            return manual_curation($feat, $best_guess, \@RETURN );
        } elsif ( $input =~ /^f/i || $input =~ /find/i ) {

            my ($query, @query_results);

            if ($input =~ /(?:^f|find)\s+?(.*?)$/) { $query = $1;
            } else {

                #do a SO search
                print "Type your search query\n:";
                while (<STDIN>) { chomp($query = $_); last }
            }

            for (keys(%$TYPE_MAP), keys(%$SYN_MAP)) {
                SO_fuzzy_match($_, $query, '', $_, \@query_results, '(?i)');
            }

            return manual_curation($feat, $best_guess, \@query_results);

        } elsif ( $input =~ /^i/i || $input =~ /input/i ) {

            #NML fast input for later
            #my $query;
            #if ($input =~ /(?:^i|input)\s+?(.*?)$/) { $query = $1 };

            #manual input
            print "Type the term you want to use\n:";
            while (<STDIN>) {
                chomp(my $input = $_);

                if (! $TYPE_MAP->{$input}) {

                    print "\"$input\" doesn't appear to be a valid SO term. Are ".
                        "you sure you want to use it? (y or n)\n:";

                    while (<STDIN>) {

                        chomp(my $choice = $_);

                        if ($choice eq 'y') {
                            print 
                                "\nWould you like to save your preference for " .
                                "future use (so you don't have to redo manual " .
                                "curation for this feature everytime you run " . 
                                "the converter)? (y or n)\n";

                            #NML: all these while loops are a mess. Really should condense it.
                            while (<STDIN>) {

                                chomp(my $choice = $_);

                                if ($choice eq 'y') {
                                    curation_save($feat, $input);
                                    return $input;
                                } elsif ($choice eq 'n') {
                                    return $input
                                } else {
                                    print "\nDidn't recognize that command. Please " . 
                                        "type y or n.\n:" 
                                }
                            }

                                
                        } elsif ($choice eq 'n') {
                            return options_cycle($start, $stop, $msg, $feat,
                                    $directions, $best_guess, @opt)
                        } else {
                            print "\nDidn't recognize that command. Please " . 
                                "type y or n.\n:" 
                        }
                    }

                } else { 
                    print 
                        "\nWould you like to save your preference for " .
                        "future use (so you don't have to redo manual " .
                        "curation for this feature everytime you run  " . 
                        "the converter)? (y or n)\n";

                    #NML: all these while loops are a mess. Really should condense it.
                    while (<STDIN>) {

                        chomp(my $choice = $_);

                        if ($choice eq 'y') {
                            curation_save($feat, $input);
                            return $input;
                        } elsif ($choice eq 'n') {
                            return $input
                        } else {
                            print "\nDidn't recognize that command. Please " . 
                                "type y or n.\n:" 
                        }
                    }

                } 

            }
        } else { 
            print "\nDidn't recognize that command. Please re-enter your choice.\n:" 
        }
    }

}

sub GenBank_entry {
    my ($f, $delimiter, $num) = @_;

    $delimiter ||= "\n";


    my $entry  = 

        ($num ? ' [1] ' : ' ' x 5) . $f->primary_tag 
        . ($num 
            ? ' ' x (12 - length $f->primary_tag ) . ' [2] '
            : ' ' x (15 - length $f->primary_tag) 
          )
        . $f->start.'..'.$f->end

        . "$delimiter";

    if ($num) {
        words_tag($f, \$entry);
    } else {
        for my $tag ($f->all_tags) {
            for my $val ( $f->each_tag_value($tag) ) {
                $entry .= ' ' x 20;
                #$entry .= "/$tag=\"$val\"$delimiter";
                $entry .= $val eq '_no_value'
                    ? "/$tag$delimiter"
                    : "/$tag=\"$val\"$delimiter";
            }
        }

    }

    return $entry;
}


sub gff_validate {
    warn "Validating GFF file\n" if $DEBUG;
    my @feat = @_;

    my (%parent2child, %all_ids, %descendants, %reserved);

    for my $f (@feat) {
        for my $aTags (['Parent', \%parent2child], ['ID', \%all_ids]) {
            map { push @{$$aTags[1]->{$_}}, $f } $f->get_tag_values($$aTags[0])
                if $f->has_tag($$aTags[0]); 
        }
    }

    if ($SO_FILE) {
        while (my ($parentID, $aChildren) = each %parent2child) {
            parent_validate($parentID, $aChildren, \%all_ids, \%descendants, \%reserved);
        }
    }

    id_validate(\%all_ids, \%reserved);        
}

sub parent_validate {
    my ($parentID, $aChildren, $hAll, $hDescendants, $hReserved) = @_;

    my $aParents = $hAll->{$parentID};

    map { 
        my $child = $_;
        $child->add_tag_value( validation_error => 
        "feature tried to add Parent tag, but no Parent found with ID $parentID"
        );
        my %parents;
        map { $parents{$_} = 1 } $child->get_tag_values('Parent');
        delete $parents{$parentID};
        my @parents = keys %parents;

        $child->remove_tag('Parent');

        unless ($child->has_tag('ID')) {
            my $id = gene_name($child);
            $child->add_tag_value('ID', $id);
            push @{$hAll->{$id}}, $child
        }

        $child->add_tag_value('Parent', @parents) if @parents;

    } @$aChildren and return unless scalar(@$aParents);

    my $par = join(',', map { $_->primary_tag } @$aParents);
    warn scalar(@$aParents)." POSSIBLE PARENT(S): $par" if $DEBUG;

    #NML manual curation needs to happen here


    my %parentsToRemove;

    CHILD:
    for my $child (@$aChildren) {
        my $childType  = $child->primary_tag;

        warn "WORKING ON $childType at ".$child->start.' to '.$child->end 
            if $DEBUG;

        for (my $i = 0; $i < scalar(@$aParents); $i++) {
            my $parent = $aParents->[$i];
            my $parentType = $parent->primary_tag;

            warn "CHECKING $childType against $parentType" if $DEBUG;

            #cache descendants so we don't have to do repeat searches
            unless ($hDescendants->{$parentType}) {

                for my $term ($ONTOLOGY->find_terms(
                        -name => $parentType
                    ) ) {
                    
                    map {
                        $hDescendants->{$parentType}{$_->name}++
                    } $ONTOLOGY->get_descendant_terms($term);

                }

                # NML: hopefully temporary fix.
                # SO doesn't consider exon/CDS to be a child of mRNA
                # even though common knowledge dictates that they are
                # This cheat fixes the false positives for now
                if ($parentType eq 'mRNA') {
                    $hDescendants->{$parentType}{'exon'} = 1;
                    $hDescendants->{$parentType}{'CDS'} = 1;
                }

            }

            warn "\tCAN $childType at " . $child->start . ' to ' . $child->end .
                " be a child of $parentType?" if $DEBUG;
            if ($hDescendants->{$parentType}{$childType}) {
                warn "\tYES, $childType can be a child of $parentType" if $DEBUG;

                #NML need to deal with multiple children matched to multiple different
                #parents. This model only assumes the first parent id that matches a child will
                #be the reserved feature. 

                $hReserved->{$parentID}{$parent}{'parent'} = $parent;
                push @{$hReserved->{$parentID}{$parent}{'children'}}, $child;

                #mark parent for later removal from all IDs 
                #so we don't accidentally change any parents

                $parentsToRemove{$i}++;

                next CHILD;
            } 
        }


        
        #NML shouldn't have to check this; somehow child can lose Parent
        #it's happening W3110
        #need to track this down
        if ( $child->has_tag('Parent') ) {

            warn "\tNO, @{[$child->primary_tag]} cannot be a child of $parentID"
                if $DEBUG;

            my %parents;

            map { $parents{$_} = 1 } $child->get_tag_values('Parent');

            delete $parents{$parentID};
            my @parents = keys %parents;

            warn 'VALIDATION ERROR '.$child->primary_tag." at ".$child->start .
                ' to ' . $child->end . " cannot be a child of ID $parentID"
                if $DEBUG;

            $child->add_tag_value( validation_error => 
                    "feature cannot be a child of $parentID");

            $child->remove_tag('Parent');

            unless ($child->has_tag('ID')) {
                my $id = gene_name($child);
                $child->add_tag_value('ID', $id);
                push @{$hAll->{$id}}, $child
            }

            $child->add_tag_value('Parent', @parents) if @parents;
        }

    }
    
    #delete $aParents->[$_] for keys %parentsToRemove;
    splice(@$aParents, $_, 1) for keys %parentsToRemove;
}

sub id_validate {
    my ($hAll, $hReserved) = @_;


    for my $id (keys %$hAll) {

        #since 1 feature can have this id, 
        #let's just shift it off and uniquify
        #the rest (unless it's reserved)

        shift @{$hAll->{$id}} unless $hReserved->{$id};
        for my $feat (@{$hAll->{$id}}) {
            id_uniquify(0, $id, $feat, $hAll);
        }
    }

    for my $parentID (keys %$hReserved) {

        my @keys = keys %{$hReserved->{$parentID}};

        shift @keys;

        for my $k (@keys) {
            my $parent = $hReserved->{$parentID}{$k}{'parent'};
            my $aChildren= $hReserved->{$parentID}{$k}{'children'};

            my $value = id_uniquify(0, $parentID, $parent, $hAll);
            for my $child (@$aChildren) {

                my %parents;
                map { $parents{$_}++ } $child->get_tag_values('Parent');
                $child->remove_tag('Parent');
                delete $parents{$parentID};
                $parents{$value}++;
                my @parents = keys %parents;
                $child->add_tag_value('Parent', @parents);
            }

        }
    }
}

sub id_uniquify {
    my ($count, $value, $feat, $hAll) = @_;

    warn "UNIQUIFYING $value" if $DEBUG;

    if (! $count) {
        $feat->add_tag_value(Alias => $value);
        $value .= ('.' . $feat->primary_tag) 
    } elsif ($count == 1) {
        $value .= ".$count" 
    } else { 
        chop $value;
        $value .= $count 
    }
    $count++;

    warn "ENDED UP WITH $value" if $DEBUG;
    if ( $hAll->{$value} ) { 
        warn "$value IS ALREADY TAKEN" if $DEBUG;
        id_uniquify($count, $value, $feat, $hAll);
    } else { 
        #warn "something's breaking ".$feat->primary_tag.' at '.$feat->start.' to '.$feat->end;
        $feat->remove_tag('ID');
        $feat->add_tag_value('ID', $value);
        push @{$hAll->{$value}}, $value;
    }

    $value;
}

sub conf_read {

    print "\nCannot read $CONF. Change file permissions and retry, " .
        "or enter another file\n" and conf_locate() unless -r $CONF;

    print "\nCannot write $CONF. Change file permissions and retry, " .
        "or enter another file\n" and conf_locate() unless -w $CONF;

    $YAML = LoadFile($CONF);

}

sub conf_create {

    my ($path, $input) = @_;

    print "Cannot write to $path. Change directory permissions and retry " .
        "or enter another save path\n" and conf_locate() unless -w $path;

    $CONF = $input;

    open(FH, '>', $CONF);
    close(FH);
    conf_read();


}

sub conf_write { DumpFile($CONF, $YAML) }

sub conf_locate {

    print "\nEnter the location of a previously saved config, or a new " .
        "path and file name to create a new config (this step is " .
        "necessary to save any preferences)";

    print "\n\nenter a command:";

    while (<STDIN>) {
        chomp(my $input = $_);
        my ($fn, $path, $suffix) = fileparse($input, qr/\.[^.]*/);

        if (-e $input && (! -d $input)) {

            print "\nReading $input...\n";
            $CONF = $input;

            conf_read(); 
            last;

        } elsif (! -d $input && $fn.$suffix) {

            print "Creating $input...\n";
            conf_create($path, $input);
            last;

        } elsif (-e $input && -d $input) {
            print "You only entered a directory. " .
                "Please enter BOTH a directory and filename\n";
        } else { 
            print "$input does not appear to be a valid path. Please enter a " .
                "valid directory and filename\n";
        }
        print "\nenter a command:";
    }
}

sub curation_save {

    my ($feat, $input) = @_;

    #my $error = "Enter the location of a previously saved config, or a new " .
    #    "path and file name to create a new config (this step is " .
    #    "necessary to save any preferences)\n";

    if (!$CONF) {
        print "\n\n"; 
        conf_locate();
    } elsif (! -e $CONF) {
        print "\n\nThe config file you have chosen doesn't exist.\n";
        conf_locate();
    } else { conf_read() }

    my $entry = GenBank_entry($feat, "\r", 1);

    my $msg = "Term entered: $input";
    my $directions = "Please select any/all tags that provide evidence for the term you
have entered. You may enter multiple tags by separating them by commas/dashes
(e.g 1,3,5-7). For tags with more than one word value (i.e 'note'), you have
the option of either selecting the entire note as evidence, or specific
keywords. If a tag has multiple keywords, they will be tagged alphabetically
for selection. To select a specific keyword in a tag field, you must enter the
tag number followed by the keyword letter (e.g 3a). Multiple keywords may be
selected by entering each letter separated by commas/dashes (e.g 3b,f,4a-c). The more tags you select, the more specific the GenBank entry will have
to be to match your curation. To match the GenBank entry exactly as it
appears, type every number (start-end), or just type 'all'. Remember, once the converter saves your
preference, you will no longer be prompted to choose a feature type for any
matching entries until you edit the curation.ini file.";
    my $msg_copy = $msg;
    my $dir_copy = $directions;

    my $format = "format STDOUT = \n" .
        '-' x 156 . "\n" . 
        '^' . '<' x 77 .  '| Directions:' . "\n" .
        '$msg_copy' . "\n" .
        '-' x 156 . "\n" . 
        ' ' x 78 . "|\n" .
        '^' . '<' x 77 . '| ^' . '<' x 75 . '~' . "\n" .
        '$entry' . ' ' x 74 . '$dir_copy,' . "\n" .
        (' ' x 15 . '^' . '<' x 62 . '| ^' . '<' x 75 . '~' . "\n" .
        ' ' x 15 . '$entry,' . ' ' x 58 . '$dir_copy,' . "\n") x 20  . ".\n";

    {
        # eval throws redefined warning that breaks formatting. 
        # Turning off warnings just for the eval to fix this.
        no warnings 'redefine';
        eval $format;
    }

    write;
    print '-' x 156 . "\nenter a command:";

    my @tags = words_tag($feat); 

    my $final = {};
    my $choices;
    while (<STDIN>) {

        chomp(my $choice = $_);

        if (scalar(keys %$final) && $choice =~ /^y/i) { last 

        } elsif (scalar(keys %$final) && $choice =~ /^n/i) { curation_save($feat, $input) 

        } elsif (scalar(keys %$final)) { print "\nInvalid selection. Please try again\n";

        } elsif ($choice eq 'all') {

            $choice = '';
            for (my $i=1; $i < scalar(@tags); $i++) {
                $choice .= "$i,";
            }
            chop $choice;
        } 
        #print "CHOICE [$choice]";

        my @selections;
        for (split(/(?<=\w)[^[:alnum:]\-]+(?=\d)/, $choice)) {
            if ($_ =~ /(\d+)(?:\D*)-(\d+)(.*)/) { 

                for ($1..$2) { push @selections, $_ }

                my $dangling_alphas = $3;
                alpha_expand($dangling_alphas, \@selections);

            } else { 
                alpha_expand($_, \@selections);
            }
        }

        foreach my $numbers (@selections) {

            my @c = split(/(?=[\w])/, $numbers);
            s/\W+//g foreach @c;
            my $num;
            
            {
                $^W = 0;
                $num = 0 + shift @c;
            }

            my $tag = $tags[$num];
            if (ref $tag && scalar(@c)) {
                my $no_value;
                foreach (@c) {
                    if (defined $tag->{$_}) {
                        $choices .= "${num}[$_] ";
                        my ($t,$v) = @{$tag->{$_}};
                        push @{${$final->{$input}}[0]{$t}}, $v;

                    } else { $no_value++ }
                    }

                if ($no_value) { 
                    _selection_add($tag,$final,$input,\$choices,$num);
                    #my ($t,$v) = @{$tag->{'all'}};
                    #unless (defined ${$final->{$input}}[0]{$t}) {
                        #$choices .= "$num, ";
                        #push @{${$final->{$input}}[0]{$t}}, $v
                    #}
                }

                $choices = substr($choices, 0, -2);
                $choices .= ', ';

            } elsif (ref $tag) { 
                _selection_add($tag,$final,$input,\$choices,$num);
                #my ($t,$v) = @{$tag->{'all'}};
                #unless (defined ${$final->{$input}}[0]{$t}) {
                    #$choices .= "$num, ";
                    #push @{${$final->{$input}}[0]{$t}}, $v
                #}
            } 
        }
        $choices = substr($choices, 0, -2) if $choices;
        if ($final) {
            print "\nYou have chosen the following tags:\n$choices\n";
            print "This will be written to the config file as:\n";
            print Dump $final;
            print "\nIs this correct? (y or n)\n";
        } else { print "\nInvalid selection. Please try again\n" }
    }
    push @{$YAML->{$input}}, $final->{$input};
    conf_write();
}

#  words_tag() splits each tag value string into multiple words so that the 
#  user can select the parts he/she wants to use for curation
#  it can tag 702 (a - zz) separate words; this should be enough

sub words_tag {

    my ($feat, $entry) = @_;

    my @tags;

    @tags[1,2] = ({'all' => ['primary_tag', $feat->primary_tag]}, {'all' => ['location', $feat->start.'..'.$feat->end]});
    my $i = 3;
    foreach my $tag ($feat->all_tags) {
        foreach my $value ($feat->each_tag_value($tag)) {

            my ($string, $tagged_string);

            my @words = split(/(?=\w+?)/, $value);

            my $pos = 0;


            foreach my $word (@words) {

                (my $sanitized_word = $word) =~ s/\W+?//g;
                $string .= $word;

                my $lead = int($pos/ALPHABET_DIVISOR);
                my $lag = $pos % ALPHABET_DIVISOR;

                my $a =  $lead ? ${(ALPHABET)}[$lead-1] : '';
                $a .= $lag ? ${(ALPHABET)}[$lag] : 'a';

                $tagged_string .= " ($a) $word";

                $tags[$i]{$a} = [$tag, $sanitized_word];
                $pos++;
            }

            $value = $tagged_string if scalar(@words) > 1;

            $$entry .= "[$i] /$tag=\"$value\"\r";

            $tags[$i]{'all'} = [$tag, $string];
        }
        $i++;
    }

    return @tags;
    
}

sub alpha_expand {

    my ($dangling_alphas, $selections) = @_;

    if (defined($dangling_alphas) && $dangling_alphas =~ /(\d)*([[:alpha:]]+)-([[:alpha:]]+)/) {

        my $digit = $1;
        push @$selections, $digit if $digit;

        my $start = $2;
        my $stop = $3;

        my @starts = split('', $start);
        my @stops = split('', $stop);

        my ($final_start, $final_stop);

        for ([\$final_start, \@starts], [\$final_stop, \@stops]) {

            my ($final, $splits) = @$_;

            my $int = ${(ALPHABET_TO_NUMBER)}{$$splits[0]};
            my $rem;


            if ($$splits[1]) {
                $rem = ${(ALPHABET_TO_NUMBER)}{$$splits[1]};
                $int++
            } else {
                $rem = $int;
                $int = 0;
            }


            $$final = $int * ALPHABET_DIVISOR;
            $$final += $rem;

        }

        my $last_number = pop @$selections;
        for my $pos ($final_start..$final_stop) {
            my $lead = int($pos/ALPHABET_DIVISOR);
            my $lag = $pos % ALPHABET_DIVISOR;
            my $alpha =  $lead ? ${(ALPHABET)}[$lead-1] : '';
            $alpha .= $lag ? ${(ALPHABET)}[$lag] : 'a';
            push @$selections, $last_number.$alpha;        
        }

    } elsif (defined($dangling_alphas)) { 
        if ($dangling_alphas =~ /^\d/) {
            push @$selections, $dangling_alphas;
        } elsif ($dangling_alphas =~ /^\D/) {
            #print "$dangling_alphas ".Dumper @$selections;
            my $last_number = pop @$selections;
            $last_number ||= '';
            push @$selections, $last_number.$dangling_alphas;
            #$$selections[-1] .= $dangling_alphas;
        }
    }

}

sub _selection_add {

    my ($tag, $final, $input, $choices, $num) = @_;
    my ($t,$v) = @{$tag->{'all'}};
    unless (defined ${$final->{$input}}[0]{$t}) {
        $$choices .= "$num, ";
        push @{${$final->{$input}}[0]{$t}}, $v
    }

}
