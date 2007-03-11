#!/usr/bin/perl -w
#$Id: genbank2gff3.PLS,v 1.10 2007-03-11 15:13:44 bosborne Exp $;

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


    Options:
        --dir     -d  path to a list of genbank flatfiles
        --outdir  -o  location to write GFF files
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

=head3 Note1:

In cases where the input files contain many GenBank records (for
example, the chromosome files for the mouse genome build), a very
large number of output files will be produced if the 'split' or
'nolump' options are selected.  If you do have lists of files E<gt> 6000,
use the --long_list option in bp_bulk_load_gff.pl or
bp_fast_load_gff.pl to load the gff and/ or fasta files.

=head3 Note2:

This script is designed for refseq genomic sequence entries.  It may
work for third party annotations but this has not been tested.

=head1 AUTHOR 

Sheldon McKay (mckays@cshl.edu)

Copyright (c) 2004 Cold Spring Harbor Laboratory.

=cut

use strict;

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
use Bio::Tools::GFF;
use Getopt::Long;

use vars qw/$split @filter $zip $outdir $help $ethresh
            $file @files $dir $summary $nolump
            $gene_id $rna_id $tnum %method %id %seen/;

$| = 1;

GetOptions( 'd|dir:s'   => \$dir,
	    'z|zip'     => \$zip, 
	    'h|help'    => \$help,
	    's|summary' => \$summary,
	    'o|outdir:s'=> \$outdir,
	    'x|filter:s'=> \@filter,
	    'y|split'   => \$split,
            "ethresh|e=s"=>\$ethresh,
            'n|nolump'  => \$nolump);

my $lump = 1 unless $nolump || $split;

# look for help request
pod2usage(2) if $help;

# initialize handlers
my $unflattener = Bio::SeqFeature::Tools::Unflattener->new;
$unflattener->error_threshold($ethresh) if $ethresh;
my $tm  = Bio::SeqFeature::Tools::TypeMapper->new;
my $idh = Bio::SeqFeature::Tools::IDHandler->new;

# stringify filter list if applicable
my $filter = join ' ', @filter  if @filter;

# determine input files
if ( $file ) {
    -e $file or die "file $file does not exist\n";

    if ( $file =~ m|/|) {
	($dir) = $file =~ m|(\S+)/\S+$|;
    } 

    $dir ||= '.';
    @files = ($file);
}
elsif ( $dir ) {
    if ( -d $dir ) {
	opendir DIR, $dir or die "could not open $dir for reading: $!";
	@files = grep { /\.gb.*/ } readdir DIR;
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

if ( $outdir && !-e $outdir ) {
    mkdir($outdir) or die "could not create directory $outdir: $!\n";        
}
elsif ( !$outdir ) {
    $outdir = $dir || '.';
}

$outdir .= '/' unless $outdir =~ m|/$|;

for my $file ( @files ) {
    chomp $file;
    die "$! $file" unless -e $file;
    print "Processing file $file...\n";

    my $lump_fh;
    if ( $lump ) {
		# this really doesn't do what you think it does.
        # ($lump) = $file =~ /^(\S+?)\./;
		# this is better, but still should use catfile
         
	  my ($vol,$dirs,$fileonly) = File::Spec->splitpath($file); 
        $lump   = $outdir . $fileonly . '.gff';
	open $lump_fh, ">$lump" or die "Could not create a lump outfile called ($lump) because ($!)\n";

    }
    
    my ($outfile, $outfa);
    
    # open input file, unzip if req'd
    if ( $file =~ /\.gz/ ) {
	open FH, "gunzip -c $file |";
    }
    else {
	open FH, "<$file";
    }

    my $in = Bio::SeqIO->new(-fh => \*FH, -format => 'GenBank');
    my $gffio = Bio::Tools::GFF->new( -noparse => 1, -gff_version => 3 );

    while ( my $seq = $in->next_seq ) {
	my $seq_name = $seq->accession;
	my $end = $seq->length;
	my @to_print;

        # arrange disposition of GFF output
        $outfile = $lump || $outdir . $seq->accession . ".gff";
	my $out;

	if ( $lump ) {
	    $outfile = $lump;
	    $out = $lump_fh;
	}
	else {
	    $outfile = $outdir . $seq->accession . ".gff";
	    open $out, ">$outfile";
	}

        # filter out unwanted features
        filter($seq);

	# abort if there are no features
        warn "$seq_name has no features, skipping\n" and next
	    if !$seq->all_SeqFeatures;

        # unflatten gene graphs, apply SO types, etc
        unflatten_seq($seq);

        # construct a GFF header
        print $out &gff_header($seq_name, $end);

	# Note that we use our own get_all_SeqFeatures function 
        # to rescue cloned exons
	for my $feature ( get_all_SeqFeatures($seq) ) {
	    
	    $feature->source_tag('GenBank');
	    my $method = $feature->primary_tag;
	    
	    # current gene name.  The unflattened gene features should be in order so any
            # exons, CDSs, etc that follow will belong to this gene
	    if ( $method eq 'gene' ) {
		gene_name($feature);
	    }

	    if ( $feature->has_tag('gene') || $method =~ /CDS|exon|RNA|UTR|gene/ ) {
		my $unique = gene_features($feature);
		push @to_print, $feature if $unique;
	    }
	    
	    # otherwise handle as generic feats with IDHandler labels 
	    else {
		my $gff = generic_features($feature,$gffio,$seq_name);
		print $out "$gff\n" if $gff;
	    }
	}

        for my $printme ( @to_print ) {
	    my $gff = $gffio->gff_string($printme);
            print $out "$gff\n";
        }

	# deal with the corresponding DNA
	my $dna = $seq->seq;
	$dna    =~ s/(\S{60})/$1\n/g;
	$dna   .= "\n";
        
        my ($fa_out,$fa_outfile);
        if ($split) {
	    $fa_outfile = $outfile;
	    $fa_outfile =~ s/gff$/fa/;
	    open $fa_out, ">$fa_outfile" or die $!; 
	    print $fa_out ">$seq_name\n$dna";
	}
	else {
	    print $out ">$seq_name\n$dna";
	}

	if ( $zip && !$lump ) {
	    system "gzip -f $outfile";
	    system "gzip -f $fa_outfile";
	    $outfile .= '.gz';
	    $fa_outfile .= '.gz' if $split;
	}

	print " GFF3 saved to $outfile";
	print $split ? "; DNA saved to $fa_outfile\n" : "\n";
	
        if ( $summary ) {
	    print "Summary:\nFeature\tCount\n-------\t-----\n";
	
	    for ( keys %method ) {
		print "$_  $method{$_}\n";
	    }
	    print "\n";
	}       
    
    }

    if ( $zip && $lump ) {
	system "gzip -f $lump";
    }
    
    close FH;
}

sub gene_features {
    my $f = shift;
    local $_ = $f->primary_tag;
    $method{$_}++;
    
    if ( /gene/ ) {
	#($gene_id)  = $f->get_tag_values('gene');
	#$gene_id    = 'gene:' . $gene_id;
	$f->add_tag_value( ID => $gene_id );
	$tnum   = 0;
    }
    elsif ( /mRNA/ ) {
        return 0 unless $gene_id;
	$rna_id    = $gene_id;
	$rna_id    =~ s/gene/mRNA/;
	$rna_id   .= '.t0' . ++$tnum;
	$f->add_tag_value( ID => $rna_id );
	$f->add_tag_value( Parent => $gene_id );
    }
    elsif ( /exon/ || /CDS/ ) {
	return 0 unless $rna_id;
	$f->add_tag_value( Parent => $rna_id );
    }
    else {
	return 0 unless $gene_id;
	$f->add_tag_value( Parent => $gene_id );
    }
    
    # now we can skip cloned exons
    return 0 if /exon/ && ++$seen{$f} > 1;

    return 1;
}

sub generic_features {
    my ($f, $io, $refseq) = @_;
    my $method = $f->primary_tag;
    $method{$method}++;

    if ( $f->has_tag($method) ) {
	my ($fname) = $f->get_tag_values($method);
	$f->add_tag_value( ID => "$method:$fname" )
	    unless $f->has_tag('ID');
    }
    else {
	$idh->generate_unique_persistent_id($f);
    }

    $io->gff_string($f);
}

sub gff_header {
    my ($name, $end) = @_;
    
    return <<END;
##gff-version 3
##sequence-region $name 1 $end
##source bp_genbank2gff3.pl
$name\tGenBank\tregion\t1\t$end\t.\t.\t.\tID=$name
END
}

sub unflatten_seq {
    my $seq = shift;

    print "working on contig ", $seq->accession, "..."; 
    my $uh_oh = "Possible gene unflattening error with" .  $seq->accession .
                ": consult STDERR\n";
    
    eval {
	$unflattener->unflatten_seq( -seq => $seq, 
				     -use_magic => 1 );
    };
    
    # deal with unflattening errors
    if ( $@ ) {
	warn $seq->accession . " Unflattening error:\n";
	warn "Details: $@\n";
	print $uh_oh;
    }

    return 0 if !$seq || !$seq->all_SeqFeatures;

    # map feature types to the sequence ontology
    $tm->map_types_to_SO( -seq => $seq );

    1;
}

sub filter {
    my $seq = shift;
    return unless $filter;
    my @feats;

    for my $f ( $seq->remove_SeqFeatures ) {
	my $m = $f->primary_tag;
	push @feats, $f unless $filter =~ /$m/i;
    }

    $seq->add_SeqFeature(@feats) if @feats;
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

    if ($g->has_tag('gene')) {
	($gene_id) = $g->get_tag_values('gene'); 
    }
    elsif ($g->has_tag('locus_tag')) {
	($gene_id) = $g->get_tag_values('locus_tag');
    }

    $gene_id;
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

