#!/usr/bin/perl -w

=head1 NAME

ncbi_2_gff.pl - Massage NCBI chromosome annotation into GFF-format suitable for Bio::DB::GFF

=head1 VERSION (CVS-info)

 $RCSfile$
 $Revision$
 $Author$
 $Date$


=head2 SYNOPSIS

   perl process_ncbi_human.pl [options] /path/to/gzipped/datafile(s)

=head2 DESCRIPTION

This script massages the chromosome annotation files located at

  ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/maps/mapview/chromosome_order/

into the GFF-format recognized by Bio::DB::GFF. If the resulting GFF-files are loaded into a Bio::DB:GFF database using the utilities described below, the annotation can be viewed in the Generic Genome Browser (http://www.gmod.org/ggb/) and interfaced with using the Bio::DB:GFF libraries.
  (NB these NCBI-datafiles are dumps from their own mapviewer database backend, according to their READMEs)

To produce the GFF-files, download all the chr*sequence.gz files from the FTP-directory above. While in that same directory, run the following example command (see also help clause by running script with no arguments):

process_ncbi_human.pl --locuslink [path to LL.out_hs.gz] chr*sequence.gz

This will unzip all the files on the fly and open an output file with
the name chrom[$chrom]_ncbiannotation.gff for each, read the LocusLink
records into an in-memory hash and then read through the NCBI feature
lines, lookup 'locus' features in the LocusLink hash for details on
'locus' features and print to the proper GFF files.  LL.out_hs.gz is
accessible here at the time of writing:

  ftp://ftp.ncbi.nih.gov/refseq/LocusLink/LL.out_hs.gz

Note that several of the NCBI features are skipped from the
reformatting, either because their nature is not fully known at this
time (TAG,GS_TRAN) or their sheer volume stands in the way of them
being accessibly in Bio::DB::GFF at this time (EST similarities). You
can easily change this by modifying the $SKIP variable to your liking
to add or remove features, but if you add then you will have to add
handling for those new features.

To bulk-import the GFF-files into a Bio::DB::GFF database, use the
bulk_load_gff.pl utility provided with Bio::DB::GFF

=head2 AUTHOR

Gudmundur Arni Thorisson E<lt>mummi@cshl.orgE<gt>

Copyright (c) 2002 Cold Spring Harbor Laboratory

       This code is free software; you can redistribute it
       and/or modify it under the same terms as Perl itself.

=cut

use strict;
use Getopt::Long;
use IO::File;
use Bio::DB::GFF::Util::Binning 'bin';
use File::Basename;

my $self = basename($0);
my ($doTSCSNP,$doLocuslink,$debug);
my $opt = &GetOptions ('locuslink=s'  => \$doLocuslink,
		       'tscsnp=s'     => \$doTSCSNP,
		       'debug=s'      => \$debug,
		       );
die <<USAGE if(!defined($opt) || @ARGV == 0);
Usage: $self [options] <GFF filename or wildcard pattern>
  Massage NCBI chromosome annotation datafiles into GFF-format suitable for importing into  Bio::DB::GFF database. Note that the program handles both unzipped datafiles and gzipped, bzipped or compressed ones, so do not bother with unzipping big downloads before running.
  See 'perldoc $self' for more info
Options:
   --locuslink Path to zipped LocusLink file, currently located at
               ftp://ftp.ncbi.nih.gov/refseq/LocusLink/LL.out_hs.gz
               used to lookup gene description and official symbols
   --tscsnp    DSN string to TSC MySQL database to use for auxiliary
               SNP feature attributes (CSHL internal use)
   --debug     Enable debugging output for the DBI database driver.
               (CSHL internal use)

  Options can be abbreviated.  For example, you can use -l for
--locuslink.
Author: Gudmundur Arni Thorisson <mummi\@cshl.org>
Copyright (c) 2002 Cold Spring Harbor Laboratory
       This library is free software; you can redistribute it
       and/or modify it under the same terms as Perl itself.

USAGE
;

#Prepare decompression streams for input files, if necessary
my %FH;
print "\nPreparing input and output streams:\n";
foreach (@ARGV) {
    my $chrom=basename($_);                       # NG 02-10-24
    ($chrom) = $chrom=~/0?([0-9,XYxy]{1,2})/;     # NG 02-10-24
    unless($chrom)
    {
	print "can't get chrom name from filename '$_', SKIPPING";
	next;
    }
    $FH{'Chr'.$chrom} = IO::File->new("chrom$chrom\_ncbiannotation.gff",">") or die $_,": $!";
    $_ = "gunzip -c $_ |" if /\.gz$/;
    $_ = "uncompress -c $_ |" if /\.Z$/;
    $_ = "bunzip2 -c $_ |" if /\.bz2$/;
}

#If TSC SNP processing is to be performed, connect to db and prepare query
my $dbh;
my $tsc_sth;
if($doTSCSNP)
{
    #Is this an argstring or file with the string in it?
    my $dbistring = -f $doTSCSNP ? `cat $doTSCSNP` : $doTSCSNP;
    $dbh = &dbConnect($dbistring);
    $dbh->trace($debug) if $debug;
    print "\nConnecting to TSC database using '$dbistring'\n";
    my $query = qq/
SELECT ta.snp_id,
       ta.variation,
       tr.locus_id,
       tr.gene_symbol,
       tr.fxn_class,
       tf.institute_code as lab,
       tf.pop_type,
       tf.outcome,
       tf.pooled_data,
       tf.num_people_typed,
       tf.allele_a_freq as A,
       tf.allele_c_freq as C,
       tf.allele_g_freq as G,
       tf.allele_t_freq as T
FROM tbl_dbsnp_2_tsc dbsnp
LEFT JOIN tbl_refsnp_gene_fxn tr on tr.refsnp_id=dbsnp.rs_id
LEFT JOIN tbl_allele_freq tf on tf.dbsnp_id=dbsnp.ss_id
LEFT JOIN TBL_SNP_ALL ta on ta.snp_id = dbsnp.tsc_id
WHERE dbsnp.rs_id = ? limit 1/;
    $tsc_sth = $dbh->prepare($query);
}

#If Locuslink-processing is to be performed, Read 
#previously cached data structure from disk
my $llData;


if($doLocuslink)
{
    $doLocuslink = "gunzip -c $doLocuslink |" if $doLocuslink =~ /\.gz$/;
    $doLocuslink = "uncompress -c $doLocuslink |" if $doLocuslink =~ /\.Z$/;
    $doLocuslink = "bunzip -c $doLocuslink |" if $doLocuslink =~ /\.bz$/;
    open LL,$doLocuslink || die $!;
    my $l = 0;
    while(<LL>)
    {
	$l++;
	print "\r--$l LocusLink records loaded" if $l % 100 ==0;
	my ($id,$osym,$isym,$mim,$chrom,$loc,$desc,$taxid,$db) = split /\t/;
	my $name = $osym || $isym;
	#print "  Loading in Locuslink id='$id',osym='$osym',isym='$isym',name='$name',desc='$desc'\n";
	$llData->{$id}->{name} = $name;
	$llData->{$id}->{isym} = $isym;
	$llData->{$id}->{mim}  = $mim;
	$llData->{$id}->{chrom}= $chrom;
	$llData->{$id}->{loc}  = $loc;
	$llData->{$id}->{desc} = $desc;
	$llData->{$id}->{taxid}= $taxid;
	$llData->{$id}->{db}   = $db;
    }
    close LL;
}


my %sources     = (snp        => 'dbSNP',
		   sts        => 'UniSTS',
		   locus      => 'LocusLink',
		   transcript => 'RefSeq',
		   transcript_mouse => 'RefSeq',
		   transcript_human => 'RefSeq',
		   component  => 'Genbank',
		   contig     => 'RefSeq',
		   tag        => 'SAGE',
		   gs_tran    => 'GenomeScan',
		   );
my %classes = (component    => 'Sequence',
	       sts          => 'STS',
	       snp          => 'SNP',
	       locus        => 'Locus',
	       transcript   => 'Transcript',
	       transcript_human   => 'Transcript',
	       transcript_mouse   => 'Transcript',
	       contig       => 'Contig',
	       clone        => 'Clone',
	       tag          => 'SAGE_tag',
	       gs_tran      => 'Transcript',
	       );

my %subcomponents = (transcript => 'exon',
		     transcript_human => 'exon',
		     transcript_mouse => 'exon',
		     gs_tran => 'exon',
		     locus      => 'exon',		     
		     component  => 'subcomponent',
		     );


#And now process all incoming data streams
my $i = 0;
my %maxCoords;
my $SKIP = q/^EST|component|clone/;
my %groups   = (); #aggregate parent features
my %density  = ();
my $binSize  = 100000;
my $max = 0;
print "\nStarting main loop:\n";
while(<>)
{
    chomp;
    next if /^\#/;
    my ($type,$objId,$name,$chrom,$start,$stop,$strand) = split "\t";
    my ($class,$source);
    my $score = '.';
    $strand = '.' unless $strand =~ m/^(\+|\-)$/;
    $type = lc $type;
    if($type eq 'gs_tran')
    {
	$type   = 'transcript';
	$source = 'GenomeScan';
    }
    elsif($type eq 'est_human' && $name =~ /^NM_/)
    {
	$type   = 'transcript'; 
	$source = 'RefSeq-human';
    }
    elsif($type eq 'est_mouse' && $name =~ /^NM_/)
    {
	$type   = 'transcript';
	$source ='RefSeq-mouse';
    }
    next if $type =~ /$SKIP/i;
    $i++;
    #my ($chrom,$ctg) = split /\|/,$chromctg;
    next if $chrom =~ /NT/; #ambigously placed NT-contig at start of chrom
    $chrom = "Chr$chrom";
    $max = $stop if $stop > $max;
    $class ||= $classes{$type};
    unless($class)
    {
	print "need class for type '$type': '$_' (OR add type to \$SKIP pattern\n";
	next;
    }
    my $method = $type;
    $source ||= $sources{$type} || die "ERROR: need source for type '$type'";
    $objId = $name if $type =~ /transcript|snp|contig|gs_tran|tag/;
    my $attributes = qq/$class $objId/;
    $attributes .= qq/; Name $name/ unless $objId eq $name;
    my $bin = &bin($start,$stop,$binSize);
    $bin =~ s/^[10]+\.[0]+//;
    $bin ||= 0;
    #print "\$bin='$bin' ($start=>$stop)\n";
    $density{$chrom}->{$method.'_dens:'.$source}->{$bin}++;
    
    #Deduce start/stop for certain parent features to be printed
    #to output file AFTER we've processed everything. This is 
    #necessary because NCBI only gives start/stop values for the child
    #features, like exons in a gene, but not the whole parent feature
    if($type =~ /transcript|locus/)
    {
	$groups{$type}->{$objId}->{$chrom}->{name} = $name;
	$groups{$type}->{$objId}->{$chrom}->{start} ||= 9999999999999;
	$groups{$type}->{$objId}->{$chrom}->{stop} ||= 0;
	$groups{$type}->{$objId}->{$chrom}->{start} = $start 
	    if  $start < $groups{$type}->{$objId}->{$chrom}->{start};
	$groups{$type}->{$objId}->{$chrom}->{stop} = $stop 
	    if $stop > $groups{$type}->{$objId}->{$chrom}->{stop}; 
	$groups{$type}->{$objId}->{$chrom}->{source} = $source;
	$groups{$type}->{$objId}->{$chrom}->{strand} = $strand;
	$groups{$type}->{$objId}->{$chrom}->{method} = $method;
	$groups{$type}->{$objId}->{$chrom}->{class} = $class;
	$method  = $subcomponents{$type};
	next if  $type eq 'locus';
    }
    #This is for internal CSHL usage
    elsif($type =~ /snp/ && $doTSCSNP)
    {
	#print "  -got refSNP ID: $name, let's do TSC lookup\n";
	if(my $tscAttributes = &queryTSCdb($dbh,$name))
	{
	    #$FH{$chrom}->print(qq/$chrom\tTSC\tsnp\t$start\t$stop\t.\t$strand\t.\t$tscAttributes\n/);
	    $attributes .= $tscAttributes;
	    #print "\$attributes='$attributes'\n";
	}
    }

    #Trying to work around the contig pile-up at the start of a chromosome 
    if($method eq 'contig' && $stop == 0)
    {
	print STDERR "SKIPPING, contig '$name' as stop = $stop and start = $start.\n";
    }

    #And finally print to the proper output stream
    $FH{$chrom}->print(qq/$chrom\t$source\t$method\t$start\t$stop\t.\t$strand\t$score\t$attributes\n/);

    #Collect max coordinates, to deduce chromosome sizes
    $maxCoords{$chrom} ||= 0;
    $maxCoords{$chrom} = $stop if $stop > $maxCoords{$chrom};
    
    #Progress indicator
    if ( $i % 1000 == 0) 
    {
	my ($chrom) = $ARGV=~ /0?([0-9,XYxy]{1,2})/;
	print STDERR "$i total features parsed. Now doing chromosome $chrom";
	print STDERR -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
    }
}#MAIN LOOP ENDS


#Print out group features like transcripts and genes that 
#were collected before and print to the proper output streams
print "\nPrinting out collected aggregate features\n";
foreach my $type(keys %groups)
{
    foreach my $objId (keys %{$groups{$type}})
    {
	#print "\$name='$name'\n";
	foreach my $chrom(keys %{$groups{$type}->{$objId}})
	{
	    my $name    = $groups{$type}->{$objId}->{$chrom}->{name};
	    my $start   = $groups{$type}->{$objId}->{$chrom}->{start};
	    my $stop    = $groups{$type}->{$objId}->{$chrom}->{stop};
	    my $strand  = $groups{$type}->{$objId}->{$chrom}->{strand};
	    my $method  = $groups{$type}->{$objId}->{$chrom}->{method};
	    my $class   = $groups{$type}->{$objId}->{$chrom}->{class};
	    my $source  = $groups{$type}->{$objId}->{$chrom}->{source};
	    if($type eq 'locus' && $doLocuslink)
	    {
		my $llInfo = '';
		my $ll    = $llData->{$objId};
		my $id    = $ll->{id};
		my $note  = $ll->{desc} ? qq/Note "$name:$ll->{desc}"/ : ' ';
		$note =~ s/;/\\;/g;
		$FH{$chrom}->print( qq/$chrom\t$source\t$method\t$start\t$stop\t.\t$strand\t.\tLocus $objId/);
		$FH{$chrom}->print(qq/; Name $name/) unless $objId eq $name;
		$FH{$chrom}->print(qq/; $note\n/);
	    }
	    else
	    {
		$FH{$chrom}->print(qq/$chrom\t$source\t$method\t$start\t$stop\t.\t$strand\t.\t$class $name; Name $name\n/);
	    }
	}
    }
}

#   	$density{$method.'_dens:'.$source}->{$bin}++;
#Print out the collected binned density stats
print "Printing out density stats\n";

foreach my $chrom(sort keys %density)
{
    foreach my $meth(sort keys %{$density{$chrom}})
    {
	my $bc = 0;
	foreach my $bin(sort {$a<=>$b}keys %{$density{$chrom}->{$meth}})
	{
	    $bc++;
	    my $count = $density{$chrom}->{$meth}->{$bin};
	    my $binstart = $bin*$binSize;
	    my $binstop  = $binstart+$binSize;
	    print "  \$bin=$bin,\$binstart=$binstart,\$count=$count\n";
	    $FH{$chrom}->print(qq/$chrom\tNCBI\t$meth\t$binstart\t$binstop\t$count\t.\t.\t\n/);
	}
	print " $bc bins for method $meth\n";
    }
}

#Print a line for the reference sequences themselves
while(my ($chrom,$max) = each %maxCoords)
{
    $FH{$chrom}->print(qq/$chrom\tassembly\tchromosome\t1\t$max\t.\t+\t.\tSequence \"$chrom\"\n/);
}

print "\nDONE. $i features parsed\n\n";

#------------------------------------------------
# Subroutines
#------------------------------------------------

#For internal CSHL use. Queries our inhouse MySQL database with 
#SNP Consortium data for various auxiliary data on some SNPs
sub queryTSCdb
{
    my $dbh   = shift;
    my $rs_id = shift;
    my $attributes;
    $rs_id  =~ s/rs//;

    #Baeat vid herna, na i classification string fra dbSNP
    my ($note,$tsc_id,$var,$lab,$dbsnp_id,$class,$gene_symbol,$locus_id,$freq);
    $tsc_sth->execute($rs_id) || die $@;
    my $tscInfo = $tsc_sth->fetchrow_hashref() || return undef;
    do{
	#while(my($k,$v) = each %$tscInfo){print "  '$k'=>'$v'\n";}
	#dbSNP stuff, may apply to more than just TSC snps
	$locus_id = $tscInfo->{locus_id};
	$class = $tscInfo->{fxn_class};
	$gene_symbol = $tscInfo->{gene_symbol};
 	$attributes .= qq/; SNPClass $class/ if $class;

	#TSC specific stuff
	$tsc_id = $tscInfo->{snp_id} || return $attributes;
	$tsc_id = sprintf("TSC%7.7d", $tsc_id);
	$var = lc $tscInfo->{variation};
	$lab = $tscInfo->{institute_code};
	$dbsnp_id = $tscInfo->{dbsnp_id};
	$freq = 1 if $tscInfo->{pop_type} && $tscInfo->{outcome} eq 'S';
	$attributes .= qq/; Alias $tsc_id/;
	#$attributes .= qq/; Variation $var/ if $var;
	$attributes .= qq/; AllFreq 1/ if $freq;
    }while($tscInfo = $tsc_sth->fetchrow_hashref());
    $note = qq/; Note "$tsc_id($var)"/;
    $rs_id = $dbh=$tsc_id=$var=$lab=$dbsnp_id=$class=$gene_symbol=$locus_id= $tscInfo=$freq= undef;
    return  $attributes.$note;
}

sub dbConnect
{
    my $dsn = shift;
    my $dbh;
    use DBI;
    eval{$dbh = DBI->connect($dsn,
			     {
				 RaiseError => 1,
				 FetchHashKeyName => 'NAME_lc',
			     }
                             )
         };
    if($@ || !$dbh)
    {
        print STDERR "ERROR, cannot connect to DB! $@\n";
        die $DBI::errstr;
    }
    return $dbh;
}

