#!/usr/bin/perl

use strict;
# use lib './blib/lib';
use DBI;
use IO::File;
use Getopt::Long;
use Bio::DB::GFF;
use Bio::DB::GFF::Util::Binning 'bin';
use Bio::DB::GFF::Adaptor::dbi::mysqlopt;

use constant MYSQL => 'mysql';

use constant FDATA      => 'fdata';
use constant FTYPE      => 'ftype';
use constant FGROUP     => 'fgroup';
use constant FDNA       => 'fdna';
use constant FATTRIBUTE => 'fattribute';
use constant FATTRIBUTE_TO_FEATURE => 'fattribute_to_feature';

=head1 NAME

bulk_load_gff.pl - Bulk-load a Bio::DB::GFF database from GFF files.

=head1 SYNOPSIS

  % bulk_load_gff.pl -d testdb dna1.fa dna2.fa features1.gff features2.gff ...

=head1 DESCRIPTION

This script loads a Bio::DB::GFF database with the features contained
in a list of GFF files and/or FASTA sequence files.  You must use the
exact variant of GFF described in L<Bio::DB::GFF>.  Various
command-line options allow you to control which database to load and
whether to allow an existing database to be overwritten.

This script differs from bp_load_gff.pl in that it is hard-coded to use
MySQL and cannot perform incremental loads.  See L<bp_load_gff.pl> for an
incremental loader that works with all databases supported by
Bio::DB::GFF, and L<bp_fast_load_gff.pl> for a MySQL loader that supports
fast incremental loads.

=head2 NOTES

If the filename is given as "-" then the input is taken from standard
input. Compressed files (.gz, .Z, .bz2) are automatically
uncompressed.

FASTA format files are distinguished from GFF files by their filename
extensions.  Files ending in .fa, .fasta, .fast, .seq, .dna and their
uppercase variants are treated as FASTA files.  Everything else is
treated as a GFF file.  If you wish to load -fasta files from STDIN,
then use the -f command-line swith with an argument of '-', as in 

    gunzip my_data.fa.gz | bp_fast_load_gff.pl -d test -f -

The nature of the bulk load requires that the database be on the local
machine and that the indicated user have the "file" privilege to load
the tables and have enough room in /usr/tmp (or whatever is specified
by the \$TMPDIR environment variable), to hold the tables transiently.

Local data may now be uploaded to a remote server via the --local option
with the database host specified in the dsn, e.g. dbi:mysql:test:db_host

The adaptor used is dbi::mysqlopt.  There is currently no way to
change this.

About maxfeature: the default value is 100,000,000 bases.  If you have
features that are close to or greater that 100Mb in length, then the
value of maxbin should be increased to 1,000,000,000.

Note that Windows users must use the --create option.

If the list of GFF or fasta files exceeds the kernel limit for the 
maximum number of command-line arguments, use the 
--long_list /path/to/files option. 


=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options.
e.g. -d instead of --database.

   --database <dsn>      Mysql database name (default dbi:mysql:test)
   --create              Reinitialize/create data tables without asking
   --user                Username to log in as
   --fasta               File or directory containing fasta files to load
   --long_list           Directory containing a very large number of 
                         GFF and/or FASTA files
   --password            Password to use for authentication
   --local               Flag to indicate that the data source is local
   --maxfeature          Set the value of the maximum feature size
   --group               A list of one or more tag names (comma or space separated)
                         to be used for grouping in the 9th column.

=head1 SEE ALSO

L<Bio::DB::GFF>, L<fast_load_gff.pl>, L<load_gff.pl>

=head1 AUTHOR

Lincoln Stein, lstein@cshl.org

Copyright (c) 2002 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

package Bio::DB::GFF::Adaptor::faux;

use Bio::DB::GFF::Adaptor::dbi::mysqlopt;
use vars '@ISA';
@ISA = 'Bio::DB::GFF::Adaptor::dbi::mysqlopt';

sub insert_sequence {
  my $self = shift;
  my ($id,$offset,$seq) = @_;
  print join("\t",$id,$offset,$seq),"\n";
}

package main;

eval "use Time::HiRes"; undef $@;
my $timer = defined &Time::HiRes::time;

my $bWINDOWS = 0;    # Boolean: is this a MSWindows operating system?
if ($^O =~ /MSWin32/i) {
    $bWINDOWS = 1;
}

my ($DSN,$FORCE,$USER,$PASSWORD,$FASTA,$LOCAL,$MAX_BIN,$GROUP_TAG,$LONG_LIST);

GetOptions ('database:s'    => \$DSN,
	    'create'        => \$FORCE,
	    'user:s'        => \$USER,
	    'password:s'    => \$PASSWORD,
	    'fasta:s'       => \$FASTA,
            'local'         => \$LOCAL,
            'maxbin:s'      => \$MAX_BIN,
            'maxfeature:s'  => \$MAX_BIN,
	    'group:s'       => \$GROUP_TAG,
	    'long_list:s'   => \$LONG_LIST
	   ) or (system('pod2text', $0), exit -1);

$DSN ||= 'dbi:mysql:test';

if ($bWINDOWS && not $FORCE) {
  die "Note that Windows users must use the --create option.\n";
}

unless ($FORCE) {
  die "This will delete all existing data in database $DSN.  If you want to do this, rerun with the --create option.\n"
    if $bWINDOWS;
  open (TTY,"/dev/tty") or die "/dev/tty: $!\n";  #TTY use removed for win compatability
  print STDERR "This operation will delete all existing data in database $DSN.  Continue? ";
  my $f = <TTY>;
  die "Aborted\n" unless $f =~ /^[yY]/;
  close TTY;
}

my (@auth,$AUTH);
if (defined $USER) {
  push @auth,(-user=>$USER);
  $AUTH .= " -u$USER";
}
if (defined $PASSWORD) {
  push @auth,(-pass=>$PASSWORD);
  $AUTH .= " -p$PASSWORD";
}

$DSN=~s/database=//i;
$DSN=~s/;host=/:/i; #cater for dsn in the form of "dbi:mysql:database=$dbname;host=$host"

my($DBI,$DBD,$DBNAME,$HOST)=split /:/,$DSN;
if (defined $HOST) {
  $AUTH .= " -h$HOST";  
}
$DBNAME=$DSN unless $DSN=~/:/;
if (defined $DBNAME) {
  $AUTH .= " -D$DBNAME";
}
if (defined $LOCAL) {
  $LOCAL='local';
  $AUTH.=' --local-infile=1';
}else {
  $LOCAL='';
}


my $db = Bio::DB::GFF->new(-adaptor=>'faux',-dsn => $DSN,@auth)
  or die "Can't open database: ",Bio::DB::GFF->error,"\n";

$MAX_BIN ? $db->initialize(-erase=>1,-MAX_BIN=>$MAX_BIN) : $db->initialize(1);

# deal with really long lists of files
if ($LONG_LIST) {
  -d $LONG_LIST or die "The --long_list argument must be a directory\n";
  opendir GFFDIR,$LONG_LIST or die "Could not open $LONG_LIST for reading: $!";
  @ARGV = map { "$LONG_LIST\/$_" } readdir GFFDIR;
  closedir GFFDIR;

  if (defined $FASTA && -d $FASTA) {
    opendir FASTA,$FASTA or die "Could not open $FASTA for reading: $!";
    push @ARGV, map { "$FASTA\/$_" } readdir FASTA;
    closedir FASTA;
  }
  elsif (defined $FASTA && -f $FASTA) {
    push @ARGV, $FASTA;
  }
}

foreach (@ARGV) {
  $_ = "gunzip -c $_ |" if /\.gz$/;
  $_ = "uncompress -c $_ |" if /\.Z$/;
  $_ = "bunzip2 -c $_ |" if /\.bz2$/;
}

my (@gff,@fasta);
foreach (@ARGV) {
  if (/\.(fa|fasta|dna|seq|fast)$/i) {
    push @fasta,$_;
  } else {
    push @gff,$_;
  }
}
@ARGV = @gff;
push @fasta,$FASTA if defined $FASTA;

# drop everything that was there before
my %FH;
my $tmpdir = $ENV{TMPDIR} || $ENV{TMP} || '/usr/tmp';
$tmpdir =~ s!\\!\\\\!g if $bWINDOWS; #eliminates backslash mis-interpretation
my @files = (FDATA,FTYPE,FGROUP,FDNA,FATTRIBUTE,FATTRIBUTE_TO_FEATURE);
foreach (@files) {
  $FH{$_} = IO::File->new(">$tmpdir/$_.$$") or die $_,": $!";
  $FH{$_}->autoflush;
}

my $FID     = 1;
my $GID     = 1;
my $FTYPEID = 1;
my $ATTRIBUTEID = 1;
my %GROUPID     = ();
my %FTYPEID     = ();
my %ATTRIBUTEID = ();
my %DONE        = ();
my $FEATURES    = 0;

my $count;
my $fasta_sequence_id;
my $gff3;

$db->preferred_groups(split (/[,\s]+/,$GROUP_TAG)) if defined $GROUP_TAG;

my $last  = Time::HiRes::time() if $timer;
my $start = $last;

while (<>) {
  chomp;
  my ($ref,$source,$method,$start,$stop,$score,$strand,$phase,$group);

  # close sequence filehandle if required
  if ( /^\#|\s+|^$/ && defined $FH{FASTA}) {
    $FH{FASTA}->close;
    delete $FH{FASTA};
  }

  # print to fasta file if the handle is open
  if ( defined $FH{FASTA} ) {
    $FH{FASTA}->print("$_\n");
    next;
  }

  elsif (/^>(\S+)/) {  # uh oh, sequence coming
    $FH{FASTA} = IO::File->new(">$tmpdir/$1\.fa") or die "FASTA: $!\n";
    $FH{FASTA}->print("$_\n");
    push @fasta, "$tmpdir/$1\.fa";
    next;
  }
  
  elsif (/^\#\#\s*gff-version\s+3/) {
    $gff3++;
    next;
  }

  elsif (/^\#\#\s*group-tags\s+(.+)/) {
    $db->preferred_groups(split(/\s+/,$1));
    next;
  }

  elsif (/^\#\#\s*sequence-region\s+(\S+)\s+(\d+)\s+(\d+)/i) { # header line
    ($ref,$source,$method,$start,$stop,$score,$strand,$phase,$group) =
	($1,'reference','Component',$2,$3,'.','.','.',$gff3 ? "ID=Sequence:$1": qq(Sequence "$1"));
  }

  elsif (/^\#/) {
    next;
  }

  else {
    ($ref,$source,$method,$start,$stop,$score,$strand,$phase,$group) = split "\t";
  }
  next unless defined $ref;
  $FEATURES++;

  warn "Feature $group is larger than $MAX_BIN. You will have trouble retrieving this feature.\nRerun script with --maxfeature set to a higher value.\n" if $stop-$start+1 > $MAX_BIN;

  $source = '\N' unless defined $source;
  $score  = '\N' if $score  eq '.';
  $strand = '\N' if $strand eq '.';
  $phase  = '\N' if $phase  eq '.';

  my ($group_class,$group_name,$target_start,$target_stop,$attributes) = $db->split_group($group,$gff3);
  $group_class  ||= '\N';
  $group_name   ||= '\N';
  $target_start ||= '\N';
  $target_stop  ||= '\N';
  $method       ||= '\N';
  $source       ||= '\N';

  my $fid     = $FID++;
  my $gid     = $GROUPID{$group_class,$group_name} ||= $GID++;
  my $ftypeid = $FTYPEID{$source,$method}          ||= $FTYPEID++;

  my $bin = bin($start,$stop,$db->min_bin);
  $FH{ FDATA()  }->print(    join("\t",$fid,$ref,$start,$stop,$bin,$ftypeid,$score,$strand,$phase,$gid,$target_start,$target_stop),"\n"   );
  $FH{ FGROUP() }->print(    join("\t",$gid,$group_class,$group_name),"\n"              ) unless $DONE{"G$gid"}++;
  $FH{ FTYPE()  }->print(    join("\t",$ftypeid,$method,$source),"\n"                   ) unless $DONE{"T$ftypeid"}++;

  foreach (@$attributes) {
    my ($key,$value) = @$_;
    my $attributeid = $ATTRIBUTEID{$key}   ||= $ATTRIBUTEID++;
    $FH{ FATTRIBUTE() }->print( join("\t",$attributeid,$key),"\n"                       ) unless $DONE{"A$attributeid"}++;
    $FH{ FATTRIBUTE_TO_FEATURE() }->print( join("\t",$fid,$attributeid,$value),"\n");
  }

  if ( $fid % 1000 == 0) {
    my $now    = Time::HiRes::time() if $timer;
    my $elapsed = $timer ? sprintf(" in %5.2fs",$now - $last) : '';
    $last = $now;
    print STDERR "$fid features parsed$elapsed...";
    print STDERR -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
  }

}

for my $file (@fasta) {
  warn "Preparing DNA file $file....\n";
  my $old = select($FH{FDNA()});
  $db->load_fasta($file);
  warn "done...\n";
  select $old;
}

$_->close foreach values %FH;

printf STDERR "Total parse time %5.2fs\n",(Time::HiRes::time() - $start) if $timer;
warn "Loading feature data.  You may see duplicate key warnings here...\n";

$start = time();

my $success = 1;
my $TERMINATEDBY = $bWINDOWS ? q( LINES TERMINATED BY '\r\n') : ''; 
foreach (@files) {
  my $command =<<END;
${\MYSQL} $AUTH
-e "lock tables $_ write; delete from $_; load data $LOCAL infile '$tmpdir/$_.$$' replace into table $_  $TERMINATEDBY; unlock tables"
END
;
  $command =~ s/\n/ /g;
  $success &&= system($command) == 0;
  unlink "$tmpdir/$_.$$";
}
printf STDERR "Total load time %5.2fs\n",(time() - $start) if $timer;
warn "done...\n";

if ($success) {
  print "$FEATURES features successfully loaded\n";
  exit 0;
} else {
  print "FAILURE: Please see standard error for details\n";
  exit -1;
}

__END__
