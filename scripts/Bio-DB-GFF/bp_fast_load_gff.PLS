#!/usr/bin/perl

use strict;
use warnings;
# use lib './blib/lib';
use DBI;
use IO::File;
use Getopt::Long;
use Bio::DB::GFF::Util::Binning 'bin';
use Bio::DB::GFF::Adaptor::dbi::mysqlopt;

use constant MYSQL => 'mysql';

use constant FDATA      => 'fdata';
use constant FTYPE      => 'ftype';
use constant FGROUP     => 'fgroup';
use constant FDNA       => 'fdna';
use constant FATTRIBUTE => 'fattribute';
use constant FATTRIBUTE_TO_FEATURE => 'fattribute_to_feature';

my $DO_FAST = eval "use POSIX 'WNOHANG'; 1;";

=head1 NAME

bp_fast_load_gff.pl - Fast-load a Bio::DB::GFF database from GFF files.

=head1 SYNOPSIS

  % bp_fast_load_gff.pl -d testdb dna1.fa dna2.fa features1.gff features2.gff ...

=head1 DESCRIPTION

This script loads a Bio::DB::GFF database with the features contained
in a list of GFF files and/or FASTA sequence files.  You must use the
exact variant of GFF described in L<Bio::DB::GFF>.  Various
command-line options allow you to control which database to load and
whether to allow an existing database to be overwritten.

This script is similar to load_gff.pl, but is much faster.  However,
it is hard-coded to use MySQL and probably only works on Unix
platforms due to its reliance on pipes.  See L<bp_load_gff.pl> for an
incremental loader that works with all databases supported by
Bio::DB::GFF, and L<bp_bulk_load_gff.pl> for a fast MySQL loader that
supports all platforms.

=head2 NOTES

If the filename is given as "-" then the input is taken from
standard input. Compressed files (.gz, .Z, .bz2) are automatically
uncompressed.

FASTA format files are distinguished from GFF files by their filename
extensions.  Files ending in .fa, .fasta, .fast, .seq, .dna and their
uppercase variants are treated as FASTA files.  Everything else is
treated as a GFF file.  If you wish to load -fasta files from STDIN,
then use the -f command-line swith with an argument of '-', as in 

    gunzip my_data.fa.gz | bp_fast_load_gff.pl -d test -f -

The nature of the load requires that the database be on the local
machine and that the indicated user have the "file" privilege to load
the tables and have enough room in /usr/tmp (or whatever is specified
by the \$TMPDIR environment variable), to hold the tables transiently.
If your MySQL is version 3.22.6 and was compiled using the "load local
file" option, then you may be able to load remote databases with local
data using the --local option.

About maxfeature: the default value is 100,000,000 bases.  If you have
features that are close to or greater that 100Mb in length, then the
value of maxfeature should be increased to 1,000,000,000. This value
must be a power of 10.

If the list of GFF or fasta files exceeds the kernel limit for the
maximum number of command-line arguments, use the
--long_list /path/to/files option.

The adaptor used is dbi::mysqlopt.  There is currently no way to
change this.

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options.
e.g. -d instead of --database.

   --database <dsn>      Mysql database name
   --create              Reinitialize/create data tables without asking
   --local               Try to load a remote database using local data.
   --user                Username to log in as
   --fasta               File or directory containing fasta files to load
   --password            Password to use for authentication
   --long_list           Directory containing a very large number of
                         GFF and/or FASTA files
   --maxfeature          Set the value of the maximum feature size (default 100Mb; must be a power of 10)
   --group               A list of one or more tag names (comma or space separated)
                         to be used for grouping in the 9th column.
   --gff3_munge          Activate GFF3 name munging (see Bio::DB::GFF)
   --summary             Generate summary statistics for drawing coverage histograms.
                           This can be run on a previously loaded database or during
                           the load.
   --Temporary           Location of a writable scratch directory

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bulk_load_gff.pl>, L<load_gff.pl>

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
  print join "\t",$id,$offset,$seq,"\n";
}

package main;

eval "use Time::HiRes"; undef $@;
my $timer = defined &Time::HiRes::time;

my ($DSN,$CREATE,$USER,$PASSWORD,$FASTA,$FAILED,$LOCAL,%PID,$MAX_BIN,$GROUP_TAG,$LONG_LIST,$MUNGE,$TMPDIR,$SUMMARY_STATS);

if ($DO_FAST) {
  $SIG{CHLD} = sub {
    while ((my $child = waitpid(-1,&WNOHANG)) > 0) {
      delete $PID{$child} or next;
      $FAILED++ if $? != 0;
    }
  }
};

$SIG{INT} = $SIG{TERM} = sub {cleanup(); exit -1};

GetOptions ('database:s'    => \$DSN,
	    'create'        => \$CREATE,
	    'user:s'        => \$USER,
            'local'         => \$LOCAL,
	    'password:s'    => \$PASSWORD,
	    'fasta:s'       => \$FASTA,
	    'group:s'       => \$GROUP_TAG,
	    'long_list:s'   => \$LONG_LIST,
            'maxbin|maxfeature:s'    => \$MAX_BIN,
	    'gff3_munge'    => \$MUNGE,
	    'summary'       => \$SUMMARY_STATS,
	    'Temporary:s'   => \$TMPDIR,
	   ) or (system('pod2text',$0), exit -1);

$DSN ||= 'test';
$MAX_BIN ||= 1_000_000_000; # to accomodate human-sized chromosomes

my (@args,$AUTH);
if (defined $USER) {
  push @args,(-user=>$USER);
  $AUTH .= " -u$USER";
}
if (defined $PASSWORD) {
  push @args,(-pass=>$PASSWORD);
  $AUTH .= " -p$PASSWORD";
}
push @args,(-preferred_groups=>[split(/[,\s+]+/,$GROUP_TAG)]) if defined $GROUP_TAG;

my $db = Bio::DB::GFF->new(-adaptor=>'faux',-dsn => $DSN,@args)
  or die "Can't open database: ",Bio::DB::GFF->error,"\n";

$db->gff3_name_munging(1) if $MUNGE;

if ($CREATE) {
  $SUMMARY_STATS++;
  $MAX_BIN ? $db->initialize(-erase=>1,-MAX_BIN=>$MAX_BIN) : $db->initialize(1);
}

$MAX_BIN ||= $db->meta('max_bin') || 100_000_000;

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
}

foreach (@ARGV) {
  $_ = "gunzip -c $_ |" if /\.gz$/;
  $_ = "uncompress -c $_ |" if /\.Z$/;
  $_ = "bunzip2 -c $_ |" if /\.bz2$/;
}
my(@fasta,@gff);
foreach (@ARGV) {
  if (/\.(fa|fasta|dna|seq|fast)(?:\.|$)/i) {
    push @fasta,$_;
  } else {
    push @gff,$_;
  }
}
@ARGV = @gff;
push @fasta,$FASTA if defined $FASTA;

# initialize state variables
my $FID     = 1;
my $GID     = 1;
my $FTYPEID = 1;
my $ATTRIBUTEID = 1;
my %GROUPID     = ();
my %FTYPEID     = ();
my %ATTRIBUTEID = ();
my %DONE        = ();
my $FEATURES    = 0;

load_tables($db->dbh) unless $CREATE;
my ($major,$minor,$sub) = split /\./,$db->dbh->get_info(18); # SQL_DBMS_VER
my $can_disable_indexes = ($major >= 4 and $minor >= 0);

# open up pipes to the database
my (%FH,%COMMAND);
my $MYSQL  = MYSQL;
my $tmpdir = $TMPDIR || $ENV{TMPDIR} || $ENV{TMP} || File::Spec->tmpdir();
-d $tmpdir or die <<END;
I could not find a suitable temporary directory to write scratch files into ($tmpdir by default).
Please select a directory and indicate its location by setting the TMP environment variable, or
by using the --Temporary switch.
END

my @fasta_files_to_be_unlinked;
my @files = (FDATA,FTYPE,FGROUP,FDNA,FATTRIBUTE,FATTRIBUTE_TO_FEATURE);
foreach (@files) {
  my $file = "$tmpdir/$_.$$";
  print STDERR "creating load file $file...";
  $DO_FAST &&= (system("mkfifo $file") == 0);  # for system(), 0 = success
  print STDERR "ok\n";
  my $delete = $CREATE ? "delete from $_" : '';
  my $local  = $LOCAL ? 'local' : '';
  my $analyze = "analyze table $_";
  my $command =<<END;
$MYSQL $AUTH
-N
-s
-e "lock tables $_ write; $delete; load data $local infile '$file' replace into table $_; unlock tables; $analyze"
$DSN
END
;
  $command =~ s/\n/ /g;
  $COMMAND{$_} = $command;

  if ($DO_FAST) {
    if (my $pid = fork) {
      $PID{$pid} = $_;
      print STDERR "pausing for 0.5 sec..." if $DO_FAST;
      select(undef,undef,undef,0.50); # work around a race condition
      print STDERR "ok\n";
    } else {  # THIS IS IN CHILD PROCESS
      die "Couldn't fork: $!" unless defined $pid;
      exec $command || die "Couldn't exec: $!";
      exit 0;
    }
  }
  print STDERR "opening load file for writing...";
  $FH{$_} = IO::File->new($file,'>') or die $_,": $!";
  print STDERR "ok\n";
  $FH{$_}->autoflush;
}

print STDERR "Fast loading enabled\n"    if $DO_FAST;

my ($count,$gff3,$last,$start,$beginning,$current_file);

$last  = Time::HiRes::time() if $timer;
$beginning = $start = $last;

# avoid hanging on standalone --fasta load
if (!@ARGV) {
    $FH{NULL} = IO::File->new(">$tmpdir/null");
    push @ARGV, "$tmpdir/null";
}

while (<>) {

  # reset GFF3 flag if new filehandle
  $current_file ||= $ARGV;
  unless ($current_file eq $ARGV) {
    undef $gff3;
    $current_file = $ARGV;
  }

  chomp;
  my ($ref,$source,$method,$start,$stop,$score,$strand,$phase,$group);

  # close sequence filehandle if required
  if ( /^\#|\s+|^$|^>|\t/ && defined $FH{FASTA}) {
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
      push @fasta_files_to_be_unlinked,"$tmpdir/$1\.fa";
      print STDERR "Processing embedded sequence $1\n";
      next;
  }

  elsif (/^\#\#\s*group-tags\s+(.+)/) {
    $db->preferred_groups(split(/\s+/,$1));
    next;
  }

  elsif (/^\#\#\s*gff-version\s+(\d+)/) {
    $gff3 = ($1 >= 3);
    $db->print_gff3_warning() if $gff3;
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

  warn "Feature $group is larger than $MAX_BIN. You will have trouble retrieving this feature.\nRerun script with --maxfeature set to a higher power of 10.\n" if $stop-$start+1 > $MAX_BIN;

  $source = '\N' unless defined $source;
  $score  = '\N' if $score  eq '.';
  $strand = '\N' if $strand eq '.';
  $phase  = '\N' if $phase  eq '.';

  my ($gclass,$gname,$target_start,$target_stop,$attributes) = $db->split_group($group,$gff3);
  # GFF2/3 transition
  $gclass = [$gclass] unless ref $gclass;
  $gname  = [$gname]  unless ref $gname;

  for (my $i=0; $i < @$gname; $i++) {
    my $group_class = $gclass->[$i];
    my $group_name  = $gname->[$i];
    $group_class  ||= '\N';
    $group_name   ||= '\N';
    $target_start ||= '\N';
    $target_stop  ||= '\N';
    $method       ||= '\N';
    $source       ||= '\N';

    my $fid     = $FID++;
    my $gid     = $GROUPID{lc join($;,$group_class,$group_name)} ||= $GID++;
    my $ftypeid = $FTYPEID{lc join($;,$source,$method)}          ||= $FTYPEID++;

    my $bin = bin($start,$stop,$db->min_bin);
    $FH{ FDATA()  }->print(    join("\t",$fid,$ref,$start,$stop,$bin,$ftypeid,$score,$strand,$phase,$gid,$target_start,$target_stop),"\n"   );
    $FH{ FGROUP() }->print(    join("\t",$gid,$group_class,$group_name),"\n"              ) unless $DONE{"fgroup$;$gid"}++;
    $FH{ FTYPE()  }->print(    join("\t",$ftypeid,$method,$source),"\n"                   ) unless $DONE{"ftype$;$ftypeid"}++;

    foreach (@$attributes) {
      my ($key,$value) = @$_;
      my $attributeid = $ATTRIBUTEID{lc $key}   ||= $ATTRIBUTEID++;
      $FH{ FATTRIBUTE() }->print( join("\t",$attributeid,$key),"\n"                       ) unless $DONE{"fattribute$;$attributeid"}++;
      $FH{ FATTRIBUTE_TO_FEATURE() }->print( join("\t",$fid,$attributeid,$value),"\n");
    }

    if ( $FEATURES % 1000 == 0) {
      my $now    = Time::HiRes::time() if $timer;
      my $elapsed = $timer ? sprintf(" in %5.2fs",$now - $last) : '';
      $last = $now;
      print STDERR "$fid features parsed$elapsed...";
      print STDERR -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
    }
  }
}

$FH{FASTA}->close if exists $FH{FASTA};

printf STDERR "Feature load time %5.2fs\n",(Time::HiRes::time() - $start) if $timer;
$start = time();

for my $fasta (@fasta) {
  warn "Loading fasta ",(-d $fasta?"directory":"file"), " $fasta\n";
  my $old = select($FH{FDNA()});
  my $loaded = $db->load_fasta($fasta);
  warn "$fasta: $loaded records loaded\n";
  select $old;
}

printf STDERR "Fasta load time %5.2fs\n",(Time::HiRes::time() - $start) if $timer;
$start = time();

my $success = 1;
if ($DO_FAST) {
  warn "Indexing and analyzing tables.  This may take some time (you may see database messages during the process)...\n";
}

$_->close foreach values %FH;

if (!$DO_FAST) {
  warn "Loading feature data and analyzing tables.  You may see database messages here...\n";
  $success &&= system($COMMAND{$_}) == 0 foreach @files;
}

# wait for children
while (%PID) {
  sleep;
}
$success &&= !$FAILED;

cleanup();

printf STDERR "Total parse & load time %5.2fs\n",(Time::HiRes::time() - $beginning) if $timer;

if ($success) {
  print "SUCCESS: $FEATURES features successfully loaded\n";
  exit 0;
} else {
  print "FAILURE: Please see standard error for details\n";
  exit -1;
}

if ($SUMMARY_STATS) {
    warn "Building summary statistics for coverage histograms...\n";
    $db->build_summary_statistics;
}

exit 0;

sub cleanup {
  foreach (@files,@fasta_files_to_be_unlinked) {
    unlink "$tmpdir/$_.$$";
  }
}

# load copies of some of the tables into memory
sub load_tables {
  my $dbh = shift;
  print STDERR "loading normalized group, type and attribute information...";
  $FID         = 1 + get_max_id($dbh,'fdata','fid');
  $GID         = 1 + get_max_id($dbh,'fgroup','gid');
  $FTYPEID     = 1 + get_max_id($dbh,'ftype','ftypeid');
  $ATTRIBUTEID = 1 + get_max_id($dbh,'fattribute','fattribute_id');
  get_ids($dbh,\%DONE,\%GROUPID,'fgroup','gid','gclass','gname');
  get_ids($dbh,\%DONE,\%FTYPEID,'ftype','ftypeid','fsource','fmethod');
  get_ids($dbh,\%DONE,\%ATTRIBUTEID,'fattribute','fattribute_id','fattribute_name');
  print STDERR "ok\n";
}

sub get_max_id {
  my $dbh = shift;
  my ($table,$id) = @_;
  my $sql = "select max($id) from $table";
  my $result = $dbh->selectcol_arrayref($sql) or die $dbh->errstr;
  $result->[0];
}

sub get_ids {
  my $dbh = shift;
  my ($done,$idhash,$table,$id,@columns) = @_;
  my $columns = join ',',$id,@columns;
  my $sql = "select $columns from $table";
  my $sth = $dbh->prepare($sql) or die $dbh->errstr;
  $sth->execute or die $dbh->errstr;
  while (my($id,@cols) = $sth->fetchrow_array) {
    my $key = lc join $;,@cols;
    $idhash->{$key} = $id;
    $done->{$table,$id}++;
  }
}

__END__
