#!/usr/bin/perl
# $Id$

use strict;
use lib './blib/lib';
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

my ($DSN,$CREATE,$USER,$PASSWORD,$FASTA,$FAILED,$LOCAL,%PID);

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
	   ) or die <<USAGE;
Usage: $0 [options] <gff file 1> <gff file 2> ...
Fast load a Bio::DB::GFF database from GFF files.

 Options:
   --database <dsn>      Mysql database name
   --create              Reinitialize/create data tables (will drop existing data)
   --local               Try to load a remote database using local data
   --user                Username to log in as
   --fasta               File or directory containing fasta files to load
   --password            Password to use for authentication

Options can be abbreviated.  For example, you can use -d for
--database.

NOTE: If no arguments are provided, then the input is taken from
standard input. Compressed files (.gz, .Z, .bz2) are automatically
uncompressed.  Fasta files must end in .fa optionally followed by a
compression suffix in order to be recognized.

The nature of the bulk load requires that the database be on the local
machine and that the indicated user have the "file" privilege to load
the tables and have enough room in /usr/tmp (or whatever is specified
by the \$TMPDIR environment variable), to hold the tables transiently.

The --remote option *may* allow you to load remote databases, for
example:

   fast_load_gff.pl -local -d 'test;host=brie3.cshl.org' test_data.gff

However, for this to work, MySQL must have been compiled with the
--enable-local-file option, thereby allowing the "load data local
infile" syntax to work.

The adaptor used is dbi::mysqlopt.  There is currently no way to
change this.

USAGE
;

$DSN ||= 'test';

my (@auth,$AUTH);
if (defined $USER) {
  push @auth,(-user=>$USER);
  $AUTH .= " -u$USER";
}
if (defined $PASSWORD) {
  push @auth,(-pass=>$PASSWORD);
  $AUTH .= " -p$PASSWORD";
}

my $db = Bio::DB::GFF->new(-adaptor=>'faux',-dsn => $DSN,@auth)
  or die "Can't open database: ",Bio::DB::GFF->error,"\n";

$db->initialize(1) if $CREATE;

foreach (@ARGV) {
  $_ = "gunzip -c $_ |" if /\.gz$/;
  $_ = "uncompress -c $_ |" if /\.Z$/;
  $_ = "bunzip2 -c $_ |" if /\.bz2$/;
}

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

# open up pipes to the database
my (%FH,%COMMAND);
my $MYSQL = MYSQL;
my $tmpdir = $ENV{TMPDIR} || $ENV{TMP} || '/usr/tmp';
my @files = (FDATA,FTYPE,FGROUP,FDNA,FATTRIBUTE,FATTRIBUTE_TO_FEATURE);
foreach (@files) {
  my $file = "$tmpdir/$_.$$";
  warn "creating $file...\n";
  $DO_FAST &&= (system("mkfifo $file") == 0);  # for system(), 0 = success
  warn "...ok\n";
  my $delete = $CREATE ? "delete from $_" : '';
  my $local  = $LOCAL ? 'local' : '';
  my $command =<<END;
$MYSQL $AUTH
-e "lock tables $_ write; $delete; load data $local infile '$file' replace into table $_; unlock tables"
$DSN
END
;
  $command =~ s/\n/ /g;
  $COMMAND{$_} = $command;

  if ($DO_FAST) {
    if (my $pid = fork) {
      $PID{$pid} = $_;
    } else {
      die "Couldn't fork: $!" unless defined $pid;
      exec $command || die "Couldn't exec: $!";
    }
  }
  select(undef,undef,undef,0.2); # work around a race condition
  warn "opening file for writing...\n";
  $FH{$_} = IO::File->new($file,'>') or die $_,": $!";
  warn "...ok\n";
  $FH{$_}->autoflush;
}

warn "Fast loading enabled\n" if $DO_FAST;

my $count;
while (<>) {
  chomp;
  next if /^\#/;
  my ($ref,$source,$method,$start,$stop,$score,$strand,$phase,$group) = split "\t";
  $FEATURES++;

  $source = '\N' unless defined $source;
  $score  = '\N' if $score  eq '.';
  $strand = '\N' if $strand eq '.';
  $phase  = '\N' if $phase  eq '.';

  # handle group parsing
  $group =~ s/\\;/$;/g;  # protect embedded semicolons in the group
  $group =~ s/( \"[^\"]*);([^\"]*\")/$1$;$2/g;
  my @groups = split(/\s*;\s*/,$group);
  foreach (@groups) { s/$;/;/g }

  my ($group_class,$group_name,$target_start,$target_stop,$attributes) = Bio::DB::GFF->_split_group(@groups);
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
  $FH{ FGROUP() }->print(    join("\t",$gid,$group_class,$group_name),"\n"              ) unless $DONE{"fgroup$;$gid"}++;
  $FH{ FTYPE()  }->print(    join("\t",$ftypeid,$method,$source),"\n"                   ) unless $DONE{"ftype$;$ftypeid"}++;

  foreach (@$attributes) {
    my ($key,$value) = @$_;
    my $attributeid = $ATTRIBUTEID{$key}   ||= $ATTRIBUTEID++;
    $FH{ FATTRIBUTE() }->print( join("\t",$attributeid,$key),"\n"                       ) unless $DONE{"fattribute$;$attributeid"}++;
    $FH{ FATTRIBUTE_TO_FEATURE() }->print( join("\t",$fid,$attributeid,$value),"\n");
  }

  if ( $FEATURES % 1000 == 0) {
    print STDERR "$FEATURES features parsed...";
    print STDERR -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
  }
}

if ($FASTA) {
  warn "Loading fasta ",(-d $FASTA?"directory":"file"), " $FASTA\n";
  my $old = select($FH{FDNA()});
  my $loaded = $db->load_fasta($FASTA);
  warn "$FASTA: $loaded records loaded\n";
  select $old;
}

my $success = 1;
$_->close foreach values %FH;

if (!$DO_FAST) {
  warn "Loading feature data.  You may see duplicate key warnings here...\n";
  $success &&= system($COMMAND{$_}) == 0 foreach @files;
}

# wait for children
while (%PID) {
  sleep;
}
$success &&= !$FAILED;

cleanup();

if ($success) {
  print "SUCCESS: $FEATURES features successfully loaded\n";
  exit 0;
} else {
  print "FAILURE: Please see standard error for details\n";
  exit -1;
}

exit 0;

sub cleanup {
  foreach (@files) {
    unlink "$tmpdir/$_.$$";
  }
}

# load copies of some of the tables into memory
sub load_tables {
  my $dbh = shift;
  $FID         = 1 + get_max_id($dbh,'fdata','fid');
  $GID         = 1 + get_max_id($dbh,'fgroup','gid');
  $FTYPEID     = 1 + get_max_id($dbh,'ftype','ftypeid');
  $ATTRIBUTEID = 1 + get_max_id($dbh,'fattribute','fattribute_id');
  get_ids($dbh,\%DONE,\%GROUPID,'fgroup','gid','gclass','gname');
  get_ids($dbh,\%DONE,\%FTYPEID,'ftype','ftypeid','fsource','fmethod');
  get_ids($dbh,\%DONE,\%ATTRIBUTEID,'fattribute','fattribute_id','fattribute_name');
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
    my $key = join $;,@cols;
    $idhash->{$key} = $id;
    $done->{$table,$id}++;
  }
}

__END__
