#!/usr/bin/perl

use strict;
use DBI;
use IO::File;
use Getopt::Long;
use Bio::DB::GFF::Util::Binning 'bin';

use constant MYSQL => 'mysql';

use constant FDATA      => 'fdata';
use constant FTYPE      => 'ftype';
use constant FGROUP     => 'fgroup';
use constant FNOTE      => 'fnote';
use constant MIN_BIN    => 1000;

my ($DSN,$FORCE,$USER,$PASSWORD);

GetOptions ('database:s'    => \$DSN,
	    'create'         => \$FORCE,
	    'user:s'        => \$USER,
	    'password:s'    => \$PASSWORD,
	   ) or die <<USAGE;
Usage: $0 [options] <gff file 1> <gff file 2> ...
Bulk-load a Bio::DB::GFF database from GFF files.

 Options:
   --database <dsn>      Mysql database name
   --adaptor  <adaptor>  Schema adaptor ('dbi::mysql','dbi::mysqlopt')
   --create              Reinitialize/create data tables without asking
   --user                Username to log in as
   --password            Password to use for authentication

NOTE: If no arguments are provided, then the input is taken from
standard input. Compressed files (.gz, .Z, .bz2) are automatically
uncompressed.

The nature of the bulk load requires that the database be on the local
machine and that the indicated user have the "file" privilege to load
the tables and have enough room in /usr/tmp (or whatever is specified
by the \$TMPDIR environment variable), to hold the tables transiently.
USAGE
;

$DSN ||= 'test';

unless ($FORCE) {
  open (TTY,"/dev/tty") or die "/dev/tty: $!\n";
  print STDERR "This operation will delete all existing data in database $DSN.  Continue? ";
  my $f = <TTY>;
  die "Aborted\n" unless $f =~ /^[yY]/;
  close TTY;
}

foreach (@ARGV) {
  $_ = "gunzip -c $_ |" if /\.gz$/;
  $_ = "uncompress -c $_ |" if /\.Z$/;
  $_ = "bunzip2 -c $_ |" if /\.bz2$/;
}

# drop everything that was there before
my %FH;
my $tmpdir = $ENV{TMPDIR} || $ENV{TMP} || '/usr/tmp';
foreach (FDATA,FTYPE,FGROUP,FNOTE) {
  $FH{$_} = IO::File->new("$tmpdir/$_",">") or die $_,": $!";
  $FH{$_}->autoflush;
}

my $FID     = 1;
my $GID     = 1;
my $FTYPEID = 1;
my %GROUPID = ();
my %FTYPEID = ();
my %DONE    = ();

my $count;
while (<>) {
  chomp;
  next if /^\#/;
  my ($ref,$source,$method,$start,$stop,$score,$strand,$phase,$group) = split "\t";

  $source = '\N' unless defined $source;
  $score  = '\N' if $score  eq '.';
  $strand = '\N' if $strand eq '.';
  $phase  = '\N' if $phase  eq '.';

  # handle group parsing
  $group =~ s/(\"[^\"]*);([^\"]*\")/$1$;$2/g;  # protect embedded semicolons in the group
  my @groups = split(/\s*;\s*/,$group);
  foreach (@groups) { s/$;/;/g }

  my ($group_class,$group_name,$target_start,$target_stop,$notes) = split_group(@groups);
  $group_class  ||= '\N';
  $group_name   ||= '\N';
  $target_start ||= '\N';
  $target_stop  ||= '\N';
  $method       ||= '\N';
  $source       ||= '\N';

  my $fid     = $FID++;
  my $gid     = $GROUPID{$group_class,$group_name} ||= $GID++;
  my $ftypeid = $FTYPEID{$source,$method}          ||= $FTYPEID++;

  my $bin = bin($start,$stop,MIN_BIN);
  $FH{ FDATA()  }->print(    join("\t",$fid,$ref,$start,$stop,$bin,$ftypeid,$score,$strand,$phase,$gid,$target_start,$target_stop),"\n"   );
  $FH{ FGROUP() }->print(    join("\t",$gid,$group_class,$group_name),"\n"              ) unless $DONE{"G$gid"}++;
  $FH{ FTYPE()  }->print(    join("\t",$ftypeid,$method,$source),"\n"                   ) unless $DONE{"T$ftypeid"}++;
  $FH{ FNOTE()  }->print(    join("\t",$fid,$_),"\n"                                    ) foreach (@$notes);

  if ( $fid % 1000 == 0) {
    print STDERR "$fid features parsed...";
    print STDERR -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
  }

}
$_->close foreach values %FH;

warn "loading...\n";

my $AUTH = '';
$AUTH .= " -u$USER"     if defined $USER;
$AUTH .= " -p$PASSWORD" if defined $PASSWORD;

foreach (FDATA,FGROUP,FTYPE,FNOTE) {
  my $command =<<END;
${\MYSQL} $AUTH
-e "lock tables $_ write; delete from $_; load data infile '$tmpdir/$_' replace into table $_; unlock tables"
$DSN
END
;
    $command =~ s/\n/ /g;
    system $command;
    unlink "$tmpdir/$_";
}
warn "done...\n";

sub split_group {
  my ($gclass,$gname,$tstart,$tstop,@notes);

  for (@_) {

    my ($tag,$value) = /^(\S+)\s*(.*)/;
    $value =~ s/\\t/\t/g;
    $value =~ s/\\r/\r/g;
    $value =~ s/^"//;
    $value =~ s/"$//;

    # if the tag is "Note", then we add this to the
    # notes array
   if ($tag eq 'Note') {  # just a note, not a group!
     push @notes,$value;
   }

    # if the tag eq 'Target' then the class name is embedded in the ID
    # (the GFF format is obviously screwed up here)
    elsif ($tag eq 'Target' && $value =~ /([^:\"]+):([^\"]+)/) {
      ($gclass,$gname) = ($1,$2);
      ($tstart,$tstop) = /(\d+) (\d+)/;
    }

    elsif (!$value) {
      push @notes,$tag;  # e.g. "Confirmed_by_EST"
    }

    # otherwise, the tag and value correspond to the
    # group class and name
    else {
      ($gclass,$gname) = ($tag,$value);
    }
  }

  return ($gclass,$gname,$tstart,$tstop,\@notes);
}
