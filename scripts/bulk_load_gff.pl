#!/usr/bin/perl

# README:
# THIS IS A BULK LOADER FOR THE BIO::DB::GFF DATABASE USING THE DBI::MYSQLOPT ADAPTOR
# CURRENTLY IT IS HARD-CODED TO USE A DATABASE NAMED ELEGANS ON THE LOCAL MACHINE
# BECAUSE IT IS PRIMARILY FOR TESTING AND DEVELOPMENT PURPOSES.
# THIS WILL BE REPLACED WITH A MORE CONFIGURABLE BULK LOADER WHEN THESE MODULES
# REACH RELEASE STATE.
# L STEIN

use strict;
use DBI;
use IO::File;
use Bio::DB::GFF::Util::Binning 'bin';

#NOTE: in test database, I have deleted the following
my %EXCLUDED_SOURCES = ();
my %EXCLUDED_METHODS = () ;

use constant DUMP => '.';
use constant MYSQL => 'mysql';
use constant AUTH  => '';
use constant DB    => 'elegans';

use constant FDATA      => 'fdata';
use constant FTYPE      => 'ftype';
use constant FGROUP     => 'fgroup';
use constant FNOTE      => 'fnote';
use constant MIN_BIN    => 1000;
use constant REF_CLASS  => 'Sequence';

@ARGV = <${\DUMP}/*.gff.gz> unless @ARGV;
foreach (@ARGV) {
  $_ = "gunzip -c $_ |" if /\.gz$/;
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
  my $class = REF_CLASS;

  next if $EXCLUDED_SOURCES{$source};
  next if $EXCLUDED_METHODS{$method};

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

  # special munging for the Sanger GFF conventions...
  if ($source eq 'Link') {
    $group_class = 'Link';
  }

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

foreach (FDATA,FGROUP,FTYPE,FNOTE) {
  my $command =<<END;
${\MYSQL} ${\AUTH}
-e "lock tables $_ write; delete from $_; load data infile '$tmpdir/$_' replace into table $_; unlock tables"
${\DB}
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

    my ($tag,$value) = /^(\S+)(?:\s+\"([^\"]+)\")?/;
    $value =~ s/\\t/\t/g;
    $value =~ s/\\r/\r/g;

    # if the tag is "Note", then we add this to the
    # notes array
   if ($tag eq 'Note') {  # just a note, not a group!
     push @notes,$value;
   }

    # if the tag eq 'Target' then the class name is embedded in the ID
    # (the GFF format is obviously screwed up here)
    elsif ($tag eq 'Target' && /^\"([^:\"]+):([^\"]+)\"/) {
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
