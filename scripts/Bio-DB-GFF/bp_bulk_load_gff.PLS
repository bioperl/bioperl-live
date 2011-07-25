#!/usr/bin/perl

use strict;
use warnings;
# use lib './blib/lib';
use DBI;
use IO::File;
use File::Spec;
use Getopt::Long;
use Bio::DB::GFF;
use Bio::DB::GFF::Util::Binning 'bin';

use constant MYSQL => 'mysql';
use constant FDATA      => 'fdata';
use constant FTYPE      => 'ftype';
use constant FGROUP     => 'fgroup';
use constant FDNA       => 'fdna';
use constant FATTRIBUTE => 'fattribute';
use constant FATTRIBUTE_TO_FEATURE => 'fattribute_to_feature';

=head1 NAME

bp_bulk_load_gff.pl - Bulk-load a Bio::DB::GFF database from GFF files.

=head1 SYNOPSIS

  % bp_bulk_load_gff.pl -d testdb dna1.fa dna2.fa features1.gff features2.gff ...

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
value of maxfeature should be increased to 1,000,000,000. This value
must be a power of 10.

Note that Windows users must use the --create option.

If the list of GFF or fasta files exceeds the kernel limit for the 
maximum number of command-line arguments, use the 
--long_list /path/to/files option. 


=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options.
e.g. -d instead of --database.

   --database <dsn>      Database name (default dbi:mysql:test)
   --adaptor             Adaptor name (default mysql)
   --create              Reinitialize/create data tables without asking
   --user                Username to log in as
   --fasta               File or directory containing fasta files to load
   --long_list           Directory containing a very large number of 
                         GFF and/or FASTA files
   --password            Password to use for authentication
                           (Does not work with Postgres, password must be
                           supplied interactively or be left empty for
                           ident authentication)
   --maxbin              Set the value of the maximum bin size
   --local               Flag to indicate that the data source is local
   --maxfeature          Set the value of the maximum feature size (power of 10)
   --group               A list of one or more tag names (comma or space separated)
                         to be used for grouping in the 9th column.
   --gff3_munge          Activate GFF3 name munging (see Bio::DB::GFF)
   --summary             Generate summary statistics for drawing coverage histograms.
                           This can be run on a previously loaded database or during
                           the load.
   --Temporary           Location of a writable scratch directory

=head1 SEE ALSO

L<Bio::DB::GFF>, L<fast_load_gff.pl>, L<load_gff.pl>

=head1 AUTHOR

Lincoln Stein, lstein@cshl.org

Copyright (c) 2002 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

package Bio::DB::GFF::Adaptor::fauxmysql;

use Bio::DB::GFF::Adaptor::dbi::mysqlopt;
use vars '@ISA';
@ISA = 'Bio::DB::GFF::Adaptor::dbi::mysqlopt';

sub insert_sequence {
  my $self = shift;
  my ($id,$offset,$seq) = @_;
  print join("\t",$id,$offset,$seq),"\n";
};

package Bio::DB::GFF::Adaptor::fauxmysqlcmap;

use Bio::DB::GFF::Adaptor::dbi::mysqlcmap;
use vars '@ISA';
@ISA = 'Bio::DB::GFF::Adaptor::dbi::mysqlcmap';

sub insert_sequence {
  my $self = shift;
  my ($id,$offset,$seq) = @_;
  print join("\t",$id,$offset,$seq),"\n";
};

package Bio::DB::GFF::Adaptor::fauxpg;

use Bio::DB::GFF::Adaptor::dbi::pg;
use vars '@ISA';
@ISA = 'Bio::DB::GFF::Adaptor::dbi::pg';

#these two subs are to separate the table creation from the
#index creation
sub do_initialize {
  my $self = shift;
  my $erase = shift;
  $self->drop_all if $erase;
                                                                                
  my $dbh = $self->features_db;
  my $schema = $self->schema;
  foreach my $table_name ($self->tables) {
    my $create_table_stmt = $schema->{$table_name}{table} ;
    $dbh->do($create_table_stmt) ||  warn $dbh->errstr;
  #  $self->create_other_schema_objects(\%{$schema->{$table_name}});
  }
  1;
}

sub _create_indexes_etc {
  my $self = shift;

  my $dbh = $self->features_db;
  my $schema = $self->schema;
  foreach my $table_name ($self->tables) {
    $self->create_other_schema_objects(\%{$schema->{$table_name}});
  }
}

sub insert_sequence {
  my $self = shift;
  my ($id,$offset,$seq) = @_;
  print "$id\t$offset\t$seq\n";
}

package main;

eval "use Time::HiRes"; undef $@;
my $timer = defined &Time::HiRes::time;

my $bWINDOWS = 0;    # Boolean: is this a MSWindows operating system?
if ($^O =~ /MSWin32/i) {
    $bWINDOWS = 1;
}

my ($DSN,$ADAPTOR,$FORCE,$USER,$PASSWORD,$FASTA,$LOCAL,$MAX_BIN,$GROUP_TAG,$LONG_LIST,$MUNGE,$TMPDIR);

GetOptions ('database:s'    => \$DSN,
	    'adaptor:s'     => \$ADAPTOR,
	    'create'        => \$FORCE,
	    'user:s'        => \$USER,
	    'password:s'    => \$PASSWORD,
	    'fasta:s'       => \$FASTA,
	    'local'         => \$LOCAL,
	    'maxbin|maxfeature:s'    => \$MAX_BIN,
	    'group:s'       => \$GROUP_TAG,
	    'long_list:s'   => \$LONG_LIST,
	    'gff3_munge'    => \$MUNGE,
	    'Temporary:s'   => \$TMPDIR,
	   ) or (system('pod2text', $0), exit -1);

# If called as pg_bulk_load_gff.pl behave as that did.
if ($0 =~/pg_bulk_load_gff.pl/){
    $ADAPTOR ||= 'Pg';
    $DSN     ||= 'test';
}
$DSN     ||= 'dbi:mysql:test';
$MAX_BIN ||= 1_000_000_000; # to accomodate human-sized chromosomes


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

# postgres DBD::Pg allows 'database', but also 'dbname', and 'db':
# and it must be Pg (not pg)
$DSN=~s/pg:database=/Pg:/i;
$DSN=~s/pg:dbname=/Pg:/i;
$DSN=~s/pg:db=/Pg:/i;

# leave these lines for mysql
$DSN=~s/database=//i;
$DSN=~s/;host=/:/i; #cater for dsn in the form of "dbi:mysql:database=$dbname;host=$host"


my($DBI,$DBD,$DBNAME,$HOST)=split /:/,$DSN;
$DBNAME=$DSN unless $DSN=~/:/;
$ADAPTOR ||= $DBD; 
$ADAPTOR ||= 'mysql';

if ($DBD eq 'Pg') {
        # rebuild DSN, DBD::Pg requires full dbname=<name> format
        $DSN = "dbi:Pg:dbname=$DBNAME";
        if ($HOST) { $DSN .= ";host=$HOST"; }
}

my ($use_mysql,$use_mysqlcmap,$use_pg) = (0,0,0);
if ( $ADAPTOR eq 'mysqlcmap' ) {
  $use_mysqlcmap = 1;
}
elsif ( $ADAPTOR =~ /^mysql/ ) {
  $use_mysql = 1;
}
elsif ( $ADAPTOR eq "Pg" ) {
  $use_pg = 1;
}
else{
    die "$ADAPTOR is not an acceptable database adaptor.";
}


my (@auth,$AUTH);
if (defined $USER) {
  push @auth,(-user=>$USER);
  if ( $use_mysql or $use_mysqlcmap ) {
    $AUTH .= " -u$USER";
  }
  elsif ( $use_pg ) {
    $AUTH .= " -U $USER ";
  }
}
if (defined $PASSWORD) {
  push @auth,(-pass=>$PASSWORD);
  if ( $use_mysql or $use_mysqlcmap ) {
    $AUTH .= " -p$PASSWORD";
  }
#  elsif ( $use_pg ) {
#    $AUTH .= " -W $PASSWORD ";
#  }
}

if (defined $HOST) {
  $AUTH .= " -h$HOST";  
}
if (defined $DBNAME) {
  if ( $use_mysql or $use_mysqlcmap ) {
    $AUTH .= " -D$DBNAME ";
  }
}
if (defined $LOCAL) {
  $LOCAL='local';
  $AUTH.=' --local-infile=1';
}else {
  $LOCAL='';
}

my $faux_adaptor;
if ( $use_mysqlcmap ) {
  $faux_adaptor = "fauxmysqlcmap";
}
elsif ( $use_mysql ) {
  $faux_adaptor = "fauxmysql";
}
elsif ( $use_pg ) {
  $faux_adaptor = "fauxpg";
}

my $db = Bio::DB::GFF->new(-adaptor=>$faux_adaptor,-dsn => $DSN,@auth)
  or die "Can't open database: ",Bio::DB::GFF->error,"\n";

$db->gff3_name_munging(1) if $MUNGE;

$MAX_BIN ? $db->initialize(-erase=>1,-MAX_BIN=>$MAX_BIN) : $db->initialize(1);
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
  if (/\.(fa|fasta|dna|seq|fast)(?:$|\.)/i) {
    push @fasta,$_;
  } else {
    push @gff,$_;
  }
}
@ARGV = @gff;
push @fasta,$FASTA if defined $FASTA;

# drop everything that was there before
my %FH;
my $tmpdir = File::Spec->tmpdir() || '/tmp';
$tmpdir =~ s!\\!\\\\!g if $bWINDOWS; #eliminates backslash mis-interpretation
-d $tmpdir or die <<END;
I could not find a suitable temporary directory to write scratch files into ($tmpdir by default).
Please select a directory and indicate its location by setting the TMP environment variable, or
by using the --Temporary switch.
END
my @fasta_files_to_be_unlinked;
my @files = (FDATA,FTYPE,FGROUP,FDNA,FATTRIBUTE,FATTRIBUTE_TO_FEATURE);
foreach (@files) {
  $FH{$_} = IO::File->new(">$tmpdir/$_.$$") or die $_,": $!";
  $FH{$_}->autoflush;
}

if ( $use_pg ) {
  $FH{FDATA()                }->print("COPY fdata (fid, fref, fstart, fstop, fbin, ftypeid, fscore, fstrand, fphase, gid, ftarget_start, ftarget_stop) FROM stdin;\n");
  $FH{FTYPE()                }->print("COPY ftype (ftypeid, fmethod, fsource) FROM stdin;\n");
  $FH{FGROUP()               }->print("COPY fgroup (gid, gclass, gname) FROM stdin;\n");
  $FH{FATTRIBUTE()           }->print("COPY fattribute (fattribute_id, fattribute_name) FROM stdin;\n");
  $FH{FATTRIBUTE_TO_FEATURE()}->print("COPY fattribute_to_feature (fid, fattribute_id, fattribute_value) FROM stdin;\n");
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

my %tmpfiles; # keep track of temporary fasta files
my $count;
my $fasta_sequence_id;
my $gff3;
my $current_file; #used to reset GFF3 flag in mix of GFF and GFF3 files

$db->preferred_groups(split (/[,\s]+/,$GROUP_TAG)) if defined $GROUP_TAG;

my $last  = Time::HiRes::time() if $timer;
my $start = $last;

  # avoid hanging on standalone --fasta load
if (!@ARGV) {
    $FH{NULL} = IO::File->new(">$tmpdir/null");
    push @ARGV, "$tmpdir/null";
}

my ($cmap_db);
if ($use_mysqlcmap){
  my $options = {
		 AutoCommit       => 1,
		 FetchHashKeyName => 'NAME_lc',
		 LongReadLen      => 3000,
		 LongTruncOk      => 1,
		 RaiseError       => 1,
		};

  $cmap_db = DBI->connect( $DSN, $USER, $PASSWORD, $options );
}
# Only load CMap::Utils if using cmap
unless (!$use_mysqlcmap or
	eval {
	  require Bio::GMOD::CMap::Utils;
	  Bio::GMOD::CMap::Utils->import('next_number');
	  1;
	} 
       ) {
  print STDERR "Error loading Bio::GMOD::CMap::Utils\n";
}


while (<>) {

  $current_file ||= $ARGV;

  # reset GFF3 flag if new filehandle
  unless($current_file eq $ARGV){
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
    print STDERR "Preparing embedded sequence $1\n";
    push @fasta, "$tmpdir/$1\.fa";
    push @fasta_files_to_be_unlinked,"$tmpdir/$1\.fa";
    $tmpfiles{"$tmpdir/$1\.fa"}++;
    next;
  }

  elsif (/^\#\#\s*gff-version\s+(\d+)/) {
    $gff3 = ($1 >= 3);
    $db->print_gff3_warning() if $gff3;
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
  if ( not defined( $ref ) or length ($ref) == 0) {
    warn "\$ref is null.  source = $source, method = $method, group = $group\n";
    next;
  }
  $FEATURES++;
  my $size = $stop-$start+1;
  warn "Feature $group ($size) is larger than $MAX_BIN. You will have trouble retrieving this feature.\nRerun script with --maxfeature set to a higher power of 10.\n" if $size > $MAX_BIN;

  $source = '\N' unless defined $source;
  $score  = '\N' if $score  eq '.';
  $strand = '\N' if $strand eq '.';
  $phase  = '\N' if $phase  eq '.';

  my ($group_class,$group_name,$target_start,$target_stop,$attributes) = $db->split_group($group,$gff3);

  # GFF2/3 transition
  $group_class = [$group_class] unless ref $group_class;
  $group_name  = [$group_name]  unless ref $group_name;

  for (my $i=0; $i < @$group_name; $i++) {
    $group_class->[$i]  ||= '\N';
    $group_name->[$i]   ||= '\N';
    $target_start ||= '\N';
    $target_stop  ||= '\N';
    $method       ||= '\N';
    $source       ||= '\N';

    my $fid     = $FID++;
    my $gid     = $GROUPID{lc join('',$group_class->[$i],$group_name->[$i])}  ||= $GID++;
    my $ftypeid = $FTYPEID{lc join('',$source,$method)}                       ||= $FTYPEID++;

    my $bin = bin($start,$stop,$db->min_bin);
    $FH{ FDATA()  }->print(    join("\t",$fid,$ref,$start,$stop,$bin,$ftypeid,$score,$strand,$phase,$gid,$target_start,$target_stop),"\n"   );
    if ($use_mysqlcmap){
      my $feature_id    = next_number(
				      db         => $cmap_db,
				      table_name => 'cmap_feature',
				      id_field   => 'feature_id',
				     )
	or die 'No feature id';
      my $direction = $strand eq '-' ? -1:1;
      $FH{ FGROUP() }->print(    
			     join("\t",$feature_id,$feature_id,'NULL',0, $group_name->[$i],0,0,'NULL',1,$direction, $group_class->[$i],)
			     ,"\n"
			    ) unless $DONE{"G$gid"}++;
    }
    else {
      $FH{ FGROUP() }->print(    join("\t",$gid,$group_class->[$i],$group_name->[$i]),"\n") unless $DONE{"G$gid"}++;
    }
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
}

$FH{FASTA}->close if exists $FH{FASTA};

for my $file (@fasta) {
  warn "Preparing DNA file $file....\n";
  if ($use_pg){
    $FH{FDNA() }->print("COPY fdna (fref, foffset, fdna) FROM stdin;\n");
  }
  my $old = select($FH{FDNA()});
  $db->load_fasta($file) or warn "Couldn't load fasta file $file: $!";
  if ($use_pg){
    $FH{FDNA() }->print("\\.\n\n");
  }
  warn "done...\n";
  select $old;
  unlink $file if $tmpfiles{$file};
}

if ($use_pg) { 
  $FH{FDATA()                }->print("\\.\n\n");
  $FH{FTYPE()                }->print("\\.\n\n");
  $FH{FGROUP()               }->print("\\.\n\n");
  $FH{FATTRIBUTE()           }->print("\\.\n\n");
  $FH{FATTRIBUTE_TO_FEATURE()}->print("\\.\n\n");
}


$_->close foreach values %FH;
printf STDERR "Total parse time %5.2fs\n",(Time::HiRes::time() - $start) if $timer;
warn "Loading feature data and analyzing tables.  You may see RDBMS messages here...\n";

if ($use_pg){
  warn "Loading feature data.  You may see Postgres comments...\n";

  foreach (@files) {
    my $file = "$tmpdir/$_.$$";

    $AUTH ? system("psql $AUTH -f $file $DBNAME")
          : system('psql','-f', $file, $DBNAME);

    unlink $file;
  }

  warn "Updating sequences ...\n";
  $db->update_sequences();

  warn "Creating indexes ...\n";
  $db->_create_indexes_etc();

  warn "done...\n";

}

elsif( $use_mysql or $use_mysqlcmap ) {
  $start = time();

  my $success = 1;
  my $TERMINATEDBY = $bWINDOWS ? q( LINES TERMINATED BY '\r\n') : ''; 
  for my $f (@files) {
    my $table = function_to_table($f,$ADAPTOR);
    my $sql = join ('; ',
		    "lock tables $table write",
		    "delete from $table",
		    "load data $LOCAL infile '$tmpdir/$f.$$' replace into table $table $TERMINATEDBY",
		    "unlock tables");
    my $command = MYSQL . qq[$AUTH -s -e "$sql"];
    $command =~ s/\n/ /g;
    $success &&= system($command) == 0;
    unlink "$tmpdir/$f.$$";
  }
  printf STDERR "Total load time %5.2fs\n",(time() - $start) if $timer;
  print STDERR "done...\n";

  print STDERR "Analyzing/optimizing tables. You will see database messages...\n";
  $start = time();
  my $sql = '';
  for my $f (@files) {
    my $table = function_to_table($f,$ADAPTOR);
    $sql       .= "analyze table $table;";
  }
  my $command = MYSQL . qq[$AUTH -N -s -e "$sql"];
  $success &&= system($command) == 0;
  printf STDERR "Optimization time time %5.2fs\n",(time() - $start);

  if ($success) {
    print "$FEATURES features successfully loaded\n";
  } else {
    print "FAILURE: Please see standard error for details\n";
    exit -1;
  }
}

foreach (@fasta_files_to_be_unlinked) {
  unlink "$tmpdir/$_.$$";
}

warn "Building summary statistics for coverage histograms...\n";
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
my $db = Bio::DB::GFF->new(-adaptor=>"dbi::$ADAPTOR",-dsn => $DSN,@args)
  or die "Can't open database: ",Bio::DB::GFF->error,"\n";
$db->build_summary_statistics;

exit 0;

sub function_to_table {
    my $function = shift;
    my $adaptor  = shift;

    if ($function eq 'fdata'){
        return 'fdata';
    }
    elsif ($function eq 'ftype'){
        return 'ftype';
    }
    elsif ($function eq 'fgroup'){
        return 'cmap_feature' if ($adaptor eq 'mysqlcmap');
        return 'fgroup';
    }
    elsif ($function eq 'fdna'){
        return 'fdna';
    }
    elsif ($function eq 'fattribute'){
        return 'fattribute';
    }
    elsif ($function eq 'fattribute_to_feature'){
        return 'fattribute_to_feature';
    }
    return '';
}

__END__
