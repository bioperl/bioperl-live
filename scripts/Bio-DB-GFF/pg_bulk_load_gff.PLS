#!/usr/bin/perl

use strict;
use lib './blib/lib';
use DBI;
use IO::File;
use Getopt::Long;
use Bio::DB::GFF::Util::Binning 'bin';

use constant FDATA      => 'fdata';
use constant FTYPE      => 'ftype';
use constant FGROUP     => 'fgroup';
use constant FDNA       => 'fdna';
use constant FATTRIBUTE => 'fattribute';
use constant FATTRIBUTE_TO_FEATURE => 'fattribute_to_feature';

=head1 NAME

bulk_load_gff.pl - Bulk-load a Bio::DB::GFF database from GFF files.

=head1 SYNOPSIS

  % bulk_load_gff.pl -d testdb -f dna.fa features1.gff features2.gff ...

=head1 DESCRIPTION

This script loads a Bio::DB::GFF database with the features contained
in a list of GFF files.  You must use the exact variant of GFF
described in L<Bio::DB::GFF>.  Various command-line options allow you
to control which database to load, whether to load DNA (in fasta
format) as well, and whether to allow an existing database to be
overwritten.  This script differs from load_gff.pl in that it is
hard-coded to use PostgreSQL and cannot perform incremental loads.  See
L<load_gff.pl> for an incremental loader that works with all databases
supported by Bio::DB::GFF.

=head2 NOTES

Note that this bulk loader is different from the mysql bulk loader,
as this script generates a Postgres dump file that it prints to
standard out.  Therefore, the most common way to use it (ie, to
load to a PostgreSQL database), pipe the output to psql:

  pg_bulk_load_gff.pl <arguments> | psql <dbname>

The dump file created is only compatible with PostgreSQL version 7.3 or
greater.  Note that the create option will NOT create the database, that
must be done beforehand.  The create option is a flag to indicate
that it is OK to wipe the database (ie, drop the tables), and if you
don't supply it, you will be prompted for that OK before tables are
dropped.  Note that the Postgres user must have table create permissions
or the load will fail.

Also note that if no arguments are provided, or the filename is given as "-"
then the input is taken from standard input. Compressed files (.gz,
.Z, .bz2) are automatically uncompressed.  Fasta files must end in .fa
optionally followed by a compression suffix in order to be recognized.

To load FASTA files without GFF data, run like this:

 load_gff.pl --database my_database --fasta my_fasta_file_path </dev/null

The nature of the bulk load requires that the database be on the local
machine and have enough room in /usr/tmp (or whatever is specified
by the \$TMPDIR environment variable), to hold the tables transiently.

Note that Windows users must use the --create option.

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options.
e.g. -d instead of --database.

   --database <dsn>      PostgreSQL database name (default 'dbi:Pg:test')
   --create              Reinitialize/create data tables without asking
   --user                Username to log in as
   --fasta               File or directory containing fasta files to load
   --password            Password to use for authentication
   --webowner		 User that Apache runs as (default 'nobody')

=head1 SEE ALSO

L<Bio::DB::GFF>, L<fast_load_gff.pl>, L<load_gff.pl>

=head1 AUTHOR

Scott Cain, cain@cshl.org

Copyright (c) 2003 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

package Bio::DB::GFF::Adaptor::faux;

use Bio::DB::GFF::Adaptor::dbi::pg;
use vars '@ISA';
@ISA = 'Bio::DB::GFF::Adaptor::dbi::pg';

sub insert_sequence {
  my $self = shift;
  my ($id,$offset,$seq) = @_;
  print "$id\t$offset\t$seq\n";
}

package main;

my $bWINDOWS = 0;    # Boolean: is this a MSWindows operating system?
if ($^O =~ /MSWin32/i) {
    $bWINDOWS = 1;
}

my ($DSN,$FORCE,$USER,$PASSWORD,$FASTA,$WEBOWNER);

GetOptions ('database:s'    => \$DSN,
	    'create'         => \$FORCE,
	    'user:s'        => \$USER,
	    'password:s'    => \$PASSWORD,
	    'fasta:s'       => \$FASTA,
	    'webowner:s'    => \$WEBOWNER
	   ) or die <<USAGE;
Usage: $0 [options] <gff file 1> <gff file 2> ...
Bulk-load a Bio::DB::GFF database from GFF files.

 Options:
   --database <dsn>      Mysql database name
   --create              Reinitialize/create data tables without asking
   --user                Username to log in as
   --fasta               File or directory containing fasta files to load
   --password            Password to use for authentication
   --webowner		 The user that owns Apache processes (default 'nobody')

Options can be abbreviated.  For example, you can use -d for
--database.

NOTE: If no arguments are provided, then the input is taken from
standard input. Compressed files (.gz, .Z, .bz2) are automatically
uncompressed.  Fasta files must end in .fa optionally followed by a
compression suffix in order to be recognized.

The nature of the bulk load requires that the database be on the local
machine and have enough room in /usr/tmp (or whatever is specified
by the \$TMPDIR environment variable), to hold the tables transiently.

Note that Windows users must use the --create option.

USAGE
;

$DSN ||= 'dbi:Pg:test';

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
  $AUTH .= " -U $USER";
}
if (defined $PASSWORD) {
  push @auth,(-pass=>$PASSWORD);
  $AUTH .= " -W $PASSWORD";
}

my $db = Bio::DB::GFF->new(-adaptor=>'faux',-dsn => $DSN,@auth)
  or die "Can't open database: ",Bio::DB::GFF->error,"\n";

#$db->initialize(1); let this script do that, to separate table and index creation

foreach (@ARGV) {
  $_ = "gunzip -c $_ |" if /\.gz$/;
  $_ = "uncompress -c $_ |" if /\.Z$/;
  $_ = "bunzip2 -c $_ |" if /\.bz2$/;
}

# drop everything that was there before
my %FH;
my $tmpdir = $ENV{TMPDIR} || $ENV{TMP} || '/usr/tmp';
$tmpdir =~ s!\\!\\\\!g if $bWINDOWS; #eliminates backslash mis-interpretation
my @files = (FDATA,FTYPE,FGROUP,FDNA,FATTRIBUTE,FATTRIBUTE_TO_FEATURE);
foreach (@files) {
  $FH{$_} = IO::File->new("$tmpdir/$_.$$",">") or die $_,": $!";
  $FH{$_}->autoflush;
}

$FH{FDATA()                }->print("COPY fdata (fid, fref, fstart, fstop, fbin, ftypeid, fscore, fstrand, fphase, gid, ftarget_start, ftarget_stop) FROM stdin;\n");
$FH{FTYPE()                }->print("COPY ftype (ftypeid, fmethod, fsource) FROM stdin;\n");
$FH{FGROUP()               }->print("COPY fgroup (gid, gclass, gname) FROM stdin;\n");
$FH{FATTRIBUTE()           }->print("COPY fattribute (fattribute_id, fattribute_name) FROM stdin;\n");
$FH{FATTRIBUTE_TO_FEATURE()}->print("COPY fattribute_to_feature (fid, fattribute_id, fattribute_value) FROM stdin;\n");

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
while (<>) {
  chomp;
  my ($ref,$source,$method,$start,$stop,$score,$strand,$phase,$group);
  if (/^\#\#\s*sequence-region\s+(\S+)\s+(\d+)\s+(\d+)/i) { # header line
    ($ref,$source,$method,$start,$stop,$score,$strand,$phase,$group) = 
      ($1,'reference','Component',$2,$3,'.','.','.',qq(Sequence "$1"));
  } elsif (/^\#/) {
    next;
  } else {
    ($ref,$source,$method,$start,$stop,$score,$strand,$phase,$group) = split "\t";
  }

  if (length ($ref) == 0) {
    warn "\$ref is null.  source = $source, method = $method, group = $group\n";
    next;
  }
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
  $FH{ FGROUP() }->print(    join("\t",$gid,$group_class,$group_name),"\n"              ) unless $DONE{"G$gid"}++;
  $FH{ FTYPE()  }->print(    join("\t",$ftypeid,$method,$source),"\n"                   ) unless $DONE{"T$ftypeid"}++;

  foreach (@$attributes) {
    my ($key,$value) = @$_;
    my $attributeid = $ATTRIBUTEID{$key}   ||= $ATTRIBUTEID++;
    $FH{ FATTRIBUTE() }->print( join("\t",$attributeid,$key),"\n"                       ) unless $DONE{"A$attributeid"}++;
    $FH{ FATTRIBUTE_TO_FEATURE() }->print( join("\t",$fid,$attributeid,$value),"\n");
  }

  if ( $fid % 1000 == 0) {
    print STDERR "$fid features parsed...";
    print STDERR -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
  }

}

$WEBOWNER ||='nobody';
table_create($WEBOWNER);

if ($FASTA) {
  warn "Preparing DNA files....\n";
  $FH{FDNA()                 }->print("COPY fdna (fref, foffset, fdna) FROM stdin;\n");

  my $old = select($FH{FDNA()});
  $db->load_fasta($FASTA);
  $FH{FDNA()                 }->print("\\.\n\n");
  warn "done...\n";
  select $old;
}

$FH{FDATA()                }->print("\\.\n\n");
$FH{FTYPE()                }->print("\\.\n\n");
$FH{FGROUP()               }->print("\\.\n\n");
$FH{FATTRIBUTE()           }->print("\\.\n\n");
$FH{FATTRIBUTE_TO_FEATURE()}->print("\\.\n\n");


$_->close foreach values %FH;

warn "Loading feature data.  You will see Postgres comments...\n";

foreach (@files) {
  my $file = "$tmpdir/$_.$$";

  system('cat', $file);

  unlink $file;
}
cleanup();
warn "done...\n";


sub table_create {
  my $webowner = shift;

  print<<END;
SET search_path = public, pg_catalog;

/* CREATE USER $webowner;  */

DROP TABLE fmeta;
CREATE TABLE fmeta (
    fname character varying(255) DEFAULT '' NOT NULL,
    fvalue character varying(255) DEFAULT '' NOT NULL
);

REVOKE ALL ON TABLE fmeta FROM PUBLIC;
GRANT SELECT ON TABLE fmeta TO $webowner;

DROP TABLE fgroup;
CREATE TABLE fgroup (
    gid serial NOT NULL,
    gclass character varying(100),
    gname character varying(100)
);

REVOKE ALL ON TABLE fgroup FROM PUBLIC;
GRANT SELECT ON TABLE fgroup TO $webowner;

DROP TABLE fdata;
CREATE TABLE fdata (
    fid serial NOT NULL,
    fref character varying(100) DEFAULT '' NOT NULL,
    fstart integer DEFAULT '0' NOT NULL,
    fstop integer DEFAULT '0' NOT NULL,
    fbin double precision DEFAULT '0.000000' NOT NULL,
    ftypeid integer DEFAULT '0' NOT NULL,
    fscore double precision,
    fstrand character varying(3),
    fphase character varying(3),
    gid integer DEFAULT '0' NOT NULL,
    ftarget_start integer,
    ftarget_stop integer,
    CONSTRAINT chk_fdata_fphase CHECK ((((fphase = '0'::character varying) OR (fphase = '1'::character varying)) OR (fphase = '2'::character varying))),
    CONSTRAINT chk_fdata_fstrand CHECK (((fstrand = '+'::character varying) OR (fstrand = '-'::character varying)))
);

REVOKE ALL ON TABLE fdata FROM PUBLIC;
GRANT SELECT ON TABLE fdata TO $webowner;

DROP TABLE fattribute_to_feature;
CREATE TABLE fattribute_to_feature (
    fid integer DEFAULT '0' NOT NULL,
    fattribute_id integer DEFAULT '0' NOT NULL,
    fattribute_value text
);

REVOKE ALL ON TABLE fattribute_to_feature FROM PUBLIC;
GRANT SELECT ON TABLE fattribute_to_feature TO $webowner;

DROP TABLE fdna;
CREATE TABLE fdna (
    fref character varying(100) DEFAULT '' NOT NULL,
    foffset integer DEFAULT '0' NOT NULL,
    fdna text
);

REVOKE ALL ON TABLE fdna FROM PUBLIC;
GRANT SELECT ON TABLE fdna TO $webowner;

DROP TABLE  fattribute;
CREATE TABLE fattribute (
    fattribute_id serial NOT NULL,
    fattribute_name character varying(255) DEFAULT '' NOT NULL
);

REVOKE ALL ON TABLE fattribute FROM PUBLIC;
GRANT SELECT ON TABLE fattribute TO $webowner;

DROP TABLE ftype;
CREATE TABLE ftype (
    ftypeid serial NOT NULL,
    fmethod character varying(100) DEFAULT '' NOT NULL,
    fsource character varying(100)
);

REVOKE ALL ON TABLE ftype FROM PUBLIC;
GRANT SELECT ON TABLE ftype TO $webowner;

END
;
}
 
sub cleanup {
  print<<END;
CREATE UNIQUE INDEX fgroup_gclass_idx ON fgroup USING btree (gclass, gname);

CREATE INDEX fdata_ftypeid_idx ON fdata USING btree (ftypeid);

CREATE INDEX fdata_gid_idx ON fdata USING btree (gid);

CREATE UNIQUE INDEX fdata_fref_idx ON fdata USING btree (fref, fbin, fstart, fstop, ftypeid, gid);

CREATE INDEX fattribute_to_feature_fid ON fattribute_to_feature USING btree (fid, fattribute_id);

CREATE INDEX ftype_fmethod_idx ON ftype USING btree (fmethod);

CREATE INDEX ftype_fsource_idx ON ftype USING btree (fsource);

CREATE UNIQUE INDEX ftype_ftype_idx ON ftype USING btree (fmethod, fsource);

ALTER TABLE ONLY fmeta
    ADD CONSTRAINT pk_fmeta PRIMARY KEY (fname);

ALTER TABLE ONLY fgroup
    ADD CONSTRAINT pk_fgroup PRIMARY KEY (gid);

ALTER TABLE ONLY fgroup
    ADD CONSTRAINT gclass_fgroup UNIQUE (gclass, gname);

ALTER TABLE ONLY fdata
    ADD CONSTRAINT pk_fdata PRIMARY KEY (fid);

ALTER TABLE ONLY fdata
    ADD CONSTRAINT fref_fdata UNIQUE (fref, fbin, fstart, fstop, ftypeid, gid);

ALTER TABLE ONLY fdna
    ADD CONSTRAINT pk_fdna PRIMARY KEY (fref, foffset);

ALTER TABLE ONLY fattribute
    ADD CONSTRAINT pk_fattribute PRIMARY KEY (fattribute_id);

ALTER TABLE ONLY ftype
    ADD CONSTRAINT pk_ftype PRIMARY KEY (ftypeid);

ALTER TABLE ONLY ftype
    ADD CONSTRAINT ftype_ftype UNIQUE (fmethod, fsource);

END
;
}

__END__
