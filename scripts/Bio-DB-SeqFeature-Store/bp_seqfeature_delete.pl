#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use File::Spec;
use Bio::DB::SeqFeature::Store;

my $DSN      = 'dbi:mysql:test';
my $USER     = '';
my $PASS     = '';
my $ADAPTOR  = 'DBI::mysql';
my $NAME     = 0;
my $TYPE     = 0;
my $ID       = 0;
my $VERBOSE  = 1;
my $TEST     = 0;
my $FAST     = 0;

GetOptions(
	   'dsn|d=s'       => \$DSN,
	   'adaptor=s'   => \$ADAPTOR,
	   'verbose!'    => \$VERBOSE,
           'dryrun|dry-run' => \$TEST,
           'name|n'      => \$NAME,
           'type|t'      => \$TYPE,
           'id'          => \$ID,
           'fast|f'          => \$FAST,
	   'user=s'      => \$USER,
	   'password=s'  => \$PASS,
	   ) || die <<END;
Usage: $0 [options] <feature1> <feature2> <feature3>
  Options:
          -d --dsn        The database name ($DSN)
          -a --adaptor    The storage adaptor to use ($ADAPTOR)
          -n --name       Delete features based on name or wildcard pattern (default)
          -t --type       Delete features based on type
          -i --id         Delete features based on primary id
          -v --verbose    Turn on verbose progress reporting (default)
             --noverbose  Turn off verbose progress reporting
          --dryrun        Dry run; report features to be deleted without actually deleting them
          -u --user       User to connect to database as
          -p --password   Password to use to connect to database
          -f --fast       Deletes each item instantly not atomic for full dataset (mainly for deleting massive datasets linked to a type)

Examples:
  
 Delete from mysql database volvox features named f08 f09 f10
     $0 -d volvox -n f08 f09 f10

 Delete features whose names start with f  
     $0 -d volvox -n 'f*'

 Delete all features of type remark, source example
     $0 -d volvox -t remark:example

 Delete all remark features, regardless of source
     $0 -d volvox -t 'remark:*'

 Delete the feature with ID 1234
     $0 -d volvox -i 1234

 Delete all features named f* from a berkeleydb database
     $0 -a berkeleydb -d /usr/local/share/db/volvox -n 'f*'

Remember to protect wildcards against shell interpretation by putting
single quotes around them!
END
    ;

if ($NAME+$TYPE+$ID > 1) {
    die "Please provide only one of the --name, --type or --id options.\nRun \"$0 --help\" for usage.\n";
}

unless (@ARGV) {
    die "Please provide a list of feature names, types or ids.\n Run \"$0 --help\" for usage.\n";
}

my $mode = $ID   ? 'id'
          :$TYPE ? 'type'
          :$NAME ? 'name'
          :'name';


my @options;
@options = ($USER,$PASS) if $USER || $PASS;

my $store = Bio::DB::SeqFeature::Store->new(
					    -dsn     => $DSN,
					    -adaptor => $ADAPTOR,
					    -user    => $USER,
					    -pass    => $PASS,
					    -write    => 1,
    )
  or die "Couldn't create connection to the database";

my @features = retrieve_features($store,$mode,\@ARGV);

if ($VERBOSE || $TEST) {
    print scalar (@features)," feature(s) match.\n\n";
    my $heading;
    foreach (@features) {
	printf "%-20s %-20s %-12s\n%-20s %-20s %-12s\n",
	       'Name','Type','Primary ID',
	       '----','----','----------'
		   unless $heading++;
	printf "%-20s %-20s %-12d\n",$_->display_name,$_->type,$_->primary_id;
    }
    print "\n";
}

if (@features && !$TEST) {
    if($FAST) {
      my $del = 0;
      foreach my $feat(@features) {
        my @tmp_feat = ($feat);
        my $deleted = $store->delete(@tmp_feat);
        $del++ if($deleted);
        if ($VERBOSE && $deleted) {
          print 'Feature ',$del," successfully deleted.\n";
        } elsif (!$deleted) {
          die "An error occurred. Some or all of the indicated features could not be deleted.";
        }
      }
    }
    else {
        my $deleted = $store->delete(@features);
        if ($VERBOSE && $deleted) {
	        print scalar(@features)," features successfully deleted.\n";
        } elsif (!$deleted) {
	        die "An error occurred. Some or all of the indicated features could not be deleted.";
        }
    }
}

exit 0;

sub retrieve_features {
    my($db,$mode,$list) = @_;
    my @features;
    if ($mode eq 'name') {
	@features = map {$db->get_features_by_alias($_)} @$list;
    }
    elsif ($mode eq 'type') {
	my $regexp = glob2regexp(@$list);
	my @types  = grep {/$regexp/} $db->types;
	@features  = $db->get_features_by_type(@types) if @types;
    }
    elsif ($mode eq 'id') {
	@features  = grep {defined $_} map {$db->get_feature_by_primary_id($_)} @$list;
    }
    return @features;
}

sub glob2regexp {
    my @globs = map {
	$_ = quotemeta($_);
	s/\\\*/.*/g;
	s/\?/./g;
	$_ } @_;
    return '^(?:'.join('|',@globs).')$';
 }
