#!/usr/bin/perl -w

#---------------------------------------------------------
# A minimal demo script for the Bio::Root::Utilities.pm module.
# Steve A. Chervitz (sac@genome.stanford.edu)
# $Id$
# Installation: 
#   Edit the use lib line to point at your Bioperl lib
#---------------------------------------------------------

use lib "/home/sac";
use Bio::Root::Utilities qw(:obj);

select(STDOUT);$|=1;

$file = $ARGV[0] || __FILE__;

printf "%-15s: %s\n", 'Default date', $Util->date_format();
printf "%-15s: %s\n", 'Full date', $Util->date_format('full');
printf "%-15s: %s\n", "File date of $file", $Util->file_date($file);

$tmp = 0;

if(-T $file) {
    if(-o $file) {
	$compressed = $Util->compress($file);
	$tmp = 1;
    } else {
	$compressed = $Util->compress($file, 'tmp');
    }
    printf "%-15s: %s\n", "Compressed", $compressed;
} else {
    $compressed = $file;
}

$uncompressed = $tmp ? $Util->uncompress($compressed, 1) : $Util->uncompress($compressed);
printf "%-15s: %s\n", "Uncompressed", $uncompressed;



