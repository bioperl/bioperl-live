#!/usr/bin/perl -w

# Author Jason Stajich <jason@bioperl.org>
# 
# Make a non-redundant database based on sequence (not on ID!)
# This script is still in progress but is intended to mimic what 
# Warren Gish's nrdb does

# It requires that Digest::MD5 is installed (for now)

use strict;
use Bio::SeqIO;
use Getopt::Long;

use Digest::MD5 qw(md5_hex);
my $output;
GetOptions(
	   'o|output:s'    => \$output,
	   
	   );

my @files;
if( @ARGV == 1 ) {
    push @files, ['', shift @ARGV];
}

while( @ARGV ) {
    my ($db,$file) = (shift @ARGV, shift @ARGV);
    if( ! $file ) { 
	die("Must provide 'id val' pairing of db name and dbfile");
    }
    push @files, [ $db,$file];
}


my $out;
if( defined $output ) {
    $out = new Bio::SeqIO(-file => ">$output");
} else { 
    $out = new Bio::SeqIO(); # use STDOUT
}
my %unique;
foreach my $pair ( @files ) {
    my ($db,$file) = @$pair;
    my $in = new Bio::SeqIO(-file => $file);
    while( my $seq = $in->next_seq ) {
	if( $db ) { 
	    $seq->display_id("$db:".$seq->display_id);
	}
	my $s = lc($seq->seq());
	my $md5sum = md5_hex($s);
	if( defined $unique{$md5sum} ) {
	    my $desc = $unique{$md5sum}->description;
	    my $id2 = sprintf("; %s:%s %s",$db,$seq->display_id,$seq->description);
	    $unique{$md5sum}->desc($desc . $id2);
	    $unique{$md5sum} = $seq;
	} else { 
	    $unique{$md5sum} = $seq;	
	}
    }
}

foreach my $seq ( values %unique ) {
    $out->write_seq($seq);
}
