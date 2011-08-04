#!perl

# Author Jason Stajich <jason-at-bioperl-dot-org>
# 
# Make a non-redundant database based on sequence (not on ID!)
# This script is still in progress but is intended to mimic what 
# Warren Gish's nrdb does

# It requires that Digest::MD5 is installed (for now)

=head1 NAME

bp_nrdb.PLS - a script to emulate Warren Gish's nrdb, make a unique sequence database from a set of input databases

=head1 SYNOPSIS


Usage: 
  bp_nrdb.PLS [options] file1 file2 file3

Alternative usage
  bp_nrdb.PLS -p [options] file1 id1 file2 id2 file3 id3

=head1 DESCRIPTION

This script will create a unique database of sequences
(quasi-nonredundant).  The options are:

   -o filename          - the filename the db is written (STDOUT by default)
   -a filename          - the filename to append the db to
   -l#                  - minimum required sequence length
   -i                   - do not check for duplicates
   -n#                  - max number of descriptions to report per seq
   -d#                  - delimiter to use between consecutive descriptions
   -p                   - use database id prefixes from command line

=head1 AUTHOR

Jason Stajich, jason-at-bioperl-dot-org

=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

use Digest::MD5 qw(md5_hex);
my ($output,$append,$min_len, 
    $no_duplicate_check,$desc_count,
    $delimiter, $expect_prefixes,$help);
$delimiter = ';';

GetOptions(
	   'o|output:s'    => \$output,
	   'a|append:s'    => \$append,
	   'n:s'           => \$desc_count,
	   'l:s'           => \$min_len,
	   'd:s'           => \$delimiter,
	   'p'             => \$expect_prefixes,
	   'i'             => \$no_duplicate_check,
	   'h'             => \$help,
	   );

die("must supply a positive integer for -d") if ( defined $desc_count &&
						  ( $desc_count !~ /^\d+$/ ||
						    $desc_count < 1) );
die("must supply a positive integer for -l") if ( defined $min_len &&
						  ( $min_len !~ /^\d+$/ ||
						    $min_len < 1) );
my @files;

if( $help || ! @ARGV ) {
    exec('perldoc',$0);
    exit(0);
}
while( @ARGV ) {
    
    my ($file, $id) = (undef,'');
    if( $expect_prefixes ) {
	($file,$id) = (shift @ARGV, shift @ARGV);
	if( ! $id ) { 
	    die("Must provide 'name id' pairing of dbfile and id");
	}
    } else { 
	$file = shift @ARGV;
    }
    push @files, [ $file,$id];
}


my $out;
if( $append ) {
    $out = new Bio::SeqIO(-file => ">>$append");
} elsif( $output ) { 
    $out = new Bio::SeqIO(-file => ">$output");
} else {
    $out = new Bio::SeqIO(); # use STDOUT
}

my %unique;
my %seqcount;
my $counter = 0;
foreach my $pair ( @files ) {
    my ($file,$id) = @$pair;
    my $in = new Bio::SeqIO(-file => $file);
    while( my $seq = $in->next_seq ) {
	next if defined $min_len && $seq->length < $min_len;
	if( $id ) { 
	    $seq->display_id("$id:".$seq->display_id);
	}
	my $s = lc($seq->seq());
	my $md5sum = md5_hex($s);
	if( $no_duplicate_check ) {
	    $md5sum = $counter++;
	}
	    
	if( defined $unique{$md5sum} ) {
	    $seqcount{$md5sum}++;
	    next if defined $desc_count && $seqcount{$md5sum++} > $desc_count;
	    my $desc = $unique{$md5sum}->description;	    
	    my $id2 = sprintf("%s %s:%s %s",$delimiter,
			      $id,$seq->display_id,$seq->description);
	    $unique{$md5sum}->desc($desc . $id2);
	} else { 
	    $unique{$md5sum} = $seq;	
	}
    }
}

foreach my $seq ( values %unique ) {
    $out->write_seq($seq);
}

__END__
