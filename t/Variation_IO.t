## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..26\n";
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::Variation::IO;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run.

sub test ($$;$) {
    my($num, $true, $msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
    $true ? return 1 : return 0;
}

my $no = 2;

sub fileformat ($) {
    my ($file) = shift;
    my $format;
    if ($file =~ /.*dat$/) {
	$format = 'flat';
    }
    elsif ($file =~ /.*xml$/ ) {
	$format = 'xml';
    } else {
	print "Wrong extension! [$file]";
	exit;
    }
    return $format;
}

sub ext ($) {
    my ($file) = @_;
    my ($name) = $file =~ /.*.(...)$/;
    return $name;
}

sub filename ($) {
    my ($file) = @_;
    my ($name) = $file =~ /(.*)....$/;
    return $name;
}

#
#  read and write test
#
sub io {
    my ($t_file, $o_file) = @_; 
    my $res;

    my ($t_ext) = ext ($t_file);
    my ($o_ext) = ext ($o_file);
    my ($t_format) = fileformat ($t_file);
    my ($o_format) = fileformat ($o_file);
    my ($t_name) = filename($t_file);
    my ($o_name) = filename($o_file);

    my( $before );
    {
        local $/ = undef;
        local *BEFORE;
        open BEFORE, "$t_name.$o_ext";
        $before = <BEFORE>;
        close BEFORE;
    }

    test $no++, $before, "Error in reading input file [$t_name.$o_ext]";

    # Test reading
    my $in = Bio::Variation::IO->new( -file => $t_file);
    my  @entries ;
    while (my $e = $in->next) {
        push @entries, $e;
    }
    my $count = scalar @entries;
    test $no++, @entries > 0, "No SeqDiff objects [$count]";
    
    # Test writing
    my $out = Bio::Variation::IO->new( -FILE => "> $o_file", -FORMAT => $o_format);
    my $out_ok = 1;
    foreach my $e (@entries) {
        $out->write($e) or $out_ok = 0;
    }
    undef($out);  # Flush to disk
    test $no++, $out_ok,  "error writing variants";

    my( $after );
    {
        local $/ = undef;
        local *AFTER;
        open AFTER, $o_file;
        $after = <AFTER>;
        close AFTER;
    }

    test $no++, $after, "Error in reading in again the output file [$o_file]";

    #Test that input and output files are identical
    $res = test $no++, $before eq $after, "test output file differs from input";
    print STDERR `diff $t_file $o_file` if ! $res;
    unlink($o_file); 
}



#
# First test flat file IO
#

# mutations from one allele to an other
io  't/mutations.dat', 't/mutations.out.dat'; #2..6

#complex sequence difference: two variations, one with multiple alleles
io  't/polymorphism.dat', 't/polymorphism.out.dat'; #7..11

#
# XML IO
#

# first check whether the necessary XML modules are installed
eval {
    require Bio::Variation::IO::xml;
};

if( $@ ) {
    print STDERR
	"\nThe XML-format conversion requires the CPAN modules ",
	"XML::Node, XML::Writer, and IO::String to be installed ",
	"on your system, which they probably aren't. Skipping these tests.\n";
    while($no <= 26) {
	test($no++, 1, "");
    }
    exit(0);
}

# mutations from one allele to an other
eval {
	io  't/mutations.xml', 't/mutations.out.xml'; #12..16
};

#complex sequence difference: two variations, one with multiple alleles
eval {
	io  't/polymorphism.xml', 't/polymorphism.out.xml'; #17..21
};

#
# flat <=> XML
#

# mutations from one allele to an other
eval { 
	io  't/mutations.dat', 't/mutations.out.xml'; #22..26
};

