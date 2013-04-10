#!/usr/bin/perl -w
use version;
our $API_VERSION = qv('1.1.7');

use strict;
use File::Path qw(mkpath rmtree);

BEGIN {
    use lib '.';    # for core package test scripts only
    use Bio::Root::Test;

    test_begin(
        -tests               => 89,
        -requires_modules    => [q(Bio::SeqIO::msout)],
        -requires_networking => 0
    );

    use_ok('Bio::SeqIO::msout');

}

# skip tests if the msout.pm module is too old.
my $api_version = $Bio::SeqIO::msout::API_VERSION;
cmp_ok( $api_version, '>=', qv('1.1.5'),
    "Bio::SeqIO::msout is at least api version 1.1.5" );

test_file_1( 0, "msout/msout_infile1" );
test_file_2( 0, "msout/msout_infile2" );
test_file_3( 0, "msout/msout_infile3" );

# tests to run for api versions >= 1.1.6
SKIP: {
    skip q($Bio::SeqIO::msout::API_VERSION < 1.1.6), 22
      unless $api_version >= qv('1.1.6');
    test_file_4( 0, q(msout/msout_infile4) );
}

# tests to run for api versions >= 1.1.7
SKIP: {
    skip q($Bio::SeqIO::msout::API_VERSION < 1.1.7), 4
      unless $api_version >= qv('1.1.7');
    bad_test_file_1( 0, q(msout/bad_msout_infile1) );
    bad_test_file_2( 0, q(msout/bad_msout_infile2) );
}

sub create_dir {

    my $dir = shift;

    $dir = test_input_file($dir);

    unless ( -d $dir ) {
        mkpath($dir);
    }
}

sub remove_dir {

    my $dir = shift;

    $dir = test_input_file($dir);

    if ( -d $dir ) {
        rmtree($dir);
    }
    else { warn "Tried to remove $dir, but it does not exist" }
}

sub test_file_1 {
##############################################################################
## Test file 1
##############################################################################

    my $gzip   = shift;
    my $infile = shift;
    $infile = test_input_file($infile);

    print_file1( $infile, $gzip );

    my $file_sequence = $infile;
    if ($gzip) {
        $file_sequence = "gzip -dc <$file_sequence |";
    }
    my $msout = Bio::SeqIO->new(
        -file   => "$file_sequence",
        -format => 'msout',
    );

    isa_ok( $msout, 'Bio::SeqIO::msout' );

    my $rh_base_conversion_table = $msout->get_base_conversion_table;

    my %attributes = (
        RUNS              => 3,
        SEGSITES          => 7,
        SEEDS             => [qw(1 1 1)],
        MS_INFO_LINE      => 'ms 6 3 -s 7 -I 3 3 2 1',
        TOT_RUN_HAPS      => 6,
        POPS              => [qw(3 2 1)],
        NEXT_RUN_NUM      => 1,
        LAST_READ_HAP_NUM => 0,
        POSITIONS => [qw(79.1001 80.1001 81.101 82.101 83.10001 84.801 85)],
        CURRENT_RUN_SEGSITES => 7
    );

    foreach my $attribute ( keys %attributes ) {
        my $func = lc($attribute);

        if ( $attribute =~ m/POPS|SEEDS|POSITIONS/ ) {
            $func = ucfirst($func);
        }

        $func = 'get_' . $func;
        my @returns = $msout->$func();
        my ( $return, $got );

        # If there were more than one return value, then compare references to
        # arrays instead of scalars
        unless ( @returns > 1 ) {
            $got = shift @returns;
        }
        else { $got = \@returns }

        my $expected = $attributes{$attribute};

        if ( defined $got && defined $expected ) {
            is_deeply( $got, $expected, "Get $attribute" );
        }
        else { is_deeply( $got, $expected, "Get $attribute" ) }
    }

    # Testing next_hap at beginning of run
    my @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_seq );
    my @data_expected = qw(1111111);
    is_deeply( \@data_got, \@data_expected,
        "Get next_hap at beginning of run" );

    # Testing next_hap after beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_seq );
    @data_expected = qw(5555555);
    is_deeply( \@data_got, \@data_expected,
        "Get next_hap after beginning of run" );

    # Surprise test! testing msout::outgroup
    my $outgroup = $msout->outgroup;
    is( $outgroup, 1, "Testing msout::outgroup" );

    # Testing next_pop after beginning of pop
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_pop );
    @data_expected = qw(4444444);
    is_deeply( \@data_got, \@data_expected,
        "Get next_pop after beginning of pop" );

    # Testing next_pop at beginning of pop
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_pop );
    @data_expected = qw(4444444 5555555);
    is_deeply( \@data_got, \@data_expected,
        "Get next_pop at beginning of pop" );

    # Testing next_run after beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_run );
    @data_expected = qw(4444444);
    is_deeply( \@data_got, \@data_expected,
        "Get next_run after beginning of run" );

    # Testing next_pop at beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_pop );
    @data_expected = qw(5555555 5555555 5555555);
    is_deeply( \@data_got, \@data_expected,
        "Get next_pop at beginning of run" );

    # Testing next_hap after pop
    @data_got      = $msout->get_next_hap;
    @data_expected = qw(1010101);
    is_deeply( \@data_got, \@data_expected, "Get next_hap after pop" );

    # Testing next_run after pop and hap
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_run );
    @data_expected = qw(1111111 1515151);
    is_deeply( \@data_got, \@data_expected, "Get next_run after pop and hap" );

    # Testing next_run at beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_run );
    @data_expected = qw(
      1414141
      1414141
      1515151
      1414141
      1515151
      1515151);
    is_deeply( \@data_got, \@data_expected,
        "Get next_run at beginning of run" );

    is( $msout->get_next_run_num, undef, 'have all lines been read?' );
}

sub test_file_2 {
##############################################################################
## Test file 2
##############################################################################

    my $gzip   = shift;
    my $infile = shift;
    $infile = test_input_file($infile);

    print_file2( $infile, $gzip );

    my $file_sequence = $infile;
    if ($gzip) {
        $file_sequence = "gzip -dc <$file_sequence |";
    }

    my $msout = Bio::SeqIO->new(
        -file   => "$file_sequence",
        -format => 'msout',
    );

    isa_ok( $msout, 'Bio::SeqIO::msout' );

    my %attributes = (
        RUNS              => 3,
        SEGSITES          => 7,
        SEEDS             => [qw(1 1 1)],
        MS_INFO_LINE      => 'ms 6 3',
        TOT_RUN_HAPS      => 6,
        POPS              => 6,
        NEXT_RUN_NUM      => 1,
        LAST_READ_HAP_NUM => 0,
        POSITIONS => [qw(79.1001 80.1001 81.101 82.101 83.10001 84.801 85)],
        CURRENT_RUN_SEGSITES => 7
    );

    foreach my $attribute ( keys %attributes ) {
        my $func = lc($attribute);

        if ( $attribute =~ m/POPS|SEEDS|POSITIONS/ ) {
            $func = ucfirst($func);
        }

        $func = 'get_' . $func;
        my @returns = $msout->$func();
        my ( $return, $got );

        # If there were more than one return value, then compare references to
        # arrays instead of scalars
        unless ( @returns > 1 ) {
            $got = shift @returns;
        }
        else { $got = \@returns }

        my $expected = $attributes{$attribute};

        if ( defined $got && defined $expected ) {
            is_deeply( $got, $expected, "Get $attribute" );
        }
        else { is_deeply( $got, $expected, "Get $attribute" ) }
    }

    my $rh_base_conversion_table = $msout->get_base_conversion_table;

    # Testing next_hap at beginning of run
    my @data_got      = $msout->get_next_hap;
    my @data_expected = '1111111';
    is_deeply( \@data_got, \@data_expected,
        "Get next_hap at beginning of run" );

    # Testing next_hap after beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_seq );
    @data_expected = '5555555';
    is_deeply( \@data_got, \@data_expected,
        "Get next_hap after beginning of run" );

    # Surprise test! testing msout::outgroup
    my $outgroup = $msout->outgroup;
    is( $outgroup, 0, "Testing msout::outgroup" );

    # Testing next_pop after beginning of pop
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_pop );
    @data_expected = qw(
      4444444
      4444444
      5555555
      4444444);
    is_deeply( \@data_got, \@data_expected,
        "Get next_pop after beginning of pop" );

    # Testing next_pop at beginning of pop/run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_pop );
    @data_expected = qw(
      5555555
      5555555
      5555555
      1010101
      1111111
      1515151);
    is_deeply( \@data_got, \@data_expected,
        "Get next_pop at beginning of pop/run" );

    # Testing next_run at beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_run );
    @data_expected = qw(
      1414141
      1414141
      1515151
      1414141
      1515151
      1515151);
    is_deeply( \@data_got, \@data_expected,
        "Get next_run at beginning of run" );

    # Testing next_hap at beginning of run 2
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_seq );
    @data_expected = '1515151';
    is_deeply( \@data_got, \@data_expected,
        "Get next_hap at beginning of run 2" );

    # Testing next_run after hap
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_run );
    @data_expected = qw(
      5050505
      5151515
      5555555
      5454545
      5454545);
    is_deeply( \@data_got, \@data_expected, "Get next_run after hap" );

    is( $msout->get_next_run_num, 5, 'next run should be 5.' );

    # getting the last hap of the file via next hap
    # Testing next_run after hap
    @data_got      = $msout->get_next_hap;
    @data_expected = qw( 5555555 );
    is_deeply( \@data_got, \@data_expected, "Get last hap through next_hap" );

}

sub test_file_3 {
##############################################################################
## Test file 3
##############################################################################

    my $gzip   = shift;
    my $infile = shift;
    $infile = test_input_file($infile);

    print_file3( $infile, $gzip );

    my $file_sequence = $infile;
    if ($gzip) {
        $file_sequence = "gzip -dc <$file_sequence |";
    }
    my $msout = Bio::SeqIO->new(
        -file   => "$file_sequence",
        -format => 'msout',
    );

    isa_ok( $msout, 'Bio::SeqIO::msout' );

    my $rh_base_conversion_table = $msout->get_base_conversion_table;

    isa_ok( $msout, 'Bio::SeqIO::msout' );

    my %attributes = (
        RUNS              => 1,
        SEGSITES          => 7,
        SEEDS             => [qw(1 1 1)],
        MS_INFO_LINE      => 'ms 3 1',
        TOT_RUN_HAPS      => 3,
        POPS              => 3,
        NEXT_RUN_NUM      => 1,
        LAST_READ_HAP_NUM => 0,
        POSITIONS => [qw(79.1001 80.1001 81.101 82.101 83.10001 84.801 85)],
        CURRENT_RUN_SEGSITES => 7
    );

    foreach my $attribute ( keys %attributes ) {
        my $func = lc($attribute);

        if ( $attribute =~ m/POPS|SEEDS|POSITIONS/ ) {
            $func = ucfirst($func);
        }

        $func = 'get_' . $func;
        my @returns = $msout->$func();
        my ( $return, $got );

        # If there were more than one return value, then compare references to
        # arrays instead of scalars
        unless ( @returns > 1 ) {
            $got = shift @returns;
        }
        else { $got = \@returns }

        my $expected = $attributes{$attribute};

        if ( defined $got && defined $expected ) {
            is_deeply( $got, $expected, "Get $attribute" );
        }
        else { is_deeply( $got, $expected, "Get $attribute" ) }
    }

    # Testing next_hap at beginning of run
    my @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_pop );
    my @data_expected = qw(1111111 5555555 4444444);
    is_deeply( \@data_got, \@data_expected, "Get next_pop at end of run" );

    is( $msout->get_next_run_num, undef, 'have all lines been read?' );

    # Testing what happens when we read from empty stream
    @data_got      = $msout->get_next_pop;
    @data_expected = ();
    is_deeply( \@data_got, \@data_expected, "Get next_pop at eof" );

    # Testing what happens when we read from empty stream
    @data_got      = $msout->get_next_run;
    @data_expected = ();
    is_deeply( \@data_got, \@data_expected, "Get next_run at eof" );

    # Testing what happens when we read from empty stream
    @data_got      = $msout->get_next_hap;
    @data_expected = undef;
    is_deeply( \@data_got, \@data_expected, "Get next_hap at eof" );

    # Testing what happens when we read from empty stream
    @data_got      = $msout->get_next_seq;
    @data_expected = ();
    is_deeply( \@data_got, \@data_expected, "Get next_seq at eof" );

}

sub test_file_4 {
##############################################################################
## Test file 4
##############################################################################

  # All this does is test to see if Bio::SeqIO::msout can handle ms output files
  # with multiple newline characters randomly interspersed in the file.

    my $gzip   = shift;
    my $infile = shift;
    $infile = test_input_file($infile);

    print_file4( $infile, $gzip );

    my $file_sequence = $infile;
    if ($gzip) {
        $file_sequence = "gzip -dc <$file_sequence |";
    }
    my $msout = Bio::SeqIO->new(
        -file   => "$file_sequence",
        -format => 'msout',
    );

    isa_ok( $msout, 'Bio::SeqIO::msout' );

    my $rh_base_conversion_table = $msout->get_base_conversion_table;

    my %attributes = (
        RUNS              => 3,
        SEGSITES          => 7,
        SEEDS             => [qw(1 1 1)],
        MS_INFO_LINE      => 'ms 6 3 -s 7 -I 3 3 2 1',
        TOT_RUN_HAPS      => 6,
        POPS              => [qw(3 2 1)],
        NEXT_RUN_NUM      => 1,
        LAST_READ_HAP_NUM => 0,
        POSITIONS => [qw(79.1001 80.1001 81.101 82.101 83.10001 84.801 85)],
        CURRENT_RUN_SEGSITES => 7
    );

    foreach my $attribute ( keys %attributes ) {
        my $func = lc($attribute);

        if ( $attribute =~ m/POPS|SEEDS|POSITIONS/ ) {
            $func = ucfirst($func);
        }

        $func = 'get_' . $func;
        my @returns = $msout->$func();
        my ( $return, $got );

        # If there were more than one return value, then compare references to
        # arrays instead of scalars
        unless ( @returns > 1 ) {
            $got = shift @returns;
        }
        else { $got = \@returns }

        my $expected = $attributes{$attribute};

        if ( defined $got && defined $expected ) {
            is_deeply( $got, $expected, "Get $attribute" );
        }
        else { is_deeply( $got, $expected, "Get $attribute" ) }
    }

    # Testing next_hap at beginning of run
    my @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_seq );
    my @data_expected = qw(1111111);
    is_deeply( \@data_got, \@data_expected,
        "Get next_hap at beginning of run" );

    # Testing next_hap after beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_seq );
    @data_expected = qw(5555555);
    is_deeply( \@data_got, \@data_expected,
        "Get next_hap after beginning of run" );

    # Surprise test! testing msout::outgroup
    my $outgroup = $msout->outgroup;
    is( $outgroup, 1, "Testing msout::outgroup" );

    # Testing next_pop after beginning of pop
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_pop );
    @data_expected = qw(4444444);
    is_deeply( \@data_got, \@data_expected,
        "Get next_pop after beginning of pop" );

    # Testing next_pop at beginning of pop
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_pop );
    @data_expected = qw(4444444 5555555);
    is_deeply( \@data_got, \@data_expected,
        "Get next_pop at beginning of pop" );

    # Testing next_run after beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_run );
    @data_expected = qw(4444444);
    is_deeply( \@data_got, \@data_expected,
        "Get next_run after beginning of run" );

    # Testing next_pop at beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_pop );
    @data_expected = qw(5555555 5555555 5555555);
    is_deeply( \@data_got, \@data_expected,
        "Get next_pop at beginning of run" );

    # Testing next_hap after pop
    @data_got      = $msout->get_next_hap;
    @data_expected = qw(1010101);
    is_deeply( \@data_got, \@data_expected, "Get next_hap after pop" );

    # Testing next_run after pop and hap
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_run );
    @data_expected = qw(1111111 1515151);
    is_deeply( \@data_got, \@data_expected, "Get next_run after pop and hap" );

    # Testing next_run at beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $msout->get_next_run );
    @data_expected = qw(
      1414141
      1414141
      1515151
      1414141
      1515151
      1515151);
    is_deeply( \@data_got, \@data_expected,
        "Get next_run at beginning of run" );

    is( $msout->get_next_run_num, undef, 'have all lines been read?' );
}

sub print_file1 {

    my $destination = shift;
    my $gzip        = shift;

    my $out = <<END
ms 6 3 -s 7 -I 3 3 2 1
1 1 1

//
segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
1111111
5555555
4444444
4444444
5555555
4444444
//
segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
5555555
5555555
5555555
1010101
1111111
1515151
//
segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
1414141
1414141
1515151
1414141
1515151
1515151
END
      ;

    if ($gzip) {
        $gzip = "| gzip";
    }
    else { $gzip = ' '; }
    open OUT, "$gzip >$destination" or die "Unable to open $destination\n";

    print OUT $out;
    close OUT;
}

sub print_file2 {

    my $destination = shift;
    my $gzip        = shift;

    my $out = <<END
ms 6 3
1 1 1

//
segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
1111111
5555555
4444444
4444444
5555555
4444444
//
segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
5555555
5555555
5555555
1010101
1111111
1515151
//
segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
1414141
1414141
1515151
1414141
1515151
1515151
//
segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
1515151
5050505
5151515
5555555
5454545
5454545
//
segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
5555555
END
      ;

    if ($gzip) {
        $gzip = "| gzip";
    }
    else { $gzip = ' '; }

    open OUT, "$gzip >$destination" or die "Unable to open $destination\n";

    print OUT $out;
    close OUT;
}

sub print_file3 {

    my $destination = shift;
    my $gzip        = shift;

    my $out = <<END ;
ms 3 1
1 1 1

//
segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
1111111
5555555
4444444
END

    if ($gzip) {
        $gzip = "| gzip";
    }
    else { $gzip = ' '; }

    open OUT, "$gzip >$destination" or die "Unable to open $destination\n";

    print OUT $out;
    close OUT;
}

sub print_file4 {

    my $destination = shift;
    my $gzip        = shift;

    my $out = <<END




ms 6 3 -s 7 -I 3 3 2 1



1 1 1

//

segsites: 7


positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85


1111111
5555555


4444444

4444444
5555555

4444444


//

segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
5555555
5555555



5555555
1010101
1111111
1515151

//
segsites: 7



positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85

1414141
1414141
1515151
1414141
1515151
1515151





END
      ;

    if ($gzip) {
        $gzip = "| gzip";
    }
    else { $gzip = ' '; }
    open OUT, "$gzip >$destination" or die "Unable to open $destination\n";

    print OUT $out;
    close OUT;
}

sub print_to_file {
    my ( $ra_in, $out ) = @_;
    unless ( open OUT, ">$out" ) {
        die "\nCould not open outfile $out!!\n\n";
    }
    print OUT ("@$ra_in");
    close OUT;
}

sub convert_bases_to_nums {

    my ( $rh_base_conversion_table, @seqs ) = @_;

    my @out_seqstrings;
    foreach my $seq (@seqs) {
        my $seqstring = $seq->seq;
        foreach my $base ( keys %{$rh_base_conversion_table} ) {
            $seqstring =~ s/($base)/$rh_base_conversion_table->{$base}/g;
        }
        push @out_seqstrings, $seqstring;
    }

    return @out_seqstrings;

}

sub bad_test_file_1 {
##############################################################################
## Bad Test file 1
##############################################################################

    # This sub tests to see if msout.pm will catch if the msinfo line's
    # advertized haps are less than are actually in the file

    my $gzip   = shift;
    my $infile = shift;
    $infile = test_input_file($infile);

    my $file_sequence = $infile;
    if ($gzip) {
        $file_sequence = "gzip -dc <$file_sequence |";
    }
    my $msout = Bio::SeqIO->new(
        -file   => "$file_sequence",
        -format => 'msout',
    );

    isa_ok( $msout, 'Bio::SeqIO::msout' );

    throws_ok { $msout->get_next_run }
qr/msout file has only 2 hap\(s\), which is less than indicated in msinfo line \( 9 \)/,
      q(Caught error in bad msout file 1);

}

sub bad_test_file_2 {
##############################################################################
## Bad Test file 2
##############################################################################

    # This sub tests to see if msout.pm will catch if the msinfo line's
    # advertized haps are more than are actually in the file

    my $gzip   = shift;
    my $infile = shift;
    $infile = test_input_file($infile);

    my $file_sequence = $infile;
    if ($gzip) {
        $file_sequence = "gzip -dc <$file_sequence |";
    }
    my $msout = Bio::SeqIO->new(
        -file   => "$file_sequence",
        -format => 'msout',
    );

    isa_ok( $msout, 'Bio::SeqIO::msout' );

    throws_ok { $msout->get_next_run }
qr/\'\/\/\' not encountered when expected. There are more haplos in one of the msOUT runs than advertised in the msinfo line/,
      q(Caught error in bad msout file 2);

}
