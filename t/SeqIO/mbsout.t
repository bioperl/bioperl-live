#!/usr/bin/perl
use version;
our $API_VERSION = qv('1.1.3');

use strict;
use File::Path qw(mkpath rmtree);
use Carp;

BEGIN {
    use lib '.';    # for core package test scripts only
    use Bio::Root::Test;

    test_begin(
        -tests               => 74,
        -requires_modules    => [q(Bio::SeqIO::mbsout)],
        -requires_networking => 0
    );

    use_ok('Bio::SeqIO::mbsout');
}

# skip tests if the msout.pm module is too old.
cmp_ok( $Bio::SeqIO::mbsout::API_VERSION,
    '>=', qv('1.1.3'), "Bio::SeqIO::mbsout is at least api version 1.1.3" );

create_dir("mbsout");
test_file_1( 0, "mbsout/mbsout_infile1" );
test_file_2( 0, "mbsout/mbsout_infile2" );
test_file_3( 0, "mbsout/mbsout_infile3" );

sub create_dir {

    my $dir = shift;

    $dir = test_input_file($dir);

    unless ( -d $dir ) {
        mkpath($dir);
    }
}

sub test_file_1 {
##############################################################################
## Test file 1
##############################################################################

    my $gzip   = shift;
    my $infile = shift;
    $infile = test_input_file($infile);

    #print_file1( $infile, $gzip );

    my $file_sequence = $infile;
    if ($gzip) {
        $file_sequence = "gunzip -c <$file_sequence |";
    }
    my $mbsout = Bio::SeqIO->new(
        -file   => "$file_sequence",
        -format => 'mbsout',
    );

    isa_ok( $mbsout, 'Bio::SeqIO::mbsout' );

    my $rh_base_conversion_table = $mbsout->get_base_conversion_table;

    isa_ok( $mbsout, 'Bio::SeqIO::mbsout' );

    my %attributes = (
        RUNS     => 3,
        SEGSITES => 7,
        MBS_INFO_LINE =>
          'command:        mbs 6 -t 0.001 -r 0.00025 -s 5000 2500 -f 3 1 traj ',
        TOT_RUN_HAPS      => 6,
        NEXT_RUN_NUM      => 1,
        LAST_READ_HAP_NUM => 0,
        POSITIONS => [qw(79.1001 80.1001 81.101 82.101 83.10001 84.801 85)],
        CURRENT_RUN_SEGSITES      => 7,
        POP_MUT_PARAM_PER_SITE    => 0.001,
        POP_RECOMB_PARAM_PER_SITE => 0.00025,
        NSITES                    => 5000,
        SELPOS                    => 2500,
        NFILES                    => 3,
        NREPS                     => 1,
        TRAJ_FILENAME             => 'traj'
    );

    foreach my $attribute ( keys %attributes ) {
        my $func = lc($attribute);

        if ( $attribute =~ m/POPS|SEEDS|POSITIONS/ ) {
            $func = ucfirst($func);
        }

        $func = 'get_' . $func;
        my @returns = $mbsout->$func();
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
      convert_bases_to_nums( $rh_base_conversion_table, $mbsout->get_next_seq );
    my @data_expected = qw(1111111);
    is_deeply( \@data_got, \@data_expected,
        "Get next_hap at beginning of run" );

    # Testing next_hap after beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $mbsout->get_next_seq );
    @data_expected = qw(5555555);
    is_deeply( \@data_got, \@data_expected,
        "Get next_hap after beginning of run" );

    # Testing next_pop after beginning of pop
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $mbsout->get_next_seq );
    @data_expected = qw(4444444);
    is_deeply( \@data_got, \@data_expected,
        "Get next_pop after beginning of pop" );

    # Testing next_hap
    @data_got      = $mbsout->get_next_hap;
    @data_expected = qw(4444444);
    is_deeply( \@data_got, \@data_expected, "Get next_hap" );

    # Testing next_hap
    @data_got      = $mbsout->get_next_hap;
    @data_expected = qw(5555555);
    is_deeply( \@data_got, \@data_expected, "Get next_hap" );

    # Testing next_run after beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $mbsout->get_next_run );
    @data_expected = qw(4444444);
    is_deeply( \@data_got, \@data_expected,
        "Get next_run after beginning of run" );

    # Testing next_run at beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $mbsout->get_next_run );
    @data_expected = qw(5555555 5555555 5555555 1010101 1111111 1515151);
    is_deeply( \@data_got, \@data_expected,
        "Get next_run at beginning of run" );

    # Testing next_run at beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $mbsout->get_next_run );
    @data_expected = qw(
      1414141
      1414141
      1515151
      1414141
      1515151
      1515151);
    is_deeply( \@data_got, \@data_expected,
        "Get next_run at beginning of run" );

    is( $mbsout->get_next_run_num, undef, 'have all lines been read?' );
}

sub test_file_2 {
##############################################################################
## Test file 2
##############################################################################

    my $gzip   = shift;
    my $infile = shift;
    $infile = test_input_file($infile);

    #print_file2( $infile, $gzip );

    my $file_sequence = $infile;
    if ($gzip) {
        $file_sequence = "gunzip -c <$file_sequence |";
    }

    my $mbsout = Bio::SeqIO->new(
        -file   => "$file_sequence",
        -format => 'mbsout',
    );

    isa_ok( $mbsout, 'Bio::SeqIO::mbsout' );

    my %attributes = (
        RUNS     => 5,
        SEGSITES => 7,
        MBS_INFO_LINE =>
          'command:        mbs 6 -t 0.001 -r 0.00025 -s 5000 2500 -f 5 1 traj ',
        TOT_RUN_HAPS      => 6,
        NEXT_RUN_NUM      => 1,
        LAST_READ_HAP_NUM => 0,
        POSITIONS => [qw(79.1001 80.1001 81.101 82.101 83.10001 84.801 85)],
        CURRENT_RUN_SEGSITES      => 7,
        POP_MUT_PARAM_PER_SITE    => 0.001,
        POP_RECOMB_PARAM_PER_SITE => 0.00025,
        NSITES                    => 5000,
        SELPOS                    => 2500,
        NFILES                    => 5,
        NREPS                     => 1,
        TRAJ_FILENAME             => 'traj'
    );

    foreach my $attribute ( keys %attributes ) {
        my $func = lc($attribute);
        if ( $attribute =~ m/POSITIONS/ ) {
            $func = ucfirst($func);
        }
        elsif ( $attribute =~ m/\_file/ ) {
            $func = q(infile);
        }

        $func = 'get_' . $func;
        my @returns = $mbsout->$func();
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

    my $rh_base_conversion_table = $mbsout->get_base_conversion_table;

    # Testing next_hap at beginning of run
    my @data_got      = $mbsout->get_next_hap;
    my @data_expected = qw(1111111);
    is_deeply( \@data_got, \@data_expected,
        "Get next_hap at beginning of run" );

    # Testing next_hap after beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $mbsout->get_next_seq );
    @data_expected = '5555555';
    is_deeply( \@data_got, \@data_expected,
        "Get next_hap after beginning of run" );

    # Surprise test! testing mbsout::outgroup
    my $outgroup = $mbsout->outgroup;
    is( $outgroup, 0, "Testing mbsout::outgroup" );

    # Testing next_run after beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $mbsout->get_next_run );
    @data_expected = qw(
      4444444
      4444444
      5555555
      4444444);
    is_deeply( \@data_got, \@data_expected,
        "Get next_run after beginning of run" );

    # Testing next_run after beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $mbsout->get_next_run );
    @data_expected = qw(
      5555555
      5555555
      5555555
      1010101
      1111111
      1515151);
    is_deeply( \@data_got, \@data_expected,
        "Get next_run after beginning of run" );

    # Testing next_run at beginning of run
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $mbsout->get_next_run );
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
      convert_bases_to_nums( $rh_base_conversion_table, $mbsout->get_next_seq );
    @data_expected = '1515151';
    is_deeply( \@data_got, \@data_expected,
        "Get next_hap at beginning of run 2" );

    # Testing next_run after hap
    @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $mbsout->get_next_run );
    @data_expected = qw(
      5050505
      5151515
      5555555
      5454545
      5454545);
    is_deeply( \@data_got, \@data_expected, "Get next_run after hap" );

    is( $mbsout->get_next_run_num, 5, 'next run should be 5.' );

}

sub test_file_3 {
##############################################################################
## Test file 3
##############################################################################

    my $gzip   = shift;
    my $infile = shift;
    $infile = test_input_file($infile);

    #print_file3( $infile, $gzip );

    my $file_sequence = $infile;
    if ($gzip) {
        $file_sequence = "gunzip -c <$file_sequence |";
    }
    my $mbsout = Bio::SeqIO->new(
        -file   => "$file_sequence",
        -format => 'mbsout',
    );

    isa_ok( $mbsout, 'Bio::SeqIO::mbsout' );

    my $rh_base_conversion_table = $mbsout->get_base_conversion_table;

    my %attributes = (
        RUNS     => 1,
        SEGSITES => 7,
        MBS_INFO_LINE =>
          'command:        mbs 3 -t 0.001 -r 0.00025 -s 5000 2500 -f 1 1 traj ',
        TOT_RUN_HAPS      => 3,
        NEXT_RUN_NUM      => 1,
        LAST_READ_HAP_NUM => 0,
        POSITIONS => [qw(79.1001 80.1001 81.101 82.101 83.10001 84.801 85)],
        CURRENT_RUN_SEGSITES      => 7,
        POP_MUT_PARAM_PER_SITE    => 0.001,
        POP_RECOMB_PARAM_PER_SITE => 0.00025,
        NSITES                    => 5000,
        SELPOS                    => 2500,
        NFILES                    => 1,
        NREPS                     => 1,
        TRAJ_FILENAME             => 'traj'
    );

    foreach my $attribute ( keys %attributes ) {
        my $func = lc($attribute);

        if ( $attribute =~ m/POPS|SEEDS|POSITIONS/ ) {
            $func = ucfirst($func);
        }

        $func = 'get_' . $func;
        my @returns = $mbsout->$func();
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

    # Testing next_run at beginning of run
    my @data_got =
      convert_bases_to_nums( $rh_base_conversion_table, $mbsout->get_next_run );
    my @data_expected = qw(1111111 5555555 4444444);
    is_deeply( \@data_got, \@data_expected,
        "Get next_run at end/beginning of run" );

    is( $mbsout->get_next_run_num, undef, 'have all lines been read?' );

    # Testing what happens when we read from empty stream
    @data_got      = $mbsout->get_next_run;
    @data_expected = ();
    is_deeply( \@data_got, \@data_expected, "Get next_run at eof" );

    # Testing what happens when we read from empty stream
    @data_got      = $mbsout->get_next_hap;
    @data_expected = undef;
    is_deeply( \@data_got, \@data_expected, "Get next_hap at eof" );

    # Testing what happens when we read from empty stream
    @data_got      = $mbsout->get_next_seq;
    @data_expected = ();
    is_deeply( \@data_got, \@data_expected, "Get next_seq at eof" );

}

sub print_file1 {

    my $destination = shift;
    my $gzip        = shift;

    my $out = <<END
command:        mbs 6 -t 0.001 -r 0.00025 -s 5000 2500 -f 3 1 traj 

//0-1   allele: a a a d d d d 
segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
1111111
5555555
4444444
4444444
5555555
4444444
//1-1   allele: a a a a d d d 
segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
5555555
5555555
5555555
1010101
1111111
1515151
//2-1   allele: a d d d d d d 
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
    open my $OUT, "$gzip >$destination" or croak "Could not write file '$destination': $!\n";

    print $OUT $out;
    close $OUT;
}

sub print_file2 {

    my $destination = shift;
    my $gzip        = shift;

    my $out = <<END
command:        mbs 6 -t 0.001 -r 0.00025 -s 5000 2500 -f 5 1 traj 

//0-1   allele: a a a d d d d 
segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
1111111
5555555
4444444
4444444
5555555
4444444
//1-1   allele: a a a d d d d 
segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
5555555
5555555
5555555
1010101
1111111
1515151
//2-1   allele: a a a d d d d 

segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
1414141
1414141
1515151
1414141
1515151
1515151
//3-1   allele: a a a a a d d 

segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
1515151
5050505
5151515
5555555
5454545
5454545
//4-1   allele: a a d d d d d 

segsites: 7
positions: 79.1001 80.1001 81.101 82.101 83.10001 84.801 85
5555555
END
      ;

    if ($gzip) {
        $gzip = "| gzip";
    }
    else { $gzip = ' '; }

    open my $OUT, "$gzip >$destination" or croak "Could not write file '$destination': $!\n";

    print $OUT $out;
    close $OUT;
}

sub print_file3 {

    my $destination = shift;
    my $gzip        = shift;

    my $out = <<END ;
command:        mbs 3 -t 0.001 -r 0.00025 -s 5000 2500 -f 1 1 traj 

//0-1   allele: a a a d d d d 
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

    open my $OUT, "$gzip >$destination" or croak "Could not write file '$destination': $!\n";

    print $OUT $out;
    close $OUT;
}

sub print_to_file {
    my ( $ra_in, $out ) = @_;
    open my $OUT, '>', $out or croak "\nCould not write outfile '$out': $!\n";
    print $OUT ("@$ra_in");
    close $OUT;
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
