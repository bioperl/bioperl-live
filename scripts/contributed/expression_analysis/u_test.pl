#!/usr/local/bin/perl

# Permutation U-test (2001)

open (PREIN, '*.csv') or die "$!";      # initial raw data file
@prefile = <PREIN>;
close (PREIN);
open (OUT, '>rank.csv') or die "$!";     # convert to rank
foreach $preline (@prefile) {
    $preline =~ s/\n//;
    ($name, $a, $b, $c, $d, $e, $f, $g, $h) = split(/,/,$preline);
    @a = ($a, $b, $c, $d, $e, $f, $g, $h);
    @rr = &rank(@a);
    print OUT "$name,$rr[0],$rr[1],$rr[2],$rr[3],$rr[4],$rr[5],$rr[6],$rr[7]\n";
}
close OUT;

open (IN, 'rank.csv') or die "$!";
@file = <IN>;
close (IN);
open (OUT, '>result.csv') or die "$!";

foreach $line (@file) {
    $line =~ s/\n//;
    ($name, $a, $b, $c, $d, $e, $f, $g, $h) = split(/,/,$line);
    @array = ($a, $b, $c, $d, $e, $f, $g, $h);
    @permute = "";

    for ($tra = 0; $tra <= 9999; ++ $tra) {          # number of permutations
        $U = $array[0] + $array[1] + $array[2] + $array[3];
       #  $U1 = 26 - ($array[0] + $array[1] + $array[2] + $array[3]);
       #  $U2 = 26 - ($array[4] + $array[5] + $array[6] + $array[7]);
       #  if ($U1 >= $U2) {
       #       $U = $U2;
       #  } else {
       #       $U = $U1;
       #  }
       #  print "$U, ";
         push (@permute, $U);
         &FYS;
    }
    @permute = splice(@permute,1,10000);            # number of permutations
    $original = $permute[0];
    @sorted = sort {$a <=> $b} @permute;
    $low = ($sorted[248] + $sorted[249]) / 2;
    $high = ($sorted [9750] + $sorted[9751]) / 2;                        # P = .05
    if (($original < $low or $original > $high) and ($original <= 6 or $original >= 22)) {
        print OUT "$name, $a, $b, $c, $d, $e, $f, $g, $h\n";
        $gene += 1;
    }
}

print $gene."genes";

sub FYS {                           # Fisher_Yates Shuffle
    my $i;
    for ($i = 7; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($array[$i],$array[$j]) = ($array[$j],$array[$i]);
    }
}

sub rank {                           # convert to rank
    my @rank = map {
        [0,$_];
    }
    @_;
    my $i = 0;
    for( sort {$a->[1] <=> $b->[1]} @rank ) {
        $_->[0] = $i++;
    }
    return map {
        $_->[0];
    }
    @rank;
}

close (OUT);


