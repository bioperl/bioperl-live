#!/usr/local/bin/perl

# Permutation t-test (2001)

open (IN, '*.csv') or die "$!";        # file name
@file = <IN>;
close (IN);
open (OUT, '>result.csv') or die "$!";

foreach $line (@file) {
    $line =~ s/\n//;
    ($name, $a, $b, $c, $d, $e, $f, $g, $h) = split(/,/,$line);
    @array = ($a, $b, $c, $d, $e, $f, $g, $h);
    @permute = "";

    for ($tra = 0; $tra <= 9999; ++ $tra) {      # number of permutations
        $m1 = ($array[0] + $array[1] + $array[2] + $array[3]) / 4;
        $m2 = ($array[4] + $array[5] + $array[6] + $array[7]) / 4;
        $S1_2 = (($array[0] - $m1)**2 + ($array[1] - $m1)**2 
            + ($array[2] - $m1)**2 + ($array[3] - $m2)**2) / 3;
        $S2_2 = (($array[4] - $m2)**2 + ($array[5] - $m2)**2 
            + ($array[6] - $m2)**2 + ($array[7] - $m2)**2) / 3;
        $Sp_2 = (3 * $S1_2 + 3 * $S2_2) / 6;
        $t = abs($m1 - $m2) / sqrt($Sp_2 * (1/4 +1/4));
        push (@permute, $t);
        &FYS;
    }
@permute = splice(@permute,1,10000);         # number of permutations
$original = $permute[0];
    @sorted = sort {$a <=> $b} @permute;
    $threshold = ($sorted[9989] + $sorted[9990]) / 2;             # P=.001
        if ($original >= $threshold and $original >= 2.447) {     # P=.05
            print OUT "$name, $a, $b, $c, $d, $e, $f, $g, $h\n";
            $gene += 1;
        }
}


sub FYS {                          # Fisher_Yates Shuffle
    my $i;
    for ($i = 7; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($array[$i],$array[$j]) = ($array[$j],$array[$i]);
    }
}

close (OUT);
print $gene."genes\n";


