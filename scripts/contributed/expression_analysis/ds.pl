#!/usr/local/bin/perl

# Discrimination score (2001)

open (IN, '*.csv') or die "$!";       # file name
@file = <IN>;
close (IN);
open (OUT, '>result.csv') or die "$!";

foreach $line (@file) {
    $line =~ s/\n//;
    ($name, $a, $b, $c, $d, $e, $f, $g, $h) = split(/,/,$line);
    @array = ($a, $b, $c, $d, $e, $f, $g, $h);
    @perwt = "";$wt = 0;

    for ($trai = 0; $trai <= 9999; ++ $trai) {          # number of permutations
        $m1 = ($array[0] + $array[1] + $array[2]+ $array[3])/3;
        $m2 = ($array[4] + $array[5] + $array[6] + $array[7])/4;
        $sd1 = sqrt (($array[0]-$m1)**2+ ($array[1]-$m1)**2+ ($array[2]-$m1)**2
               + ($array[3]-$m1)**2)/3;
        $sd2 = sqrt (($array[4]-$m2)**2+ ($array[5]-$m2)**2+ ($array[6]-$m2)**2
               + ($array[7]-$m2)**2)/3;
        $wt = abs ($m1 - $m2)/($sd1 + $sd2);
        push (@perwt, $wt);
        &FYS;
    }

    @perwt = splice(@perwt,1,10000);                    # number of permutations
    $origwt = $perwt[0];
    @sortwt = sort {$a <=> $b} @perwt;
    $threwt = ($sortwt[9989] + $sortwt[9990]) / 2;                # P=.001
    if ($origwt > $threwt) {
        print OUT "$name,$a,$b,$c,$d,$e,$f,$g,$h\n"; $gene += 1;
    }
}

print $gene."genes\n";

close (OUT);

sub FYS {                       #  Fisher Yates Shuffle
    my $i;
    for ($i = 7; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($array[$i], $array[$j]) = ($array[$j], $array[$i]);
    }
}
