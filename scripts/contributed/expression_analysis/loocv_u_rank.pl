#!/usr/local/bin/perl

# LOOCV-U-test, predict simple rank (2001)

open (IN, '*.csv') or die "$!";              # file name
@file = <IN>;
close (IN);
@miss = "";

for ($trao = 0; $trao <= 999; ++ $trao) {                 # outer
    $miss = 0;
    $dec = 0;
    print "No.".$trao.": ";
    for ($minus =0; $minus <=7; ++$minus)  {          # LOOCV
        $gene = 0;
        foreach $line (@file) {
            $line =~ s/\n//;
            ($name, $a, $b, $c, $d, $e, $f, $g, $h) = split(/,/,$line);
            @arrayminus = @array = ($a, $b, $c, $d, $e, $f, $g, $h);
            $sample = splice(@arrayminus,$minus,1);
            @original = @minusone = @arrayminus;
            @permute = "";
            for ($p = 7 - $minus; $p > 0; --$p) {
                $minusone[$minus + $p] = $minusone[$minus + $p] - 1;
            }
            if ($minus <= 3) {
                for ($trai = 0; $trai <= 9999; ++ $trai) {        # inner
                    $U = $rankminus[0] + $rankminus[1] + $rankminus[2];
                   # $U1 = 18 - ($minusone[0] + $minusone[1] + $minusone[2]);
                   # $U2 = 22 - ($minusone[3] + $minusone[4] + $minusone[5] + $minusone[6]);
                   # if ($U1 >= $U2) {
                   #    $U = $U2;
                   # } else {
                   #    $U = $U1;
                   # }
                   push (@permute, $U);
                   &FYS7;
                }
                # $rank1 = ($original[0] + $original[1] + $original[2]) / 3;
                # $rank2 = ($original[3] + $original[4] + $original[5] + $original[6]) / 4;
            }
            if ($minus >= 4) {
                for ($trai = 0; $trai <= 9999; ++ $trai) {
                    $U = $rankminus[4] + $rankminus[5] + $rankminus[6];
                   # $U1 = 22 - ($minusone[0] + $minusone[1] + $minusone[2] + $minusone[3]);
                   # $U2 = 18 - ($minusone[4] + $minusone[5] + $minusone[6]);
                   # if ($U1 >= $U2) {
                   #     $U = $U2;
                   # } else {
                   #     $U = $U1;
                   # }
                    push (@permute, $U);
                    &FYS7;
                }
               # $rank1 = ($original[0] + $original[1] + $original[2] + $original[3]) / 4;
               # $rank2 = ($original[4] + $original[5] + $original[6]) / 3;
            }
            @permute = splice(@permute,1,10000);
            $original = $permute[0];
            @sorted = sort {$a <=> $b} @permute;
            $low = ($sorted[249] + $sorted[250]) / 2;
            $high = ($sorted [9749] + $sorted[9750]) / 2;                  # P=.05
            if (($original < $low or $original > $high) and ($original <= 6 or $original >=18)) {
                $gene += 1;
                if ($minus <=3 and $rank1 < $rank2 and $sample < 4.5) {
                    $dec += 1
                }
                if ($minus <=3 and $rank1 > $rank2 and $sample > 4.5) {
                    $dec += 1
                }
                if ($minus <=3 and $rank1 < $rank2 and $sample > 4.5) {
                    $dec -= 1
                }
                if ($minus <=3 and $rank1 > $rank2 and $sample < 4.5) {
                    $dec -= 1
                }
                if ($minus >=4 and $rank1 < $rank2 and $sample > 4.5) {
                    $dec += 1
                }
                if ($minus >=4 and $rank1 > $rank2 and $sample < 4.5) {
                    $dec += 1
                }
                if ($minus >=4 and $rank1 < $rank2 and $sample < 4.5) {
                    $dec -= 1
                }
                if ($minus >=4 and $rank1 > $rank2 and $sample > 4.5) {
                    $dec -= 1
                }
            }
        }

        print $gene."g-";
        if ($dec > 0) {
            print "O, ";
        } else {
            print "X, ";
            $miss += 1;
        }
    }                                              # LOOCV end

    print "miss:$miss\n";
    push (@miss, $miss);

    open (IN, '*.csv') or die "$!";                  # permute 8 cells
    @file = <IN>;
    close (IN);
    open (OUT, '>permuted.csv') or die "$!";
    foreach $line (@file) {
        $line =~ s/\n//;
        ($name, $a, $b, $c, $d, $e, $f, $g, $h) = split(/,/,$line);
        @array = ($a, $b, $c, $d, $e, $f, $g, $h);
        &FYS8;
         print OUT
"$name,$array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$array[6],$array[7]\n";
    }
    close (OUT);
    open (IN, 'permuted.csv') or die "$!";
    @file = <IN>;
    close (IN);
}

@miss = splice(@miss,1,1000);
print "miss:@miss\n";
$missoriginal = $miss[0];
@missorted = sort {$a <=> $b} @miss;
print "misorrt:@missorted\n";
$missthr = ($missorted[49] + $missorted[50])/2;          # P=.05
print "missoriginal:$missoriginal\n";
print "missthreshold:$missthr\n";

sub FYS7 {                       # Fisher Yates Shuffle
    my $i;
    for ($i = 6; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($minusone[$i],$minusone[$j]) = ($minusone[$j], $minusone[$i]);
    }
}

sub FYS8 {                        # Fisher Yates Shuffle
    my $i;
    for ($i = 7; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($array[$i], $array[$j]) = ($array[$j], $array[$i]);
    }
}
