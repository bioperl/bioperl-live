#!/user/local/bin/perl

# LOOCV t-test weighted vote (2001)

open (IN, '*.csv') or die "$!";           # file name
@file = <IN>;
close (IN);
@miss = "";

for ($trao = 0; $trao <= 999; ++ $trao) {             # outer
    $miss = 0;
    print "No.".$trao.": ";
    for ($minus =0; $minus <=7; ++$minus)  {        # LOOCV
        $gene = 0;$vote = 0;
        foreach $line (@file) {
            $line =~ s/\n//;
            ($name, $a, $b, $c, $d, $e, $f, $g, $h) = split(/,/,$line);
            @valueminus = ($a, $b, $c, $d, $e, $f, $g, $h);
            $valuesampl = splice(@valueminus,$minus,1);
            @original = @minusone = @valueminus;
            @permute = "";
            if ($minus <= 3) {
                for ($trai = 0; $trai <= 9999; ++ $trai) {     # inner
                    $m1 = ($minusone[0] + $minusone[1] + $minusone[2]) / 3;
                    $m2 = ($minusone[3] + $minusone[4] + $minusone[5] + $minusone[6]) / 4;
                    $S1_2 = (($minusone[0] - $m1)**2 + ($minusone[1] - $m1)**2 
                        + ($minusone[2] - $m2)**2) / 2;
                    $S2_2 = (($minusone[3] - $m2)**2 + ($minusone[4] - $m2)**2 
                        + ($minusone[5] - $m2)**2 + ($minusone[6] - $m2)**2) / 3;
                    $Sp_2 = (2 * $S1_2 + 3 * $S2_2) / 5;
                    $t = abs($m1 - $m2)/ sqrt($Sp_2 * (1/3 +1/4));
                    push (@permute, $t);
                    &FYS7;
                }
            }
            if ($minus >= 4) {  
                for ($trai = 0; $trai <= 9999; ++ $trai) {      # inner
                    $m1 = ($minusone[0] + $minusone[1] + $minusone[2] + $minusone[3]) / 4;
                    $m2 = ($minusone[4] + $minusone[5] + $minusone[6]) / 3;
                    $S1_2 = (($minusone[0] - $m1)**2 + ($minusone[1] - $m1)**2 
                        + ($minusone[2] - $m1)**2 + ($minusone[3] - $m2)**2) / 3;
                    $S2_2 = (($minusone[4] - $m2)**2 + ($minusone[5] - $m2)**2 
                        + ($minusone[6] - $m2)**2) / 2;
                    $Sp_2 = (3 * $S1_2 + 2 * $S2_2) / 5;
                    $t = abs($m1 - $m2)/ sqrt($Sp_2 * (1/4 +1/3));
                    push (@permute, $t);
                    &FYS7;
                }
            }
            @permute = splice(@permute,1,10000);
            $original = $permute[0];
            @sorted = sort {$a <=> $b} @permute;
            $threshold = ($sorted [9989] + $sorted[9990]) / 2;        # P = .001
            if ($original >= $threshold and $original > 2.571) {             # p = .05
                $gene += 1;
                if ($minus <= 3) {
                    $m1 = ($valueminus[0] + $valueminus[1] + $valueminus[2])/3;
                    $m2 = ($valueminus[3] + $valueminus[4] + $valueminus[5] + $valueminus[6])/4;
                    $m = ($m1 + $m2)/2;
                    $sd1 = sqrt ((($valueminus[0]-$m1)**2+ ($valueminus[1]-$m1)**2
                        + ($valueminus[2]-$m1)**2)/2);
                    $sd2 = sqrt ((($valueminus[3]-$m2)**2+ ($valueminus[4]-$m2)**2
                        + ($valueminus[5]-$m2)**2 + ($valueminus[6]-$m2)**2)/3);
                    $wt = abs ($m1 - $m2)/($sd1 + $sd2);
                    $delta1 = abs($m1 - $valuesampl);
                    $delta2 = abs($m2 - $valuesampl);
                    if ($delta1 < $delta2) {
                        $vote += (abs($m - $valuesampl) * $wt);
                    } elsif ($delta1 > $delta2) {
                        $vote += (-1)*(abs($valuesampl - $m) * $wt);
                    }
                }
                if ($minus >= 4) {
                    $m1 = ($valueminus[0] + $valueminus[1] + $valueminus[2] + $valueminus[3])/4;
                    $m2 = ($valueminus[4] + $valueminus[5] + $valueminus[6])/3;
                    $m = ($m1 + $m2)/2;
                    $sd1 = sqrt ((($valueminus[0]-$m1)**2+ ($valueminus[1]-$m1)**2
                        + ($valueminus[2]-$m1)**2 + ($valueminus[3]-$m1)**2)/3);
                    $sd2 = sqrt ((($valueminus[4]-$m2)**2+ ($valueminus[5]-$m2)**2
                        + ($valueminus[6]-$m2)**2)/2);
                    $wt = abs ($m1 - $m2)/($sd1 + $sd2);
                    $delta1 = abs($m1 - $valuesampl);
                    $delta2 = abs($m2 - $valuesampl);
                    if ($delta1 > $delta2) {
                        $vote += (abs($m - $valuesampl) * $wt);
                    } elsif ($delta1 < $delta2) {
                        $vote += (-1)*(abs($valuesampl - $m) * $wt);
                    }
                }
            }
        }

        print $gene."g-";
        if ($vote < 0) {
            $miss += 1;
            print "X, ";
        } else {
            print "O, "
        }
    }                                         # LOOCV end

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
$missthr = ($missorted[48] + $missorted[49])/2;           # P = .05
print "Number of miss in original set: $missoriginal\n";
print "Threshold number of miss from random permutations (P=.05): $missthr\n";
for ($l = 0; $l <= 999; ++$l) {
    if ($missoriginal >= $missorted[$l] ) {
        $count += 1;
    }
}
print "Number of total random permutations :".($trao+1)."\n";
print "Number of random permutations with fewer miss than original set :$count\n";

sub FYS7 {                # Fisher Yates Shuffle
    my $i;
    for ($i = 6; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($minusone[$i],$minusone[$j]) = ($minusone[$j], $minusone[$i]);
    }
}

sub FYS8 {                 # Fisher Yates Shuffle
    my $i;
    for ($i = 7; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($array[$i], $array[$j]) = ($array[$j], $array[$i]);
    }
}
