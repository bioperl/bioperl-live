#!/user/local/bin/perl

# LOOCV-discrimination score (2001)

open (IN, '*.csv') or die "$!";             # file name
@file = <IN>;
close (IN);
@miss = "";

for ($trao = 0; $trao <= 999; ++ $trao) {     # number of permutations (outer)
    $miss = 0;
    print "No.".$trao.": ";
    for ($minus =0; $minus <=7; ++$minus) {     # LOOCV
        $vote = 0;
        $gene = 0;
        foreach $line (@file) {
            $line =~ s/\n//;
            ($name, $a, $b, $c, $d, $e, $f, $g, $h) = split(/,/,$line);
            @valueminus = ($a, $b, $c, $d, $e, $f, $g, $h);
            $valuesampl = splice(@valueminus,$minus,1);
            @valueorig = @valueminus;
            @perwt = "";
            $wt = 0;
            if ($minus <= 3) {
                for ($trai = 0; $trai <= 9999; ++ $trai) {     # number of permutations (inner)
                    $m1 = ($valueminus[0] + $valueminus[1] + $valueminus[2])/3;
                    $m2 = ($valueminus[3] + $valueminus[4] + $valueminus[5] + $valueminus[6])/4;
                    $m = ($m1 + $m2)/2;
                    $sd1 = sqrt (($valueminus[0]-$m1)**2+ ($valueminus[1]-$m1)**2
                        + ($valueminus[2]-$m1)**2)/2;
                    $sd2 = sqrt (($valueminus[3]-$m2)**2+ ($valueminus[4]-$m2)**2
                        + ($valueminus[5]-$m2)**2 + ($valueminus[6]-$m2)**2)/3;
                    $wt = abs ($m1 - $m2)/($sd1 + $sd2);
                    push (@perwt, $wt);
                    &FYS7;
                }
            }

            if ($minus >= 4) {
                for ($trai = 0; $trai <= 9999; ++ $trai) {      # number of permutations (inner)
                    $m1 = ($valueminus[0] + $valueminus[1] + $valueminus[2]+$valueminus[3])/4;
                    $m2 = (+ $valueminus[4] + $valueminus[5] + $valueminus[6])/3;
                    $m = ($m1 + $m2)/2;
                    $sd1 = sqrt (($valueminus[0]-$m1)**2+ ($valueminus[1]-$m1)**2
                        + ($valueminus[2]-$m1)**2 + ($valueminus[3]-$m1)**2)/3;
                    $sd2 = sqrt (($valueminus[4]-$m2)**2+ ($valueminus[5]-$m2)**2
                        + ($valueminus[6]-$m2)**2)/2;
                    $wt = abs ($m1 - $m2)/($sd1 + $sd2);
                    push (@perwt, $wt);
                    &FYS7;
                }
            }

            @perwt = splice(@perwt,1,10000);
            $origwt = $perwt[0];
            @sortwt = sort {$a <=> $b} @perwt;
            $threwt = ($sortwt[9989] + $sortwt[9990]) / 2;             # P=.001
            if ($origwt >= $threwt) {
                if ($minus <= 3) {
                    $m1 = ($valueorig[0] + $valueorig[1] + $valueorig[2])/3;
                    $m2 = ($valueorig[3] + $valueorig[4] + $valueorig[5] + $valueorig[6])/4;
                    $m = ($m1 + $m2)/2;
                    $sd1 = sqrt ((($valueorig[0]-$m1)**2+ ($valueorig[1]-$m1)**2
                        + ($valueorig[2]-$m1)**2)/2);
                    $sd2 = sqrt ((($valueorig[3]-$m2)**2+ ($valueorig[4]-$m2)**2
                        + ($valueorig[5]-$m2)**2 + ($valueminus[6]-$m2)**2)/3);
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
                    $m1 = ($valueorig[0] + $valueorig[1] + $valueorig[2] + $valueorig[3])/4;
                    $m2 = ($valueorig[4] + $valueorig[5] + $valueorig[6])/3;
                    $m = ($m1 + $m2)/2;
                    $sd1 = sqrt ((($valueorig[0]-$m1)**2+ ($valueorig[1]-$m1)**2
                        + ($valueorig[2]-$m1)**2 + ($valueorig[3]-$m1)**2)/3);
                    $sd2 = sqrt ((($valueorig[4]-$m2)**2+ ($valueorig[5]-$m2)**2
                        + ($valueorig[6]-$m2)**2)/2);
                    $wt = abs ($m1 - $m2)/($sd1 + $sd2);
                    $delta1 = abs($m1 - $valuesampl);
                    $delta2 = abs($m2 - $valuesampl);
                    if ($delta1 > $delta2) {
                        $vote += (abs($m - $valuesampl) * $wt);
                    } elsif ($delta1 < $delta2) {
                        $vote += (-1)*(abs($valuesampl - $m) * $wt);
                    }
                }

                $gene += 1;
            }
        }
        print $gene."g-";

        if ($vote < 0) {
            $miss += 1;
            print "X, ";
        } else {
            print "O, ";
        }
    }                                     # LOOCV end

    print "miss:$miss\n";
    push (@miss, $miss);

    open (IN, '*.csv') or die "$!";                 # permute 8 cells
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
$missthr = ($missorted[48] + $missorted[49])/2;                           # P=.05
print "Number of miss in original set: $missoriginal\n";
print "Threshold number of miss from random permutations (P=.05): $missthr\n";
for ($l = 0; $l <= 999; ++$l) {
    if ($missoriginal >= $missorted[$l] ) {
        $count += 1;
    }
}
print "Number of total random permutations :".($trao+1)."\n";
print "Number of random permutations with fewer miss than original set :$count\n";

sub FYS7 {                          # Fisher Yates Shuffle
    my $i;
    for ($i = 6; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($valueminus[$i],$valueminus[$j]) = ($valueminus[$j], $valueminus[$i]);
    }
}

                      
sub FYS8 {                           # Fisher Yates Shuffle
    my $i;
    for ($i = 7; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($array[$i], $array[$j]) = ($array[$j], $array[$i]);
    }
}
