#!/user/local/bin/perl

# LOOCV t-test, compound covariate (2001)

open (IN, '*.csv') or die "$!";             # file name
@file = <IN>;
close (IN);
@miss = "";

for ($trao = 0; $trao <= 999; ++ $trao) {    # outer
    $miss = 0;
    print "No.".$trao.": ";
    for ($minus =0; $minus <=7; ++$minus)  {       # LOOCV
        @CC = "";
        $gene = 0;
        foreach $line (@file) {
            $line =~ s/\n//;
            ($name, $a, $b, $c, $d, $e, $f, $g, $h) = split(/,/,$line);
            @arrayminus = @array = ($a, $b, $c, $d, $e, $f, $g, $h);
            $sample = splice(@arrayminus,$minus,1);
            @original = @minusone = @arrayminus;
            @permute = "";
            if ($minus <= 3) {
                for ($trai = 0; $trai <= 9999; ++ $trai) {        # inner
                    $m1 = ($minusone[0] + $minusone[1] + $minusone[2]) / 3;
                    $m2 = ($minusone[3] + $minusone[4] + $minusone[5] + $minusone[6]) / 4;
                    $S1_2 = (($minusone[0] - $m1)**2 + ($minusone[1] - $m1)**2 
                        + ($minusone[2] - $m2)**2) / 2;
                    $S2_2 = (+ ($minusone[3] - $m2)**2 + ($minusone[4] - $m2)**2 
                        + ($minusone[5] - $m2)**2 + ($minusone[6] - $m2)**2) / 3;
                    $Sp_2 = (2 * $S1_2 + 3 * $S2_2) / 5;
                    $t = abs($m1 - $m2)/ sqrt($Sp_2 * (1/3 +1/4));
                    push (@permute, $t);
                    &FYS7;
                }
            }
            if ($minus >= 4) {  
                for ($trai = 0; $trai <= 9999; ++ $trai) {        # inner
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
            $threshold = ($sorted [94] + $sorted[95]) / 2;           # P = .05
            if ($original >= $threshold and $original > 2.571) {      # P = .05
                $sampleCC += $sample * $permute[0];
                $CC[0] += $original[0] * $permute[0]; 
                $CC[1] += $original[1] * $permute[0]; 
                $CC[2] += $original[2] * $permute[0]; 
                $CC[3] += $original[3] * $permute[0]; 
                $CC[4] += $original[4] * $permute[0]; 
                $CC[5] += $original[5] * $permute[0]; 
                $CC[6] += $original[6] * $permute[0]; 
                $gene += 1;
            }
        }
        print $gene."g-";

        if ($minus <= 3) {                        # count miss
            $meanCC1 = ($CC[0] + $CC[1] + $CC[2]) / 3;
            $meanCC2 = ($CC[3] + $CC[4] + $CC[5] + $CC[6]) / 4;
            $delta1 = abs ($meanCC1 - $sampleCC);
            $delta2 = abs ($meanCC2 - $sampleCC);            
            if ($delta1 > $delta2) {
                $miss += 1; print "X, ";
            } else {
                print "O, ";
            }
        }
        if ($minus >= 4) {
            $meanCC1 = ($CC[0] + $CC[1] + $CC[2] + $CC[3]) /4;
            $meanCC2 = ($CC[4] + $CC[5] + $CC[6]) / 3;
            $delta1 = abs ($meanCC1 - $sampleCC);
            $delta2 = abs ($meanCC2 - $sampleCC);
            if ($delta1 < $delta2) {
                $miss += 1;
                print "X, ";
            } else {
                print "O, ";
            }
        }      
    }                             # LOOCV end
    print "miss:$miss\n";
    push (@miss, $miss);

    open (IN, '*.csv') or die "$!";       # permute 8 cells
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
$missthr = ($missorted[48] + $missorted[49])/2;             # P = .05
print "Number of miss in original set: $missoriginal\n";
print "Threshold number of miss from random permutations (P=.05): $missthr\n";
for ($l = 0; $l <= 999; ++$l) {
    if ($missoriginal >= $missorted[$l] ) {
        $count += 1;
    }
}
print "Number of total random permutations :".($trao+1)."\n";
print "Number of random permutations with fewer miss than original set :$count\n";

sub FYS7 {                   # Fisher Yates Shuffle
    my $i;
    for ($i = 6; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($minusone[$i],$minusone[$j]) = ($minusone[$j], $minusone[$i]);
    }
}

sub FYS8 {                   # Fisher Yates Shuffle
    my $i;
    for ($i = 7; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($array[$i], $array[$j]) = ($array[$j], $array[$i]);
    }
}
