#!/usr/local/bin/perl

# LOOCV U-test, weighted vote (2001)

open (IN, '*.csv') or die "$!";          # file name
@file = <IN>;
close (IN);
@miss = "";

for ($trao = 0; $trao <= 999; ++ $trao) {        # outer
    $miss = 0;
    print "No.".$trao.": \n";
    for ($minus =0; $minus <=7; ++$minus)  {        # LOOCV
        $vote = 0;
        $gene = 0;
        foreach $line (@file) {
            $line =~ s/\n//;
            ($name, $a, $b, $c, $d, $e, $f, $g, $h) = split(/,/,$line);
            @valueminus = ($a, $b, $c, $d, $e, $f, $g, $h);
            $valuesampl = splice(@valueminus,$minus,1);
            @rankminus = &rank(@valueminus);
            @permute = "";

            if ($minus <= 3) {
                for ($trai = 0; $trai <= 9999; ++ $trai) {     # inner
                    $U = $rankminus[0] + $rankminus[1] + $rankminus[2];
                    push (@permute, $U);
                    &FYS7;
                }
            }
            if ($minus >= 4) {
                for ($trai = 0; $trai <= 9999; ++ $trai) {      # inner
                    $U = $rankminus[4] + $rankminus[5] + $rankminus[6];
                    push (@permute, $U);
                    &FYS7;
                }
            }
            @permute = splice(@permute,1,10000);
            $original = $permute[0];
            @sorted = sort {$a <=> $b} @permute;
            $low = ($sorted[249] + $sorted[250]) / 2;
            $high = ($sorted[9749] + $sorted[9750]) / 2;         # P=.05
            if (($original < $low and $original<=3) or ($original > $high and $original >=15)) {
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
        print "vote: ".$vote.", ";
        if ($vote < 0) {
            $miss += 1;
            print "X, \n";
        } else {
            print "O, \n";
        }
    }                                    # LOOCV end

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
}                                                        # outer end

@miss = splice(@miss,1,1000);                            # count miss
print "miss:@miss\n";
$missoriginal = $miss[0];
@missorted = sort {$a <=> $b} @miss;
print "misorrt:@missorted\n";
$missthr = ($missorted[49] + $missorted[50])/2;                           # P=.05
print "Number of miss in original set: $missoriginal\n";
print "Threshold number of miss from random permutations (P=.05): $missthr\n";
for ($l = 0; $l <= 999; ++$l) {
    if ($missoriginal >= $missorted[$l] ) {
        $count += 1;
    }
}
print "Number of total random permutations :".($trao+1)."\n";
print "Number of random permutations with fewer miss than original set :$count\n";

sub FYS7 {                         # Fisher Yates Shuffle
    my $i;
    for ($i = 6; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($rankminus[$i],$rankminus[$j]) = ($rankminus[$j], $rankminus[$i]);
    }
}

sub FYS8 {                          # Fisher Yates Shuffle
    my $i;
    for ($i = 7; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($array[$i], $array[$j]) = ($array[$j], $array[$i]);
    }
}

sub rank {          # convert raw data to rank
    my @rank = map { [0,$_]; } @_;
    my $i = 0; for( sort {$a->[1] <=> $b->[1]} @rank ) {
        $_->[0] = $i++;
    }
    return map { $_->[0]; } @rank;
}
