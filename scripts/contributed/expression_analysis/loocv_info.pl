#!/user/local/bin/perl

# LOOCV-infoscore (TNoM) (2001)

open (IN, '*.csv') or die "$!";             # file name
@file = <IN>;
close (IN);
@miss = "";

for ($trao = 0; $trao <= 999; ++ $trao) {            # outer
    $miss = 0; 
    print "No.".$trao.": ";
    for ($minus =0; $minus <=7; ++$minus)  {          # LOOCV
        $vote = 0;
        $gene = 0;
        foreach $line (@file) {
            $line =~ s/\n//;
            ($name, $a, $b, $c, $d, $e, $f, $g, $h) = split(/,/,$line);
            @valueminus = ($a, $b, $c, $d, $e, $f, $g, $h);
            @rank = &rank(@valueminus);
            @groupminus[$rank[0], $rank[1], $rank[2], $rank[3], $rank[4], $rank[5], $rank[6], $rank[7]] = (S,S,S,S,N,N,N,N);

            $valuesampl = splice(@valueminus,$minus,1);
            $groupsampl = splice(@groupminus,$minus,1); 

            @pertnom = "";
            @perinfo = "";
            if ($minus <= 3) {
                for ($trai = 0; $trai <= 9999; ++ $trai) {        # inner
                    $seqlt = $groupminus[0].$groupminus[1].$groupminus[2];
                    $seqrt = $groupminus[3].$groupminus[4].$groupminus[5].$groupminus[6];
                    if ($seqlt eq "SSS" or ($seqlt eq "NNN" and $groupminus[4] eq "N")) {
                        $TNoM = 0;
                        $info = 0;
                    } elsif (($seqlt eq "NNN" or $seqlt eq "SSN" or $seqlt eq "SNS" or $seqlt eq "NSS") and ($groupminus[4] eq "S")) {
                        $TNoM = 1;
                        $info = 0.2442;
                    } elsif (($seqlt eq "SSN") and ($groupminus[4] eq "N")) {
                        $TNoM = 1;
                        $info = 0.2173;
                    } elsif (($seqlt eq "SNS" or $seqlt = "NSS") and ($groupminus[4] eq "N")) {
                        $TNoM = 2;
                        $info = 0.5206;
                    } elsif (($seqlt eq "SNN" or $seqlt = "NSN" or $seqlt = "NNS") and ($seqrt eq "SSNN")) {
                        $TNoM = 2;
                        $info = 0.2923;
                    } elsif (($seqlt eq "SNN" or $seqlt = "NSN" or $seqlt = "NNS") and ($seqrt eq "NNSS")) {
                        $TNoM = 1;
                        $info = 0.2173;
                    } elsif (($seqlt eq "SNN" or $seqlt = "NSN") and ($groupminus[4].$groupminus[5] eq "SN")) {
                        $TNoM = 3;
                        $info = 0.5774;
                    } elsif (($seqlt eq "SNN" or $seqlt = "NSN") and ($groupminus[4].$groupminus[5] eq "NS")) {
                        $TNoM = 2;
                        $info = 0.5206;
                    } elsif ($seqlt eq "NNS" and $groupminus[4].$groupminus[5] eq "SN") {
                        $TNoM = 2;
                        $info = 0.2923;
                    } elsif ($seqlt eq "NNS" and $seqrt eq "NSSN") {
                        $TNoM = 2;
                        $info = 0.5206;
                    } elsif ($seqlt eq "NNS" and $seqrt eq "NSNS") {
                        $TNoM = 2;
                        $info = 0.2923;
                    }

                    push (@pertnom, $TNoM);
                    push (@perinfo, $info);
                    &FYS7;
                }
            }

            if ($minus >= 4) {
                for ($trai = 0; $trai <= 9999; ++ $trai) {          # inner
                    $seqlt = $groupminus[0].$groupminus[1].$groupminus[2].$groupminus[3]; 
                    $seqrt = $groupminus[4].$groupminus[5].$groupminus[6];
                    if ($seqlt eq "SSSS" and $seqlt eq "NNNS") {
                        $TNoM = 0;
                        $info = 0;
                    } elsif (($seqlt eq "SSSN" or $seqlt eq "SSNS" or $seqlt eq "SNSS" or $seqlt eq "NSSS") and ($groupminus[4] eq "S")) {
                        $TNoM = 1;
                        $info = 0.2173;
                    } elsif ($seqlt eq "SSSN" and $groupminus[4] eq "N") {
                        $TNoM = 1;
                        $info = 0.2442;
                    } elsif ($seqlt eq "SSNS" and $groupminus[4] eq "N") {
                        $TNoM = 2;
                        $info = 0.2923;
                    } elsif (($seqlt eq "SNSS" or $seqlt eq "NSSS") and $groupminus[4] eq "N") {
                        $TNoM = 2;
                        $info = 0.5206;
                    } elsif ($seqlt.$seqrt eq "SSNNSSN" or $seqlt.$seqrt eq "SSNNNSS" 
                        or $seqlt.$seqrt eq "SNSNNSS" or $seqlt.$seqrt eq "SNNSNSS" 
                        or $seqlt.$seqrt eq "NSNSNSS" or $seqlt.$seqrt eq "NSSNNSS") {
                        $TNoM = 2;
                        $info = 0.2923;
                    } elsif ($seqlt.$seqrt eq "SSNNSNS" or $seqlt.$seqrt eq "SNSNSSN" 
                        or $seqlt.$seqrt eq "SNSNSNS" or $seqlt.$seqrt eq "NSSNSSN" 
                        or $seqlt.$seqrt eq "NSSNSNS") {
                        $TNoM = 3;
                        $info = 0.5774;
                    } elsif ($seqlt.$seqrt eq "SNNSSSN" or $seqlt.$seqrt eq "SNNSSNS" 
                        or $seqlt.$seqrt eq "NSNSSSN" or $seqlt.$seqrt eq "NSNSSNS") {
                        $TNoM = 2;
                        $info = 0.5206;
                    } elsif ($seqlt.$seqrt eq "NNSSSSN" or $seqlt.$seqrt eq "NNSSSNS" 
                        or $seqlt.$seqrt eq "NNSSNSS" or $seqlt.$seqrt eq "NNSNSSS") {
                        $TNoM = 1;
                        $info = 0.2173;
                    } elsif ($seqlt.$seqrt eq "SNNSSSS" or $seqlt.$seqrt eq "NSNNSSS") {
                        $TNoM = 1;
                        $info = 0.2442;
                    }
                    push (@pertnom, $TNoM);
                    push (@perinfo, $info);
                    &FYS7;    
	        }
	    }

            @perinfo = splice(@perinfo,1,10000);
            $originfo = $perinfo[0];
            @sortinfo = sort {$a <=> $b} @perinfo;
            $threinfo = ($sortinfo [8] + $sortinfo[9]) / 2;           # P=.001
            if ($originfo <= $threinfo) {
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

                $gene += 1;
            }
        }

        print $gene."g-";

        if ($vote < 0) {
            $miss += 1;
            print "X, ";
        } else {
            print "O, "
        }
    }                                                     # LOOCV end

    print "miss:$miss\n";
    push (@miss, $miss);

    open (IN, '*.csv') or die "$!";                     # permute 8 cells
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

@miss = splice(@miss,1,1000);                                    # number of miss count
print "miss:@miss\n";
$missoriginal = $miss[0];
@missorted = sort {$a <=> $b} @miss;
print "misorrt:@missorted\n";
$missthr = ($missorted[48] + $missorted[49])/2;                      # P=.05
print "Number of miss in original set: $missoriginal\n";
print "Threshold number of miss from random permutations (P=.05): $missthr\n";
for ($l = 0; $l <= 999; ++$l) {
    if ($missoriginal >= $missorted[$l] ) {
        $count += 1;
    }
}

print "Number of total random permutations :".($trao+1)."\n";
print "Number of random permutations with fewer miss than original set :$count\n";

sub FYS7 {                    # Fisher Yates Shuffle
    my $i;
    for ($i = 6; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($groupminus[$i],$groupminus[$j]) = ($groupminus[$j], $groupminus[$i]);
     }
}

sub FYS8 {                    # Fisher Yates Shuffle
    my $i;
    for ($i = 7; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($array[$i], $array[$j]) = ($array[$j], $array[$i]);
    }
}

sub rank {                    # covert raw data to rank
    my @rank = map { [0,$_]; } @_;
    my $i = 0;
    for( sort {$a->[1] <=> $b->[1]} @rank ) {
        $_->[0] = $i++;
    }
    return map { $_->[0]; } @rank;
}


