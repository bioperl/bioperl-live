#!/usr/local/bin/perl

# Categorization (2001)

open (IN, '*.csv') or die "$!";          # file name
@file = <IN>;
close (IN);
open (OUT, '>result.csv') or die "$!";

foreach $line (@file) { 
    $line =~ s/\n//;
    ($name, $a, $b, $c, $d, $e, $f, $g, $h) = split(/,/,$line);
    @valueminus = @value = ($a, $b, $c, $d, $e, $f, $g, $h);

    for ($trai = 0; $trai <= 9999; ++ $trai) {            # number of permutations
        $uplt=0;$downlt=0;$nclt=0;$uprt=0;$downrt=0;$ncrt=0;$cat=0;
        for ($s=0;$s<=3;++$s) {
            if ($valueminus[$s]>=1) {
                $uplt += 1;
            }elsif($valueminus[$s]<=-1) {
                $downlt += 1;
            }else{
                $nclt += 1;
            }
        }
        for ($s=4;$s<=7;++$s) {
            if ($valueminus[$s]>=1) {
                $uprt += 1;
            }elsif($valueminus[$s]<=-1) {
                $downrt += 1;
            }else{
                $ncrt += 1;
            }
        }
        if ($uplt<$uprt) {
            $cat+=$uplt;
        }else{
            $cat+=$uprt;
        }
        if ($downlt<$downrt) {
            $cat+=$downlt;
        }else{
            $cat+=$downrt;
        }
        if ($nclt<$ncrt) {
            $cat+=$nclt;
        }else{
            $cat+=$ncrt;
        }
        push (@percat, $cat);
        &FYS;
    }

    @percat = splice(@percat,1,10000);              # number of permutations
    $origcat = $percat[0];
    @sortcat = sort {$a <=> $b} @percat;
    $threcat = ($sortcat[8] + $sortcat[9]) / 2;         # P=.001
    if ($origcat <= $threcat) {
        print OUT "$name,$value[0],$value[1],$value[2],$value[3],$value[4],$value[5],$value[6],$value[7]\n";
        $gene += 1;
    }
}

prinit $gene."genes\n";

sub FYS {                               # Fisher Yates Shuffle
    my $i;
    for ($i = 7; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($valueminus[$i],$valueminus[$j]) = ($valueminus[$j], $valueminus[$i]);
    }
}

