#!/usr/local/bin/perl

# Relation (2001)

open (IN1,'relation.csv') or die "$!"; @file1 = <IN1>;close (IN1);
open (IN2,'relation.csv') or die "$!";@file2 = <IN2>;close (IN2);
open (OUT,'>sampler2.csv');

# threshold value derived from rel_rand.pl

# $low = -0.978977408494646; $high = 0.996330902712367;       # Lowest,Highest #
# $low = -0.890099507457944; $high = 0.946794155163645;       # P = 0.0001 #
 $low = -0.827761659646113; $high = 0.881589918836945;        # P = 0.001 #
# $low = -0.67945868942256; $high = 0.69981828413013;         # P = 0.01 #


foreach $line1 (@file1) { $line1 =~ s/\n//;
    ($name1,$a1,$b1,$c1,$d1,$e1,$f1,$g1,$h1) = split(/,/,$line1);

    foreach $line2 (@file2) { $line2 =~ s/\n//;
        ($name2,$a2,$b2,$c2,$d2,$e2,$f2,$g2,$h2) = split(/,/,$line2);
        unless ($name1 eq $name2) {
            $m1 = ($a1+$b1+$c1+$d1+$e1+$f1+$g1+$h1)/8;
            $m2 = ($a2+$b2+$c2+$d2+$e2+$f2+$g2+$h2)/8;
            $sd1 = sqrt((($a1-$m1)**2+($b1-$m1)**2 +($c1-$m1)**2+($d1-$m1)**2+($e1-$m1)**2
                  +($f1-$m1)**2+($g1-$m1)**2+($h1-$m1)**2)/7);
            $sd2= sqrt((($a2-$m2)**2+($b2-$m2)**2 +($c2-$m2)**2+($d2-$m2)**2+($e2-$m2)**2
                  +($f2-$m2)**2+($g2-$m2)**2+($h2-$m2)**2)/7);
            $r = (($a1-$m1)*($a2-$m2)+($b1-$m1)*($b2-$m2) +($c1-$m1)*($c2-$m2)
                +($d1-$m1)*($d2-$m2)+($e1-$m1)*($e2-$m2)+($f1-$m1)*($f2-$m2)
                +($g1-$m1)*($g2-$m2)+($h1-$m1)*($h2-$m2))/(7*$sd1*$sd2);
            $r2 = $r**3 / abs($r);
            if ($r2 < $low or $r2 > $high) {
                print OUT "$name1,$name2,$r2\n";
                $relation += 1;
            }
        }
    }
}

print $relation." relations\n";
close (OUT);

