#!/usr/local/bin/perl

# Random relation (2001)

for ($tra = 0; $tra <= 99; ++$tra) {             # 100 permutations
    open (IN,'*.csv');                          # file name
    @file = <IN>;
    close (IN);
    open (PER,'>permuted.csv');

    foreach $line (@file) {
        $line =~ s/\n//; ($name,$a,$b,$c,$d,$e,$f,$g,$h) = split(/,/,$line);
        @per = ($a,$b,$c,$d,$e,$f,$g,$h);
        &FYS8;
        print PER "$name,$per[0],$per[1],$per[2],$per[3],$per[4],$per[5],$per[6],$per[7]\n";
    }
    close (PER);

    open (INF,'permuted.csv') or die "$!";
    @file1 = <INF>; close (INF);
    open (INS,'permuted.csv') or die "$!";
    @file2 = <INS>; close (INS);

    foreach $line1 (@file1) {
        $line1 =~ s/\n//;
        ($name1,$a1,$b1,$c1,$d1,$e1,$f1,$g1,$h1) = split(/,/,$line1);
        @inr2 =""; $trai = 1;
        foreach $line2 (@file2) {
            $line2 =~ s/\n//;
            ($name2,$a2,$b2,$c2,$d2,$e2,$f2,$g2,$h2) = split(/,/,$line2);
            unless ($name1 eq $name2) {
                $m1 = ($a1+$b1+$c1+$d1+$e1+$f1+$g1+$h1)/8;
                $m2 = ($a2+$b2+$c2+$d2+$e2+$f2+$g2+$h2)/8;
                $sd1 = sqrt((($a1-$m1)**2+($b1-$m1)**2 +($c1-$m1)**2+($d1-$m1)**2+($e1-$m1)**2
                    +($f1-$m1)**2+($g1-$m1)**2+($h1-$m1)**2)/7);
                $sd2= sqrt((($a2-$m2)**2+($b2-$m2)**2 +($c2-$m2)**2+($d2-$m2)**2+($e2-$m2)**2
                    +($f2-$m2)**2+($g2-$m2)**2+($h2-$m2)**2)/7);
                $r = (($a1-$m1)*($a2-$m2)+($b1-$m1)*($b2-$m2) +($c1-$m1)*($c2-$m2)+($d1-$m1)*($d2-$m2)
                    +($e1-$m1)*($e2-$m2)+($f1-$m1)*($f2-$m2) +($g1-$m1)*($g2-$m2)+($h1-$m1)*($h2-$m2))
                    /(7*$sd1*$sd2);
                $r2 = $r**3 / abs($r);
                push (@inr2, $r2);
            }
        }
        @midr2 = splice(@inr2,1,*);                       # *:(number-1) of features (genes)
        @sortmidr2 = sort {$a<=>$b} @midr2;
        for ($mabiki = 1;$mabiki <= **;$mabiki += 2) {    # **:(number-2) of features (genes)
            push (@r2, $sortmidr2[$mabiki]);
        }
        @r2 = splice(@r2,1,***);                          # ***:(number-1)/2 of features (genes)
        push (@low, $r2[0]);push (@low, $r2[1]);push (@low, $r2[2]);push (@low, $r2[3]);
        push (@high, $r2[-4]);push (@high, $r2[-3]);push (@high, $r2[-2]);push (@high, $r2[-1]);
        @r2 = "";
    }
}

print scalar(@low).", ".scalar(@high)."\n";

@sortlow = sort {$a<=>$b} @low;
@sorthigh = sort {$a<=>$b} @high;
$lowest = $sortlow[0]; $highest = $sorthigh[-1];
$low0001 = ($sortlow[] + $sortlow[])/2; $high0001 = ($sorthigh[] + $sorthigh[])/2;
$low001 = ($sortlow[] + $sortlow[])/2;  $high001 = ($sorthigh[] + $sorthigh[])/2;
$low01 = ($sortlow[] + $sortlow[])/2; $high01 = ($sorthigh[] + $sorthigh[])/2;

print "Lowest :".$lowest."\n"."Highest : ".$highest."\n\n";
print "Low0001: ".$low0001."\n". "High0001: ".$high0001."\n\n";
print "Low001: ".$low001."\n". "High001: ".$high001."\n\n";
print "Low01: ".$low01."\n". "high01: ".$high01."\n";

sub FYS8 {
    my $i;
    for ($i = 7; $i >=0; --$i) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($per[$i],$per[$j]) = ($per[$j],$per[$i]);
    }
}
