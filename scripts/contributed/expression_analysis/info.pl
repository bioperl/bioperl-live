#!/usr/local/bin/perl

# TNoM, Infoscore (2001)

open (PREIN, '*.csv') or die "$!";      # initial raw data file
@prefile = <PREIN>;
close (PREIN);
open (OUT, '>rank.csv') or die "$!";      # convert to rank
foreach $preline (@prefile) {
    $preline =~ s/\n//;
    ($name, $a, $b, $c, $d, $e, $f, $g, $h) = split(/,/,$preline);
    @a = ($a, $b, $c, $d, $e, $f, $g, $h);
    @rr = &rank(@a);
    print OUT "$name,$rr[0],$rr[1],$rr[2],$rr[3],$rr[4],$rr[5],$rr[6],$rr[7]\n";
}
close OUT;

open (IN, 'rank.csv') or die "$!";
@file = <IN>;
close (IN);

# open (TNOM, '>result.csv') or die "$!";           # for TNoM
open (INFO, '>result.csv') or die "$!";            # for Infoscore

foreach $line (@file) {
    $line =~ s/\n//;
    ($name, $a, $b, $c, $d, $e, $f, $g, $h) = split(/,/,$line);
    @group[$a,$b,$c,$d,$e,$f,$g,$h] = (S,S,S,S,N,N,N,N);    # indicate "S" and "N" group
   # @pertnom ="";                       # for TNoM
    @perinfo = "";                       # for Infoscore

    for ($tra = 0; $tra <= 9999; ++ $tra) {              # number of permutations
        $seqlt = $group[0].$group[1].$group[2].$group[3];
        $seqrt = $group[4].$group[5].$group[6].$group[7];
        if ($seqlt eq ("SSSS" or "NNNN")) {
            $TNoM = 0;$info = 0;
        } elsif  (($seqlt eq "SSSN") or ($seqlt eq "NNNS")) {
            $TNoM = 1; $info = 0.1358;
        } elsif  (($seqlt eq "SSNS") or ($seqlt eq "SNSS") or ($seqlt eq "NSSS")) {
            if ($group[4] eq "S") {
                $TNoM = 1; $info = 0.1358;
            } elsif ($group[4] eq "N") {
                $TNoM = 2; $info = 0.2442;
            }
        } elsif  (($seqlt eq "NNSN") or ($seqlt eq "NSNN") or ($seqlt eq "SNNN")) {
            if ($group[4] eq "N") {
                $TNoM = 1; $info = 0.1358;
            } elsif ($group[4] eq "S") {
                $TNoM = 2; $info = 0.2442;
            }
        } elsif (($seqlt eq "SSNN") or ($seqlt eq "NNSS")) {
            $TNoM = 2; $info = 0.2073;
        } elsif (($seqrt eq "SSNN") or ($seqlt eq "NNSS")) {
            $TNoM = 2; $info = 0.2073;
        } else {
            $TNoM =3; $info = .2863;
        }

   # push (@pertnom, $TNoM);                # for TNoM
    push (@perinfo, $info);                # for Infoscore
    &FYS;
    }

 #   @mpertnom = splice(@pertnom,1,10000);        # for TNoM
 #   $oritnom = $mpertnom[0];
 #   @sortedt = sort {$a <=> $b} @mpertnom;
 #   $thret = ($sortedt[8] + $sortedt[9]) / 2;     # P=.001
 #   if ($oritnom <= $thret) {
 #       print TNOM "$name, $a, $b, $c, $d, $e, $f, $g, $h\n";
 #       $tnom += 1;
 #   }

    @mperinfo = splice(@perinfo,1,1000);
    $oriinfo = $mperinfo[0];
    @sortedi = sort {$a <=> $b} @mperinfo;
    $threi = ($sortedi[8] + $sortedi[9]) / 2;
        if ($oriinfo <= $threi) {
              print INFO "$name, $a, $b, $c, $d, $e, $f, $g, $h\n";
              $infonum += 1;
        }
}

# print "Selected genes by TNoM :".$tnom."\n";     # for TNoM
print "Selected genes by Info :".$infonum."\n";   # Infoscore

sub FYS {                         #  Fisher Yates Shuffle
    my $i;
    for ($i = 7; $i>=0; --$i ) {
        my $j = int rand($i + 1);
        next if $i == $j;
        ($group[$i],$group[$j]) = ($group[$j],$group[$i]);
    }
}

sub rank {
    my @rank = map {
        [0,$_];
    }
    @_;
    my $i = 0;
    for( sort {$a->[1] <=> $b->[1]} @rank ) {
        $_->[0] = $i++;
    }
    return map {
        $_->[0];
    }
    @rank;
}

# close (TNOM);
close (INFO);
