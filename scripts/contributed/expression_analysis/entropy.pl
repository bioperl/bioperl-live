#!/usr/local/bin/perl

# Entropy (2001)

open (ENT,'*.csv') or die "$!";     # file name
@allfeat = <ENT>;
close (ENT);

open (ENTOUT,'>entropy.csv');

foreach $feat (@allfeat) {
    $feat =~ s/\n//;
    ($name,$a,$b,$c,$d,$e,$f,$g,$h) = split(/,/,$feat);
    @feature = ($a,$b,$c,$d,$e,$f,$g,$h);
    @sortfeat = sort {$a <=> $b} @feature;
    $interval = ($sortfeat[7] - $sortfeat[0])/9;
    @compart = (0,0,0,0,0,0,0,0,0,0);
    for ($fe = 0;$fe <= 7;++$fe) {
        $cell = int (($feature[$fe]-$sortfeat[0])/$interval);
        $compart[$cell] += 0.125;
    }
    $entropy = -($compart[0]*log($compart[0]+0.00001)/log(2)
        +$compart[1]*log($compart[1]+0.00001)/log(2)
        +$compart[2]*log($compart[2]+0.00001)/log(2)
        +$compart[3]*log($compart[3]+0.00001)/log(2)
        +$compart[4]*log($compart[4]+0.00001)/log(2)
        +$compart[5]*log($compart[5]+0.00001)/log(2)
        +$compart[6]*log($compart[6]+0.00001)/log(2)
        +$compart[7]*log($compart[7]+0.00001)/log(2)
        +$compart[8]*log($compart[8]+0.00001)/log(2)
        +$compart[9]*log($compart[9]+0.00001)/log(2));
    push (@entropy,$entropy);
    print ENTOUT "$name,$a,$b,$c,$d,$e,$f,$g,$h,$entropy\n";
}

close (ENTOUT);

# Cut features genes with lowest 5% entropy

@sortent = sort{$a <=> $b} @entropy;
$entther = ($sortent[*] + $sortent[*])/2;   # *:lowest 5% number of features (genes)

open (CENT,'entropy.csv') or die "$!";
@allfeat2 = <CENT>;
close (CENT);

open (CENTOUT,'>relation.csv');
open (DROPOUT,'>dropout.csv');

foreach $feat2 (@allfeat2) {
    $feat2 =~ s/\n//;
    ($name,$a,$b,$c,$d,$e,$f,$g,$h,$entropy) = split(/,/,$feat2);
    # print "$entther\n";
    if ($entropy >= $entther) {
        print CENTOUT "$name,$a,$b,$c,$d,$e,$f,$g,$h\n";
        $entgene +=1;
    } else {
        print DROPOUT "$name,$a,$b,$c,$d,$e,$f,$g,$h\n";
        $dropgene +=1;
    }
}

print "Gene with sufficient entoropy:$entgene\n";
print "Dropout gene:$dropgene\n";

close (CENTOUT);
close (DROPOUT);
