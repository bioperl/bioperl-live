use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 6);

    use_ok 'Bio::Seq';
    use_ok 'Bio::SeqFeature::Primer';
    use_ok 'Bio::SeqFeature::Amplicon';
}

my ($amplicon, $fwd_primer, $rev_primer, $template);




$amplicon = Bio::SeqFeature::Amplicon->new();
isa_ok $amplicon, 'Bio::SeqFeature::Amplicon';
isa_ok $amplicon, 'Bio::SeqFeature::Generic';

$template = Bio::Seq->new( -seq => 'AAAAACCCCCGGGGGTTTTT' );
$fwd_primer = Bio::SeqFeature::Primer->new( -seq => 'ATGG' );
$fwd_primer = Bio::SeqFeature::Primer->new( -seq => 'ATTA' );
ok $amplicon = Bio::SeqFeature::Amplicon->new(
    -seq        => $template,
###    -fwd_primer => $fwd_primer,
);

use Data::Dumper;
print Dumper($amplicon);

##print $amplicon->seq->seq."\n";

#my $ans = $amplicon->fwd_primer;



#ok $amplicon->rev_primer($rev_primer);



#is_deeply $amplicon->fwd_primer(), $fwd_primer;
#is_deeply $amplicon->rev_primer(), $rev_primer;

