use strict;
use warnings;

BEGIN {
    use Bio::Root::Test;
    test_begin( -tests => 8, -requires_module => 'List::MoreUtils');
}

use_ok 'Bio::Tools::SeqPattern::Backtranslate';

can_ok 'Bio::Tools::SeqPattern::Backtranslate', '_reverse_translate_motif';

my @data;
while (<DATA>) {
    chomp;
    push @data, [split];
}

foreach my $line (@data) {
    if ( $line->[1] ) {
        is _reverse_translate_motif( $line->[0] ), $line->[1];
    } else {
        throws_ok { _reverse_translate_motif( $line->[0] ) }
        qr/syntax token/;
    }
}

__DATA__
MAEELKAVAP ATGGCNGARGARYTNAARGCNGTNGCNCCN
LKGHB[WhYq]Q YTNAARGGNCAYRAYYRNCAR
(LK){2,3}[^GHB][WHYQ]Q (YTNAAR){2,3}HBNYRNCAR
LK[^GHB][WHYQ]Q YTNAARHBNYRNCAR
(LK){2,3}[^GHB][WHYQ]QX.X (YTNAAR){2,3}HBNYRNCARNNNNNNNNN
/&%/&%(/)(/%
