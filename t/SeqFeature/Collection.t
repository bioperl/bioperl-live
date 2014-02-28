# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(
        -tests => 24,
        -requires_module => 'DB_File'
    );

    use_ok('Bio::SeqFeature::Collection');
    use_ok('Bio::Location::Simple');
    use_ok('Bio::Tools::GFF');
    use_ok('Bio::SeqIO');
}

my $verbose = test_debug();

#First of all we need to create an flat db
my $simple = Bio::SeqIO->new(
    -format => 'genbank',
    -file   =>  test_input_file('AB077698.gb')
);

my @features;
my $seq = $simple->next_seq();
@features = $seq->top_SeqFeatures();
is(scalar @features, 11);

ok my $col = Bio::SeqFeature::Collection->new(-verbose => $verbose);

is($col->add_features( \@features), 11);
my @feat = $col->features_in_range(
    -range => (
        Bio::Location::Simple->new(
            -start  => 100,
            -end    => 300,
            -strand => 1,
        )
    ),
    -contain => 0,
);
is(scalar @feat, 5);
if( $verbose ) {    
    for my $f ( @feat ) {
        print "location: ", $f->location->to_FTstring(), "\n";
    }
}

is(scalar $col->features_in_range(
    -range => (
        Bio::Location::Simple->new(
            -start => 100,
            -end   => 300,
            -strand => -1,
        )
    ),
    -strandmatch => 'ignore',
    -contain => 1,
), 2);

@feat = $col->features_in_range(
    -start => 79,
    -end   => 1145,
    -strand => 1,
    -strandmatch => 'strong',
    -contain => 1
);
is(scalar @feat, 5);
if( $verbose ) {    
    for my $f ( sort { $a->start <=> $b->start} @feat ) {
        print $f->primary_tag, " ", $f->location->to_FTstring(), "\n";
    }
}

is($feat[0]->primary_tag, 'CDS');
ok($feat[0]->has_tag('gene'));

$verbose = 0;
# specify input via -fh or -file
my $gffio = Bio::Tools::GFF->new(
    -file => test_input_file('myco_sites.gff'), 
    -gff_version => 2,
);
@features = ();
# loop over the input stream
while(my $feature = $gffio->next_feature()) {
    # do something with feature
    push @features, $feature;
}
$gffio->close();

is(scalar @features, 412);
$col = Bio::SeqFeature::Collection->new(
    -verbose => $verbose,
    -usefile => 1,
);

ok($col);

is($col->add_features( \@features), 412);

my $r = Bio::Location::Simple->new(
    -start => 67700,
    -end   => 150000,
    -strand => 1,
);

@feat = $col->features_in_range(
    -range => $r,
    -strandmatch => 'ignore',
    -contain => 0,
);

is(scalar @feat, 56);
is($col->feature_count, 412);
my $count = $col->feature_count;
$col->remove_features( [$features[58], $features[60]]);

is( $col->feature_count, 410);
@feat = $col->features_in_range(
    -range => $r,
    -strandmatch => 'ignore',
    -contain => 0,
);
is( scalar @feat, 54);
# add the removed features back in in order to get the collection back to size 

$col->add_features([$features[58], $features[60]]);

# let's randomize so we aren't removing and adding in the same order
# and hopefully randomly deal with a bin's expiration
fy_shuffle(\@features);

for my $f ( @features ) {
    $count--, next unless defined $f;
    $col->remove_features([$f]);
#    ok( $col->feature_count, --$count);
}
is($col->feature_count, 0);

# explicitly destroy old instances above (should clear out any open filehandles
# w/o -keep flag set)
undef $col; 

my $filename = test_output_file();
my $newcollection = Bio::SeqFeature::Collection->new(
    -verbose => $verbose,
    -keep    => 1,
    -file    => $filename,
);
$newcollection->add_features(\@feat);
is($newcollection->feature_count, 54);
undef $newcollection;
ok(-s $filename);
$newcollection = Bio::SeqFeature::Collection->new(
    -verbose => $verbose,
    -file    => $filename,
);
is($newcollection->feature_count, 54);
undef $newcollection;
ok( ! -e $filename);
# without -keep => 1, $filename was deleted as expected.
# to stop Bio::Root::Test complaining that the temp file was already deleted,
# we'll just create it again
open my $TMP, '>', $filename or die "Could not write file '$filename': $!\n";
print $TMP "temp\n";
close $TMP;

if( $verbose ) {
    my @fts =  sort { $a->start <=> $b->start}  
    grep { $r->overlaps($_,'ignore') } @features;
    
    if( $verbose ) {
        for my $f ( @fts ) {
            print $f->primary_tag, "    ", $f->location->to_FTstring(), "\n";
        }
        print "\n";
    }

    my %G = map { ($_,1) } @feat; 
    my $c = 0;
    for my $A ( @fts ) {
        if( ! $G{$A} ) {
            print "missing ", $A->primary_tag, " ", $A->location->to_FTstring(), "\n";
        } else { 
            $c++;
        }
    }
    print "Number of features correctly retrieved $c\n";
    for my $f ( sort { $a->start <=> $b->start} @feat ) {
        print $f->primary_tag, "    ", $f->location->to_FTstring(), "\n";
    }
}

sub fy_shuffle { 
    my $array = shift;
    my $i;
    for( $i = @$array; $i--; ) { 
        my $j = int rand($i+1);
        next if $i==$j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}
