use strict;
use vars qw($DEBUG $TESTCOUNT);
my $error;

BEGIN {     
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    $TESTCOUNT = 10;
    plan tests => $TESTCOUNT;

    eval { require IO::String; require Bio::DB::GenBank; };
    if( $@ ) {
	print STDERR "Bio::DB::GenBank unable to be loaded. This means the Bio::DB::* modules are not usable. so can't test remote locations.\n";
	for( 1..$TESTCOUNT ) {
	    skip("Not possible to test without DB access",1);
	}
       $error = 1; 
    }
};

if( $error ==  1 ) {
    exit(0);
}


use Bio::Seq;
use Bio::SeqIO;


my $str = Bio::SeqIO->new( '-file'=> Bio::Root::IO->catfile("t","data","test.genbank"), 
			'-format' => 'GenBank');

ok $str;
my $seq;

ok ( $seq = $str->next_seq() );

foreach my $ft ( $seq->top_SeqFeatures ) {
	my $t = $ft->spliced_seq();
}

ok(1);

my $db = Bio::DB::GenBank->new();

foreach my $ft ( $seq->top_SeqFeatures ) {
	my $t = $ft->spliced_seq();
	print "Got ",$t->seq,"\n";
}




