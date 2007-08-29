# -*-Perl-*- Test Harness script for Bioperl
# $Id: SeqFeatAnnotated.t,v 1.50 2007/06/27 10:16:37 sendu Exp $

use strict;

BEGIN { 
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 26,
	-requires_module => 'URI::Escape');
	
	use_ok('Bio::SeqFeature::Generic');
	use_ok('Bio::SeqFeature::Annotated');
}

SKIP: {
    my $sfa = Bio::SeqFeature::Annotated->new(-start => 1,
					      -end => 5,
					      -strand => "+",
					      -frame => 2,
					      -phase => 2,
					      -score => 12,
					      -display_name => 'test.annot',
					      -seq_id => 'test.displayname' );
    
    ok (defined $sfa);
    my $loc = $sfa->location;
    ok $loc->isa("Bio::Location::Simple");    
    ok $sfa->display_name eq 'test.annot';

	#test bsfa::from_feature
    my $sfg = Bio::SeqFeature::Generic->new ( -start => 400,
					      -end => 440,
					      -strand => 1,
					      -primary => 'nucleotide_motif',
					      -source  => 'program_a',
					      -tag => {
						  silly => 20,
						  new => 1
						  }
					      );
    my $sfa2;
    $sfa2 = Bio::SeqFeature::Annotated->new(-feature => $sfg);
    
    is $sfa2->type->name,'nucleotide_motif';
    is $sfa2->primary_tag, 'nucleotide_motif';
    is $sfa2->source->display_text,'program_a';
    is $sfa2->source_tag,'program_a';
    is $sfa2->strand,1;
    is $sfa2->start,400;
    is $sfa2->end,440;
    is $sfa2->get_Annotations('silly')->value,20;
    is $sfa2->get_Annotations('new')->value,1;
    
    my $sfa3 = Bio::SeqFeature::Annotated->new( -start => 1,
						-end => 5,
						-strand => "+",
						-frame => 2,
						-phase => 2,
						-score => 12,
						-display_name => 'test.annot',
						-seq_id => 'test.displayname' );
    $sfa3->from_feature($sfg);
    
    is $sfa3->type->name,'nucleotide_motif', 'type->name';
    is $sfa3->primary_tag, 'nucleotide_motif', 'primary_tag';
    is $sfa3->source->display_text,'program_a';
    is $sfa3->source_tag,'program_a';
    is $sfa3->strand,1;
    is $sfa3->start,400;
    is $sfa3->end,440;
    is $sfa3->get_Annotations('silly')->value,20;
    is $sfa3->get_Annotations('new')->value,1;
    # Note there is an API conflict with SeqFeature::Generic, where score is a
    # simple scalar, and here it is a Bio::Annotation::SimpleValue
    # By popular vote there is no operator overloading, so this needs to be
    # resolved
    is $sfa3->score(), 12; 
    $sfa3->score(11);
    is $sfa3->score(), 11; 
    $sfa3->score(0);
    is $sfa3->score(), 0;	# test that setting to 0 no longer is overriddent to set score to '.' (fixed in Bio::SeqFeature::Annotated version 1.3.7)
}