use strict;

BEGIN {
   use lib '.';
   use Bio::Root::Test;
   
   test_begin(-tests => 14,
              -requires_modules => [ qw(XML::SAX) ]
              );
   
   use_ok('Bio::SeqIO::gbxml');
}

my $verbose = test_debug();

my $in = Bio::SeqIO->new(-format  => 'gbxml',
                          -verbose => $verbose,
                          # This class can't parse dbEST data yet...
                          # -file    => test_input_file('roa1.gbxml'));
                          # So let's try a <GBSeq> file:
                          -file    => test_input_file('EG352462.gbxml'));
isa_ok($in, 'Bio::SeqIO');
my $seq = $in->next_seq();
is($seq->molecule,   'mRNA',                                            'molecule');
is($seq->alphabet,   'dna',                                             'alphabet');
is($seq->primary_id,  116038450,                                        'primary_id');
is($seq->display_id, 'EG352462',                                        'display_id');
is($seq->version,     1,                                                'version');
is($seq->is_circular, 0,                                                'is_circular');

is(substr($seq->description, 0, 10), 'SAAH-aad23',                      'description');
is(substr($seq->seq,         0, 10), 'aataaaatta',                      'sequence');

my @class = $seq->species->classification;
is($class[$#class],'Eukaryota',                                         'classification');

my ($feat) = $seq->get_SeqFeatures;
is_deeply([ $feat->get_tag_values('clone_lib') ], [ 'Agen 0058' ],      'feat - clone_lib');
is_deeply([ $feat->get_tag_values('db_xref') ],   [ 'taxon:79327' ],    'feat - db_xref');
is_deeply([ $feat->get_tag_values('lab_host') ],  [ 'DH10B cells' ],    'feat - lab_host');



