# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 22);
	
    use_ok('Bio::Tools::Lucy');
}

my @params = (adv_stderr => 1, seqfile => test_input_file('lucy.seq'), rev_desig => 'R'); 
# Bio::Tools::Lucy will find .qual, .info, and .stderr files in this folder 

my $lucyObj = Bio::Tools::Lucy->new(@params);
isa_ok $lucyObj,'Bio::Tools::Lucy';
ok $lucyObj->seqfile();
$lucyObj->adv_stderr(1);
my $stderr = $lucyObj->adv_stderr();
ok $stderr;
my $names =$lucyObj->get_sequence_names();
ok $names;
my $seq = shift @$names;
is $seq, 'TaLr1010B10R';
is $lucyObj->length_raw("$seq"), 1060;
is $lucyObj->length_clear("$seq"), 420;
is $lucyObj->start_clear("$seq"), 86;
is $lucyObj->end_clear("$seq"), 505;
ok $lucyObj->avg_quality("$seq");
ok $lucyObj->full_length("$seq");
ok $lucyObj->polyA("$seq");
is $lucyObj->direction("$seq"), 'R';
ok $lucyObj->per_GC("$seq");
ok $lucyObj->sequence("$seq");
ok $lucyObj->quality("$seq");
my $seqObj = $lucyObj->get_Seq_Obj("$seq");
ok $seqObj;
my $seqObjs = $lucyObj->get_Seq_Objs();
ok $seqObjs;

my $rejects = $lucyObj->get_rejects();
ok $rejects;
my ($key) = (sort keys %$rejects);
is $key, 'TaLr1011A07R';
is $rejects->{$key}, 'Q'; 
