# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

use strict;
BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;

    plan tests => 22;
}


use Bio::Tools::Lucy;
use Bio::Root::IO;

ok(1);

my @params = (adv_stderr => 1, seqfile => Bio::Root::IO->catfile("t","data","lucy.seq"), rev_desig => 'R'); 
# Bio::Tools::Lucy will find .qual, .info, and .stderr files in this folder 

my $lucyObj = Bio::Tools::Lucy->new(@params);
ok $lucyObj->isa('Bio::Tools::Lucy');
ok $lucyObj->seqfile();
$lucyObj->adv_stderr(1);
my $stderr = $lucyObj->adv_stderr();
ok $stderr;
my $names =$lucyObj->get_sequence_names();
ok $names;
my $seq = shift @$names;
ok $seq, 'TaLr1010B10R';
ok $lucyObj->length_raw("$seq"), 1060;
ok $lucyObj->length_clear("$seq"), 420;
ok $lucyObj->start_clear("$seq"), 86;
ok $lucyObj->end_clear("$seq"), 505;
ok $lucyObj->avg_quality("$seq");
ok $lucyObj->full_length("$seq");
ok $lucyObj->polyA("$seq");
ok $lucyObj->direction("$seq"), 'R';
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
ok $key, 'TaLr1011A07R';
ok $rejects->{$key}, 'Q'; 
