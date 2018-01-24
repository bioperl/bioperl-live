use strict;
use warnings;
BEGIN {
  use Bio::Root::Test;
  test_begin(
    -tests => 4,
   );
  use_ok 'Bio::DB::HIV';
  use_ok 'Bio::DB::Query::HIVQuery';
}

SKIP: {
  test_skip(-tests => 2,
	    -requires_networking => 1);
  my $db = new Bio::DB::HIV;
  my $seq = $db->get_Seq_by_id('94284');                                 # LANL sequence id
  $seq = $db->get_Seq_by_acc('EF432710');                             # GenBank accession
  my $q = new Bio::DB::Query::HIVQuery( " (C D)[subtype] SI[phenotype] (symptomatic AIDS)[health_status] " );
  my $seqio = $db->get_Stream_by_query($q);
  $seq = $seqio->next_seq();
  like $seq->annotation->get_value('Virus','subtype'), qr/C|D/, "returns 'C or D'";
  is $seq->annotation->get_value('Patient','patient_health'), 'symptomatic', "returns 'symptomatic'";

}
