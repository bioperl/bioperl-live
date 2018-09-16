#!/usr/bin/perl
# Author: Damien Mattei C.N.R.S / U.N.S.A - UMR 6549
# example: ./idfetch.pl AP001266

use Bio::DB::GenBank;

$gb = new Bio::DB::GenBank();

# this returns a Seq object :
$seq1 = $gb->get_Seq_by_acc($ARGV[0]);
print $seq1->display_id() . "\n" ;


foreach $feat ($seq1->all_SeqFeatures()) {

  #print $feat->primary_tag . " " . $feat->source_tag() . "\n" ;

  print "Feature from ", $feat->start, " to ",
   $feat->end, " Primary tag  ", $feat->primary_tag,
   ", produced by ", $feat->source_tag(), "\n";

  if( $feat->strand == 0 ) {
    print "Feature applicable to either strand\n";
  } else {
    print "Feature on strand ", $feat->strand,"\n"; # -1,1
  }

  foreach $tag ( $feat->all_tags() ) {
    print "Feature has tag ", $tag, " with values, ",
    join(' ',$feat->each_tag_value($tag)), "\n";
  }

  print "new feature\n" if $feat->has_tag('new');

}


exit;

__END__

It will display something like that:

[dmattei@pclgmch2 gmap]$ ./idfetch.pl AP001266
AP001266
Feature from 1 to 168978 Primary tag  source, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag chromosome with values, 11
Feature has tag map with values, 11q13
Feature has tag clone with values, RP11-770G2
Feature has tag organism with values, Homo sapiens
Feature has tag db_xref with values, taxon:9606
Feature from 1 to 31550 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 31651 to 48510 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 48611 to 64044 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 64145 to 78208 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 78309 to 89008 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 89109 to 99704 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 99805 to 107965 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 108066 to 116032 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 116133 to 124010 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 124111 to 130494 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 130595 to 136072 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 136173 to 139649 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 139750 to 144590 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 144691 to 148482 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 148583 to 152279 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 152380 to 153632 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment clone_end:T7 
vector_side:left
Feature from 153733 to 155746 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 155847 to 156405 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment clone_end:SP6 
vector_side:right
Feature from 156506 to 158398 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 158499 to 161333 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 161434 to 163304 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 163405 to 164604 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 164705 to 166693 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment
Feature from 166794 to 168978 Primary tag  misc_feature, produced by 
EMBL/GenBank/SwissProt
Feature on strand 1
Feature has tag note with values, assembly_fragment

