# -*-Perl-*- Test Harness script for Bioperl
# $Id: epost.t 15112 2008-12-08 18:12:38Z sendu $

use strict;
use warnings;
use Data::Dumper;
BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 83,
			   -requires_module => 'XML::Simple');
	
    use_ok('Bio::Tools::EUtilities');
    use_ok('Bio::Tools::EUtilities::EUtilParameters');
}

my $eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'esummary',
    -file       => test_input_file('eutils','esummary1.xml'));

isa_ok($eutil, 'Bio::Tools::EUtilities::Summary');

# note that XML output does not contain the database; in order to retrieve this
# (and similar missing params) one should pass in the EUtilParameters object

is(join(',',$eutil->get_databases), '');

# we'll add in some parameters (normally passed in via Bio::DB::EUtilities)
my $p = Bio::Tools::EUtilities::EUtilParameters->new(
                        -eutil => 'esummary',
                        -db => 'protein',
                        -id => [1621261,89318838,68536103,20807972,730439]);

isa_ok($eutil, 'Bio::Tools::EUtilities::Summary');

$eutil->parameter_base($p);
is(join(',',$eutil->get_databases), 'protein');

# for esummary, DocSums contain the IDs, but we glob them together when called from the parser
is(join(',',$eutil->get_ids), '1621261,89318838,68536103,20807972,730439', 'get_ids');

my @ds = $eutil->get_DocSums;
is(scalar(@ds), 5);

isa_ok($ds[0], 'Bio::Tools::EUtilities::Summary::DocSum');

# One ID per DocSum (get_ids is implemented JIC)
is($ds[0]->get_id, '1621261');
is(join(',',$ds[0]->get_ids), '1621261');

# test two DocSums: get Items
my @items = $ds[0]->get_all_Items;
is(scalar(@items), 12);

isa_ok($items[0], 'Bio::Tools::EUtilities::Summary::Item');
# each Item has four possible pieces of data: ID, name, type, and content
# also, an Item may have sub-Items (up to 3 hierarchal layers: Item, ListItem, StructureItem)
is($items[0]->get_id,1621261);
is($items[0]->get_name,'Caption');
is($items[0]->get_type,'String');
is($items[0]->get_content,'CAB02640');
is(scalar($items[0]->get_ListItems), 0);

is($items[3]->get_id,1621261);
is($items[3]->get_name,'Gi');
is($items[3]->get_type,'Integer');
is($items[3]->get_content,1621261);
is(scalar($items[3]->get_ListItems), 0);

is($items[7]->get_id,1621261);
is($items[7]->get_name,'TaxId');
is($items[7]->get_type,'Integer');
is($items[7]->get_content,83332);
is(scalar($items[7]->get_ListItems), 0);

@items = $ds[2]->get_all_Items;
is(scalar(@items), 12);

isa_ok($items[0], 'Bio::Tools::EUtilities::Summary::Item');
# each Item has four possible pieces of data: ID, name, type, and content
# also, an Item may have sub-Items (up to 3 hierarchal layers: Item, ListItem, StructureItem)
is($items[0]->get_id,68536103);
is($items[0]->get_name,'Caption');
is($items[0]->get_type,'String');
is($items[0]->get_content,'YP_250808');
is(scalar($items[0]->get_ListItems), 0);

is($items[3]->get_id,68536103);
is($items[3]->get_name,'Gi');
is($items[3]->get_type,'Integer');
is($items[3]->get_content,68536103);
is(scalar($items[3]->get_ListItems), 0);

is($items[7]->get_id,68536103);
is($items[7]->get_name,'TaxId');
is($items[7]->get_type,'Integer');
is($items[7]->get_content,306537);
is(scalar($items[7]->get_ListItems), 0);

# getting data directly from DocSum

is($ds[0]->get_type_by_name('Gi'), 'Integer');
is(join(',',$ds[0]->get_contents_by_name('CreateDate')), '2003/11/21');

is($ds[1]->get_type_by_name('Status'), 'String');
is(join(',',$ds[1]->get_contents_by_name('Extra')), 'gi|89318838|gb|EAS10332.1|[89318838]');

is($ds[3]->get_type_by_name('TaxId'), 'Integer');
is(join(',',$ds[3]->get_contents_by_name('Title')), 'pyrimidine regulatory protein PyrR [Thermoanaerobacter tengcongensis MB4]');


$eutil = Bio::Tools::EUtilities->new(
    -eutil      => 'esummary',
    -file       => test_input_file('eutils','esummary2.xml'));

isa_ok($eutil, 'Bio::Tools::EUtilities::Summary');

# note that XML output does not contain the database; in order to retrieve this
# (and similar missing params) one should pass in the EUtilParameters object

is(join(',',$eutil->get_databases), '');

# we'll add in some parameters (normally passed in via Bio::DB::EUtilities)
$p = Bio::Tools::EUtilities::EUtilParameters->new(
                        -eutil => 'esummary',
                        -db => 'homologene');

isa_ok($eutil, 'Bio::Tools::EUtilities::Summary');

$eutil->parameter_base($p);
is(join(',',$eutil->get_databases), 'homologene');

# for esummary, DocSums contain the IDs, but we glob them together when called from the parser
is(join(',',$eutil->get_ids), '32049,45614', 'get_ids');

@ds = $eutil->get_DocSums;
is(scalar(@ds), 2);

isa_ok($ds[0], 'Bio::Tools::EUtilities::Summary::DocSum');

# One ID per DocSum (get_ids is implemented JIC)
is($ds[0]->get_id, '32049');
is(join(',',$ds[0]->get_ids), '32049');

# flattened list
@items = $ds[0]->get_all_Items;
is(scalar(@items), 62);

# Items are layered when caling get_Items; flattened list is not the same as
# normal list
@items = $ds[0]->get_Items;
is(scalar(@items), 2);

isa_ok($items[0], 'Bio::Tools::EUtilities::Summary::Item');
# each Item has four possible pieces of data: ID, name, type, and content
# also, an Item may have sub-Items (up to 3 hierarchal layers: Item, ListItem, StructureItem)
is($items[0]->get_id,32049);
is($items[0]->get_name,'HomoloGeneDataList');
is($items[0]->get_type,'List');
is($items[0]->get_content, undef); # List contents are other Items

# access List layer from top Item
my @li = $items[0]->get_ListItems;
is(scalar(@li), 10);
@li = $items[0]->get_Items;
is(scalar(@li), 10);

# access Structure Layer from List Item
my @si = $li[1]->get_StructureItems;
is(scalar(@si), 5);
@si = $li[1]->get_StructureItems;
is(scalar(@si), 5);

# test List Item
is($li[1]->get_id,32049);
is($li[1]->get_name,'HomoloGeneData');
is($li[1]->get_type,'Structure');
is($li[1]->get_content,undef); # Structure contents are other Items

# test Structure Item
is($si[2]->get_id,32049);
is($si[2]->get_name,'Symbol');
is($si[2]->get_type,'String');
is($si[2]->get_content,'NOTCH1');

# getting data directly from DocSum

is($ds[0]->get_type_by_name('HomoloGeneData'), 'Structure');
is(join(',',$ds[0]->get_contents_by_name('Symbol')), 'NOTCH1,NOTCH1,NOTCH1,NOTCH1,Notch1,Notch1,NOTCH1,notch1b,N,AgaP_AGAP001015');

is($ds[1]->get_type_by_name('HomoloGeneDataList'), 'List');
is(join(',',$ds[1]->get_contents_by_name('TaxId')), '9606,9913,10090,10116,9031,7955');

