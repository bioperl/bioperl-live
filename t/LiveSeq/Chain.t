# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN { 
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 45);
    
	use_ok('Bio::LiveSeq::Chain');
}

# this script tests Bio::LiveSeq::Chain methods directly since it isn't OO:
# Bio::LiveSeq::Chain->new() isn't expected here.

my $chain = Bio::LiveSeq::Chain::string2chain("abcdefghijklmnopqrstuvwxyz");
ok defined $chain;
is( Bio::LiveSeq::Chain::down_chain2string($chain), 
    "abcdefghijklmnopqrstuvwxyz");
is( Bio::LiveSeq::Chain::down_chain2string($chain,undef,4),
    "abcd"); # default start=1
my ($warning,$output);
eval {
  local $SIG{__WARN__}=sub{ $warning=$_[0]};
  $output=Bio::LiveSeq::Chain::down_chain2string($chain,1,4,6);
};
is (((index($warning,"Warning chain2string: argument LAST:6 overriding LEN:4!")==0)&&($output eq "abcdef")),1);
my $arrayref=Bio::LiveSeq::Chain::down_labels($chain,1,4);
is $arrayref->[1], 2;
$arrayref=Bio::LiveSeq::Chain::up_labels($chain,4,1);
is $arrayref->[1], 3;
$arrayref=Bio::LiveSeq::Chain::up_labels($chain);
is scalar(@{$arrayref}), 26; # total number of labels should be 26
is Bio::LiveSeq::Chain::start($chain), '1';
is Bio::LiveSeq::Chain::end($chain), '26';
ok Bio::LiveSeq::Chain::label_exists($chain,'4');
is Bio::LiveSeq::Chain::label_exists($chain,'28'), '0';
is Bio::LiveSeq::Chain::down_get_pos_of_label($chain,4), '4';
is Bio::LiveSeq::Chain::down_get_pos_of_label($chain,4,4), '1';
is Bio::LiveSeq::Chain::up_get_pos_of_label($chain,26,1), '1';
is Bio::LiveSeq::Chain::down_subchain_length($chain,1,4), '4';
is Bio::LiveSeq::Chain::up_subchain_length($chain,4,1), '4';
ok Bio::LiveSeq::Chain::invert_chain($chain);
ok Bio::LiveSeq::Chain::invert_chain($chain);
is Bio::LiveSeq::Chain::down_get_value_at_pos($chain,4), 'd';
is Bio::LiveSeq::Chain::down_get_value_at_pos($chain,1,4), 'd';
is Bio::LiveSeq::Chain::up_get_value_at_pos($chain,4), 'w';

ok Bio::LiveSeq::Chain::up_set_value_at_pos($chain,'W',4);
is Bio::LiveSeq::Chain::up_get_value_at_pos($chain,4), 'W';

ok Bio::LiveSeq::Chain::down_set_value_at_pos($chain,'D',4); 
is Bio::LiveSeq::Chain::down_get_value_at_pos($chain,4), 'D';

ok Bio::LiveSeq::Chain::set_value_at_label($chain,'d',4);
is Bio::LiveSeq::Chain::get_value_at_label($chain,4), 'd';

is Bio::LiveSeq::Chain::down_get_label_at_pos($chain,1,4), '4';
is Bio::LiveSeq::Chain::up_get_label_at_pos($chain,4), '23';
ok Bio::LiveSeq::Chain::is_downstream($chain,3,4);
is Bio::LiveSeq::Chain::is_downstream($chain,4,3), '0';
ok Bio::LiveSeq::Chain::is_upstream($chain,4,3);
is Bio::LiveSeq::Chain::is_upstream($chain,3,4), '0';
is Bio::LiveSeq::Chain::splice_chain($chain,4,2), 'de';
is Bio::LiveSeq::Chain::splice_chain($chain,7,undef,9), 'ghi';

my @array=Bio::LiveSeq::Chain::praeinsert_string($chain,"ghi",10);
is $array[0],27;
is $array[1],29;

@array=Bio::LiveSeq::Chain::postinsert_string($chain,"de",3);
is $array[0], 30;
is $array[1], 31;
is Bio::LiveSeq::Chain::up_chain2string($chain), "zyxWvutsrqponmlkjihgfedcba";

@array=Bio::LiveSeq::Chain::check_chain($chain);
is $array[0], 1;
is $array[1], 1;
is $array[2], 1;
is $array[3], 1;
