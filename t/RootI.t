# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

use strict;

BEGIN {
    eval {require Test::More;};
	if ($@) {
		use lib 't/lib';
	}
	use Test::More;
	
    plan tests => 23;
	
	use_ok('Bio::Root::Root');
}

ok my $obj = new Bio::Root::Root();
isa_ok($obj, 'Bio::Root::RootI');

eval { $obj->throw('Testing throw') };
ok $@ =~ /Testing throw/;# 'throw failed';

# doesn't work in perl 5.00405
#my $val;
#eval {
#    my ($tfh,$tfile) = $obj->tempfile();
#    local * STDERR = $tfh;
#    $obj->warn('Testing warn');
#    close $tfh;    
#    open(IN, $tfile) or die("cannot open $tfile");    
#    $val = join("", <IN>) ;
#    close IN;
#    unlink $tfile;
#};
#ok $val =~ /Testing warn/;
#'verbose(0) warn did not work properly' . $val;

$obj->verbose(-1);
eval { $obj->throw('Testing throw') };
ok $@=~ /Testing throw/;# 'verbose(-1) throw did not work properly' . $@;

eval { $obj->warn('Testing warn') };
ok !$@;

$obj->verbose(1);
eval { $obj->throw('Testing throw') };
ok $@ =~ /Testing throw/;# 'verbose(1) throw did not work properly' . $@;

# doesn't work in perl 5.00405
#undef $val;
#eval {
#    my ($tfh,$tfile) = $obj->tempfile();
#    local * STDERR = $tfh;
#    $obj->warn('Testing warn');
#    close $tfh;
#    open(IN, $tfile) or die("cannot open $tfile");    
#    $val = join("", <IN>);
#    close IN;
#    unlink $tfile;
#};
#ok $val =~ /Testing warn/;# 'verbose(1) warn did not work properly' . $val;

my @stack = $obj->stack_trace();
is scalar @stack, 2;

my $verbobj = new Bio::Root::Root(-verbose=>1,-strict=>1);
is $verbobj->verbose(), 1;

$Bio::Root::Root::DEBUG = 1;
require Bio::Seq;
my $seq = new Bio::Seq;
is $seq->verbose, 1;

# test for bug #1343
my @vals = Bio::Root::RootI->_rearrange([qw(apples pears)], 
					-apples => 'up the',
					-pears  => 'stairs');
eval { $obj->throw_not_implemented() };
ok $@ =~ /Bio::Root::NotImplemented/;

is shift @vals, 'up the';
is shift @vals, 'stairs';

# tests for _set_from_args
{
	no warnings 'redefine';
	
	# simplest form
	local *Bio::Root::Root::new = sub {
		my $class = shift;
		my $self = {};
		bless $self, ref($class) || $class;
	
		$self->_set_from_args(\@_);
		
		return $self;
	};
	
	$obj = Bio::Root::Root->new(-verbose => 1, t1 => 1, '--Test-2' => 2);
	ok ! $obj->can('t1'), 'arg not callable';
	
	# with method creation
	local *Bio::Root::Root::new = sub {
		my $class = shift;
		my $self = {};
		bless $self, ref($class) || $class;
	
		$self->_set_from_args(\@_, -create => 1);
		
		return $self;
	};
	
	$obj = Bio::Root::Root->new(-verbose => 1, t3 => 1, '--Test-4' => 2);
	ok $obj->can('t3'), 'arg callable since method was created';
	ok $obj->can('test_4'), 'mal-formed arg callable since method was created with good name';
	
	# with method creation, but limited to certain methods
	local *Bio::Root::Root::new = sub {
		my $class = shift;
		my $self = {};
		bless $self, ref($class) || $class;
	
		$self->_set_from_args(\@_, -methods => ['verbose', 't5'], -create => 1);
		
		return $self;
	};
	
	$obj = Bio::Root::Root->new(-verbose => 1, t5 => 1, '--Test-6' => 2);
	ok $obj->can('t5'), 'arg callable since method was created';
	ok ! $obj->can('test_6'), 'arg not in method list not created';
	
	# with synonyms
	local *Bio::Root::Root::new = sub {
		my $class = shift;
		my $self = {};
		bless $self, ref($class) || $class;
		
		my %args = @_;
		
		$self->_set_from_args(\%args, -methods => {(verbose => 'v',
									         		test7 => 't7',
													test_8 => 't8')},
								      -create => 1);
		
		return $self;
	};
	
	$obj = Bio::Root::Root->new(-verbose => 1, t7 => 1, '--Test-8' => 2);
	is $obj->verbose, 1, 'verbose was set correctly';
	is $obj->t7, 1, 'synonym was set correctly';
	is $obj->test7, 1, 'real method of synonym was set correctly';
	is $obj->test_8, 2, 'mal-formed arg correctly resolved to created method';
	is $obj->t8, 2, 'synonym of set method was set correctly';
}

1;

