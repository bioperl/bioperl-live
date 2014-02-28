# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 27);

    use_ok('Bio::Tools::Run::WrapperBase');
}

my @params = qw(test1 test_2);
my @switches = qw(Test3 test_4);
*Bio::Tools::Run::WrapperBase::new = sub {
    my ($class, @args) = @_;
    my $self = $class->Bio::Tools::Run::WrapperBase::SUPER::new(@args);
    
    $self->_set_from_args(\@args, -methods => [@params, @switches],
                                  -create => 1);
    
    return $self;
};
my $new = *Bio::Tools::Run::WrapperBase::new; # just to avoid warning
my $obj = Bio::Tools::Run::WrapperBase->new(-test_2 => 2, -test3 => 1, -test_4 => 0);
isa_ok($obj, 'Bio::Tools::Run::WrapperBase');

# it is interface-like with throw_not_implemented methods; check their
# existance
foreach my $method (qw(run program_dir program_name version)) {
    ok $obj->can($method), "$method() exists";
}

## most methods are defined; check their function

# simple get/setters
foreach my $method (qw(error_string arguments no_param_checks save_tempfiles
                       outfile_name quiet)) {
    $obj->$method(1);
    is $obj->$method(), 1, "$method could be set";
}

# tempdir

$obj->save_tempfiles(0);
my $tmpdir = $obj->tempdir();
ok -d $tmpdir, 'tempdir created a directory';
ok open(my $test, '>', File::Spec->catfile($tmpdir, 'test')), 'could create file in tempdir';
print $test "test\n";
close $test;

# cleanup

$obj->cleanup();
ok ! -d $tmpdir, 'following cleanup() with save_tempfiles unset, tempdir was deleted';

# io
my $io1 = $obj->io;
my $io2 = $obj->io;
isa_ok($io1, 'Bio::Root::IO');
is $io1, $io2, 'io() always returns the same instance of IO';

# program_dir and program_name need to be defined for program_path and
# executable to work
{
    no warnings 'redefine';
    *Bio::Tools::Run::WrapperBase::program_dir = sub {
        my $self = shift;
        if (@_) { $self->{pdir} = shift }
        return $self->{pdir} || '';
    };
    *Bio::Tools::Run::WrapperBase::program_name = sub {
        my $self = shift;
        if (@_) { $self->{pname} = shift }
        return $self->{pname} || '';
    };
}
$obj->program_dir('test_dir');
$obj->program_name('test_name');

# program_path
is $obj->program_path, File::Spec->catfile('test_dir', 'test_name'.($^O =~ /mswin/i ?'.exe':'')), 'program_path was correct';

# executable
throws_ok { $obj->executable } qr/Cannot find executable/, 'pretend program name not found as executable';
$obj->program_name('perl');
ok $obj->executable, 'perl found as executable';

# _setparams
my $params = $obj->_setparams(-params => \@params,
                              -switches => \@switches);
is $params, ' test_2 2 Test3', 'params string correct';
$params = $obj->_setparams(-params => \@params,
                           -switches => \@switches,
                           -join => '=');
is $params, ' test_2=2 Test3', 'params string correct';
$params = $obj->_setparams(-params => \@params,
                           -switches => \@switches,
                           -join => '=',
                           -lc => 1);
is $params, ' test_2=2 test3', 'params string correct';
$params = $obj->_setparams(-params => \@params,
                           -switches => \@switches,
                           -join => '=',
                           -lc => 1,
                           -dash => 1);
is $params, ' -test_2=2 -test3', 'params string correct';
$params = $obj->_setparams(-params => \@params,
                           -switches => \@switches,
                           -join => '=',
                           -lc => 1,
                           -double_dash => 1);
is $params, ' --test_2=2 --test3', 'params string correct';
$params = $obj->_setparams(-params => \@params,
                           -switches => \@switches,
                           -join => '=',
                           -lc => 1,
                           -double_dash => 1,
                           -underscore_to_dash => 1);
is $params, ' --test-2=2 --test3', 'params string correct';
$params = $obj->_setparams(-params => {(test1 => 't1', test_2 => 't_2')},
                           -switches => {(Test3 => 'T3', test_4 => 't4')},
                           -join => '=',
                           -lc => 1,
                           -double_dash => 1,
                           -underscore_to_dash => 1);
is $params, ' --t-2=2 --t3', 'params string correct';
