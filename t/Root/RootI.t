# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 62);

    use_ok 'Bio::Root::Root';
}

ok my $obj = Bio::Root::Root->new();
isa_ok $obj, 'Bio::Root::RootI';

throws_ok { $obj->throw('Testing throw') } qr/Testing throw/;# 'throw failed';

# test throw_not_implemented()
throws_ok { $obj->throw_not_implemented() } qr/EXCEPTION: Bio::Root::NotImplemented/;

{
    package Bio::FooI;
    use base qw(Bio::Root::RootI);
    sub new {
            my $class = shift;
            my $self = {};
            bless $self, ref($class) || $class;
            return $self;
    };
}
$obj = Bio::FooI->new();
throws_ok { $obj->throw_not_implemented() } qr/EXCEPTION /;
$obj = Bio::Root::Root->new();

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
throws_ok { $obj->throw('Testing throw') } qr/Testing throw/;# 'verbose(-1) throw did not work properly' . $@;

lives_ok { $obj->warn('Testing warn') };

$obj->verbose(1);
throws_ok { $obj->throw('Testing throw') } qr/Testing throw/;# 'verbose(1) throw did not work properly' . $@;

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

my $verbobj = Bio::Root::Root->new(-verbose=>1,-strict=>1);
is $verbobj->verbose(), 1;

$Bio::Root::Root::DEBUG = 1;
my $seq = Bio::Root::Root->new();
is $seq->verbose, 1;

# test for bug #1343
my @vals = Bio::Root::RootI->_rearrange([qw(apples pears)],
                                        -apples => 'up the',
                                        -pears  => 'stairs');
is shift @vals, 'up the';
is shift @vals, 'stairs';

# test deprecated()

# class method
{
    local $Bio::Root::Root::VERSION = 8.9;
    warning_like{ Bio::Root::Root->deprecated('Test1') } qr/Test1/, 'simple';
    warning_like{ Bio::Root::Root->deprecated(-message => 'Test2') } qr/Test2/;
    warning_like{ Bio::Root::Root->deprecated('Test3', 999.999) } qr/Test3/,
        'warns for versions below current version';
    warning_like{ Bio::Root::Root->deprecated(-message => 'Test4',
                                              -version => 999.999) } qr/Test4/,
        'warns for versions below current version';
    throws_ok{ Bio::Root::Root->deprecated('Test5', 0.001) } qr/Test5/,
        'throws for versions above current version';
    throws_ok{ Bio::Root::Root->deprecated(-message => 'Test6',
                                       -version => 0.001) } qr/Test6/,
        'throws for versions above current version';

    throws_ok{ Bio::Root::Root->deprecated(-message => 'Test6',
                                           -version => $Bio::Root::Root::VERSION) } qr/Test6/,
        'throws for versions equal to current version';

    # object method
    my $root = Bio::Root::Root->new();
    warning_like{ $root->deprecated('Test1') } qr/Test1/, 'simple';
    warning_like{ $root->deprecated(-message => 'Test2') } qr/Test2/, 'simple';
    warning_like{ $root->deprecated('Test3', 999.999) } qr/Test3/,
        'warns for versions below current version';
    warning_like{ $root->deprecated(-message => 'Test4',
                                    -version => 999.999) } qr/Test4/,
                                'warns for versions below current version';
    throws_ok{ $root->deprecated('Test5', 0.001) } qr/Test5/,
      'throws for versions above current version';
    throws_ok{ $root->deprecated(-message => 'Test6',
                                 -version => 0.001) } qr/Test6/,
                             'throws for versions above current version';

}

# tests for _set_from_args()
# Let's not pollute Bio::Root::Root namespace if possible
# Create temp classes instead which inherit Bio::Root::Root, then test

{

    package Bio::Foo1;
    use base qw(Bio::Root::Root);
    our $VERSION = '2.00';
    sub new {
        my $class = shift;
        my $self = {};
        bless $self, ref($class) || $class;

        $self->_set_from_args(\@_);

        return $self;
    };
}

$obj = Bio::Foo1->new(-verbose => 1, t1 => 1, '--Test-2' => 2);
#ok ! $obj->can('t1'), 'arg not callable';

{

    package Bio::Foo2;
    use base qw(Bio::Root::Root);
    sub new {
        my $class = shift;
        my $self = {};
        bless $self, ref($class) || $class;

        $self->_set_from_args(\@_, -create => 1);

        return $self;
    };

}

$obj = Bio::Foo2->new(-verbose => 1, t3 => 1, '--Test-4' => 2);
ok $obj->can('t3'), 'arg callable since method was created';
ok $obj->can('test_4'), 'mal-formed arg callable since method was created with good name';
for my $m (qw(t3 test_4)) {
    can_ok('Bio::Foo2',$m);
    ok ! Bio::Root::Root->can($m), "Methods don't pollute original Bio::Root::Root namespace";
}

{
    package Bio::Foo3;
    use base qw(Bio::Root::Root);
    sub new {
        my $class = shift;
        my $self = {};
        bless $self, ref($class) || $class;

        $self->_set_from_args(\@_, -methods => ['verbose', 't5'], -create => 1);

        return $self;
    };
}

$obj = Bio::Foo3->new(-verbose => 1, t5 => 1, '--Test-6' => 2);
can_ok $obj, 't5';
ok ! $obj->can('test_6'), 'arg not in method list not created';

can_ok ('Bio::Foo3','t5');
ok ! UNIVERSAL::can('Bio::Root::Root','t5'), "Methods don't pollute original Bio::Root::Root namespace";

{
    package Bio::Foo4;
    use base qw(Bio::Root::Root);
    sub new {
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
}

# with synonyms

$obj = Bio::Foo4->new(-verbose => 1, t7 => 1, '--Test-8' => 2);
is $obj->verbose, 1, 'verbose was set correctly';
is $obj->t7, 1, 'synonym was set correctly';
is $obj->test7, 1, 'real method of synonym was set correctly';
is $obj->test_8, 2, 'mal-formed arg correctly resolved to created method';
is $obj->t8, 2, 'synonym of set method was set correctly';

for my $m (qw(t7 test7 test_8 t8)) {
    can_ok 'Bio::Foo4', $m;
    ok ! UNIVERSAL::can('Bio::Root::Root','t7'), "Methods don't pollute original Bio::Root::Root namespace";
}

# test basic Root::clone()

my $clone = $obj->clone;

is $clone->t7, $obj->t7, 'clone';
is $clone->test7, $obj->test7, 'clone';
is $clone->test_8, $obj->test_8, 'clone';
$clone->test_8('xyz');
isnt $clone->test_8, $obj->test_8, 'clone changed, original didn\'t';

# test Root::clone() with parameter passing, only works for methods
# (introspection via can())

my $clone2 = $obj->clone(-t7 => 'foo');

is $clone2->t7, 'foo', 'parameters passed to clone() modify object';
is $obj->t7, 1, 'original is not modified';



# test deprecations using start_version
{
    package Bio::Foo5;
    use base qw(Bio::Root::Root);

    our $v = '18.001';
    our $VERSION = $v;

    sub not_good {
        my $self = shift;
        $self->deprecated(-message => 'This is not good',
                          -warn_version  => $v,
                          -throw_version => $v + 0.001);
    }

    sub not_good2 {
        my $self = shift;
        # note, due to _rearrange, ordering is throw version, then warn version
        $self->deprecated('This is not good',$v + 0.001,$v);
    }

    sub really_not_good {
        my $self = shift;
        $self->deprecated(-message => 'This is really not good',
                          -warn_version  => $v - 0.001,
                          -throw_version => $v,);
    }

    # version is the same as throw_version (and vice versa)
    sub still_very_bad {
        my $self = shift;
        $self->deprecated(-message => 'This is still very bad',
                          -warn_version  => $v - 0.001,
                          -version => $v);
    }

    sub okay_for_now {
        my $self = shift;
        $self->deprecated(-message => 'This is okay for now',
                          -warn_version  => $v + 0.001,
                          -throw_version => $v + 0.002);
    }

    sub plain_incorrect {
        my $self = shift;
        $self->deprecated(-message => 'This is not going to work',
                          -warn_version  => '1.2.3.4',
                          -throw_version => 'a.b.c.d');
    }
}

my $foo = Bio::Foo5->new();

throws_ok { $foo->plain_incorrect } qr/Version must be numerical/,
    'must use proper versioning scheme';

warning_like{ $foo->not_good } qr/This is not good/,
    'warns for versions >= current version';
# this tests the three-arg (non-named) form just to make sure it works, even
# though we probably won't support it
warning_like{ $foo->not_good2 } qr/This is not good/,
    'warns for versions >= current version';

throws_ok { $foo->really_not_good } qr/This is really not good/,
    'throws for versions >= current version';
throws_ok { $foo->still_very_bad } qr/This is still very bad/,
    'throws for versions >= current version';
lives_ok { $foo->okay_for_now } 'No warnings/exceptions below current version';


