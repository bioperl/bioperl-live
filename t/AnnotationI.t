#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 4;

BEGIN {
    use_ok('Bio::AnnotationI');
}

{
    package MyAnnotation;
    use base qw(Bio::AnnotationI);

    sub new {
        my $class = shift;
        return bless {}, $class;
    }

    sub tagname  { return "mock_tag"; }
    sub as_text  { return "Mock annotation as text"; }
    sub hash_tree { return {}; }
}

my $a = MyAnnotation->new();

isa_ok($a, 'Bio::AnnotationI');
is($a->tagname, 'mock_tag', 'tagname returns expected value');
is($a->as_text, 'Mock annotation as text', 'as_text returns expected value');

