#!/usr/bin/perl

use lib ".";
use strict;

require Bio::Root::Version;

foreach ( @ARGV ) {
    `perl -pi -e 's/\\\@\\\@VERSION\\\@\\\@/$Bio::Root::Version::VERSION/g;' "$_"`;
}

__END__
