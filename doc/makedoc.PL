#!/usr/bin/perl

use lib ".";
use strict;

require Bio::Root::Version;
require File::Find;

File::Find::find(\&wanted, '.');

sub wanted { 
    next unless ( /^.*\.pod\z/os ); # add .pm to be able to use @@VERSION@@
                                  # within module docs as well
    `perl -pi -e 's/==VERSION==/$Bio::Root::Version::VERSION/g;' $_`;
}

__END__
