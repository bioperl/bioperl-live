# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'
use Test;
use strict;
BEGIN { plan tests => 1 }
use Bio::Tools::Blast::Sbjct;

ok(1);
