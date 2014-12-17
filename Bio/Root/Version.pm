package Bio::Root::Version;
use strict;

# ABSTRACT: provide global, distribution-level versioning
# AUTHOR:   Aaron Mackey <amackey@virginia.edu>
# OWNER:    Aaron Mackey
# LICENSE:  Perl_5

=head1 SYNOPSIS

  package Bio::Tools::NiftyFeature;
  require Bio::Root::RootI;


  # later, in client code:
  package main;
  use Bio::Tools::NiftyFeature 3.14;


  ## alternative usage: NiftyFeature defines own $VERSION:
  package Bio::Tools::NiftyFeature;
  my $VERSION = 9.8;

  # later in client code:
  package main;

  # ensure we're using an up-to-date BioPerl distribution
  use Bio::Perl 3.14;

  # NiftyFeature has its own versioning scheme:
  use Bio::Tools::NiftyFeature 9.8;

=head1 DESCRIPTION

This module provides a mechanism by which all other BioPerl modules
can share the same $VERSION, without manually synchronizing each file.

Bio::Root::RootI itself uses this module, so any module that directly
(or indirectly) uses Bio::Root::RootI will get a global $VERSION
variable set if it's not already.

=cut

our $VERSION = '1.006924'; # pre-1.7
$VERSION = eval $VERSION;

sub import {
    # try to handle multiple levels of inheritance:
    my $i = 0;
    my $pkg = caller($i);
    no strict 'refs';
    while ($pkg) {
        if (    $pkg =~ m/^Bio::/o
            and not defined ${$pkg . "::VERSION"}
            ) {
            ${$pkg . "::VERSION"} = $VERSION;
        }
        $pkg = caller(++$i);
    }
}

1;
__END__
