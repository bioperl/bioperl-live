[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.16344.svg)](http://dx.doi.org/10.5281/zenodo.16344)
[![Build Status](https://travis-ci.org/bioperl/bioperl-live.svg?branch=master)](https://travis-ci.org/bioperl/bioperl-live)
[![Coverage Status](https://coveralls.io/repos/bioperl/bioperl-live/badge.svg?branch=master)](https://coveralls.io/r/bioperl/bioperl-live?branch=master)
[![Documentation Status](https://readthedocs.org/projects/bioperl/badge/?version=latest)](https://readthedocs.org/projects/bioperl/?badge=latest)

# About BioPerl

BioPerl is a project for development of free and open source Perl
tools for computational molecular biology.  For example, it includes
classes for biological sequences, readers of multiple formats,
sequence alignments, database searching objects, and interfaces to
multiple programs such as EMBOSS, ClustalW, and BLAST.

The BioPerl project has developed multiple module distributions for
different purposes.  The one named BioPerl (named after the project)
provides the foundation for all others distributions.

This is the repository for the BioPerl distribution only.  Other
distributions have [their own
repositories](https://github.com/bioperl/).

# Installation

BioPerl distribution has the same name as the BioPerl.  However, the
BioPerl distribution only includes a subset of the project modules.
Because of this, the meaning of "installing BioPerl" is rarely clear.
Instead of "install BioPerl", the aim must be "install module X".

[CPAN.org](https://www.cpan.org/modules/INSTALL.html) provides an
overview on how to install and manage Perl modules but the bottom-line
is:

1. find the module you need, for example `Bio::DB::EUtilities`
2. install it with `cpanm`, for example `cpanm Bio::DB::EUtilities`

Alternatively, some Linux distributions have packaged BioPerl and have
it available through their package manager.

# Documentation and Support

Documentation for individual modules is in POD and can also be read
online at [metacpan](https://metacpan.org/pod/BioPerl).  Useful
documentation in the form of example code can also be found in the
`examples/` and `bin/` directories.

Additional resources and information about the project is available on
the [project website](https://bioperl.org), with discussion happening
on the [bioperl-l@bioperl.org](mailto:bioperl-l@bioperl.org) mailing
list, and on the `#bioperl` channel of the freenode IRC server.

Bug reports are handle on the distribution github page.

# Development

See the [`HACKING.md`](HACKING.md) file for details on the project
structure, such as building from source and running the test suite.
