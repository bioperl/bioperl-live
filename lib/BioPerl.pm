package BioPerl;

use utf8;
use strict;
use warnings;

# ABSTRACT: Perl modules for biology
# AUTHOR: See the individual modules for their authors.
# OWNER: many people (see the individual modules for their copyright holders)
# LICENSE: Perl_5

=head1 SYNOPSIS

If you're new to BioPerl, you should start reading the BioPerl HOWTO's:

L<http://bioperl.org/howtos/index.html>

=head2 References for Individual Modules

For ease of maintenance and coordination amongst contributors, BioPerl
code is maintained in a modular form, as is the documentation.  Refer to
the documentation for individual modules by using perldoc, i.e.

C<perldoc Bio::Seq>

to get documentation for the Bio::Seq object.

=head1 DESCRIPTION

BioPerl is the product of a community effort to produce Perl code which is
useful in biology. Examples include Sequence objects, Alignment objects and
database searching objects. These objects not only do what they are advertised
to do in the documentation, but they also interact - Alignment objects are made
from the Sequence objects, Sequence objects have access to Annotation and
SeqFeature objects and databases, Blast objects can be converted to Alignment
objects, and so on. This means that the objects provide a coordinated and
extensible framework to do computational biology.

BioPerl development focuses on Perl classes, or code that is used to create
objects representing biological entities. There are scripts provided in the
scripts/ and examples/ directories but scripts are not the main focus of the
BioPerl developers. Of course, as the objects do most of the hard work for you,
all you have to do is combine a number of objects together sensibly to make
useful scripts.

The intent of the BioPerl development effort is to make reusable tools that aid
people in creating their own sites or job-specific applications.

The BioPerl website at L<http://bioperl.org> also attempts to maintain links
and archives of standalone bio-related Perl tools that are not affiliated or
related to the core BioPerl effort. Check the site for useful code ideas and
contribute your own if possible.

=head1 INSTALLATION

The BioPerl modules are distributed as a tar file that expands into a standard
perl CPAN distribution. Detailed installation directions can be found in the
distribution INSTALL file. Installing on windows using ActiveState Perl is
covered in the INSTALL.WIN file.  We highly suggest reading the installation
instructions on the BioPerl website:

L<http://bioperl.org/INSTALL.html>

Note that only the following are supported at this time with the current API:

=over 3

=item BioPerl-db

=item BioPerl-network

=item BioPerl-run

=item BioPerl-pedigree

=item Bio::Graphics

=back

=head1 GETTING STARTED

The distribution I<scripts/> directory has working scripts for use with BioPerl,
check the self-described I<examples/> directory as well. You are more than
welcome to contribute your script!

If you have installed BioPerl in the standard way, as detailed in the INSTALL in
the distribution, these scripts should work by just running them.

=head1 GETTING INVOLVED

BioPerl is a completely open community of developers. We are not funded and we
don't have a mission statement. We encourage collaborative code, in particular
in Perl. You can help us in many different ways, from just a simple statement
about how you have used BioPerl to doing something interesting to contributing a
whole new object hierarchy. See L<http://bioperl.org> for more information. Here
are some ways of helping us:

=head2 Asking questions and telling us you used it

We are very interested to hear how you experienced using BioPerl. Did it install
cleanly? Did you understand the documentation? Could you get the objects to do
what you wanted them to do? If BioPerl was useless we want to know why, and if
it was great - that too. Post a message to B<bioperl-l@bioperl.org>, the BioPerl
mailing list, where all the developers are.

Only by getting people's feedback do we know whether we are providing anything
useful.

=head2 Writing a script that uses it

By writing a good script that uses BioPerl you both show that BioPerl is useful
and probably save someone elsewhere writing it. If you contribute it to the
'script central' at L<http://bioperl.org> then other people can view and use it.
Don't be nervous if you've never done this sort of work, advice is freely given
and all are welcome!

=head2 Find bugs!

We know that there are bugs in this code. If you find something which you are pretty
sure is a problem, post a bug report using our Bugzilla tracking system:

L<https://github.com/bioperl/bioperl-live/issues>

Please read the main bug tracking (L<http://www.bioperl.org/articles/Bugs.html>) for an
overview of what we expect in a bug report. Specifically, having a code and
data example where appropriate helps tremendously.  We gladly accept all patches
after a quick code review.

=head2 Suggest new functionality

You can suggest areas where the objects are not ideally written and could be
done better. The best way is to find the main developer of the module (each
module was written principally by one person, except for Seq.pm). Talk to him or
her and suggest changes.

=head2 Make your own objects

If you can make a useful object we will happily include it into the core.
Probably you will want to read a lot of the documentation in L<Bio::Root::Root>
and talk to people on the BioPerl mailing list, B<bioperl-l@bioperl.org>.

=head2 Writing documentation

We appreciate good documentation. It's what tells the world what's in BioPerl,
it's what instructs the user, it's what describes the rationale and inner
workings of the package. Feel free to contribute.

=head1 ACKNOWLEDGEMENTS

For a more detailed history of the BioPerl project, we recommend the History of
BioPerl:

L<http://bioperl.org/articles/History_of_BioPerl.html>

=cut


1;
