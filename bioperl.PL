#!/usr/bin/perl
use strict;
use LocalConfig;
use File::Basename qw(basename dirname);
my $file;
($file = basename($0)) =~ s/\.PL$/\.pod/;
open OUT,">$file" or die "Can't create $file: $!";

print OUT<<"!GROK!THIS!";
=head1 NAME

Bioperl v$Local{VERSION} - Perl Modules for Biology

=head1 SYNOPSIS

If you're new to Bioperl, you should start reading F<bptutorial.pl>, an
overview of the Bioperl toolkit.

In the style of the perl documentation, core Bioperl documentation has
been split up into the following sections:

=head2 Overview

	bioperl		Bioperl overview (this document)
	biodatabases	How to use databases with Bioperl
	biodesign	A guide for authoring a Bioperl module

=head2 Tutorials

	bptutorial.pl	Bioperl tutorial for beginners

=head2 References for Individual Modules

For ease of maintenance and coordination amongst contributors, bioperl is
code is maintained in a modular form, as is the documentation.  Refer to
the documentation for individual modules by using perldoc, i.e.

C<perldoc Bio::Seq>

to get documentation for the Bio::Seq object.

=head1 DESCRIPTION

Bioperl is the product of a community effort to produce Perl code which is
useful in biology.  Examples include Sequence objects, Alignment objects
and database searching objects. These objects not only do what they are
advertised to do in the documentation, but they also interact - Alignment
objects are made from the Sequence objects, Sequence objects have access
to Annotation and SeqFeature objects and databases, Blast objects can be
converted to Alignment objects, and so on. This means that the objects
provide a coordinated and extensible framework to do computational
biology.

Bioperl development focuses on Perl classes, or code that is used to
create objects representing biological entities. There are scripts
provided in the scripts/ directory but scripts are not the main focus of
the Bioperl developers. Of course, as the objects do most of the hard work
for you, all you have to do is combine a number of objects together
sensibly to make useful scripts.

The intent of the Bioperl development effort is to make reusable tools
that aid people in creating their own sites or job-specific applications.

The bioperl.org website at http://bioperl.org also attempts to maintain 
links and archives of standalone bio-related Perl tools that are not
affiliated or related to the core Bioperl effort. Check the site for
useful code ideas and contribute your own if possible.


=head1 DOCUMENTATION

The Bio::Perl module is designed to flatten the learning curver for
newcomers to Perl/Bioperl. This is a good place to start if you want
some simple functionality. We have a cookbook tutorial in
bptutorial.pl which has embedded documentation. Start there if
learning-by-example suits you most, or examine the Bioperl online
course at
http://www.pasteur.fr/recherche/unites/sis/formation/bioperl. Make
sure to check the documentation in the modules as well - there are
over 500 modules in Bioperl, and counting, and there's detail in the
modules' documentation that will not appear in the general
documentation.


=head1 INSTALLATION

The Bioperl modules are distributed as a tar file that expands into a
standard perl CPAN distribution.  Detailed installation directions 
can be found in the distribution INSTALL file. Installing on windows
using ActiveState Perl is covered in the INSTALL.WIN file.

The Bioperl modules can now interact with local flat file and relational
databases. To learn how to set this up, look at the biodatabases.pod
documentation ('perldoc biodatabases.pod' should work once Bioperl has
been installed).

The bioperl-db, bioperl-run, bioperl-gui, corba-server,
bioperl-pipeline and corba-client packages are installed separately
from Bioperl. Please refer to their respective documentation for more
information.


=head1 GETTING STARTED

A good place to start is by reading and running the cookbook script,
bptutorial.pl.

The distribution I<scripts/> directory has working scripts for use
with Bioperl. A list and brief description of these scripts is found
in bioscripts.pod. You are more than welcome to contribute your script!

If you have installed Bioperl in the standard way, as detailed in the
README in the distribution, these scripts should work by just running
them. If you have not installed it in a standard way you will
have to change the 'use lib' to point to your installation (see INSTALL
for details).


=head1 GETTING INVOLVED

Bioperl is a completely open community of developers. We are not
funded and we don't have a mission statement. We encourage
collaborative code, in particular in Perl. You can help us in many
different ways, from just a simple statement about how you have used
Bioperl to doing something interesting to contributing a whole new object
hierarchy. See http://bioperl.org for more information. Here are
some ways of helping us:

=head2 Asking questions and telling us you used it

We are very interested to hear how you experienced using Bioperl. Did
it install cleanly? Did you understand the documentation? Could you
get the objects to do what you wanted them to do? If Bioperl was useless
we want to know why, and if it was great - that too. Post a message to
B<bioperl-l\@bioperl.org>, the Bioperl mailing list, where all the developers
are.

Only by getting people's feedback do we know whether we are providing
anything useful.

=head2 Writing a script that uses it

By writing a good script that uses Bioperl you both show that Bioperl
is useful and probably save someone elsewhere writing it. If you
contribute it to the 'script central' at http://bioperl.org then other
people can view and use it. Don't be nervous if you've never done this
sort of work, advice is freely given and all are welcome!

=head2 Find bugs!

We know that there are bugs in there. If you find something which you are
pretty sure is a problem, post a note to bioperl-bugs\@bioperl.org and
we will get on it as soon as possible. You can also access the bug
system through the web pages.

=head2 Suggest new functionality

You can suggest areas where the objects are not ideally written and
could be done better. The best way is to find the main developer
of the module (each module was written principally by one person,
except for Seq.pm). Talk to him or her and suggest changes.

=head2 Make your own objects

If you can make a useful object we will happily include it into the
core. Probably you will want to read a lot of the documentation 
in the L<Bio::Root>, talk to people on the Bioperl mailing list,
B<bioperl-l\@bioperl.org>, and read biodesign.pod. biodesign.pod provides
documentation on the conventions and ideas used in Bioperl, it's definitely
worth a read if you would like to be a Bioperl developer.


=head1 ACKNOWLEDGEMENTS

Bioperl owes its early organizational support to its association with
the award-winning VSNS-BCD BioComputing Courses; some students of the
1996 course (Chris Dagdigian, Richard Resnick, Lew Gramer, Alessandro
Guffanti, and others) have contributed code and commentary. Georg
Fuellen, the VSNS-BCD chief organizer was one of the early driving forces
behind Bioperl. Steven Brenner, who was an early adopter of Perl for
bioinformatics provided some of the early work on Bioperl. Lincoln Stein
has long provided guidance and code.

Bioperl was then taken up by people developing code at the large
genome centres. In particular Steve Chervitz at Stanford, Ian Korf at
the Genome Sequencing Centre (St. Louis) and Ewan Birney at the Sanger
Centre (Cambridge UK).  All of the C code XS extensions were
provided by Ewan Birney. Bioperl is used in anger at these sites,
indicating that is both useful and that it works.

Jason Stajich and Hilmar Lapp joined Bioperl for the drive towards a
0.7 release over 2000 and the first part of 2001, which includes a
revised feature location model, richer feature objects (in particular
genes) and more and better tools. Peter Schattner and Lorenz Pollak
contributed serious chunks of code, being the AlignIO and bptutorial
scripts and the BPLite port to Bioperl respectively. At this time
Bioperl was being used in absolute earnest by the Ensembl group which
shook out a number of problems in the code base. Additional
compatibility with the Sequence Workbench (Bioperl-gui, Mark
Wilkinson and David Block) and Biocorba (Jason Stajich, Brad Chapman
and Alan Robinson) and finally Game-XML (Brad Marshall) provided more
interoperability.

From 2001-2002 we saw the developer ranks swell and contributions
increase with folks volunteering everything from bibliographic modules
(Martin Senger) to new SCF chromatogram file parsers (Chad Matsalla
and Aaron Mackey) to a new relational database and Bioperl layer for
fast access to feature data (Lincoln Stein) as well graphical
rendering of sequences with features (Lincoln Stein).  In early 2002
Open-Bioinformatics foundation was formally incorperated as a
non-profit organization (see http://www.open-bio.org) in the early
part of 2002. We also participated in the 1st Bio-Hackathon that
February which drew developers from several Open-Bio projects to
lovely Cape Town, South Africa and the hospitality of the co-sponsor
Electric Genetics (O'Reilly and Associates also co-sponsored the
event, and AstraZeneca contributed to the sponorship).  In early May
we released bioperl 1.0 as an effort to put a stamp of completeness on
the toolkit and submitted a paper detailing the efforts thus far (see
the FAQ).  In late 2002 we released a developer release 1.1 for
testing purposes and released bioperl 1.2 at the beginning of 2003.


Current server hardware for bioperl.org (and other open-bio.org hosted
projects) is provided by Sun. A previous server was provided by Compaq
Computer Corporation.  The compaq donation was facilitated by both the
Pharmaceutical Sales and High Performance Technical Computing (HPTC)
groups. 

The Bioperl servers reside in Cambridge, Massachusetts USA with
colocation facilities and Internet bandwidth donated by Genetics
Institute. In particular Dr. Steven Howes, Kenny Grant & 
Rich DiNunno have made significant efforts on our behalf.

Finally, this document is maintained by Allen Day.  Contact
<allenday\@ucla.edu> for comments/corrections.

$Local{AUTHORS}

=head1 COPYRIGHT

Copyright (c) 1996-2003 Georg Fuellen, Richard Resnick, Steven E. Brenner,
Chris Dagdigian, Steve Chervitz, Ewan Birney, James Gilbert, Elia Stupka,
and others. All Rights Reserved. This module is free software;  you can
redistribute it and/or modify it under the same terms as Perl itself.

=cut
!GROK!THIS!
