#!/usr/bin/perl -w

# Deob_detail.cgi
# part of the Deobfuscator package
# by Laura Kavanaugh and Dave Messina
#
# cared for by Dave Messina <dave-pause@davemessina.net>
#
# POD documentation - main docs before the code

=head1 NAME

deob_detail.cgi - displays a web page of detailed information about a BioPerl method

=head1 VERSION

This document describes deob_detail.cgi version 0.0.3


=head1 SYNOPSIS

This program is designed to be called by deob_interface.cgi. See
L</"DESCRIPTION"> for details.

To install deob_detail.cgi and the rest of the Deobfuscator package, see the
README.


=head1 DESCRIPTION

Deob_detail.cgi is called by deob_interface.cgi when a user clicks on a
method name. This program extracts the documentation about that method from
the Deobfuscator Berkeley DBs and returns it in some simple HTML formatting.


=head1 DIAGNOSTICS

None.


=head1 CONFIGURATION AND ENVIRONMENT

This program expects to have the 'methods.db' and 'packages.db' files in the
same directory as itself. These two files are automatically generated when
L<deob_index.pl> is run. If your installation requires that they be in a
different location, change the $BerkeleyDB_packages and $BerkeleyDB_methods
variables below to be fully qualified paths to the db files.


=head1 DEPENDENCIES

L<version>, L<CGI>, L<Deobfuscator>


=head1 INCOMPATIBILITIES

None reported.


=head1 BUGS AND LIMITATIONS

No bugs have been reported.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                       - General discussion
  http://www.bioperl.org/wiki/Mailing_lists   - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.bioperl.org/


=head1 SEE ALSO

L<Deobfuscator>, L<deob_interface.cgi>, L<deob_index.pl>


=head1 AUTHOR

Laura Kavanaugh


=head1 CONTRIBUTORS

=over

=item Dave Messina C<< <dave-pause@davemessina.net> >>

=item David Curiel

=back


=head1 ACKNOWLEDGMENTS

This software was developed originally at the Cold Spring Harbor Laboratory's
Advanced Bioinformatics Course between Oct 12-25, 2005. Many thanks to David
Curiel, who provided much-needed guidance and assistance on this project.


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2005-6 Laura Kavanaugh and Dave Messina. All Rights Reserved.

This module is free software; you may redistribute it and/or modify it under the same terms as Perl itself. See L<perlartistic>.


=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut


# Let the code begin...

## HARDCODED VALUES ##
# Change these to fit your installation.
use lib './lib';
my $BerkeleyDB_packages = './packages.db';
my $BerkeleyDB_methods  = './methods.db';

## You shouldn't need to change anything below here ##

use version; $VERSION = qv('0.0.2');
use warnings;
use strict;
use CGI ':standard';
use Deobfuscator;

# Open BerkeleyDBs
my $packages_ref = Deobfuscator::open_db($BerkeleyDB_packages);
my $methods_ref  = Deobfuscator::open_db($BerkeleyDB_methods);

# 'method' is the name of the method passed in from deob_interface.cgi
my $class_method = param('method');

# Get all of the documentation fields out of the db
my $title
    = Deobfuscator::get_method_docs( $methods_ref, $class_method, "title" );
if ( $title eq "0" ) { $title = "not documented"; }

my $usage
    = Deobfuscator::get_method_docs( $methods_ref, $class_method, "usage" );
if ( $usage eq "0" ) { $usage = "not documented"; }

my $function = Deobfuscator::get_method_docs( $methods_ref, $class_method,
    "function" );
if ( $function eq "0" ) { $function = "not documented"; }

my $returns
    = Deobfuscator::get_method_docs( $methods_ref, $class_method, "returns" );
if ( $returns eq "0" ) { $returns = "not documented"; }

my $args
    = Deobfuscator::get_method_docs( $methods_ref, $class_method, "args" );
if ( $args eq "0" ) { $args = "not documented"; }

### Make the output page

# Start the page
print header;
print start_html($class_method);

# Define some styles
my $style1
    = qq{style="border-collapse:collapse;border:solid black 1px;font-family:verdana;font-size:10px;background-color:lightgrey"};
my $style2
    = qq{style="border-collapse:collapse;border:solid black 1px;font-family:verdana;font-size:10px"};
my $style3
    = qq{style="border-collapse:collapse;border:solid black 1px;font-family:verdana;font-size:14px"};

# open the table
print '<div style="border:solid black 1px; width:100%; height:200; overflow:auto">';
print '<table width="100%" $style3>';
print "<tr><td colspan=4><center>$class_method</center></td></tr>";

my @sections = ('Usage', 'Function', 'Returns', 'Args');
my $sec_ndx = 0;

foreach my $section ($usage, $function, $returns, $args) {

	my $section_html = Deobfuscator::htmlify($section);
	print "<tr><td $style1>$sections[$sec_ndx++]</td><td $style2>$section_html</td></tr>\n";
}

# close the table
print "</table></div>";

# finish the page
print end_html;

# close BerkeleyDB
Deobfuscator::close_db($BerkeleyDB_packages);
Deobfuscator::close_db($BerkeleyDB_methods);

__END__
