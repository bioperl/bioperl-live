package Deobfuscator;

# module for retrieving method-specific documentation from a
# Berkeley database
#
# first version by Dave Messina (dmessina@watson.wustl.edu) at the
# Cold Spring Harbor Laboratory Advanced Bioinformatics Course
# Oct 12-25, 2005

# part of the Deobfuscator package
# by Laura Kavanaugh and Dave Messina
#
# cared for by Dave Messina <dave-pause@davemessina.com>
#
# POD documentation - main docs before the code

=head1 NAME

Deobfuscator - get BioPerl method and package information from a Berkeley DB


=head1 VERSION

This document describes Deobfuscator version 0.0.3


=head1 SYNOPSIS

    use Deobfuscator;

    # get all the methods available to objects belonging to a class
    # (including those inherited from parent classes)
    my $hashref = Deobfuscator::return_methods('Bio::SeqIO', 'Bio::AlignIO');

    # retrieve the return values for a method
    my $method_db_ref = Deobfuscator::open_db('methods.db');
    my $ret_vals = Deobfuscator::get_method_docs( $method_db_ref,
                                                  'Bio::SeqIO::next_seq',
                                                  'returns' );
    close_db($method_db_ref);

    # retrieve the synopsis documentation for a class
    my $pkg_db_ref = Deobfuscator::open_db('packages.db');
    my $synopsis = Deobfuscator::get_pkg_docs( $pkg_db_ref,
                                              'Bio::SeqIO',
                                              'synopsis' );
    close_db($pkg_db_ref);


=head1 DESCRIPTION

The Deobfuscator module contains functions which relate to retrieving
specific types of documentation about BioPerl packages and methods.

The deob_index.pl script reads through all of the BioPerl files, extracts
the documentation, and stores it in two BerkeleyDB databases. This module
is then used to query those databases for information about a given method
or package. (see the deob_index.pl documentation for more info.)

The types of information available for individual methods include: the
usage statement, the return values, the arguments to give to the method, the
description of the function, and an example of how to use the method.

The Deobfuscator module can be used also to retrieve the synopsis and
description documentation for a given class.


=head1 DIAGNOSTICS

=over

=item C<< error: couldn't eval $module >>

A package couldn't be loaded (eval'd), which would prevent us from determining
what its methods are.

=item C<< error: couldn't open $filename >>

One of the Berkeley databases couldn't be opened. Possible causes are:
deob_index.pl wasn't run and so the databases weren't created, or the database
files aren't in the correct place.

=item C<< error: couldn't close database >>

One of the Berkeley databases couldn't be closed. This might just be a 
transient filesystem error.

=back

=item C<< error: couldn't load [module] >>

The BioPerl modules aren't in the Perl lib (PERL5LIB) and so can't be searched
(the Deobfuscator uses I<Class::Inspector> for this. Check that the value of
your PERL5LIB includes BioPerl's modules. If need be, you can set a use lub directive
at the beginning of deob_interface.cgi.

=back


=head1 CONFIGURATION AND ENVIRONMENT

This software requires:

=over

=item A working installation of the Berkeley DB

The Berkeley DB comes standard with most UNIX distributions, so you may 
already have it installed. See L<http://www.sleepycat.com> for more information.

=item BioPerl

Deobfuscator.pm recursively navigates a directory of BioPerl modules. Note
that the BioPerl module directory need not be "installed"; any old location
will do. See L<http://www.bioperl.org> for the latest version.

=back


=head1 DEPENDENCIES

L<version>, L<Class::Inspector>, L<DB_File>


=head1 INCOMPATIBILITIES

None reported.


=head1 BUGS AND LIMITATIONS

In the current implementation, Deobfuscator does not show internal or private
methods (i.e. those whose name begins with an underscore). This is simply an
option in the Class::Inspector->methods call, and so could be presented as an
option to the user (patches welcome).

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

L<deob_index.pl>


=head1 AUTHOR

Dave Messina C<< <dave-pause@davemessina.net> >>


=head1 CONTRIBUTORS

=over

=item Laura Kavanaugh

=item David Curiel

=back


=head1 ACKNOWLEDGMENTS

This software was developed originally at the Cold Spring Harbor Laboratory's
Advanced Bioinformatics Course between Oct 12-25, 2005. Many thanks to David
Curiel, who provided much-needed guidance and assistance on this project.


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2005-6 Laura Kavanaugh and Dave Messina. All Rights Reserved.

This module is free software; you may redistribute it and/or modify it under the
same terms as Perl itself. See L<perlartistic>.


=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.


=head1 APPENDIX

The rest of the documentation details each of the functions.
Internal methods are preceded with a "_".

=cut


use version; $VERSION = qv('0.0.3');
use warnings;
use strict;
use Class::Inspector;
use DB_File;

use lib './lib';


=head2 return_methods

Title   : return_methods
Usage   : $methods_hashref = Deobfuscator::return_methods('Bio::AlignIO',
          'Bio::SeqIO');
Function: traverses the inheritance tree for a given class to determine
          the methods available to objects belonging to that class

Returns : a reference to a hash. The hash keys are fully-qualified class
          names, such as 'Bio::SeqIO'. The hash values are references to
          an array of hashes, where each array element is a reference to 
          a hash containing two key-value pairs, 'method' and 'class';

Args    : a list of fully-qualified class names

=cut

sub return_methods {

    my @input = @_;

    # key: full class name
    # value: a reference to an array of hashes
    #    where each array element is a pointer to a hash
    #    which contains two key: 'method' and 'class'
    my %methods_of;


    foreach my $class (@input) {

        # fancy eval so that we can loop through different modules
        my $retval = _load_module($class);
	if ($retval) { die "error: couldn't load $class: $retval\n"; }

        # methods returned from Class::Inspector as:
        # [
        #   [ 'Class::method1',   'Class',   'method1', \&Class::method1   ],
        #   [ 'Another::method2', 'Another', 'method2', \&Another::method2 ],
        #   [ 'Foo::bar',         'Foo',     'bar',     \&Foo::bar         ],
        # ]
        my $methods_aryref3
            = Class::Inspector->methods( $class, 'expanded', 'public' );


        for ( my $i = 0; $i < scalar @{$methods_aryref3}; $i++ ) {
            foreach my $meth ( $methods_aryref3->[$i] ) {
                my $method_name    = $meth->[2];
                my $inherited_from = $meth->[1];
                push @{$methods_of{$class}}, [$method_name, $inherited_from];
            }

        }

    }
    return \%methods_of;
}

=head2 print_methods

Title   : print_methods
Usage   : print_methods('Bio::AlignIO','Bio::SeqIO');
Function: traverses the inheritance tree for a given class to determine
           the methods available to objects belonging to that class, then
           pretty-prints the resulting information.
Returns : nothing. But it does print to the current filehandle (usually
           STDOUT).
Args    : a list of fully-qualified class names

=cut

sub print_methods {

    my @input = @_;

    foreach my $class (@input) {

        # fancy eval so that we can loop through different modules
        my $retval = _load_module($class);
	if ($retval) { die "error: couldn't load $class: $retval\n"; }

        # methods returned as
        # [
        #   [ 'Class::method1',   'Class',   'method1', \&Class::method1   ],
        #   [ 'Another::method2', 'Another', 'method2', \&Another::method2 ],
        #   [ 'Foo::bar',         'Foo',     'bar',     \&Foo::bar         ],
        # ]
        my $methods_aryref3
            = Class::Inspector->methods( $class, 'expanded', 'public' );

        print "methods for $class\n";
        print "=========================================\n";

        for ( my $i = 0; $i < scalar @{$methods_aryref3}; $i++ ) {
            print "method $i\n";
            foreach my $meth ( $methods_aryref3->[$i] ) {
                print "\t           class: $meth->[1]\n";
                print "\t          method: $meth->[2]\n";
            }
            print "--------------------------------------\n";
        }

    }

}

=head2 _load_module

Title   : _load_module
Usage   : * INTERNAL USE ONLY *
Function: attempts to load a module
Returns : nothing. But it does die upon failure to load.
Args    : a module name

=cut

sub _load_module {
    my $module = shift;
    eval "require $module";
    my $err = $@ || 'eval returned undef';
    
    if ($@) { return $@ }
    else { return }
}

=head2 open_db

Title   : open_db
Usage   : open_db($filename)
Function: opens a Berkeley DB
Returns : a hashref tied to the DB
Args    : a filename as a scalar

=cut

sub open_db {
    my ($filename) = @_;
    
    my %hash;
    my $hashref = \%hash;
    
    tie %hash, "DB_File", $filename or die "error: couldn't open $filename: $!\n";
    
    return $hashref;
}

=head2 close_db

Title   : close_db
Usage   : closes a Berkeley DB
Function: closes a database
Returns : nothing.
Args    : a hashref to a tied Berkeley DB

=cut

sub close_db {
    my ($hashref) = @_;
    
    untie $hashref or die "error: couldn't close database: $!\n";
}

=head2 get_pkg_docs

Title   : get_pkg_docs
Usage   : get_pkg_docs($db_hashref, 'Class name', 'documentation type');
Function: returns a specified part of the documentation for a class
Returns : a string containing the desired documentation or ' ' if the
          documentation doesn't exist
Args    : - $db_hashref is the ref to the hash tied to the DB
          - Class name is of the form 'Bio::SeqIO'
          - documentation type is the subfield of the method's POD.
          The possible values of documentation type are:
          short_desc, synopsis, desc

=cut

sub get_pkg_docs {
    my ($db_hashref, $pkg_name, $info_type) = @_;

    # hash to store our hash value, now split out into its constituent parts
    my %record;
        
    my $rec_sep = 'DaVe-ReC-sEp';

    # if the method isn't in our db
    if ( ! exists($db_hashref->{$pkg_name}) ) {
        return 0;
    }

    # grab the constituent parts of the pkg record
    ( $record{'short_desc'}, $record{'synopsis'}, $record{'desc'} ) = 
        ( split $rec_sep, $db_hashref->{$pkg_name} );

    # return just the part that was asked for
    if ( exists($record{$info_type}) ) {
        return $record{$info_type};
    }
    else { return ' '; }
}

=head2 get_method_docs

Title   : get_method_docs
Usage   : get_method_docs($db_hashref, 'Class+method name', 'documentation type');
Example : get_method_docs($db_hashref, 'Bio::SeqIO::next_aln', 'args');
Function: returns a specified part of the documentation for a class's method
Returns : a string containing the desired documentation, or 0 if the
         desired documentation doesn't exist
Args    : - $db_hashref is the ref to the hash tied to the DB
          - Class+method name is of the form 'Bio::SeqIO::next_aln',
            where Bio::SeqIO is the class and next_aln is the method.
          - documentation type is the subfield of the method's POD.
            The possible values of documentation type are:
            title, usage, function, returns, args

=cut

sub get_method_docs {
    my ($db_hashref, $meth_name, $info_type) = @_;

    my %record;
    my $whole_record;

    my $rec_sep = 'DaVe-ReC-sEp';

    # if the method isn't in our db
    if ( !exists( $db_hashref->{$meth_name} ) ) {
        return 0;
    }

    # separate the sub-records using the record separator and field tag
    my @parts = split $rec_sep, $db_hashref->{$meth_name};

    # put individual info types into separate hash entries...
    foreach my $part (@parts) {
        if ($part =~ /^-(\w+)\|(.*)/) { $record{$1} = $2; }
    
    # ... and put the whole thing into one big string
        $whole_record .= "$part\n";
    }

    # return a specific part if that was asked for
    if ($info_type) {
        # return just the part that was asked for      
        if ( exists( $record{$info_type} ) ) {

			# if there's really nothing in there, say so.
			if ( ( $record{$info_type} =~ /^[\s\n]*$/)
         	|| ( $record{$info_type} eq '') ) { return 0; }
			else { 
				return $record{$info_type};
			}
        }
        # or return everything
        else { return $whole_record; }
    }
    # otherwise return whole record
    else {
        return $whole_record;
    }
}

=head2 htmlify

Title   : htmlify
Usage   : htmlify($string);
Example : htmlify('this is a : doc);
Function: does some crude reformatting of POD method documentation by swapping
          isolated colons (':') into HTML <br> tags
Returns : a string
Args    : a string

=cut

sub htmlify {
	my ($string) = @_;

	# change isolated colons into <br> tags
	$string =~ s/\s:\s/ <br> /g;

    # change L<> POD link into HTML link
    if ( $string =~ /L<(.+)>/ ) {
        $string = urlify_pkg($1);
    }

	return $string;
}

=head2 urlify_pkg

Title   : urlify_pkg
Usage   : urlify_pkg($string);
Example : urlify('this is a : doc);
Function: wraps a package name in an HTML href pointing to the bioperl.org
          pdoc docs for that package
Returns : a string (an href in HTML)
Args    : a string

=cut

sub urlify_pkg {
    my ($pkg_name) = @_;
    my $bioperl_doc_url = q{http://doc.bioperl.org/bioperl-live/};

    my $pkg_as_path = $pkg_name;

    # convert Bio::DB::RefSeq to Bio/DB/RefSeq
    $pkg_as_path =~ s/::/\//g;
    my $url  = $bioperl_doc_url . $pkg_as_path . '.html';
    my $href = qq{<a href="$url">$pkg_name</a>};

    return $href;
}


1;
__END__
