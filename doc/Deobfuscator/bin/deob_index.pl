#!/usr/bin/perl

# deob_index.pl
# part of the Deobfuscator package
# by Laura Kavanaugh and Dave Messina
#
# cared for by Dave Messina <dave-pause@davemessina.net>
#
# POD documentation - main docs before the code

=head1 NAME

deob_index.pl - extracts BioPerl documentation and indexes it in a database for easy retrieval

=head1 VERSION

This document describes deob_index.pl version 0.0.3


=head1 SYNOPSIS

deob_index.pl <path to BioPerl lib> <output path>

=over

=item <path to BioPerl lib>

a directory path pointing to the root of the BioPerl lib tree. e.g. /export/share/lib/perl5/site_perl/5.8.7/Bio/

=item <output path>

where you would like deob_index.pl to put its output files.

=back


=head1 DESCRIPTION

deob_index.pl goes through the entire BioPerl library tree looking for
.pm and .pl files. For each one it finds, it tries to extract module-level
POD documentation (e.g. SYNOPSIS, DESCRIPTION) and store it in a BerkeleyDB.
It also tries to extract documentation for each method in the module and
store that in a separate BerkeleyDB.

Specific parts of the documentation for a module or method may be retrieved
individually using the functions available in Deobfuscator.pm. See that module
for details.

While going through and trying to parse each module, deob_index.pl also
reports what pieces of the documentation it can't find. For example, if
a method's documentation doesn't describe the data type it returns, this 
script logs that information to a file. This type of automated documentation-
checking could be used to standardize and improve the documentation in 
BioPerl.

deob_index.pl creates four files:

=over

=item C<< package_list.txt >>

A plaintext file listing each package found in the BioPerl directory that was
searched. Packages are listed by their module names, such as 'Bio::SeqIO'.
This file is used by L<deob_interface.cgi>.

=item C<< packages.db >>

A Berkeley DB, which stores package-level documentation, such as
the synopsis and the description. Each key is a package name,
e.g. "Bio::SeqIO", and each value string is composed of the 
individual pieces of the documentation kept separate by 
unique string record separators. The individual pieces of 
documentation are pulled out of the string using the 
get_pkg_docs function in Deobfuscator.pm. See that package
for details.

=item C<< methods.db >>

Like packages.db, methods.db is also a Berkeley DB, except it 
stores various pieces of information about individual methods
available to a class. Each method might have documentation
about its usage, its arguments, its return values, an example,
and a description of its function. 

Each key is the fully-qualified method name, e.g.
"Bio::SeqIO::next_seq". Each value is a string containing all
of the pieces of documentation concatenated together and
separated by unique strings serving as record separators. The
extraction of the actual documentation in these strings is
handled by the get_method_docs subroutine in the Deobfuscator.pm
module. See that package for details.

Not all methods will have all of these types of documentation,
and some methods will not have the different pieces of
information clearly labeled and separated. For the latter type,
deob_index.pl will try to store whatever free-form
documentation that does exist, and the get_method_docs function
in Deobfuscator.pm, if called without arguments, will return
that documentation.

=item C<< deob_index.log >>

This file contains detailed information about errors
encountered while trying to extract documentation during
the indexing process.

Each line in deob_index.log is a key-value pair describing
a single parsing error.

=back


=head1 DIAGNOSTICS

These are the parsing error codes reported in 'deob_index.log'.

=head2 Package errors

=over

=item C<< PKG_NAME >>

couldn't find the name of the package

=item C<< SYNOPSIS >>

couldn't find the synopsis

=item C<< DESC >>

couldn't find the description

=item C<< METHODS >>

couldn't find any methods

=item C<< PKG_DUP >>

This package name occurs more than once

=back

=head2 Method errors

=over

=item C<< FUNCTION >>

couldn't find the function description

=item C<< EXAMPLE >>

couldn't find the example

=item C<< ARGS >>

couldn't find the method's arguments

=item C<< USAGE >>

couldn't find the usage statement

=item C<< RETURNS >>

couldn't find the return values

=item C<< FREEFORM >>

This method's documentation doesn't conform to the BioPerl standard of having
clearly-labeled fields for title, function, example, args, usage, and returns.

=item C<< METH_DUP >>

This method name occurs more than once

=back


=head1 CONFIGURATION AND ENVIRONMENT

This software requires:

=over

=item A working installation of the Berkeley DB

The Berkeley DB comes standard with most UNIX distributions, so you may 
already have it installed. See L<http://www.sleepycat.com> for more information.

=item BioPerl

deob_index.pl recursively navigates a directory of BioPerl modules. Note
that the BioPerl module directory need not be "installed"; any old location
will do. See L<http://www.bioperl.org> for the latest version.

=back


=head1 DEPENDENCIES

L<version>, L<File::Find>, L<DB_File>


=head1 INCOMPATIBILITIES

None reported.


=head1 BUGS AND LIMITATIONS

No bugs have been reported.

deob_index.pl currently expects the sections of POD in a BioPerl module to
be in a particular order, namely: NAME, SYNOPSIS, DESCRIPTION, CONSTRUCTORS,
... , APPENDIX. Those sections are expected to be marked with =head1 POD tags,
and the documentation for each method is expected to be in =head2 sections
in the APPENDIX. The order of SYNOPSIS and DESCRIPTION can be flipped, but
this behavior should not be taken as encouragement to do so.

Most, but not all BioPerl modules conform to this standard. Those that do not
will cause deob_index.pl to report them as errors. Although the consistency
of this standard is desirable for end-users of the documentation, this code
probably needs to be a little bit more flexible (patches welcome!).

This software has only been tested in a UNIX environment. 


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

  https://github.com/bioperl/bioperl-live/issues


=head1 SEE ALSO

L<Deobfuscator>, L<deob_interface.cgi>, L<deob_detail.cgi>


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

=cut

use version; $VERSION = qv('0.0.2');
use warnings;
use strict;
use File::Find;
use DB_File;
use IO::File;
use Getopt::Std;
use File::Spec;

# GetOpt::Std-related settings
$Getopt::Std::STANDARD_HELP_VERSION = 1;
getopts('s:x:');

my $DEBUG = 0;

my $usage = "
deob_index.pl - extracts and parses BioPerl POD
and stores the info in a database.

USAGE: deob_index.pl [-s bioperl-version] [-x exclude_file] <BioPerl lib dir> <output dir>

where 

<BioPerl lib dir> is the BioPerl distribution you'd like to index

    e.g. /export/share/lib/perl5/site_perl/5.8.7/Bio/
    
and

<output dir> is where the output files should be placed

OPTIONS:
-s    user-supplied string to declare BioPerl's version
      (which will be displayed by deob_interface.cgi)
-x    excluded modules file (a module paths to skip; see POD for details)
";

unless ( @ARGV == 2 ) { die $usage; }

my ( $source_dir, $dest_dir ) = @ARGV;

# check source_dir for full path and repair if it's a relative path
unless ( File::Spec->file_name_is_absolute( $source_dir ) ) {
    $source_dir = File::Spec->rel2abs( $source_dir ) ;
}

# check dest_dir for full path and repair if it's a relative path
unless ( File::Spec->file_name_is_absolute( $dest_dir ) ) {
    $dest_dir = File::Spec->rel2abs( $dest_dir ) ;
}

# NOTE: we're allowing only one source directory, but File::Find supports
# passing an array of dirs.

# read in an optional list of modules to exclude from indexing
# - this is aimed at modules with external dependencies that are often not
# - present and thus will prevent deob_interface.cgi from loading them
our ($opt_s, $opt_x);
my %exclude;
if (defined $opt_x) {
    my $exclude_fh = IO::File->new($opt_x, "r")
        or die "couldn't open $opt_x\n";
    while (<$exclude_fh>) {
        chomp;
        next if ( /^\#/ || /^\s*$/ ); # ignore comments and blank lines
        $exclude{$_} = 1;
    }
    print STDERR "Found ", scalar keys %exclude, " modules to be excluded.\n";
}


# save a list of the BioPerl modules to a file
my $list; # filehandle
my $list_file = $dest_dir . "/package_list.txt";
if ( -e $list_file) { unlink($list_file); }
open $list, ">$list_file" or die "deob_index.pl: couldn't open $list_file:$!\n";
my @list_holder; # hold all package names so we can sort them before writing.

# record misbehaving BioPerl docs to a file
my $log;    # filehandle
my $logfile = $dest_dir . "/deob_index.log";
open $log, ">$logfile" or die "deob_index.pl: couldn't open $logfile:$!\n";

# create databases
my $meth_file = $dest_dir . '/methods.db';
if ( -e $meth_file ) { unlink($meth_file); }    # remove for production?
my $meth_db = create_db($meth_file) or die "deob_index.pl: couldn't create $meth_file: $!\n";
my $pkg_file = $dest_dir . '/packages.db';
if ( -e $pkg_file ) { unlink($pkg_file); }      # remove for production?
my $pkg_db = create_db($pkg_file) or die "deob_index.pl: couldn't create $pkg_file: $!\n";

# used to make sure we're parsing in the right order
my %FLAG;

# store version string in packages.db
$pkg_db->{'__BioPerl_Version'} = $opt_s ? $opt_s : 'unknown';

# keep stats on our indexing
my %stats = ( 
              'files'    => 0,
              'pkg_name' => 0,
              'desc'     => 0,
              'synopsis' => 0,
              'methods'  => 0,
            );

# wanted points to the subroutine which is run on each found file
# ( in this program, that subroutine is &extract_pod )
# no_chdir prevents find from chdir'ing into each subsequent directory
my %FIND_OPTIONS = ( wanted => \&extract_pod);#, no_chdir => 1 );

# This is the important line - Find::File actually doing the
# traversal of the directory tree.
find( \%FIND_OPTIONS, $source_dir );

# sort and write out package list
foreach my $sorted_pkg (sort @list_holder) {
    print $list $sorted_pkg, "\n";
}

# store user-supplied BioPerl version number

# output stats
print STDOUT "\nThis indexing run found:\n";
print $log "\nThis indexing run found:\n";
foreach my $stat ( 'files', 'pkg_name', 'desc', 'synopsis', 'methods' ) {
    printf STDOUT "%5d %s\n", $stats{$stat}, $stat;
    printf $log "%5d %s\n", $stats{$stat}, $stat;
}

# close files and DBs
untie $meth_db or die "deob_index.pl: couldn't close $meth_file: $!\n";
untie $pkg_db  or die "deob_index.pl: couldn't close $pkg_file: $!\n";
close $list    or die "deob_index.pl: couldn't close $list: $!\n";
close $log     or die "deob_index.pl: couldn't close $log: $!\n";
my $mode = 0666;
chmod($mode, $pkg_file, $meth_file, $list_file);

### Parsing subroutines ###
sub extract_pod {
    my ($file) = $_;
    my $long_file = $File::Find::name;
    
    # skip if it's on our exclude list
    foreach my $one (keys %exclude) {
        if ($File::Find::name =~ /$one$/) {
            print STDERR "Excluding $file\n";
            print $log "Excluding $file\n";
            return;
        }
    }

    # skip unless it's a perl file that exists
    return unless ( $file =~ /\.PLS$/ ) or ( $file =~ /\.p[ml]$/ );
    return unless -e $file;

    $stats{'files'}++;

    open my $fh, '<', $File::Find::name or die "deob_index.pl: could not read file '$file': $!\n";

    # these have to be done in order
    my ( $pkg_name, $short_desc ) = get_pkg_name($fh);
    my ($synopsis, $desc);
    LOOP: while (my ($type, $section) = get_generic($fh) ) {
        if    ($type eq 'synopsis')    { $synopsis = $section; }
        elsif ($type eq 'description') { $desc     = $section; }
        else { last LOOP; }
    }

    my $constructors = get_constructors($fh);
    my $methods      = get_methods($fh);

    # record package name to our package list file
    if ($pkg_name) { push @list_holder, $pkg_name; }

    # store valid package data here
    my @pkg_data;

    # error reporting
    if ($pkg_name) {
        $stats{'pkg_name'}++;
        print $pkg_name, "\n" if $DEBUG == 1;
    }
    else {
        print $log " PKG_NAME: $long_file\n";
    }
    if ($short_desc) {
        $stats{'short_desc'}++;
        push @pkg_data, $short_desc;
        print $short_desc, "\n" if $DEBUG == 1;
    }
    else {
		push @pkg_data, 'no short description available'; # store something
        print $log "SHORT_DESC: $long_file\n";
    }
    if ($synopsis) {
        $stats{'synopsis'}++;
        print $synopsis, "\n" if $DEBUG == 1;
        push @pkg_data, $synopsis;
    }
    else {
		push @pkg_data, 'no synopsis available'; # store something
        print $log " SYNOPSIS: $long_file\n";
    }
    if ($desc) {
        $stats{'desc'}++;
        print $desc, "\n" if $DEBUG == 1;
        push @pkg_data, $desc;
    }
    else {
		push @pkg_data, 'no description available'; # store something
        print $log "     DESC: $long_file\n";
    }
    if ($methods) {
        my $method_count = scalar keys %$methods;
        print "**** Found $method_count methods in $pkg_name\n"
            if $DEBUG == 2;
        foreach my $method ( keys %$methods ) {
            $stats{'methods'}++;
            print $method, "\n//\n" if $DEBUG == 2;
        }
    }
    else {
        print $log "  METHODS: $long_file\n";
    }

    # prepare data for databases
    my $pkg_record   = pkg_prep(@pkg_data);
    my $meth_records = meth_prep( $pkg_name, $methods );

    # load data in databases
    if ($pkg_name) {
        pkg_load( $pkg_db, $pkg_name, $pkg_record );
        meth_load( $meth_db, $meth_records );
    }
}

sub slurp_until_next {
    my ($fh) = @_;

    my @lines;
    my $prev_line = $_;


    LINE: while (<$fh>) {
        next LINE if $_ eq $prev_line;

        # if it's a POD directive
        if (/^\=/) {

            # reset our position to the beginning of the line
            # so it is seen as part of the next POD section
            seek $fh, -length($_), 1;
            last LINE;
        }
        else {
            push @lines, $_;
        }
    }
    return join q{}, @lines;
}

sub get_pkg_name {
    my ($fh) = @_;

    my $pkg_name;
    my $short_desc;

    LINE: while (<$fh>) {
        chomp;
        print "**", $_, "\n" if $DEBUG == 2;

        # grab package name
        # - "short desc" is the one-line description of the package
        if ( $_ =~ /^\=head1\s+NAME/ ) {
            <$fh>;
            my $next_line = <$fh>;
            ( $pkg_name, $short_desc ) = split /\s+/, $next_line, 2;
			$short_desc .= slurp_until_next($fh);

            # strip off leading dash
            $short_desc =~ s/^(\-)+\s+//;

			# strip off trailing spaces
			$short_desc =~ s/\s+$//;

			# strip any newlines
			$short_desc =~ s/\n/ /;

            print $pkg_name, "\n" if $DEBUG == 1;

            last LINE;
        }

        # we've hit a =head1, but it's the wrong one
        elsif ( $_ =~ /^\=head1\s+/ ) {
            last LINE;
        }
    }
    if ($pkg_name) {
        $FLAG{'pkg_name'} = 1;
        return $pkg_name, $short_desc;
    }
}

sub get_generic {
    my ($fh) = @_;

    my $section;

    LINE: while (<$fh>) {
        chomp;
        print "**", $_, "\n" if $DEBUG == 2;

        if ( $_ =~ /^\=head1\s+SYNOPSIS/ ) {
            $section = slurp_until_next($fh);
            if ($section) {
                $FLAG{'synopsis'} = 1;
                return ('synopsis', $section);
            }
            else { last LINE; }
        }
        elsif ( $_ =~ /^\=head1\s+DESCRIPTION/ ) {
            $section = slurp_until_next($fh);
            if ($section) {
                $FLAG{'description'} = 1;
                return ('description', $section);
            }
            else { last LINE; }
        }

        # if we hit the APPENDIX, time to stop
        elsif (/^\=head1\s+APPENDIX/) {

            # reset our position to the beginning of the line
            # so it is seen by the next parser
            seek $fh, -length($_)*2, 1;
            last LINE;
        }
    }
}

sub get_synopsis {
    my ($fh) = @_;

    my $synopsis;

    LINE: while (<$fh>) {
        chomp;
        print "**", $_, "\n" if $DEBUG == 2;

        if ( $_ =~ /^\=head1\s+SYNOPSIS/ ) {
            $synopsis = slurp_until_next($fh);
            last LINE;
        }

        # we've hit a =head1, but it's the wrong one
        elsif ( $_ =~ /^\=head1\s+/ ) {
            last LINE;
        }
    }
    if ($synopsis) {
        $FLAG{'synopsis'} = 1;
        return $synopsis;
    }
}

sub get_desc {
    my ($fh) = @_;

    my $desc;

    LINE: while (<$fh>) {
        chomp;
        print "**", $_, "\n" if $DEBUG == 2;

        if ($_ =~ /^=head1\s+VERSION/ ) {
            slurp_until_next($fh);
        }

        if ( $_ =~ /^\=head1\s+DESCRIPTION/ ) {
            $desc = slurp_until_next($fh);
            last LINE;
        }

        # we've hit a =head1, but it's the wrong one
        elsif ( $_ =~ /^\=head1\s+/ ) {
            last LINE;
        }
    }
    if ($desc) {
        $FLAG{'description'} = 1;
        return $desc;
    }
}

sub get_constructors {

    # not implemented

    # should return a hashref
}

sub get_methods {
    my ($fh) = @_;
    my %methods;

    # we shouldn't see any methods until after the APPENDIX
    my $seen_appendix = 0;

    # there's an '=cut' after we enter the APPENDIX
    # we know the method '=head2' tags will come after it
    my $seen_first_cut = 0;

    LINE: while (<$fh>) {
        if ( $_ =~ /^\=head1\s+APPENDIX/ ) {
            $seen_appendix = 1;
        }

        # this should be the first tag after the APPENDIX
        if ( $seen_appendix && $_ =~ /^\=cut/ ) {
            $seen_first_cut = 1;
        }

        # this should be a method
        if ( $seen_first_cut && $_ =~ /^\=head2\s+(\S+)/ ) {
            $methods{$1} = slurp_until_next($fh);
        }
    }

    # returns a hashref
    return \%methods;
}

### Database subroutines ###
sub create_db {
    my ($filename) = @_;

    my %hash;
    my $hashref = \%hash;

    tie %hash, "DB_File", $filename
        or die "ERROR: couldn't open $filename:$!\n";

    return $hashref;
}

sub pkg_prep {

    # unique string on which to split our sub-records
    my $rec_sep = 'DaVe-ReC-sEp';

    my $record = join $rec_sep, @_;

    return $record;
}

sub meth_prep {
    my ( $pkg_name, $methods ) = @_;
    my %records;

    foreach my $entry ( keys %$methods ) {
        my $key = $pkg_name . '::' . $entry;
        my $record;    # what will be stored in the db
        my $rec_sep = 'DaVe-ReC-sEp';

        # if the method conforms to the BioPerl doc spec,
        # we will split it into constituent pieces before storing
        # it in the db. If not, we store the whole thing as one lump.

        my $last;      # for grabbing multi-line entries
        my %fields = (
            'title'    => '',
            'usage'    => '',
            'function' => '',
            'example'  => '',
            'returns'  => '',
            'args'     => '',
        );


        my @lines = split "\n", $methods->{$entry};
        foreach my $line (@lines) {
            if ( $line =~ /^\s+Title\s+:(.*)/ ) {
                next if $1 =~ /^\s+$/;
                $fields{'title'} = $1;
                $last = \$fields{'title'};
            }
            elsif ( $line =~ /^\s+Usage\s+:(.*)/ ) {
                next if $1 =~ /^\s+$/;
                $fields{'usage'} = $1;
                $last = \$fields{'usage'};
            }
            elsif ( $line =~ /^\s+Function\s?:(.*)/ ) {
                next if $1 =~ /^\s+$/;
                $fields{'function'} = $1;
                $last = \$fields{'function'};
            }
            elsif ( $line =~ /^\s+Example\s+:(.*)/ ) {
                next if $1 =~ /^\s+$/;
                $fields{'example'} = $1;
                $last = \$fields{'example'};
            }
            elsif ( $line =~ /^\s+Returns\s+:(.*)/ ) {
                next if $1 =~ /^\s+$/;
                $fields{'returns'} = $1;
                $last = \$fields{'returns'};
            }
            elsif ( $line =~ /^\s+Args\s+:(.*)/ ) {
                next if $1 =~ /^\s+$/;
                $fields{'args'} = $1;
                $last = \$fields{'args'};
            }

            # grab multi-line entries
            elsif ( $line =~ /^\s{8,}(\s.*)/ ) { $$last .= $1; }
        }

        # debugging
        if ( $DEBUG == 2 ) {
            print "** $entry **\n";
            foreach my $field ( keys %fields ) {
                print STDOUT $field, "\t", $fields{$field}, "\n";
            }
            print "\n";
        }

        # if any of our fields have a value, store subrecords
        my $filled_fields = grep /\w+/, values %fields;
        print STDERR $key, "\t", $filled_fields, "\n" if $DEBUG == 3;
        if ( $filled_fields > 0 ) {
            if ( !$fields{'title'} ) { print $log '    TITLE: ', $key, "\n"; }
            if ( !$fields{'usage'} ) { print $log '    USAGE: ', $key, "\n"; }
            if ( !$fields{'function'} ) {
                print $log ' FUNCTION: ', $key, "\n";
            }
            if ( !$fields{'example'} ) {
                print $log '  EXAMPLE: ', $key, "\n";
            }
            if ( !$fields{'returns'} ) {
                print $log '  RETURNS: ', $key, "\n";
            }
            if ( !$fields{'args'} ) { print $log '     ARGS: ', $key, "\n"; }

            # create the records to be stored in the db
            foreach my $field ( keys %fields ) {
                my $subrecord
                    = $rec_sep . '-' . $field . '|' . $fields{$field};
                $record .= $subrecord;
            }

            # store the records
            $records{$key} = $record;
        }

        # if no subfields, store whatever docs we do have for the method
        else {
            $record = $methods->{$entry};
            print $log ' FREEFORM: ', $key, "\n";
        }
    }
    return \%records;
}

sub pkg_load {
    my ( $pkg_db, $pkg_name, $record ) = @_;

    if ( exists $pkg_db->{$pkg_name} ) {
        print $log '  PKG_DUP: ', $pkg_name, "\n";
        warn(
            "$pkg_name already exists in package db!\n",
            "existing record:\n$pkg_db->{$pkg_name}\n",
            "attempted to add:\n$record\n",
            )
            if $DEBUG == 2;
    }
    else {
        $pkg_db->{$pkg_name} = $record;
    }
}

sub meth_load {
    my ( $meth_db, $records ) = @_;

    foreach my $method ( keys %$records ) {
        if ( exists( $meth_db->{$method} ) ) {
            print $log ' METH_DUP: ', $method, "\n";
            warn(
                "$method already exists in method db!\n",
                "existing record:\n$meth_db->{$method}\n",
                "attempted to add:\n$records->{$method}\n",
                )
                if $DEBUG == 2;
        }
        else {
            $meth_db->{$method} = $records->{$method};
        }
    }
}

__END__
