#!/usr/bin/perl -w

=head1 NAME

modules.pl - information about modules in BioPerl core

=head1 SYNOPSIS

B<modules.pl> [B<-c|--count>] | [B<-l|--list>] | [B<-u|--untested>] |
  [B<-?|-h|--help>]

=head1 DESCRIPTION

This script counts, lists and provides other information about bioperl
modules. It is mainly meant to be ran by bioperl maintainers.

=cut

use Getopt::Long;
use strict;
#use Data::Dumper;

# declare subroutines
sub modules;
sub count;
sub list_all;
sub untested;

# command line options
my ($count,$list, $untested, $help);
GetOptions(
	   'count'    => \$count,
	   'list'     => \$list,
           'untested' => \$untested,
	   'h|help|?' => sub{ exec('perldoc',$0); exit(0) }
	   );

our %MODULES;

modules; # find modules in Bio directory

# call subroutines
if ($list) {
    list_all;
}
elsif ($untested) {
    untested;
} else {
    count;
}

#
# subroutines;
#

sub modules {
    foreach (`find ../Bio  -name "*.pm" -print`) {
        s/.pm\n//;
        s/^...//;
        s|/|::|g;
        if (/.*:[a-z]/) {
            $MODULES{$_} = 'component';
        } elsif (/[^A-Z]I$/) {
            $MODULES{$_} = 'interface';
        } else {
            $MODULES{$_} = 'instance';
        }
    }
}

=head1 OPTIONS

Only one option is processed on each run of the script.

=over 4

=item B<-c | --count>

The default action if no other option is given. Gives the count of
modules broken to B<intantance> ("usable"), B<interface> (the "I"
files) and B<component> (used from instantiable parent) modules, in
addition to total number of modules.

Nota that we are unable to separate abstact superclasses from instance
classes based on name only.

=cut

sub count {
    printf "Instance : %3d\n", scalar (grep /instance/ , values %MODULES);
    printf "Interface: %3d\n", scalar (grep /interface/ , values %MODULES);
    printf "Component: %3d\n", scalar (grep /component/ , values %MODULES);
    print  "--------------\n";
    printf "Total    : %3d\n", scalar (keys %MODULES);

}

=item B<-l | --list>

Prints all the module names in alphabetical order. The output is a tab
separated list of category (see above) and module name per line. The
output can be processed with standard UNIX command line tools.

=cut

sub list_all {
    foreach ( sort keys %MODULES) {
        print "$MODULES{$_}\t$_\n";
    }
}

=item B<-u | --untested>

Prints a list of instance modules which are I<not> explicitely used by
test files in t directory.

=cut

sub untested {
    foreach (`find ../t -name "*.t" -print | xargs grep -hs "use "`) {
        s/^ *?use +//;
        next unless /^Bio/;
        s/[\W;]+$//;
        #    print "$_\n";
        delete $MODULES{$_} if $MODULES{$_};
    }
    foreach ( sort keys %MODULES) {
        print "$_\n" if $MODULES{$_} eq 'instance';
    }

}

__END__

=item B<-? | -h  | --help>

This help text.

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Heikki Lehvaslaiho

Email heikki@ebi.ac.uk

=cut



