#!/usr/bin/perl -w

=head1 NAME

modules.pl - information about modules in BioPerl core

=head1 SYNOPSIS

B<modules.pl> [B<-V|--verbose>] [B<-c|--count>] | [B<-l|--list>] |
  [B<-u|--untested>] | [B<-i|--info> class] | [B<-i|--inherit> |
  [B<-v|--version> | [B<-?|-h|--help>]

=head1 DESCRIPTION

This script counts, lists and provides other information about bioperl
modules. It is mainly meant to be run by bioperl maintainers.

=cut

#
# The helper class to store class status;
#
package BioClass;

sub new {
    my $class = shift;
    my $name = shift;
    die unless $name;

    my $self = {};
    $self->{'name'} = $name;
    $self->{'tested'} = 0;
    $self->{'type'} = '';
    $self->{'path'} = '';

    bless $self, $class;
}


sub name {
    my $self = shift;
    return $self->{'name'};
}

sub tested {
    my $self = shift;
    my $value = shift;
    $self->{'tested'} = 1 if defined $value && $value;
    return $self->{'tested'};
}

sub type {
    my $self = shift;
    my $value = shift;
    $self->{'type'} = $value if defined $value;
    return $self->{'type'};
}

sub path {
    my $self = shift;
    my $value = shift;
    $self->{'path'} = $value if defined $value;
    return $self->{'path'};
}

sub add_superclass {
    my $self = shift;
    my $superclass = shift;
    return unless $superclass;
    $self->{'superclasses'}->{$superclass} = 1 ;
}

sub each_superclass {
    my $self = shift;
    return  keys %{$self->{'superclasses'}};
}

sub add_used_class {
    my $self = shift;
    my $used_class = shift;
    return unless $used_class;
    $self->{'used_classes'}->{$used_class} = 1 ;
}

sub each_used_class {
    my $self = shift;
    return  keys %{$self->{'used_classes'}};
}

package main;

use File::Find;
use Getopt::Long;
use Data::Dumper;
use strict;


# declare subroutines
sub dir;
sub modules;
sub count;
sub list_all;
sub untested;
sub info;
sub inherit;
sub synopsis;
sub version;

# command line options
my ($dir, $count,$list, $verbose,$info,$untested, $inherit, $synopsis,
    $version);
GetOptions(
	   'dir:s'      => \$dir,
	   'count'    => \$count,
	   'list'     => \$list,
           'test_BioClass' => \&_test_BioClass,
           'V|verbose'  => \$verbose,
           'untested' => \$untested,
           'info:s' =>  \$info,
           'inherit' => \$inherit,
           'synopsis' => \$synopsis,
           'version' => \$version,
	   'h|help|?' => sub{ exec('perldoc',$0); exit(0) }
	   );


our %MODULES; # storage structure

# find modules
my $pwd = $ENV{PWD};
my $seachdir = "$pwd/../Bio"; #default
my %FIND_OPTIONS = ( wanted => \&modules );

$seachdir = "$pwd/$dir" if $dir;
find \%FIND_OPTIONS, $seachdir;


# call subroutines
if    ($list)     { list_all }
elsif ($untested) { untested }
elsif ($info)     { info($info) }
elsif ($inherit)  { inherit }
elsif ($synopsis) { synopsis }
elsif ($version)   { version }
else              { count }


################# end main ####################


#
# subroutines;
#

sub _test_BioClass {
    $a = new BioClass('Bio::Test');
    print "Class name: ", $a->name(), "\n";
    $a->add_superclass('Bio::Super');
    $a->add_superclass('Bio::Super2');
    $a->tested(1);
    $a->type('instance');
    print Dumper [$a->each_superclass] if $a->tested;
    print Dumper $a;
    exit;
}

sub modules {
    return unless /\.pm$/ ;
    #return unless -e $_;
    #print "file: $_\n" if $verbose;
    open (F, $_) or warn "can't open file $_: $!" && return;
    my $class;
    while (<F>) {
        if (/^package\s+([\w:]+)\s*;/) {
            #print $1, "\n" if $verbose;
            $_ = $1;
            $class = new BioClass($_);
            $MODULES{$_} = $class;
            if (/.*:[a-z]/) {
                $class->type('component');
            } elsif (/:Base/ | /Base$/) {
                $class->type('base');
            } elsif (/[^A-Z]I$/) {
                $class->type('interface');
            } else {
                $class->type('instance');
            }
            $class->path($File::Find::name);
        }
        if (/^\w*use/ && /(Bio[\w:]+)\W*;/) {
	    next unless $class;
            #print "\t$1\n" if $verbose;
            $class->add_used_class($1);
        }
        if ((/\@ISA/ || /use base/) && /Bio/) {
            next unless $class;
            my $line = $_;
            while ( $line =~ /(Bio[\w:]+)/g) {
                #print "\t$1\n" if $verbose;
                $class->add_superclass($1);
            }
        }
        if (/\@ISA/ && /Bio/) {
            next unless $class;
            my $line = $_;
            while ( $line =~ /(Bio[\w:]+)/g) {
                #print "\t$1\n" if $verbose;
                $class->add_superclass($1);
            }
        }
    }

    close F;
}

=head1 OPTIONS

Only one option is processed on each run of the script. The --verbose
is an exception, it modifies the amount of output.

=over 4

=item B<-V | --verbose>

B<INACTIVE>

Set this option if you want to see more verbose output. Often that
will mean seeing warnings normally going into STDERR.

=cut

=item B<-c | --count>

The default action if no other option is given. Gives the count of
modules broken to B<instance> ("usable"), B<base> ( (abstract)?
superclass) , B<interface> (the "I" files) and B<component> (used from
instantiable parent) modules, in addition to total number of modules.

Note that abstract superclass in bioperl is not an enforced concept and
they are not clearly indicateded in the class name.

=cut

sub count {
    printf "Instance : %3d\n",
        scalar (grep $MODULES{$_}->type =~ /instance/ , keys %MODULES);
    printf "Base     : %3d\n",
        scalar (grep $MODULES{$_}->type =~ /base/ , keys %MODULES);
    printf "Interface: %3d\n",
        scalar (grep $MODULES{$_}->type =~ /interface/ , keys %MODULES);
    printf "Component: %3d\n",
        scalar (grep $MODULES{$_}->type =~ /component/ , keys %MODULES);
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
        print $MODULES{$_}->type. "\t$_\n";
    }
}

=item B<-u | --untested>

Prints a list of instance modules which are I<not> explicitly used by
test files in t directory. Superclasess or any classes used by others
are not reported, either, since their methods are assumed to be tested
by subclass tests.

This method can not be improved much without running the tests!

=cut

sub untested {
    foreach (`find ../t -name "*.t" -print | xargs grep -hs "use "`) {
        s/^ *?use +//;
        next unless /^Bio/;
        s/[\W;]+$//;
        next unless $MODULES{$_};
        $MODULES{$_}->tested(1) 
            unless defined $MODULES{$_} and $MODULES{$_}->tested;
        foreach ($MODULES{$_}->each_superclass) {
            $MODULES{$_}->tested(1)
                unless defined $MODULES{$_} or $MODULES{$_}->tested;
        }
        foreach ($MODULES{$_}->each_used_class) {
            $MODULES{$_}->tested(1)
                unless defined $MODULES{$_} and $MODULES{$_}->tested;
        }

    }

    foreach ( sort keys %MODULES) {
        print "$_\n" if
            $MODULES{$_}->type eq 'instance' and ($MODULES{$_}->tested == 0) ;
    }

}

=item B<-i | --info> class

Dumps information about a class given as an argument.

=cut

sub info {
    my $class = shift;
    die "" unless $class;
    #print Dumper $MODULES{$class};
    my $c = $MODULES{$class};
    print $c->name, "\n";
    printf "  Type:\n\t%s\n", $c->type;
    print "  Superclasses:\n";
    foreach (sort $c->each_superclass) {
        print "\t$_\n";
    }
    print "  Used classes:\n";
    foreach (sort $c->each_used_class) {
        print "\t$_\n";
    }
}


=item B<-i | --inherit>

Finds interface modules which inherit from an instantiable class.

Could be extended to check other bad inheritance patterns.

=cut

sub inherit {
    foreach ( sort keys %MODULES) {
        my $c=$MODULES{$_};
        next unless $c->type =~ /interface/;
        foreach my $super ($c->each_superclass) {
            next if $super =~ /I$/;
            print "Check this inheritance: ", $c->name, " <-- $super\n";
        }
    }
}

=item B<-s | --synopsis>

Test SYNOPSIS section of bioperl modules for runnability

=cut

sub synopsis {
    foreach ( sort keys %MODULES) {
        my $c=$MODULES{$_};

        next unless $c->type eq "instance";
        next if $c->name eq 'Bio::Root::Version';

        my $synopsis = '';
        open (F, $c->path) or warn "can't open file ".$c->name.": $!" && return;

        my $flag = 0;
        while (<F>) {
            last if $flag && /^=/;
            $synopsis .= $_ if $flag;
            $flag = 1 if /^=head1 +SYNOPSIS/;
        }

        # remove comments
        $synopsis =~ s/[^\$]#[^\n]*//g;
        # allow linking to an other Bio module, e.g.: See L<Bio::DB::GFF>.
        $synopsis =~ s/[^\n]*L<Bio[^\n]*//g;
        # protect single quotes
        $synopsis =~ s/'/"/g;

        my $res = `perl -ce '$synopsis' 2>&1 `;
        next if $res =~ /syntax OK/;
        print $c->path, "\n";
        print $synopsis;
        print $res;
        print "-" x 70, "\n"; 
        # print "SYNOPSIS not runnable\n";
        close F;
    }
}

=item B<-v | --version>

Test the VERSION of the module against the global one set in
Bio::Root::Variation. Print out the different ones.

=cut

sub version {

    use Bio::Root::Version;
    my $version =  $Bio::Root::Version::VERSION;

    my %skip = ( # these are defined together with an other module
                 # and can not be use independently
                'Bio::AnalysisI::JobI' => 1,
                'Bio::PrimarySeq::Fasta' => 1,
                'Bio::DB::Fasta::Stream' => 1,
                'Bio::DB::GFF::ID_Iterator' => 1,
                'Bio::DB::GFF::Adaptor::dbi::faux_dbh' =>1,
                'Bio::LiveSeq::IO::SRS' =>1 # tries to call an external module
               );

    foreach ( sort keys %MODULES) {
        my $n=$MODULES{$_}->name;
        next if $skip{$n};
        my $vv= "\$${n}::VERSION";
        my $v = `perl -we 'use $n; print $vv;'`;
        printf "%50s %-3s\n", $n, $v unless $version eq $v;
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



