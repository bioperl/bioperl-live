#!/usr/bin/perl

=head1 NAME

find_mod_deps.pl - inspect B<only hard-coded> dependencies of sets of perl files

=head1 DESCRIPTION

Inspects the hard-coded dependencies of a set of perl files and prints
a summary of which modules they use (by default not including
inter-dependencies between the modules being inspected).

=head1 USAGE

find_mod_deps.pl [options]  [ path ... ]

If given any paths, inspects only the files in those paths.  Defaults
to inspecting all perl files in the current directory.

=head2 Options

=over 4

=item -i

If set, also print internal dependencies, i.e. the inter-dependencies
between the files we are inspecting.

=item -B

If set, print the dependencies in a format suitable for cutting and
pasting directly into a Build.PL (i.e. Module::Build)

=item -M

If set, print the dependencies in a format suitable for cutting and
pasting directly into a Makefile.PL (i.e. Module::Install)

=item -Z

If set, print the dependencies in a format suitable for cutting and
pasting directly into a dist.ini (i.e. Dist::Zilla).  Although, if
you're using Dist::Zilla, you probably have it configured to be
auto-discovering deps, and it will find the same deps as this script.

=back

=head1 AUTHOR

Robert Buels, rbuels@cpan.org

=cut

use strict;
use warnings;

use File::Find;
use Getopt::Std;
use IO::String;
use List::MoreUtils qw/ first_value all /;
use Module::CoreList;
use Pod::Strip;
use Pod::Usage;

use Data::Dump 'dump';
use Hash::Merge;

my %opt;
getopts('iBMZ', \%opt) or pod2usage();

-d './lib' or -d './bin' or -d './scripts' or die "run this script from the root dir of a distribution\n";

my @paths = @ARGV;

@paths = qw( t lib scripts bin cgi-bin Bio )
   unless @paths;

# expand any dirs into the perl files they contain
my @perl_files = map {
    if( -d ) {
        my @f;
        find( sub { push @f, $File::Find::name if is_perl_file($_) },
              $_,
            );
        @f
    } elsif( -e ) {
        if( is_perl_file($_) ) {
            $_
        } else {
            warn "WARNING: skipping user-specified file $_, since it is not a perl file.\n";
            ()
        }
    } else {
        ()
    }
} @paths;

my %perl_files = map { $_ => 1 } @perl_files;

my %deps;
my $merger = Hash::Merge->new('RETAINMENT_PRECEDENT');
for my $file ( @perl_files ) {
    my $deps = find_deps( $file );
    %deps = %{ $merger->merge( \%deps, $deps ) };
}

# classify the deps
my %classified;
for my $modname ( keys %deps ) {
    if( all { m|^(./)?t/| } @{$deps{$modname}} ) {
        $classified{build_requires}{$modname} = $deps{$modname};
    }
    else {
        $classified{requires}{$modname} = $deps{$modname};
    }
}

# decide which format to print in
if( $opt{B} ) {
    for ( values %classified ) {
        $_ = 0 for values %$_;
    }
    print dump \%classified;
} elsif( $opt{M} ) {
    print "requires      '$_' => 0;\n"
      for sort { lc $a cmp lc $b } keys %{$classified{requires}};
    print "test_requires '$_' => 0;\n"
      for sort { lc $a cmp lc $b } keys %{$classified{build_requires}};
} elsif( $opt{Z} ) {
    print "[Prereqs]\n";
    print "$_  =  0\n"
      for sort { lc $a cmp lc $b } keys %{$classified{requires}};

    print "\n[Prereqs / TestRequires]\n";
    print "$_  =  0\n"
      for sort { lc $a cmp lc $b } keys %{$classified{build_requires}};

} else {
    print dump \%classified;
}

exit;

################## helpers #####################

sub modfile {
    my $modname = shift;
    my $modfile = "$modname.pm";
    $modfile =~ s|::|/|g;
    return first_value {
        $_ =~ /$modfile$/;
    } @perl_files;
}

sub namespace_parent {
    my $modname = shift;
    $modname =~ s/(?:::)?[^:]+$//;
    return $modname;
}

sub find_deps {
    my ( $file ) = @_;

    my $nopod;
    { open my $p, '<', $file or die "Could not read file '$file': $!\n";
      local $/;
      my $code = <$p>;
      my $strip = Pod::Strip->new;
      $strip->output_string(\$nopod);
      $strip->parse_string_document( $code );
    }
    my $f = IO::String->new( \$nopod );

    my %deps;
    while( my $depline = <$f> ) {
        $depline =~ s/#.+//; #remove comments
        next unless $depline =~ /^\s*(use|require|extends|with)\s+.+;/;
        next unless $depline && $depline =~ /\S/;

        my @toks = $depline =~ /([\w:]{3,})/ig
            or die 'cannot parse: '.$depline;

        #warn "    adding to $k->{name}\n";
        shift @toks;
        if( @toks ) {
            if ( $toks[0] eq 'base' ) {
                shift @toks;
                shift @toks if $toks[0] eq 'qw';
            } else {
                @toks = ($toks[0]);
            }
        }

      MODNAME:
        foreach my $modname (@toks) {

            chomp $depline;
            #warn "'$depline' goes to $modname\n";

            #skip if the module is in the distribution
            my $modfile = modfile($modname);
            next if !$opt{i} && $modfile && -f $modfile;

            #skip if the module is in core before 5.6
            my $rl = Module::CoreList->first_release($modname);
            next if $rl && $rl <= 5.006;

            #skip if the module is actually defined in a parent file
            my $p = $modname;
            while( $p = namespace_parent($p) ) {
                my $p_modfile = modfile($p);
                #warn  "checking $p / $p_modfile\n";

                next unless $p_modfile && -f $p_modfile;

                open my $p, '<', $p_modfile or die "Could not read file '$p_modfile': $!\n";
                while( <$p> ) {
                    next MODNAME if /^\s*package\s+$p\b/;
                }
            }

            push @{$deps{$modname} ||= []}, $file;
        }
    }

    return \%deps;
}

sub is_perl_file {
    local $_ = shift;
    return -f && ( -x || /\.(pm|t|pl)$/ );
}
