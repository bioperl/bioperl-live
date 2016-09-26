package Bio::Root::Build;
use Bio::Root::Version;
use strict;
use warnings;

=head1 SYNOPSIS

  ...TO BE ADDED

=head1 DESCRIPTION

This is a subclass of Module::Build so we can override certain methods and do
fancy stuff

It was first written against Module::Build::Base v0.2805. Many of the methods
here are copy/pasted from there in their entirety just to change one or two
minor things, since for the most part Module::Build::Base code is hard to
cleanly override.

B<Note>: per bug 3196, the majority of the code in this module has been revised
or commented out to bring it in line with the Module::Build API. In particular,
'requires/recommends' tags in the Build.PL file were not of the same format as
those for Module::Build, and so caused serious issues with newer versions
(including giving incorrect meta data). Other problematic methods involving
automatic installation of prereq modules via CPAN were also removed as they do
not work with more modern perl tools such as perlbrew and cpanm.

=head1 AUTHOR Sendu Bala

=cut

BEGIN {
    # we really need Module::Build to be installed
    eval "use base 'Module::Build'; 1" or die "This package requires Module::Build v0.42 or greater to install itself.\n$@";

    # ensure we'll be able to reload this module later by adding its path to inc
    use Cwd;
    use lib Cwd::cwd();
}

our @extra_types = qw(options excludes_os feature_requires test); # test must always be last in the list!
our $checking_types = "requires|conflicts|".join("|", @extra_types);

our $VERSION = $Bio::Root::Version::VERSION;

=head2 find_pm_files

Our modules are in Bio, not lib
=cut

sub find_pm_files {
    my $self = shift;
    foreach my $pm (@{$self->rscan_dir('Bio', qr/\.pm$/)}) {
        $self->{properties}{pm_files}->{$pm} = File::Spec->catfile('lib', $pm);
    }

    $self->_find_file_by_type('pm', 'lib');
}

=head2 choose_scripts

Ask what scripts to install (this method is unique to bioperl)
=cut

sub choose_scripts {
    my $self = shift;
    my $accept = shift;

    # we can offer interactive installation by groups only if we have subdirs
    # in scripts and no .PLS files there
    opendir(my $scripts_dir, 'scripts') or die "Can't open directory 'scripts': $!\n";
    my $int_ok = 0;
    my @group_dirs;

    # only retain top-level script directories (the 'categories')
    while (my $thing = readdir($scripts_dir)) {
        next if $thing =~ /^\./;
        $thing = File::Spec->catfile('scripts', $thing);
        if (-d $thing) {
            $int_ok = 1;
            push(@group_dirs, $thing);
        }
    }
    closedir($scripts_dir);
    my $question = $int_ok ? "Install [a]ll BioPerl scripts, [n]one, ".
        "or choose groups [i]nteractively?" : "Install [a]ll BioPerl scripts ".
        "or [n]one?";

    my $prompt = $accept ? 'a' : $self->prompt($question, 'a');

    if ($prompt =~ /^[aA]/) {
        $self->log_info("  - will install all scripts\n");
        $self->notes(chosen_scripts => 'all');
    }
    elsif ($prompt =~ /^[iI]/) {
        $self->log_info("  - will install interactively:\n");

        my @chosen_scripts;
        foreach my $group_dir (@group_dirs) {
            my $group = File::Basename::basename($group_dir);
            print "    * group '$group' has:\n";

            my @script_files = @{$self->rscan_dir($group_dir, qr/\.PLS$|\.pl$/)};
            foreach my $script_file (@script_files) {
                my $script = File::Basename::basename($script_file);
                print "      $script\n";
            }

            my $result = $self->prompt("    Install scripts for group '$group'? [y]es [n]o [q]uit", 'n');
            die if $result =~ /^[qQ]/;
            if ($result =~ /^[yY]/) {
                $self->log_info("      + will install group '$group'\n");
                push(@chosen_scripts, @script_files);
            }
            else {
                $self->log_info("      - will not install group '$group'\n");
            }
        }

        my $chosen_scripts = @chosen_scripts ? join("|", @chosen_scripts) : 'none';

        $self->notes(chosen_scripts => $chosen_scripts);
    }
    else {
        $self->log_info("  - won't install any scripts\n");
        $self->notes(chosen_scripts => 'none');
    }

    print "\n";
}

=head2 script_files

Our version of script_files doesn't take args but just installs those scripts
requested by the user after choose_scripts() is called. If it wasn't called,
installs all scripts in scripts directory
=cut

sub script_files {
    my $self = shift;

    unless (-d 'scripts') {
        return {};
    }

    my $chosen_scripts = $self->notes('chosen_scripts');
    if ($chosen_scripts) {
        return if $chosen_scripts eq 'none';
        return { map {$_, 1} split(/\|/, $chosen_scripts) } unless $chosen_scripts eq 'all';
    }

    return $_ = { map {$_,1} @{$self->rscan_dir('scripts', qr/\.PLS$|\.pl$/)} };
}

=head2 prompt

Overridden simply to not print the default answer if chosen by hitting return
=cut

sub prompt {
    my $self = shift;
    my $mess = shift or die "prompt() called without a prompt message";

    my $def;
    if ( $self->_is_unattended && !@_ ) {
        die <<EOF;
ERROR: This build seems to be unattended, but there is no default value
for this question.  Aborting.
EOF
    }
    $def = shift if @_;
    ($def, my $dispdef) = defined $def ? ($def, "[$def] ") : ('', ' ');

    local $|=1;
    print "$mess $dispdef";

    my $ans = $self->_readline();

    if ( !defined($ans)        # Ctrl-D or unattended
         or !length($ans) ) {  # User hit return
        #print "$def\n"; didn't like this!
        $ans = $def;
    }

    return $ans;
}

=head2 ACTION_manifest

We always generate a new MANIFEST instead of allowing existing files to remain
MANIFEST.SKIP is left alone
=cut

sub ACTION_manifest {
    my ($self) = @_;
    if ( -e 'MANIFEST' || -e 'MANIFEST.SKIP' ) {
        $self->log_warn("MANIFEST files already exist, will overwrite them\n");
        unlink('MANIFEST');
    }
    require ExtUtils::Manifest;  # ExtUtils::Manifest is not warnings clean.
    local ($^W, $ExtUtils::Manifest::Quiet) = (0,1);
    ExtUtils::Manifest::mkmanifest();
}

=head2 ACTION_install

Extended to run scripts post-installation
=cut

sub ACTION_install {
    my ($self) = @_;
    require ExtUtils::Install;
    $self->depends_on('build');
    ExtUtils::Install::install($self->install_map,
                               !$self->quiet,
                               0,
                               $self->{args}{uninst} || 0);
    #$self->run_post_install_scripts;
}

=head2 test_internet

For use with auto_features, which should require LWP::UserAgent as one of
its reqs

Note: as of 4-11-11, this is no longer called - if someone wants to run
network tests (off by default) w/o a network, then they are hanging themselves
by their own shoelaces.
=cut

sub test_internet {
    eval {require LWP::UserAgent;};
    if ($@) {
        # ideally this won't happen because auto_feature already specified
        # LWP::UserAgent, so this sub wouldn't get called if LWP not installed
        return "LWP::UserAgent not installed";
    }
    my $ua = LWP::UserAgent->new;
    $ua->timeout(10);
    $ua->env_proxy;
    my $response = $ua->get('http://search.cpan.org/');
    unless ($response->is_success) {
        return "Could not connect to the internet (http://search.cpan.org/)";
    }
    return;
}

=head2 ACTION_ppmdist

Don't copy across man3 docs since they're of little use under Windows and
have bad filenames
=cut

sub ACTION_ppmdist {
    my $self = shift;
    my @types = $self->install_types(1);
    $self->SUPER::ACTION_ppmdist(@_);
    $self->install_types(0);
}

=head2 install_types

When supplied a true value, pretends libdoc doesn't exist (preventing man3
installation for ppmdist). when supplied false, they exist again
=cut

sub install_types {
    my ($self, $no_libdoc) = @_;
    $self->{no_libdoc} = $no_libdoc if defined $no_libdoc;
    my @types = $self->SUPER::install_types;
    if ($self->{no_libdoc}) {
        my @altered_types;
        foreach my $type (@types) {
            push(@altered_types, $type) unless $type eq 'libdoc';
        }
        return @altered_types;
    }
    return @types;
}

=head2 ACTION_dist

We make all archive formats we want, not just .tar.gz
we also auto-run manifest action, since we always want to re-create
MANIFEST and MANIFEST.SKIP just-in-time
=cut

sub ACTION_dist {
    my ($self) = @_;

    $self->depends_on('manifest');
    $self->depends_on('distdir');

    my $dist_dir = $self->dist_dir;

    $self->make_zip($dist_dir);
    $self->make_tarball($dist_dir);
    $self->delete_filetree($dist_dir);
}

=head2 ACTION_clean

Define custom clean/realclean actions to rearrange config file cleanup
=cut

sub ACTION_clean {
    my ($self) = @_;
    $self->log_info("Cleaning up build files\n");
    foreach my $item (map glob($_), $self->cleanup) {
        $self->delete_filetree($item);
    }
    $self->log_info("Cleaning up configuration files\n");
    $self->delete_filetree($self->config_dir);
}

=head2 ACTION_realclean

Define custom clean/realclean actions to rearrange config file cleanup
=cut

sub ACTION_realclean {
    my ($self) = @_;
    $self->depends_on('clean');
    for my $method (qw(mymetafile mymetafile2 build_script)) {
        if ($self->can($method)) {
            $self->delete_filetree($self->$method);
            $self->log_info("Cleaning up $method data\n");
        }
    }
}

=head2 get_metadata

This wraps the base metafile method to add in version information from
Bio::Root::Version to META.json and META.yml if it isn't already present. Note
this should be compliant with meta_add and meta_merge, but occurs after those
steps. If a version is already set and dist_version differs from the set one, a
warning is printed.

=cut

sub get_metadata {
    my ($self, %args) = @_;
    my $metadata = $self->SUPER::get_metadata(%args);
    
    if (exists $metadata->{provides}) {
        my $ver = $self->dist_version;
        my $pkgs = $metadata->{provides};
        for my $p (keys %{$pkgs}) {
            if (!exists($pkgs->{$p}->{'version'})) {
                $pkgs->{$p}->{'version'} = $ver;
            } else {
                $self->log_warn("Note: Module $p has a set version: ".$pkgs->{$p}->{'version'}."\n")
                    if $pkgs->{$p}->{'version'} ne $ver;
            }
        }
    }
    return $metadata;
}

=head2 make_zip

Makes zip file for windows users and bzip2 files as well
=cut

sub make_zip {
    my ($self, $dir, $file) = @_;
    $file ||= $dir;

    $self->log_info("Creating $file.zip\n");
    my $zip_flags = $self->verbose ? '-r' : '-rq';
    $self->do_system($self->split_like_shell("zip"), $zip_flags, "$file.zip", $dir);

    $self->log_info("Creating $file.bz2\n");
    require Archive::Tar;
    # Archive::Tar versions >= 1.09 use the following to enable a compatibility
    # hack so that the resulting archive is compatible with older clients.
    $Archive::Tar::DO_NOT_USE_PREFIX = 0;
    my $files = $self->rscan_dir($dir);
    Archive::Tar->create_archive("$file.tar", 0, @$files);
    $self->do_system($self->split_like_shell("bzip2"), "-k", "$file.tar");
}

=head2 prompt_for_network

A method that can be called in a Build.PL script to ask the user if they want
internet tests.
Should only be called if you have tested for yourself that
$build->feature('Network Tests') is true
=cut

sub prompt_for_network {
    my ($self, $accept) = @_;

    my $proceed = $accept ? 0 : $self->y_n(  "Do you want to run tests that require connection to servers across the internet\n"
                                           . "(likely to cause some failures)? y/n", 'n');

    if ($proceed) {
        $self->notes('network' => 1);
        $self->log_info("  - will run internet-requiring tests\n");
        my $use_email = $self->y_n("Do you want to run tests requiring a valid email address? y/n",'n');
        if ($use_email) {
            my $address = $self->prompt("Enter email address:");
            $self->notes(email => $address);
        }
    }
    else {
        $self->notes(network => 0);
        $self->log_info("  - will not run internet-requiring tests\n");
    }
}

=head2 print_build_script

Override the build script warnings flag
=cut

sub print_build_script {
    my ($self, $fh) = @_;

    my $build_package = $self->build_class;

    my $closedata="";

    my $config_requires;
    if ( -f $self->metafile ) {
        my $meta = eval { $self->read_metafile( $self->metafile ) };
        $config_requires = $meta && $meta->{configure_requires}{'Module::Build'};
    }
    $config_requires ||= 0;

    my %q = map {$_, $self->$_()} qw(config_dir base_dir);

    $q{base_dir} = Win32::GetShortPathName($q{base_dir}) if $self->is_windowsish;

    $q{magic_numfile} = $self->config_file('magicnum');

    my @myINC = $self->_added_to_INC;
    @myINC = map { $_ = File::Spec->canonpath( $_ );
                   $_ =~ s/([\\\'])/\\$1/g;
                   $_;
                  } @myINC;
    # Remove duplicates
    @myINC = sort {$a cmp $b}
             keys %{ { map { $_ => 1 } @myINC } };

    foreach my $key (keys %q) {
        $q{$key} = File::Spec->canonpath( $q{$key} );
        $q{$key} =~ s/([\\\'])/\\$1/g;
    }

    my $quoted_INC = join ",\n", map "         '$_'", @myINC;
    my $shebang = $self->_startperl;
    my $magic_number = $self->magic_number;

    # unique to bioperl, shut off overly verbose warnings on windows, bug 3215
    my $w = $^O =~ /win/i ? '# no warnings (win)' : '$^W = 1;  # Use warnings';

    print $fh <<EOF;
$shebang

use strict;
use Cwd;
use File::Basename;
use File::Spec;

sub magic_number_matches {
    return 0 unless -e '$q{magic_numfile}';
    open my \$FH, '<', '$q{magic_numfile}' or return 0;
    my \$filenum = <\$FH>;
    close \$FH;
    return \$filenum == $magic_number;
}

my \$progname;
my \$orig_dir;
BEGIN {
    $w
    \$progname = basename(\$0);
    \$orig_dir = Cwd::cwd();
    my \$base_dir = '$q{base_dir}';
    if (!magic_number_matches()) {
        unless (chdir(\$base_dir)) {
            die ("Could not chdir '\$base_dir', aborting\\n");
        }
        unless (magic_number_matches()) {
            die ("Configuration seems to be out of date, please re-run 'perl Build.PL' again.\\n");
        }
    }
    unshift \@INC,
        (
$quoted_INC
        );
}

close(*DATA) unless eof(*DATA); # ensure no open handles to this script

use $build_package;
Module::Build->VERSION(q{$config_requires});

# Some platforms have problems setting \$^X in shebang contexts, fix it up here
\$^X = Module::Build->find_perl_interpreter;

if (-e 'Build.PL' and not $build_package->up_to_date('Build.PL', \$progname)) {
    warn "Warning: Build.PL has been altered.  You may need to run 'perl Build.PL' again.\\n";
}

# This should have just enough arguments to be able to bootstrap the rest.
my \$build =
    $build_package->resume( properties => { config_dir => '$q{config_dir}',
                                              orig_dir   => \$orig_dir, },
);

\$build->dispatch;
EOF
}

1;
