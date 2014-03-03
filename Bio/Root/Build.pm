package Bio::Root::Build;
use strict;
use warnings;

# ABSTRACT: a common Module::Build subclass base for BioPerl distributions
# AUTHOR:   Sendu Bala <bix@sendu.me.uk>
# OWNER:    Sendu Bala
# LICENSE:  Perl_5

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

=cut

BEGIN {
    # we really need Module::Build to be installed
    eval "use base 'Module::Build'; 1" or die "This package requires Module::Build v0.2805 or greater to install itself.\n$@";

    # ensure we'll be able to reload this module later by adding its path to inc
    use Cwd;
    use lib Cwd::cwd();
}

our $VERSION = '1.006924'; # pre-1.7
our @extra_types = qw(options excludes_os feature_requires test); # test must always be last in the list!
our $checking_types = "requires|conflicts|".join("|", @extra_types);

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

# extended to handle extra checking types
#sub features {
#    my $self = shift;
#    my $ph = $self->{phash};
#
#    if (@_) {
#        my $key = shift;
#        if ($ph->{features}->exists($key)) {
#            return $ph->{features}->access($key, @_);
#        }
#
#        if (my $info = $ph->{auto_features}->access($key)) {
#            my $failures = $self->prereq_failures($info);
#            my $disabled = grep( /^(?:\w+_)?(?:$checking_types)$/, keys %$failures ) ? 1 : 0;
#            return !$disabled;
#        }
#
#        return $ph->{features}->access($key, @_);
#    }
#
#    # No args - get the auto_features & overlay the regular features
#    my %features;
#    my %auto_features = $ph->{auto_features}->access();
#    while (my ($name, $info) = each %auto_features) {
#        my $failures = $self->prereq_failures($info);
#        my $disabled = grep( /^(?:\w+_)?(?:$checking_types)$/, keys %$failures ) ? 1 : 0;
#        $features{$name} = $disabled ? 0 : 1;
#    }
#    %features = (%features, $ph->{features}->access());
#
#    return wantarray ? %features : \%features;
#}
#*feature = \&features;

# overridden to fix a stupid bug in Module::Build and extended to handle extra
# checking types
#sub check_autofeatures {
#    my ($self) = @_;
#    my $features = $self->auto_features;
#
#    return unless %$features;
#
#    $self->log_info("Checking features:\n");
#
#    my $max_name_len = 0; # this wasn't set to 0 in Module::Build, causing warning in next line
#    $max_name_len = ( length($_) > $max_name_len ) ? length($_) : $max_name_len for keys %$features;
#
#    while (my ($name, $info) = each %$features) {
#        $self->log_info("  $name" . '.' x ($max_name_len - length($name) + 4));
#        if ($name eq 'PL_files') {
#            print "got $name => $info\n";
#            print "info has:\n";
#            while (my ($key, $val) = each %$info) {
#                print "  $key => $val\n";
#            }
#        }
#
#        if ( my $failures = $self->prereq_failures($info) ) {
#            my $disabled = grep( /^(?:\w+_)?(?:$checking_types)$/, keys %$failures ) ? 1 : 0;
#            $self->log_info( $disabled ? "disabled\n" : "enabled\n" );
#
#            my $log_text;
#            while (my ($type, $prereqs) = each %$failures) {
#                while (my ($module, $status) = each %$prereqs) {
#                    my $required = ($type =~ /^(?:\w+_)?(?:requires|conflicts)$/) ? 1 : 0;
#                    my $prefix = ($required) ? '-' : '*';
#                    $log_text .= "    $prefix $status->{message}\n";
#                }
#            }
#            $self->log_warn($log_text) if $log_text && ! $self->quiet;
#        }
#        else {
#            $self->log_info("enabled\n");
#        }
#    }
#
#    $self->log_info("\n");
#}

# TODO: STDERR output redirect is causing some installations to fail, commenting
# out until a fix is in place

# overriden just to hide pointless ugly warnings
#sub check_installed_status {
#    my $self = shift;
#
#    open (my $olderr, ">&". fileno(STDERR));
#    my $null = ($^O =~ m/mswin/i) ? 'NUL' : '/dev/null';
#    open(STDERR, $null);
#    my $return = $self->SUPER::check_installed_status(@_);
#    open(STDERR, ">&". fileno($olderr));
#    return $return;
#}

# extend to handle option checking (which takes an array ref) and code test
# checking (which takes a code ref and must return a message only on failure)
# and excludes_os (which takes an array ref of regexps).
# also handles more informative output of recommends section

#sub prereq_failures {
#    my ($self, $info) = @_;
#
#    my @types = (@{ $self->prereq_action_types }, @extra_types);
#    $info ||= {map {$_, $self->$_()} @types};
#
#    my $out = {};
#    foreach my $type (@types) {
#        my $prereqs = $info->{$type} || next;
#
#        my $status = {};
#        if ($type eq 'test') {
#            unless (keys %$out) {
#                if (ref($prereqs) eq 'CODE') {
#                    $status->{message} = &{$prereqs};
#
#                    # drop the code-ref to avoid Module::Build trying to store
#                    # it with Data::Dumper, generating warnings. (And also, may
#                    # be expensive to run the sub multiple times.)
#                    $info->{$type} = $status->{message};
#                }
#                else {
#                    $status->{message} = $prereqs;
#                }
#                $out->{$type}{'test'} = $status if $status->{message};
#            }
#        }
#        elsif ($type eq 'options') {
#            my @not_ok;
#            foreach my $wanted_option (@{$prereqs}) {
#                unless ($self->args($wanted_option)) {
#                    push(@not_ok, $wanted_option);
#                }
#            }
#
#            if (@not_ok > 0) {
#                $status->{message} = "Command line option(s) '@not_ok' not supplied";
#                $out->{$type}{'options'} = $status;
#            }
#        }
#        elsif ($type eq 'excludes_os') {
#            foreach my $os (@{$prereqs}) {
#                if ($^O =~ /$os/i) {
#                    $status->{message} = "This feature isn't supported under your OS ($os)";
#                    $out->{$type}{'excludes_os'} = $status;
#                    last;
#                }
#            }
#        }
#        else {
#            while ( my ($modname, $spec) = each %$prereqs ) {
#                $status = $self->check_installed_status($modname, $spec);
#                next if $status->{ok};
#
#                if ($type =~ /^(?:\w+_)?conflicts$/) {
#                    $status->{conflicts} = delete $status->{need};
#                    $status->{message} = "$modname ($status->{have}) conflicts with this distribution";
#                }
#                elsif ($type =~ /^(?:\w+_)?recommends$/) {
#                    my ($preferred_version, $why, $by_what) = split("/", $spec);
#                    $by_what = join(", ", split(",", $by_what));
#                    $by_what =~ s/, (\S+)$/ and $1/;
#
#                    $status->{message} = (!ref($status->{have}) && $status->{have} eq '<none>'
#                                  ? "Optional prerequisite $modname is not installed"
#                                  : "$modname ($status->{have}) is installed, but we prefer to have $preferred_version");
#
#                    $status->{message} .= "\n   (wanted for $why, used by $by_what)";
#
#                    if ($by_what =~ /\[circular dependency!\]/) {
#                        $preferred_version = -1;
#                    }
#
#                    #my $installed = $self->install_optional($modname, $preferred_version, $status->{message});
#                    #next if $installed eq 'ok';
#                    #$status->{message} = $installed unless $installed eq 'skip';
#                }
#                elsif ($type =~ /^feature_requires/) {
#                    # if there is a test code-ref, drop it to avoid
#                    # Module::Build trying to store it with Data::Dumper,
#                    # generating warnings.
#                    delete $info->{test};
#                }
#                else {
#                    my $installed = $self->install_required($modname, $spec, $status->{message});
#                    next if $installed eq 'ok';
#                    $status->{message} = $installed;
#                }
#
#                $out->{$type}{$modname} = $status;
#            }
#        }
#    }
#
#    return keys %{$out} ? $out : return;
#}

# install an external module using CPAN prior to testing and installation
# should only be called by install_required or install_optional
#sub install_prereq {
#    my ($self, $desired, $version, $required) = @_;
#
#    if ($self->under_cpan) {
#        # Just add to the required hash, which CPAN >= 1.81 will check prior
#        # to install
#        $self->{properties}{requires}->{$desired} = $version;
#        $self->log_info("   I'll get CPAN to prepend the installation of this\n");
#        return 'ok';
#    }
#    else {
#        my $question = $required ?  "$desired is absolutely required prior to installation: shall I install it now using a CPAN shell?" :
#                                    "To install $desired I'll need to open a CPAN shell right now; is that OK?";
#        my $do_install = $self->y_n($question.' y/n', 'y');
#
#        if ($do_install) {
#            # Here we use CPAN to actually install the desired module, the benefit
#            # being we continue even if installation fails, and that this works
#            # even when not using CPAN to install.
#            require Cwd;
#            require CPAN;
#
#            # Save this because CPAN will chdir all over the place.
#            my $cwd = Cwd::cwd();
#
#            CPAN::Shell->install($desired);
#            my $msg;
#            my $expanded = CPAN::Shell->expand("Module", $desired);
#            if ($expanded && $expanded->uptodate) {
#                $self->log_info("\n\n*** (back in BioPerl Build.PL) ***\n * You chose to install $desired and it installed fine\n");
#                $msg = 'ok';
#            }
#            else {
#                $self->log_info("\n\n*** (back in BioPerl Build.PL) ***\n");
#                $msg = "You chose to install $desired but it failed to install";
#            }
#
#            chdir $cwd or die "Cannot chdir() back to $cwd: $!";
#            return $msg;
#        }
#        else {
#            return $required ? "You chose not to install the REQUIRED module $desired: you'd better install it yourself manually!" :
#                               "Even though you wanted the optional module $desired, you chose not to actually install it: do it yourself manually.";
#        }
#    }
#}

# install required modules listed in 'requires' or 'build_requires' arg to
# new that weren't already installed. Should only be called by prereq_failures
#sub install_required {
#    my ($self, $desired, $version, $msg) = @_;
#
#    $self->log_info(" - ERROR: $msg\n");
#
#    return $self->install_prereq($desired, $version, 1);
#}

# install optional modules listed in 'recommends' arg to new that weren't
# already installed. Should only be called by prereq_failures
#sub install_optional {
#    my ($self, $desired, $version, $msg) = @_;
#
#    unless (defined $self->{ask_optional}) {
#        $self->{ask_optional} = $self->args->{accept}
#                              ? 'n' : $self->prompt("Install [a]ll optional external modules, [n]one, or choose [i]nteractively?", 'n');
#    }
#    return 'skip' if $self->{ask_optional} =~ /^n/i;
#
#    my $install;
#    if ($self->{ask_optional} =~ /^a/i) {
#        $self->log_info(" * $msg\n");
#        $install = 1;
#    }
#    else {
#        $install = $self->y_n(" * $msg\n   Do you want to install it? y/n", 'n');
#    }
#
#    my $orig_version = $version;
#    $version = 0 if $version == -1;
#    if ($install && ! ($self->{ask_optional} =~ /^a/i && $orig_version == -1)) {
#        return $self->install_prereq($desired, $version);
#    }
#    else {
#        my $circular = ($self->{ask_optional} =~ /^a/i && $orig_version == -1) ? " - this is a circular dependency so doesn't get installed when installing 'all' modules. If you really want it, choose modules interactively." : '';
#        $self->log_info(" * You chose not to install $desired$circular\n");
#        return 'ok';
#    }
#}

# there's no official way to discover if being run by CPAN, we take an approach
# similar to that of Module::AutoInstall
#sub under_cpan {
#    my $self = shift;
#
#    unless (defined $self->{under_cpan}) {
#        ## modified from Module::AutoInstall
#
#        my $cpan_env = $ENV{PERL5_CPAN_IS_RUNNING};
#        if ($ENV{PERL5_CPANPLUS_IS_RUNNING}) {
#            $self->{under_cpan} = $cpan_env ? 'CPAN' : 'CPANPLUS';
#        }
#
#        require CPAN;
#
#        unless (defined $self->{under_cpan}) {
#            if ($CPAN::VERSION > '1.89') {
#                if ($cpan_env) {
#                    $self->{under_cpan} = 'CPAN';
#                }
#                else {
#                    $self->{under_cpan} = 0;
#                }
#            }
#        }
#
#        unless (defined $self->{under_cpan}) {
#            # load cpan config
#            if ($CPAN::HandleConfig::VERSION) {
#                # Newer versions of CPAN have a HandleConfig module
#                CPAN::HandleConfig->load;
#            }
#            else {
#                # Older versions had the load method in Config directly
#                CPAN::Config->load;
#            }
#
#            # Find the CPAN lock-file
#            my $lock = File::Spec->catfile($CPAN::Config->{cpan_home}, '.lock');
#            if (-f $lock) {
#                # Module::AutoInstall now goes on to open the lock file and compare
#                # its pid to ours, but we're not in a situation where we expect
#                # the pids to match, so we take the windows approach for all OSes:
#                # find out if we're in cpan_home
#                my $cwd  = File::Spec->canonpath(Cwd::cwd());
#                my $cpan = File::Spec->canonpath($CPAN::Config->{cpan_home});
#
#                $self->{under_cpan} = index($cwd, $cpan) > -1;
#            }
#        }
#
#        if ($self->{under_cpan}) {
#            $self->log_info("(I think I'm being run by CPAN/CPANPLUS, so will rely on it to handle prerequisite installation)\n");
#        }
#        else {
#            $self->log_info("(I think you ran Build.PL directly, so will use CPAN to install prerequisites on demand)\n");
#            $self->{under_cpan} = 0;
#        }
#    }
#
#    return $self->{under_cpan};
#}

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

=head2 find_dist_packages

Like the Module::Build version, except that we always get version from
dist_version
=cut

sub find_dist_packages {
    my $self = shift;

    # Only packages in .pm files are candidates for inclusion here.
    # Only include things in the MANIFEST, not things in developer's
    # private stock.

    my $manifest = $self->_read_manifest('MANIFEST') or die "Can't find dist packages without a MANIFEST file - run 'manifest' action first";

    # Localize
    my %dist_files = map { $self->localize_file_path($_) => $_ } keys %$manifest;

    my @pm_files = grep {exists $dist_files{$_}} keys %{ $self->find_pm_files };

    my $actual_version = $self->dist_version;

    # First, we enumerate all packages & versions,
    # seperating into primary & alternative candidates
    my( %prime, %alt );
    foreach my $file (@pm_files) {
        next if $dist_files{$file} =~ m{^t/};  # Skip things in t/

        my @path = split( /\//, $dist_files{$file} );
        (my $prime_package = join( '::', @path[1..$#path] )) =~ s/\.pm$//;

        my $pm_info = Module::Build::ModuleInfo->new_from_file( $file );

        foreach my $package ( $pm_info->packages_inside ) {
            next if $package eq 'main';  # main can appear numerous times, ignore
            next if grep /^_/, split( /::/, $package ); # private package, ignore

            my $version = $pm_info->version( $package );
            if ($version && $version != $actual_version) {
                $self->log_warn("Package $package had version $version!\n");
            }
            $version = $actual_version;

            if ( $package eq $prime_package ) {
                if ( exists( $prime{$package} ) ) {
                    # M::B::ModuleInfo will handle this conflict
                    die "Unexpected conflict in '$package'; multiple versions found.\n";
                }
                else {
                    $prime{$package}{file} = $dist_files{$file};
                    $prime{$package}{version} = $version if defined( $version );
                }
            }
            else {
                push( @{$alt{$package}}, { file => $dist_files{$file}, version => $version } );
            }
        }
    }

    # Then we iterate over all the packages found above, identifying conflicts
    # and selecting the "best" candidate for recording the file & version
    # for each package.
    foreach my $package ( keys( %alt ) ) {
        my $result = $self->_resolve_module_versions( $alt{$package} );

        if ( exists( $prime{$package} ) ) { # primary package selected
            if ( $result->{err} ) {
                # Use the selected primary package, but there are conflicting
                 # errors amoung multiple alternative packages that need to be
                 # reported
                 $self->log_warn("Found conflicting versions for package '$package'\n" .
                                 "  $prime{$package}{file} ($prime{$package}{version})\n" . $result->{err});
            }
            elsif ( defined( $result->{version} ) ) {
                # There is a primary package selected, and exactly one
                # alternative package

                if ( exists( $prime{$package}{version} ) && defined( $prime{$package}{version} ) ) {
                    # Unless the version of the primary package agrees with the
                    # version of the alternative package, report a conflict
                    if ( $self->compare_versions( $prime{$package}{version}, '!=', $result->{version} ) ) {
                        $self->log_warn("Found conflicting versions for package '$package'\n" .
                                        "  $prime{$package}{file} ($prime{$package}{version})\n" .
                                        "  $result->{file} ($result->{version})\n");
                    }
                }
                else {
                  # The prime package selected has no version so, we choose to
                  # use any alternative package that does have a version
                  $prime{$package}{file}    = $result->{file};
                  $prime{$package}{version} = $result->{version};
                }
            }
            else {
                # no alt package found with a version, but we have a prime
                # package so we use it whether it has a version or not
            }
        }
        else { # No primary package was selected, use the best alternative
            if ( $result->{err} ) {
                $self->log_warn("Found conflicting versions for package '$package'\n" . $result->{err});
            }

            # Despite possible conflicting versions, we choose to record
            # something rather than nothing
            $prime{$package}{file}    = $result->{file};
            $prime{$package}{version} = $result->{version} if defined( $result->{version} );
        }
    }

    # Stringify versions
    for my $key ( grep { exists $prime{$_}->{version} }
                  keys %prime ) {
        $prime{$key}->{version}
            = $prime{$key}->{version}->stringify if ref($prime{$key}->{version});
    }

    return \%prime;
}

# our recommends syntax contains extra info that needs to be ignored at this
# stage
#sub _parse_conditions {
#    my ($self, $spec) = @_;
#
#    ($spec) = split("/", $spec);
#
#    if ($spec =~ /^\s*([\w.]+)\s*$/) { # A plain number, maybe with dots, letters, and underscores
#        return (">= $spec");
#    }
#    else {
#        return split /\s*,\s*/, $spec;
#    }
#}

# when generating META.yml, we output optional_features syntax (instead of
# recommends syntax). Note that as of CPAN v1.9402 nothing useful is done
# with this information, which is why we implement our own request to install
# the optional modules in install_optional().
# Also note that CPAN PLUS complains with an [ERROR] when it sees this META.yml,
# but it isn't fatal and installation continues fine.

# 'recommends' groups broken up now into separate modules and grouping the
# 'requires' instead of lumping modules together (quotes were choking YAML
# parsing). Now passes Parse::CPAN::Meta w/o errors.
# -cjfields 9-17-09

# let us store extra things persistently in _build
#sub _construct {
#    my $self = shift;
#
#    # calling SUPER::_construct will dump some of the input to this sub out
#    # with Data::Dumper, which will complain about code refs. So we replace
#    # any code refs with dummies first, then put them back afterwards
#    my %in_hash = @_;
#    my $auto_features = $in_hash{auto_features} if defined $in_hash{auto_features};
#    my %code_refs;
#    if ($auto_features) {
#        while (my ($key, $hash) = each %{$auto_features}) {
#            while (my ($sub_key, $val) = each %{$hash}) {
#                if (ref($val) && ref($val) eq 'CODE') {
#                    $hash->{$sub_key} = 'CODE_ref';
#                    $code_refs{$key}->{$sub_key} = $val;
#                }
#            }
#        }
#    }
#
#    $self = $self->SUPER::_construct(@_);
#
#    my ($p, $ph) = ($self->{properties}, $self->{phash});
#
#    if (keys %code_refs) {
#        while (my ($key, $hash) = each %{$auto_features}) {
#            if (defined $code_refs{$key}) {
#                while (my ($sub_key, $code_ref) = each %{$code_refs{$key}}) {
#                    $hash->{$sub_key} = $code_ref;
#                }
#                $ph->{auto_features}->{$key} = $hash;
#            }
#        }
#    }
#
#    foreach my $piece (qw(manifest_skip post_install_scripts)) {
#        my $file = File::Spec->catfile($self->config_dir, $piece);
#        $ph->{$piece} = Module::Build::Notes->new(file => $file);
#        $ph->{$piece}->restore if -e $file;
#    }
#
#    return $self;
#}

#sub write_config {
#    my $self = shift;
#    $self->SUPER::write_config;
#
#    # write extra things
#    $self->{phash}{$_}->write() foreach qw(manifest_skip post_install_scripts);
#
#    # be even more certain we can reload ourselves during a resume by copying
#    # ourselves to _build\lib
#    # this is only possible for the core distribution where we are actually
#    # present in the distribution
#    my $self_filename = File::Spec->catfile('Bio', 'Root', 'Build.pm');
#    -e $self_filename || return;
#
#    my $filename = File::Spec->catfile($self->{properties}{config_dir}, 'lib', 'Bio', 'Root', 'Build.pm');
#    my $filedir  = File::Basename::dirname($filename);
#
#    File::Path::mkpath($filedir);
#    warn "Could not create directory '$filedir': $!\n" unless -d $filedir;
#
#    File::Copy::copy($self_filename, $filename);
#    warn "Unable to copy 'Bio/Root/Build.pm' to '$filename'\n" unless -e $filename;
#}

# add a file to the default MANIFEST.SKIP
#sub add_to_manifest_skip {
#    my $self = shift;
#    my %files = map {$self->localize_file_path($_), 1} @_;
#    $self->{phash}{manifest_skip}->write(\%files);
#}

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

# extended to add extra things to the default MANIFEST.SKIP
#sub _write_default_maniskip {
#    my $self = shift;
#    $self->SUPER::_write_default_maniskip;
#
#    my @extra = keys %{$self->{phash}{manifest_skip}->read};
#    if (@extra) {
#        open(my $fh, '>>', 'MANIFEST.SKIP') or die "Could not append MANIFEST.SKIP file\n";
#        print $fh "\n# Avoid additional run-time generated things\n";
#        foreach my $line (@extra) {
#            print $fh $line, "\n";
#        }
#        close($fh);
#    }
#}


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

#sub add_post_install_script {
#    my $self = shift;
#    my %files = map {$self->localize_file_path($_), 1} @_;
#    $self->{phash}{post_install_scripts}->write(\%files);
#}
#
#sub run_post_install_scripts {
#    my $self = shift;
#    my @scripts = keys %{$self->{phash}{post_install_scripts}->read};
#    foreach my $script (@scripts) {
#        $self->run_perl_script($script);
#    }
#}

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

=head2 dist_dir

Nice directory names for dist-related actions
=cut

sub dist_dir {
    my ($self) = @_;
    my $version = $self->dist_version;
    if ($version =~ /^\d\.\d{6}\d$/) {
        # 1.x.x.100 returned as 1.x.x.1
        $version .= '00';
    }
    $version =~ s/00(\d)/$1./g;
    $version =~ s/\.$//;

    if (my ($minor, $rev) = $version =~ /^\d\.(\d)\.\d\.(\d+)$/) {
        my $dev = ! ($minor % 2 == 0);
        if ($rev == 100) {
            my $replace = $dev ? "_$rev" : '';
            $version =~ s/\.\d+$/$replace/;
        }
        elsif ($rev < 100) {
            $rev = sprintf("%03d", $rev);
            $version =~ s/\.\d+$/_$rev-RC/;
        }
        else {
            $rev -= 100 unless $dev;
            my $replace = $dev ? "_$rev" : ".$rev";
            $version =~ s/\.\d+$/$replace/;
        }
    }

    return "$self->{properties}{dist_name}-$version";
}

# try to be as consistent as possible with Module::Build API
#sub ppm_name {
#    my $self = shift;
#    return $self->dist_dir.'-ppm';
#}

# generate complete ppd4 version file
#sub ACTION_ppd {
#    my $self = shift;
#
#    my $file = $self->make_ppd(%{$self->{args}});
#    $self->add_to_cleanup($file);
#    #$self->add_to_manifest_skip($file);
#}

# add pod2htm temp files to MANIFEST.SKIP, generated during ppmdist most likely
#sub htmlify_pods {
#    my $self = shift;
#    $self->SUPER::htmlify_pods(@_);
#    #$self->add_to_manifest_skip('pod2htm*');
#}

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

# overridden from Module::Build::PPMMaker for ppd4 compatability

# note: no longer needed with more recent versions of Module::Build

#sub make_ppd {
#    my ($self, %args) = @_;
#
#    require Module::Build::PPMMaker;
#    my $mbp = Module::Build::PPMMaker->new();
#
#    my %dist;
#    foreach my $info (qw(name author abstract version)) {
#        my $method = "dist_$info";
#        $dist{$info} = $self->$method() or die "Can't determine distribution's $info\n";
#    }
#    $dist{codebase} = $self->ppm_name.'.tar.gz';
#    $mbp->_simple_xml_escape($_) foreach $dist{abstract}, $dist{codebase}, @{$dist{author}};
#
#    my (undef, undef, undef, $mday, $mon, $year) = localtime();
#    $year += 1900;
#    $mon++;
#    my $date = "$year-$mon-$mday";
#
#    my $softpkg_version = $self->dist_dir;
#    $softpkg_version =~ s/^$dist{name}-//;
#
#    # to avoid a ppm bug, instead of including the requires in the softpackage
#    # for the distribution we're making, we'll make a seperate Bundle::
#    # softpackage that contains all the requires, and require only the Bundle in
#    # the real softpackage
#    my ($bundle_name) = $dist{name} =~ /^.+-(.+)/;
#    $bundle_name ||= 'core';
#    $bundle_name =~ s/^(\w)/\U$1/;
#    my $bundle_dir = "Bundle-BioPerl-$bundle_name-$softpkg_version-ppm";
#    my $bundle_file = "$bundle_dir.tar.gz";
#    my $bundle_softpkg_name = "Bundle-BioPerl-$bundle_name";
#    $bundle_name = "Bundle::BioPerl::$bundle_name";
#
#    # header
#    my $ppd = <<"PPD";
#    <SOFTPKG NAME=\"$dist{name}\" VERSION=\"$softpkg_version\" DATE=\"$date\">
#        <TITLE>$dist{name}</TITLE>
#        <ABSTRACT>$dist{abstract}</ABSTRACT>
#@{[ join "\n", map "        <AUTHOR>$_</AUTHOR>", @{$dist{author}} ]}
#        <PROVIDE NAME=\"$dist{name}::\" VERSION=\"$dist{version}\"/>
#PPD
#
#    # provide section
#    foreach my $pm (@{$self->rscan_dir('Bio', qr/\.pm$/)}) {
#        # convert these filepaths to Module names
#        $pm =~ s/\//::/g;
#        $pm =~ s/\.pm//;
#
#        $ppd .= sprintf(<<'EOF', $pm, $dist{version});
#        <PROVIDE NAME="%s" VERSION="%s"/>
#EOF
#    }
#
#    # rest of softpkg
#    $ppd .= <<"PPD";
#        <IMPLEMENTATION>
#            <ARCHITECTURE NAME=\"MSWin32-x86-multi-thread-5.8\"/>
#            <CODEBASE HREF=\"$dist{codebase}\"/>
#            <REQUIRE NAME=\"$bundle_name\" VERSION=\"$dist{version}\"/>
#        </IMPLEMENTATION>
#    </SOFTPKG>
#PPD
#
#    # now a new softpkg for the bundle
#    $ppd .= <<"PPD";
#
#    <SOFTPKG NAME=\"$bundle_softpkg_name\" VERSION=\"$softpkg_version\" DATE=\"$date\">
#        <TITLE>$bundle_name</TITLE>
#        <ABSTRACT>Bundle of pre-requisites for $dist{name}</ABSTRACT>
#@{[ join "\n", map "        <AUTHOR>$_</AUTHOR>", @{$dist{author}} ]}
#        <PROVIDE NAME=\"$bundle_name\" VERSION=\"$dist{version}\"/>
#        <IMPLEMENTATION>
#            <ARCHITECTURE NAME=\"MSWin32-x86-multi-thread-5.8\"/>
#            <CODEBASE HREF=\"$bundle_file\"/>
#PPD
#
#    # required section
#    # we do both requires and recommends to make installation on Windows as
#    # easy (mindless) as possible
#    for my $type ('requires', 'recommends') {
#        my $prereq = $self->$type;
#        while (my ($modname, $version) = each %$prereq) {
#            next if $modname eq 'perl';
#            ($version) = split("/", $version) if $version =~ /\//;
#
#            # Module names must have at least one ::
#            unless ($modname =~ /::/) {
#                $modname .= '::';
#            }
#
#            # Bio::Root::Version number comes out as triplet number like 1.5.2;
#            # convert to our own version
#            if ($modname eq 'Bio::Root::Version') {
#                $version = $dist{version};
#            }
#
#            $ppd .= sprintf(<<'EOF', $modname, $version || '');
#            <REQUIRE NAME="%s" VERSION="%s"/>
#EOF
#        }
#    }
#
#    # footer
#    $ppd .= <<'EOF';
#        </IMPLEMENTATION>
#    </SOFTPKG>
#EOF
#
#    my $ppd_file = "$dist{name}.ppd";
#    my $fh = IO::File->new(">$ppd_file") or die "Cannot write to $ppd_file: $!";
#    print $fh $ppd;
#    close $fh;
#
#    $self->delete_filetree($bundle_dir);
#    mkdir($bundle_dir) or die "Cannot create '$bundle_dir': $!";
#    $self->make_tarball($bundle_dir);
#    $self->delete_filetree($bundle_dir);
#    $self->add_to_cleanup($bundle_file);
#    #$self->add_to_manifest_skip($bundle_file);
#
#    return $ppd_file;
#}

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
