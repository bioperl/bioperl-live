#!/usr/bin/perl -w

# This is a subclass of Module::Build so we can override certain methods and do
# fancy stuff

# It was first written against Module::Build::Base v0.2805. Many of the methods
# here are copy/pasted from there in their entirety just to change one or two
# minor things, since for the most part Module::Build::Base code is hard to
# cleanly override.

package ModuleBuildBioperl;

BEGIN {
    # we really need Module::Build to be installed
    unless (eval "use Module::Build; 1") {
        print "This package requires Module::Build to install itself.\n";
        
        require ExtUtils::MakeMaker;
        my $yn = ExtUtils::MakeMaker::prompt('  Install Module::Build now from CPAN?', 'y');
        
        unless ($yn =~ /^y/i) {
            die " *** Cannot install without Module::Build.  Exiting ...\n";
        }
        
        require Cwd;
        require File::Spec;
        require File::Copy;
        require CPAN;
        
        # Save this because CPAN will chdir all over the place.
        my $cwd = Cwd::cwd();
        
        my $build_pl = File::Spec->catfile($cwd, "Build.PL");
        
        File::Copy::move($build_pl, $build_pl."hidden"); # avoid bizarre bug with Module::Build tests using the wrong Build.PL if it happens to be in PERL5LIB
        CPAN::Shell->install('Module::Build');
        File::Copy::move($build_pl."hidden", $build_pl);
        CPAN::Shell->expand("Module", "Module::Build")->uptodate or die "Couldn't install Module::Build, giving up.\n";
        
        chdir $cwd or die "Cannot chdir() back to $cwd: $!";
    }
    
    eval "use base Module::Build; 1" or die $@;
    
    # ensure we'll be able to reload this module later by adding its path to inc
    use Cwd;
    use lib Cwd::cwd();
}

use strict;
use warnings;

our $VERSION = 1.005002005;
our @extra_types = qw(options test excludes_os);
our $checking_types = "requires|conflicts|".join("|", @extra_types);


# our modules are in Bio, not lib
sub find_pm_files {
    my $self = shift;
    foreach my $pm (@{$self->rscan_dir('Bio', qr/\.pm$/)}) {
        $self->{properties}{pm_files}->{$pm} = File::Spec->catfile('lib', $pm);;
    }
    
    $self->_find_file_by_type('pm', 'lib');
}

# ask what scripts to install (this method is unique to bioperl)
sub choose_scripts {
    my $self = shift;
    
    # we can offer interactive installation by groups only if we have subdirs
    # in scripts and no .PLS files there
    opendir(my $scripts_dir, 'scripts') or die "Can't open directory 'scripts': $!\n";
    my $int_ok = 0;
    my @group_dirs;
    while (my $thing = readdir($scripts_dir)) {
        next if $thing =~ /^\./;
        next if $thing eq 'CVS';
        if ($thing =~ /PLS$|pl$/) {
            $int_ok = 0;
            last;
        }
        $thing = File::Spec->catfile('scripts', $thing);
        if (-d $thing) {
            $int_ok = 1;
            push(@group_dirs, $thing);
        }
    }
    closedir($scripts_dir);
    my $question = $int_ok ? "Install [a]ll Bioperl scripts, [n]one, or choose groups [i]nteractively?" : "Install [a]ll Bioperl scripts or [n]one?";
    
    my $prompt = $self->prompt($question, 'a');
    
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

# our version of script_files doesn't take args but just installs those scripts
# requested by the user after choose_scripts() is called. If it wasn't called,
# installs all scripts in scripts directory
sub script_files {
    my $self = shift;
    
    my $chosen_scripts = $self->notes('chosen_scripts');
    if ($chosen_scripts) {
        return if $chosen_scripts eq 'none';
        return { map {$_, 1} split(/\|/, $chosen_scripts) } unless $chosen_scripts eq 'all';
    }
    
    return $_ = { map {$_,1} @{$self->rscan_dir('scripts', qr/\.PLS$|\.pl$/)} };
}

# process scripts normally, except that we change name from *.PLS to bp_*.pl
sub process_script_files {
    my $self = shift;
    my $files = $self->find_script_files;
    return unless keys %$files;
  
    my $script_dir = File::Spec->catdir($self->blib, 'script');
    File::Path::mkpath( $script_dir );
    
    foreach my $file (keys %$files) {
        my $result = $self->copy_if_modified($file, $script_dir, 'flatten') or next;
        $self->fix_shebang_line($result) unless $self->os_type eq 'VMS';
        $self->make_executable($result);
        
        my $final = File::Basename::basename($result);
        $final =~ s/\.PLS$/\.pl/;                  # change from .PLS to .pl
        $final =~ s/^/bp_/ unless $final =~ /^bp/; # add the "bp" prefix
        $final = File::Spec->catfile($script_dir, $final);
        $self->log_info("$result -> $final\n");
        File::Copy::move($result, $final) or die "Can't rename '$result' to '$final': $!";
    }
}

# extended to handle extra checking types
sub features {
    my $self = shift;
    my $ph = $self->{phash};
    
    if (@_) {
        my $key = shift;
        if ($ph->{features}->exists($key)) {
            return $ph->{features}->access($key, @_);
        }
        
        if (my $info = $ph->{auto_features}->access($key)) {
            my $failures = $self->prereq_failures($info);
            my $disabled = grep( /^(?:\w+_)?(?:$checking_types)$/, keys %$failures ) ? 1 : 0;
            return !$disabled;
        }
        
        return $ph->{features}->access($key, @_);
    }
  
    # No args - get the auto_features & overlay the regular features
    my %features;
    my %auto_features = $ph->{auto_features}->access();
    while (my ($name, $info) = each %auto_features) {
        my $failures = $self->prereq_failures($info);
        my $disabled = grep( /^(?:\w+_)?(?:$checking_types)$/, keys %$failures ) ? 1 : 0;
        $features{$name} = $disabled ? 0 : 1;
    }
    %features = (%features, $ph->{features}->access());
  
    return wantarray ? %features : \%features;
}
*feature = \&features;

# overridden to fix a stupid bug in Module::Build and extended to handle extra
# checking types
sub check_autofeatures {
    my ($self) = @_;
    my $features = $self->auto_features;
    
    return unless %$features;
    
    $self->log_info("Checking features:\n");
    
    my $max_name_len = 0; # this wasn't set to 0 in Module::Build, causing warning in next line
    $max_name_len = ( length($_) > $max_name_len ) ? length($_) : $max_name_len for keys %$features;
    
    while (my ($name, $info) = each %$features) {
        $self->log_info("  $name" . '.' x ($max_name_len - length($name) + 4));
        if ($name eq 'PL_files') {
            print "got $name => $info\n";
            print "info has:\n";
            while (my ($key, $val) = each %$info) {
                print "  $key => $val\n";
            }
        }
        
        if ( my $failures = $self->prereq_failures($info) ) {
            my $disabled = grep( /^(?:\w+_)?(?:$checking_types)$/, keys %$failures ) ? 1 : 0;
            $self->log_info( $disabled ? "disabled\n" : "enabled\n" );
            
            my $log_text;
            while (my ($type, $prereqs) = each %$failures) {
                while (my ($module, $status) = each %$prereqs) {
                    my $required = ($type =~ /^(?:\w+_)?(?:requires|conflicts)$/) ? 1 : 0;
                    my $prefix = ($required) ? '-' : '*';
                    $log_text .= "    $prefix $status->{message}\n";
                }
            }
            $self->log_warn($log_text) if $log_text && ! $self->quiet;
        }
        else {
            $self->log_info("enabled\n");
        }
    }
    
    $self->log_info("\n");
}

# extend to handle option checking (which takes an array ref) and code test
# checking (which takes a code ref and must return a message only on failure)
# and excludes_os (which takes an array ref of regexps).
# also handles more informative output of recommends section
sub prereq_failures {
    my ($self, $info) = @_;
    
    my @types = (@{ $self->prereq_action_types }, @extra_types);
    $info ||= {map {$_, $self->$_()} @types};
    
    my $out = {};
    foreach my $type (@types) {
        my $prereqs = $info->{$type} || next;
        
        my $status = {};
        if ($type eq 'test') {
            unless (keys %$out) {
                $status->{message} = &{$prereqs};
                $out->{$type}{'test'} = $status if $status->{message};
            }
        }
        elsif ($type eq 'options') {
            my @not_ok;
            foreach my $wanted_option (@{$prereqs}) {
                unless ($self->args($wanted_option)) {
                    push(@not_ok, $wanted_option);
                }
            }
            
            if (@not_ok > 0) {
                $status->{message} = "Command line option(s) '@not_ok' not supplied";
                $out->{$type}{'options'} = $status;
            }
        }
        elsif ($type eq 'excludes_os') {
            foreach my $os (@{$prereqs}) {
                if ($^O =~ /$os/i) {
                    $status->{message} = "This feature isn't supported under your OS ($os)";
                    $out->{$type}{'excludes_os'} = $status;
                    last;
                }
            }
        }
        else {
            while ( my ($modname, $spec) = each %$prereqs ) {
                $status = $self->check_installed_status($modname, $spec);
                
                if ($type =~ /^(?:\w+_)?conflicts$/) {
                    next if !$status->{ok};
                    $status->{conflicts} = delete $status->{need};
                    $status->{message} = "$modname ($status->{have}) conflicts with this distribution";
                }
                elsif ($type =~ /^(?:\w+_)?recommends$/) {
                    next if $status->{ok};
                    
                    my ($preferred_version, $why, $by_what) = split("/", $spec);
                    $by_what = join(", ", split(",", $by_what));
                    $by_what =~ s/, (\S+)$/ and $1/;
                    
                    $status->{message} = (!ref($status->{have}) && $status->{have} eq '<none>'
                                  ? "Optional prerequisite $modname is not installed"
                                  : "$modname ($status->{have}) is installed, but we prefer to have $preferred_version");
                    
                    $status->{message} .= "\n   (wanted for $why, used by $by_what)";
                    
                    my $installed = $self->install_optional($modname, $preferred_version, $status->{message});
                    next if $installed eq 'ok';
                    $status->{message} = $installed unless $installed eq 'skip';
                }
                else {
                    next if $status->{ok};
                }
                
                $out->{$type}{$modname} = $status;
            }
        }
    }
    
    return keys %{$out} ? $out : return;
}

# install optional modules listed in 'recommends' arg to new that weren't
# already installed. Should only be called by prereq_failures
sub install_optional {
    my ($self, $desired, $version, $msg) = @_;
    
    unless (defined $self->{ask_optional}) {
        $self->{ask_optional} = $self->prompt("Install [a]ll optional external modules, [n]one, or choose [i]nteractively?", 'n');
    }
    return 'skip' if $self->{ask_optional} =~ /^n/i;
    
    my $install;
    if ($self->{ask_optional} =~ /^a/i) {
        $self->log_info(" * $msg\n");
        $install = 1;
    }
    else {
        $install = $self->y_n(" * $msg\n   Do you want to install it? y/n", 'n');
    }
    
    if ($install) {
        # Here we use CPAN to actually install the desired module, the benefit
        # being we continue even if installation fails, and that this works
        # even when not using CPAN to install.
        #
        # The alternative would be to simply append the module to 'requires'
        # and let CPAN deal with required modules in the normal way, but
        # older CPANs don't do that (look only at META.yml), and we get
        # total failure for something that was only optional
        require Cwd;
        require CPAN;
        
        # Save this because CPAN will chdir all over the place.
        my $cwd = Cwd::cwd();
        
        CPAN::Shell->install($desired);
        my $msg;
        if (CPAN::Shell->expand("Module", $desired)->uptodate) {
            $self->log_info("\n\n*** (back in Bioperl Build.PL) ***\n * You chose to install $desired and it installed fine\n");
            $msg = 'ok';
        }
        else {
            $self->log_info("\n\n*** (back in Bioperl Build.PL) ***\n");
            $msg = "You chose to install $desired but it failed to install";
        }
        
        chdir $cwd or die "Cannot chdir() back to $cwd: $!";
        return $msg;
    }
    else {
        $self->log_info(" * You chose not to install $desired\n");
        return 'ok';
    }
}

# overridden simply to not print the default answer if chosen by hitting return
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

# like the Module::Build version, except that we always get version from
# dist_version
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
    for (grep exists $_->{version}, values %prime) {
        $_->{version} = $_->{version}->stringify if ref($_->{version});
    }
  
    return \%prime;
}

# our recommends syntax contains extra info that needs to be ignored at this
# stage
sub _parse_conditions {
    my ($self, $spec) = @_;
    
    ($spec) = split("/", $spec);
    
    if ($spec =~ /^\s*([\w.]+)\s*$/) { # A plain number, maybe with dots, letters, and underscores
        return (">= $spec");
    }
    else {
        return split /\s*,\s*/, $spec;
    }
}

# when generating META.yml, we output optional_features syntax (instead of
# recommends syntax). Note that as of CPAN v1.8802 nothing useful is done
# with this information, which is why we implement our own request to install
# the optional modules in install_optional()
sub prepare_metadata {
    my ($self, $node, $keys) = @_;
    my $p = $self->{properties};
    
    # A little helper sub
    my $add_node = sub {
        my ($name, $val) = @_;
        $node->{$name} = $val;
        push @$keys, $name if $keys;
    };
    
    foreach (qw(dist_name dist_version dist_author dist_abstract license)) {
        (my $name = $_) =~ s/^dist_//;
        $add_node->($name, $self->$_());
        die "ERROR: Missing required field '$_' for META.yml\n" unless defined($node->{$name}) && length($node->{$name});
    }
    $node->{version} = '' . $node->{version}; # Stringify version objects
    
    if (defined( $self->license ) && defined( my $url = $self->valid_licenses->{ $self->license } )) {
        $node->{resources}{license} = $url;
    }
    
    foreach ( @{$self->prereq_action_types} ) {
        if (exists $p->{$_} and keys %{ $p->{$_} }) {
            if ($_ eq 'recommends') {
                my $hash;
                while (my ($req, $val) = each %{ $p->{$_} }) {
                    my ($ver, $why, $used_by) = split("/", $val);
                    my $info = {};
                    $info->{description} = $why;
                    $info->{requires} = { $req => $ver };
                    $hash->{$used_by} = $info;
                }
                $add_node->('optional_features', $hash);
            }
            else {
                $add_node->($_, $p->{$_});
            }
        }
    }
    
    if (exists $p->{dynamic_config}) {
        $add_node->('dynamic_config', $p->{dynamic_config});
    }
    my $pkgs = eval { $self->find_dist_packages };
    if ($@) {
        $self->log_warn("$@\nWARNING: Possible missing or corrupt 'MANIFEST' file.\n" . "Nothing to enter for 'provides' field in META.yml\n");
    }
    else {
        $node->{provides} = $pkgs if %$pkgs;
    };
    
    if (exists $p->{no_index}) {
        $add_node->('no_index', $p->{no_index});
    }
    
    $add_node->('generated_by', "Module::Build version $Module::Build::VERSION");
    
    $add_node->('meta-spec', 
            {version => '1.2',
             url     => 'http://module-build.sourceforge.net/META-spec-v1.2.html',
            });
    
    while (my($k, $v) = each %{$self->meta_add}) {
        $add_node->($k, $v);
    }
    
    while (my($k, $v) = each %{$self->meta_merge}) {
        $self->_hash_merge($node, $k, $v);
    }
    
    return $node;
}

# let us store extra things persistently in _build
sub _construct {
    my $self = shift;
    $self = $self->SUPER::_construct(@_);
    
    my ($p, $ph) = ($self->{properties}, $self->{phash});
    
    foreach (qw(manifest_skip post_install_scripts)) {
        my $file = File::Spec->catfile($self->config_dir, $_);
        $ph->{$_} = Module::Build::Notes->new(file => $file);
        $ph->{$_}->restore if -e $file;
    }
    
    return $self;
}
sub write_config {
    my $self = shift;
    $self->SUPER::write_config;
    $self->{phash}{$_}->write() foreach qw(manifest_skip post_install_scripts);
}

# add a file to the default MANIFEST.SKIP
sub add_to_manifest_skip {
    my $self = shift;
    my %files = map {$self->localize_file_path($_), 1} @_;
    $self->{phash}{manifest_skip}->write(\%files);
}

# we always generate a new MANIFEST and MANIFEST.SKIP here, instead of allowing
# existing files to remain
sub ACTION_manifest {
    my ($self) = @_;
    
    my $maniskip = 'MANIFEST.SKIP';
    if ( -e 'MANIFEST' || -e $maniskip ) {
        $self->log_warn("MANIFEST files already exist, will overwrite them\n");
        unlink('MANIFEST');
        unlink($maniskip);
    }
    $self->_write_default_maniskip($maniskip);
    
    require ExtUtils::Manifest;  # ExtUtils::Manifest is not warnings clean.
    local ($^W, $ExtUtils::Manifest::Quiet) = (0,1);
    ExtUtils::Manifest::mkmanifest();
}

# extended to add extra things to the default MANIFEST.SKIP
sub _write_default_maniskip {
    my $self = shift;
    $self->SUPER::_write_default_maniskip;
    
    my @extra = keys %{$self->{phash}{manifest_skip}->read};
    if (@extra) {
        open(my $fh, '>>', 'MANIFEST.SKIP') or die "Could not open MANIFEST.SKIP file\n";
        print $fh "\n# Avoid additional run-time generated things\n";
        foreach my $line (@extra) {
            print $fh $line, "\n";
        }
        close($fh);
    }
}

# extended to run scripts post-installation
sub ACTION_install {
  my ($self) = @_;
  require ExtUtils::Install;
  $self->depends_on('build');
  ExtUtils::Install::install($self->install_map, !$self->quiet, 0, $self->{args}{uninst}||0);
  $self->run_post_install_scripts;
}
sub add_post_install_script {
    my $self = shift;
    my %files = map {$self->localize_file_path($_), 1} @_;
    $self->{phash}{post_install_scripts}->write(\%files);
}
sub run_post_install_scripts {
    my $self = shift;
    my @scripts = keys %{$self->{phash}{post_install_scripts}->read};
    foreach my $script (@scripts) {
        $self->run_perl_script($script);
    }
}

# for use with auto_features, which needs to require LWP::UserAgent as one of
# its reqs
sub test_internet {
    eval {require LWP::UserAgent;}; # if not installed, this sub won't actually be called
    my $ua = LWP::UserAgent->new;
    $ua->timeout(10);
    $ua->env_proxy;
    my $response = $ua->get('http://search.cpan.org/');
    unless ($response->is_success) {
        return "Could not connect to the internet (http://search.cpan.org/)";
    }
    return;
}

# nice directory names for dist-related actions
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
sub ppm_name {
    my $self = shift;
    return $self->dist_dir.'-ppm';
}

# we make all archive formats we want, not just .tar.gz
# we also auto-run manifest action, since we always want to re-create
# MANIFEST and MANIFEST.SKIP just-in-time
sub ACTION_dist {
    my ($self) = @_;
    
    $self->depends_on('manifest');
    $self->depends_on('distdir');
    
    my $dist_dir = $self->dist_dir;
    
    $self->make_zip($dist_dir);
    $self->make_tarball($dist_dir);
    $self->delete_filetree($dist_dir);
}

# makes zip file for windows users and bzip2 files as well
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

1;
