#!/usr/bin/perl

=head1 NAME

bp_netinstall.pl

=head1 SYNOPSIS

  bp_netinstall.pl -b|--build_param_str BUILD_STRING [options]

options: 

 -h|--help                Show this message 
 -d|--dev                 Use the development version of bioperl from git 
 --build_param_str=<args> Parameters that are passed in at 'perl Build.PL'
 --install_param_str=<args>
                          Use this string to predefine './Build install' 
                            parameters such as 'install_base' for
                            bioperl installation
 --bioperl_path           Path to BioPerl tarball (will not download BioPerl)
 --skip_start             Don't wait for 'Enter' at program start

=head1 DESCRIPTION

Net-based installer of BioPerl; this is based on the GBrowse netinstaller
and hopefully all references to GBrowse have been removed.  Let me know if not.

Save this to disk as "bp_netinstall.pl" and run:

   [sudo] perl bp_netinstall.pl

=head1 AUTHOR

Scott Cain scain@cpan.org

=head1 COPYRIGHT

2010. This script may be distributed under the same license as perl.

=cut



# Universal Net-based installer
# Save this to disk as "bp_netinstall.pl" and run:
#   perl bp_netinstall.pl

use warnings;
use strict;
use CPAN '!get';
use Config;
use Getopt::Long;
use Pod::Usage;
use File::Copy qw( cp move );
use File::Temp qw(tempdir);
use LWP::Simple;
use Cwd;

use constant NMAKE => 'http://download.microsoft.com/download/vc15/patch/1.52/w95/en-us/nmake15.exe';

my ( $show_help, $get_from_cvs, $build_param_string, $working_dir,
     $get_bioperl_svn, $is_cygwin, $windows,
     $binaries, $make, $tmpdir, $bioperl_path,
     $skip_start, $install_param_string, $perl_path);

BEGIN {

  GetOptions(
        'h|help'              => \$show_help,             # Show help and exit
        'd|dev'               => \$get_from_cvs,          # Use the dev svn
        'build_param_str=s'   => \$build_param_string,    # Build parameters
        'bioperl_dev'         => \$get_bioperl_svn,
        'bioperl_path=s'      => \$bioperl_path,
        'install_param_str=s' => \$install_param_string,
        'skip_start'          => \$skip_start,
        )
        or pod2usage(2);
  pod2usage(2) if $show_help;

  $perl_path = $Config{perlpath}; 

  print STDERR "\nAbout to install BioPerl and all its prerequisites.\n";
  print STDERR "\nYou will be asked various questions during this process. You can almost always";
  print STDERR "\naccept the default answer.\n";
  print STDERR "The whole process will take several minutes and will generate lots of messages.\n";
  print STDERR "\nPress return when you are ready to start!\n";
  my $h = <> unless $skip_start;
  print STDERR "*** Installing Perl files needed for a net-based install ***\n";

  $windows = $Config{osname} =~ /mswin/i;


# MAY not be necessary--we'll have to see.
#  if ($windows and $] == 5.010) {
#     print STDERR "\n\nActiveState Perl 5.10 is not compatible with GBrowse due to problems\n";
#     print STDERR "with the AS implementation.  Please remove it and install Perl 5.8 instead.\n\n\n";
#     exit(0);
#  }

# Also MAY not be necessary
#  if ($windows) {
#     print STDERR "\n\nInstalling Win32 perl module\n\n";
#     system("ppm install Win32");
#  }

  eval "CPAN::Config->load";
  eval "CPAN::Config->commit";

  $working_dir = getcwd;

  $tmpdir = tempdir(CLEANUP=>1) 
    or die "Could not create temporary directory: $!";

  $binaries = $Config{'binexp'};
  $make     = $Config{'make'};

  if ($windows) {
    system("ppm install YAML");
  }
  else {
    CPAN::Shell->install('YAML');
  }
  CPAN::Shell->install('Archive::Zip');
  CPAN::Shell->install('HTML::Tagset');
  CPAN::Shell->install('LWP::Simple');
  eval "use Archive::Zip ':ERROR_CODES',':CONSTANTS'";

  if ($windows && !-e "$binaries/${make}.exe") {

    print STDERR "Installing make utility...\n";
    -w $binaries or die "$binaries directory is not writeable. Please re-login as Admin.\n";
    chdir $tmpdir;

    my $rc = mirror(NMAKE,"nmake.zip");
    die "Could not download nmake executable from Microsoft web site."
      unless $rc == RC_OK() or $rc == RC_NOT_MODIFIED();

    my $zip = Archive::Zip->new('nmake.zip') or die "Couldn't open nmake zip file for decompression: $!";
    $zip->extractTree == AZ_OK() or die "Couldn't unzip file: $!";
    -e 'NMAKE.EXE' or die "Couldn't extract nmake.exe";

    cp('NMAKE.EXE',"$binaries/${make}.EXE") or die "Couldn't install nmake.exe: $!";
    cp('NMAKE.ERR',"$binaries/${make}.ERR"); # or die "Couldn't install nmake.err: $!"; # not fatal
  }

  CPAN::Shell->install('Archive::Tar');
  #print STDERR $@;
  #print STDERR "at end of BEGIN{}\n";
  1;
};

#print STDERR "here i am\n";
#print STDERR $@;

use Archive::Tar;
#use CPAN '!get';

$is_cygwin = 1 if ( $^O eq 'cygwin' );

if ($get_from_cvs) {
    $get_bioperl_svn = 1;
}

#if ($wincvs or ($windows and $get_from_cvs)) {
#    die "\n\nGBrowse is now in svn and fetching from svn on Windows\nis not currently supported\n ";
#}

#if ($windows and !$wincvs and $get_gbrowse_cvs ) {
#    die "\n\nThe development/cvs tags are not supported on Windows when\n"
#        ."WinCVS is not installed; exiting...\n";
#}

$build_param_string ||="";
$install_param_string ||="";

use constant BIOPERL_VERSION      => 'BioPerl-1.6.1';
use constant BIOPERL_REQUIRES     => '1.006001';  # sorry for the redundancy
use constant BIOPERL_LIVE_URL     => 'http://github.com/bioperl/bioperl-live/tarball/master';
use constant BIOPERL              => 'http://bioperl.org/DIST/'.BIOPERL_VERSION.'.tar.gz';

my %REPOSITORIES = (
                    #'BioPerl-Release-Candidates' => 'http://bioperl.org/DIST/RC',
		    'BioPerl-Regular-Releases'   => 'http://bioperl.org/DIST',
	            'Kobes'                      => 'http://theoryx5.uwinnipeg.ca/ppms',
                    'Bribes'                     => 'http://www.Bribes.org/perl/ppm',
                     'tcool'                     => 'http://ppm.tcool.org/archives/',
                    );


# this is so that ppm can be called in a pipe
$ENV{COLUMNS} = 80; # why do we have to do this?
$ENV{LINES}   = 24;

setup_ppm() if $windows;

unless ( eval "use GD 2.31; 1" ) {
   if ($windows) {
     print STDERR "Installing GD via ppm.\n";
     print STDERR "(This may take a while...\n";
     system("ppm install GD");
  }
  else {
     print STDERR "Installing GD via CPAN...\n";
     CPAN::Shell->install('GD') unless eval "use GD 2.31; 1";
  }
}

print STDERR "\n*** Installing prerequisites for BioPerl ***\n";

if ($windows and !eval "use DB_File; 1") {
  print STDERR "Installing DB_File for BioPerl.\n";

  # GBrowse doesn't like DB_File 1.820, so we explicitly get DB_File by url
  system("ppm install http://ppm.tcool.org/archives/DB_File.ppd");
}
#system("ppm install SVG") if $windows;
#CPAN::Shell->install('GD::SVG');

#needed?
CPAN::Shell->install('IO::String');
#CPAN::Shell->install('Text::Shellwords');
#if ($windows) {
#    #CGI::Session and Digest::MD5 both fail to install via cpan on windows
#    system("ppm install CGI-Session");
#    system("ppm install Digest-MD5");
#}
#else {
#    CPAN::Shell->install('CGI::Session');
#    CPAN::Shell->install('Digest::MD5');
#}
#CPAN::Shell->install('File::Temp');
#CPAN::Shell->install('Class::Base');
#CPAN::Shell->install('Statistics::Descriptive');
#CPAN::Shell->install('Data::Stag');

my $version = BIOPERL_REQUIRES;
if (!(eval "use Bio::Perl $version; 1") or $get_bioperl_svn or $bioperl_path) {
    print STDERR "\n*** Installing BioPerl ***\n";

    #would like to use ppm, but ppm won't install 1.6
    #if ($windows and !$get_bioperl_svn and !$bioperl_path) {
    #  my $bioperl_index = find_bioperl_ppm();
    #  system("ppm install --force $bioperl_index");
    #} else {
        # recent versions of Module::Build fail to install without force!
        CPAN::Shell->force('Module::Build') unless eval "require Module::Build; 1";
        do_install(BIOPERL,
                   'bioperl.tgz',
                   BIOPERL_VERSION,
                   'Build',
                   $get_bioperl_svn ? 'svn' : '',
                   $build_param_string,
                   $bioperl_path,
                   $install_param_string,
                   $perl_path);
    #}
}
else {
    print STDERR "BioPerl is up to date.\n";
}

print STDERR "\n *** Installing Bio::Graphics ***\n";


#install biographics?
CPAN::Shell->install('Bio::Graphics');

exit 0;

END {
  my $null = ($^O =~ m/mswin/i) ? 'NUL' : '/dev/null';
  open STDERR,">$null"; # windows has an annoying message when cleaning up temp file
}

sub do_install {
  my ($download,$local_name,$distribution,$method,
         $from_cvs,$build_param_string,$file_path,$install_param_string,
         $perl_path) = @_;

  $install_param_string ||= '';
  chdir $tmpdir;

  do_get_distro($download,$local_name,$distribution,$from_cvs,$file_path);

  my $build_str = $windows ? "Build" : "./Build";

  if ($method eq 'make') {
      system("$perl_path Makefile.PL $build_param_string") == 0
            or die "Couldn't run perl Makefile.PL command\n";
      system("$make install UNINST=1 $install_param_string")    == 0 ;
  }
  elsif ($method eq 'Build') {
      system("$perl_path $build_str.PL --yes=1 $build_param_string")   == 0
            or die "Couldn't run perl Build.PL command\n";
      system("$build_str install --uninst 1 $install_param_string") == 0;
  }
}

sub do_get_distro {
    my ($download,$local_name,$distribution,$distribution_method,$file_path) = @_;

    if ($file_path) {
        chdir $working_dir;
        if (-e $file_path) { #must be an absolute path
            cp($file_path, "$tmpdir/$local_name");
        }
        elsif (-e "$working_dir/$file_path") { #assume it's a rel path from the original directory
            cp("$working_dir/$file_path", "$tmpdir/$local_name");
        }
        else {
            print "Couldn't find $file_path; nothing to do so quitting...\n";
            exit(-1);
        }
        $distribution = ($local_name =~ /gbrowse/)
                      ? "Generic-Genome-Browser" : "bioperl-live"; 
        chdir $tmpdir;
        extract_tarball($local_name,$distribution);
    }
    elsif ($distribution_method) {
        my $distribution_dir;
        print STDERR "Downloading bioperl-live...\n";
        $distribution_dir = 'bioperl-live';

        my $filename = 'bioperl-live.tar.gz'; # =determine_filename();
        my $url = BIOPERL_LIVE_URL."/$filename";

        my $rc = mirror($url, $filename); 
        unless ($rc == RC_OK or $rc == RC_NOT_MODIFIED){
            print STDERR "Failed to get nightly bioperl-live file: $rc\n";
            return undef;
        }
        extract_tarball($filename,$distribution_dir);
        return 1;
        chdir $distribution_dir
            or die "Couldn't enter $distribution_dir directory: $@";
    }
    else {
        print STDERR "Downloading $download...\n";
        my $rc = mirror($download,$local_name);
        die "Could not download $distribution distribution from $download."
            unless $rc == RC_OK or $rc == RC_NOT_MODIFIED;

        extract_tarball($local_name,$distribution);
    }
    return 1;
}

#this is probably not going to be needed again, as the nightly
#bioperl build names have been simplified
sub determine_filename {
  my $listing = "dirlisting.html";
  my $rc = mirror(BIOPERL_LIVE_URL, $listing);
  die "Could not get directory listing of bioperl nightly build url: $rc\n"
      unless ($rc == RC_OK or $rc == RC_NOT_MODIFIED);

  my $filename; 
  open my $LIST, '<', $listing or die "Could not read file '$listing': $!\n";
  while (my $line = <$LIST>) {
    if ($line =~ /href="(bioperl-live.*?\.tar\.gz)"/) {
      $filename = $1;
      last;
    }
  }
  close $LIST;
  unlink $listing; 
  return $filename;
}

sub extract_tarball {
  my ($local_name,$distribution) = @_;

  print STDERR "Unpacking $local_name...\n";
  my $z = Archive::Tar->new($local_name,1)
        or die "Couldn't open $distribution archive: $@";
  my @extracted = $z->extract()
        or die "Couldn't extract $distribution archive: $@";

  if ($extracted[0]->{'name'} =~ /^(bioperl.*?)\//) {
    my $bioperl_dir = $1;
    move($bioperl_dir, $distribution) or die "couldn't move bioperl dir: $@"; 
  }

  chdir $distribution
        or die "Couldn't enter $distribution directory: $@";
  return;
}

# make sure ppm repositories are correct!
sub setup_ppm {
  open my $S, "ppm repo list --csv|" or die "Could not open ppm for listing: $!\n";
  my %repository;
  while (my $line = <$S>) {
     chomp $line;
     my ($index, $package_count, $name) = split /,/, $line;
     $repository{$name} = $index;
  }
  close $S;
  print STDERR "Adding needed PPM repositories. This may take a while....\n";
  for my $name (keys %REPOSITORIES) {
     next if $repository{$name};
     system("ppm rep add $name $REPOSITORIES{$name}");
  }
}

sub find_bioperl_ppm {
  print STDERR "Finding most recent bioperl...";
  open my $S,"ppm search bioperl |" or die "Could not open ppm for listing: $!\n";
  local $/ = ''; # paragraph mode
  my ($blessed_one, $blessed_version);
  my $best = 0;
  while (my $line = <$S>) {
    chomp $line;
    my ($number)     = ($line =~ /^(\d+): bioperl/m);
    my ($version)    = ($line =~ /^\s+Version: (.+)/m);
    my ($repository) = ($line =~ /^\s+Repo: (.+)/m);
    my $multiplier = 10000000;
    my $magnitude  = 0;
    # this dumb thing converts 1.5.1 into a real number
    foreach my $piece (split /[._]/, $version) {
      $magnitude += $piece * ($multiplier/=10);
    }
    ($blessed_one,$best,$blessed_version) = ($number,$magnitude,$version) if $best < $magnitude;
  }
  close $S;
  print STDERR $blessed_version ? "found $blessed_version\n" : "not found\n";
  return $blessed_one;
}
