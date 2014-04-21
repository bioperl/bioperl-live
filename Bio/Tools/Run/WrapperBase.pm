#
# BioPerl module for Bio::Tools::Run::WrapperBase
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::WrapperBase - A Base object for wrappers around executables

=head1 SYNOPSIS

  # do not use this object directly, it provides the following methods
  # for its subclasses

  my $errstr = $obj->error_string();
  my $exe    = $obj->executable();
  $obj->save_tempfiles($booleanflag)
  my $outfile= $obj->outfile_name();
  my $tempdir= $obj->tempdir(); # get a temporary dir for executing
  my $io     = $obj->io;  # Bio::Root::IO object
  my $cleanup= $obj->cleanup(); # remove tempfiles

  $obj->run({-arg1 => $value});

=head1 DESCRIPTION

This is a basic module from which to build executable wrapper modules.
It has some basic methods to help when implementing new modules.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track of
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 CONTRIBUTORS

Sendu Bala, bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Run::WrapperBase;
use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::Root::Root);

use File::Spec;
use File::Path qw(); # don't import anything

=head2 run

 Title   : run
 Usage   : $wrapper->run({ARGS HERE});
 Function: Support generic running with args passed in
           as a hashref
 Returns : Depends on the implementation, status OR data
 Args    : hashref of named arguments


=cut

sub run {
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}


=head2 error_string

 Title   : error_string
 Usage   : $obj->error_string($newval)
 Function: Where the output from the last analysis run is stored.
 Returns : value of error_string
 Args    : newvalue (optional)


=cut

sub error_string{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_error_string'} = $value;
    }
    return $self->{'_error_string'} || '';
}

=head2 arguments

 Title   : arguments
 Usage   : $obj->arguments($newval)
 Function: Commandline parameters
 Returns : value of arguments
 Args    : newvalue (optional)


=cut

sub arguments {
  my ($self,$value) = @_;
  if(defined $value) {
    $self->{'_arguments'} = $value;
  }
  return $self->{'_arguments'} || '';
}


=head2 no_param_checks

 Title   : no_param_checks
 Usage   : $obj->no_param_checks($newval)
 Function: Boolean flag as to whether or not we should
           trust the sanity checks for parameter values
 Returns : value of no_param_checks
 Args    : newvalue (optional)


=cut

sub no_param_checks{
   my ($self,$value) = @_;
   if( defined $value || ! defined $self->{'no_param_checks'} ) {
       $value = 0 unless defined $value;
      $self->{'no_param_checks'} = $value;
    }
    return $self->{'no_param_checks'};
}

=head2 save_tempfiles

 Title   : save_tempfiles
 Usage   : $obj->save_tempfiles($newval)
 Function: Get/set the choice of if tempfiles in the temp dir (see tempdir())
           are kept or cleaned up. Default is '0', ie. delete temp files.
           NB: This must be set to the desired value PRIOR to first creating
           a temp dir with tempdir(). Any attempt to set this after tempdir creation will get a warning.
 Returns : boolean
 Args    : none to get, boolean to set

=cut

sub save_tempfiles{
    my $self = shift;
    my @args = @_;
    if (($args[0]) && (exists ($self->{'_tmpdir'}))) {
        $self->warn ("Tempdir already created; setting save_tempfiles will not affect cleanup behavior.");
    }
    return $self->io->save_tempfiles(@_);
}

=head2 outfile_name

 Title   : outfile_name
 Usage   : my $outfile = $wrapper->outfile_name();
 Function: Get/Set the name of the output file for this run
           (if you wanted to do something special)
 Returns : string
 Args    : [optional] string to set value to


=cut

sub outfile_name{
   my ($self,$nm) = @_;
   if( defined $nm || ! defined $self->{'_outfilename'} ) {
       $nm = 'mlc' unless defined $nm;
       $self->{'_outfilename'} = $nm;
   }
   return $self->{'_outfilename'};
}


=head2 tempdir

 Title   : tempdir
 Usage   : my $tmpdir = $self->tempdir();
 Function: Retrieve a temporary directory name (which is created)
 Returns : string which is the name of the temporary directory
 Args    : none


=cut

sub tempdir{
   my ($self) = shift;

   $self->{'_tmpdir'} = shift if @_;
   unless( $self->{'_tmpdir'} ) {
       $self->{'_tmpdir'} = $self->io->tempdir(CLEANUP => ! $self->save_tempfiles );
   }
   unless( -d $self->{'_tmpdir'} ) {
       mkdir($self->{'_tmpdir'},0777);
   }
   return $self->{'_tmpdir'};
}

=head2 cleanup

 Title   : cleanup
 Usage   : $wrapper->cleanup();
 Function: Will cleanup the tempdir directory
 Returns : none
 Args    : none


=cut

sub cleanup{
   my ($self) = @_;
   $self->io->_io_cleanup();
   if( defined $self->{'_tmpdir'} && -d $self->{'_tmpdir'} ) {
      my $verbose = ($self->verbose >= 1) ? 1 : 0;
      File::Path::rmtree( $self->{'_tmpdir'}, $verbose);
   }
}

=head2 io

 Title   : io
 Usage   : $obj->io($newval)
 Function: Gets a Bio::Root::IO object
 Returns : Bio::Root::IO object
 Args    : none


=cut

sub io{
   my ($self) = @_;
   unless( defined $self->{'io'} ) {
       $self->{'io'} = Bio::Root::IO->new(-verbose => $self->verbose);
   }
    return $self->{'io'};
}

=head2 version

 Title   : version
 Usage   : $version = $wrapper->version()
 Function: Returns the program version (if available)
 Returns : string representing version of the program
 Args    : [Optional] value to (re)set version string


=cut

sub version{
   my ($self,@args) = @_;
   return;
}

=head2 executable

 Title   : executable
 Usage   : my $exe = $factory->executable();
 Function: Finds the full path to the executable
 Returns : string representing the full path to the exe
 Args    : [optional] name of executable to set path to
           [optional] boolean flag whether or not warn when exe is not found

=cut

sub executable {
    my ($self, $exe, $warn) = @_;

    if (defined $exe) {
        $self->{'_pathtoexe'} = $exe;
    }

    unless( defined $self->{'_pathtoexe'} ) {
        my $prog_path = $self->program_path;

        if ($prog_path) {
            if (-f $prog_path && -x $prog_path) {
                $self->{'_pathtoexe'} = $prog_path;
            }
            elsif ($self->program_dir) {
                $self->warn("executable not found in $prog_path, trying system path...") if $warn;
            }
        }
        unless ($self->{'_pathtoexe'}) {
            my $exe;
            if ( $exe = $self->io->exists_exe($self->program_name) ) {
                $self->{'_pathtoexe'} = $exe;
            }
            else {
                $self->warn("Cannot find executable for ".$self->program_name) if $warn;
                $self->{'_pathtoexe'} = undef;
            }
        }
    }

    # bail if we never found the executable
    unless ( defined $self->{'_pathtoexe'}) {
        $self->throw("Cannot find executable for ".$self->program_name .
            ". path=\"".$self->program_path."\"");
    }
    return $self->{'_pathtoexe'};
}

=head2 program_path

 Title   : program_path
 Usage   : my $path = $factory->program_path();
 Function: Builds path for executable
 Returns : string representing the full path to the exe
 Args    : none

=cut

sub program_path {
   my ($self) = @_;
   my @path;
   push @path, $self->program_dir if $self->program_dir;
   push @path, $self->program_name.($^O =~ /mswin/i ? '.exe' : '') if $self->program_name;
   return File::Spec->catfile(@path);
}

=head2 program_dir

 Title   : program_dir
 Usage   : my $dir = $factory->program_dir();
 Function: Abstract get method for dir of program. To be implemented
           by wrapper.
 Returns : string representing program directory
 Args    : none

=cut

sub program_dir {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 program_name

 Title   : program_name
 Usage   : my $name = $factory->program_name();
 Function: Abstract get method for name of program. To be implemented
           by wrapper.
 Returns : string representing program name
 Args    : none

=cut

sub program_name {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 quiet

 Title   : quiet
 Usage   : $factory->quiet(1);
           if ($factory->quiet()) { ... }
 Function: Get/set the quiet state. Can be used by wrappers to control if
           program output is printed to the console or not.
 Returns : boolean
 Args    : none to get, boolean to set

=cut

sub quiet {
    my $self = shift;
    if (@_) { $self->{quiet} = shift }
    return $self->{quiet} || 0;
}

=head2  _setparams()

 Title   : _setparams
 Usage   : $params = $self->_setparams(-params => [qw(window evalue_cutoff)])
 Function: For internal use by wrapper modules to build parameter strings
           suitable for sending to the program being wrapped. For each method
           name supplied, calls the method and adds the method name (as modified
           by optional things) along with its value (unless a switch) to the
           parameter string
 Example : $params = $self->_setparams(-params => [qw(window evalue_cutoff)],
                                       -switches => [qw(simple large all)],
                                       -double_dash => 1,
                                       -underscore_to_dash => 1);
           If window() and simple() had not been previously called, but
           evalue_cutoff(0.5), large(1) and all(0) had been called, $params
           would be ' --evalue-cutoff 0.5 --large'
 Returns : parameter string
 Args    : -params => [] or {}  # array ref of method names to call,
                                  or hash ref where keys are method names and
                                  values are how those names should be output
                                  in the params string
           -switches => [] or {}# as for -params, but no value is printed for
                                  these methods
           -join => string      # define how parameters and their values are
                                  joined, default ' '. (eg. could be '=' for
                                  param=value)
           -lc => boolean       # lc() method names prior to output in string
           -dash => boolean     # prefix all method names with a single dash
           -double_dash => bool # prefix all method names with a double dash
           -mixed_dash => bool  # prefix single-character method names with a
                                # single dash, and multi-character method names
                                # with a double-dash
           -underscore_to_dash => boolean # convert all underscores in method
                                            names to dashes

=cut

sub _setparams {
    my ($self, @args) = @_;

    my ($params, $switches, $join, $lc, $d, $dd, $md, $utd) =
        $self->_rearrange([qw(PARAMS
                              SWITCHES
                              JOIN
                              LC
                              DASH
                              DOUBLE_DASH
                              MIXED_DASH
                              UNDERSCORE_TO_DASH)], @args);
    $self->throw('at least one of -params or -switches is required') unless ($params || $switches);
    $self->throw("-dash, -double_dash and -mixed_dash are mutually exclusive") if (defined($d) + defined($dd) + defined($md) > 1);
    $join ||= ' ';

    my %params = ref($params) eq 'HASH' ? %{$params} : map { $_ => $_ } @{$params};
    my %switches = ref($switches) eq 'HASH' ? %{$switches} : map { $_ => $_ } @{$switches};

    my $param_string = '';
    for my $hash_ref (\%params, \%switches) {
        while (my ($method, $method_out) = each %{$hash_ref}) {
            my $value = $self->$method();
            next unless (defined $value);
            next if (exists $switches{$method} && ! $value);

            $method_out = lc($method_out) if $lc;
            my $method_length = length($method_out) if $md;
            $method_out = '-'.$method_out if ($d || ($md && ($method_length == 1)));
            $method_out = '--'.$method_out if ($dd || ($md && ($method_length > 1)));
            $method_out =~ s/_/-/g if $utd;

            if ( exists $params{$method} ) {
              # if value are quoted with " or ', re-quote it
              if ( $value =~ m{^[\'\"]+(.+)[\'\"]+$} ) {
                $value = '"'. $1 . '"';
              }
              # quote values that contain spaces
              elsif ( $value =~ m{\s+} ) {
                $value = '"'. $value . '"';
              }
            }

            $param_string .= ' '.$method_out.(exists $switches{$method} ? '' : $join.$value);
        }
    }

    return $param_string;
}

sub DESTROY {
    my $self= shift;
    unless ( $self->save_tempfiles ) {
	$self->cleanup();
    }
    $self->SUPER::DESTROY();
}


1;
