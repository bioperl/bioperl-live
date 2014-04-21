#
# BioPerl module for Bio::Tools::Run::WrapperBase::CommandExts
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Mark A. Jensen <maj -at- fortinbras -dot- us>
#
# Copyright Mark A. Jensen
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::WrapperBase::CommandExts - Extensions to WrapperBase for handling programs with commands *ALPHA*

=head1 SYNOPSIS

Devs, see L</DEVELOPER INTERFACE>.
Users, see L</USER INTERFACE>.

=head1 DESCRIPTION

This is a developer-focused experimental module. The main idea is to
extend L<Bio::Tools::Run::WrapperBase> to make it relatively easy to
create run wrappers around I<suites> of related programs, like
C<samtools> or C<blast+>.

Some definitions:

=over

=item * program

The program is the command-line frontend application. C<samtools>, for example, is run from the command line as follows:

 $ samtools view -bS in.bam > out.sam
 $ samtools faidx

=item * command

The command is the specific component of a suite run by executing the
program. In the example above, C<view> and C<faidx> are commands.

=item * command prefix

The command prefix is an abbreviation of the command name used
internally by C<CommandExts> method, and sometimes by the user of the
factory for specifying command line parameters to subcommands of
composite commands.

=item * composite command

A composite command is a pipeline or script representing a series of
separate executions of different commands. Composite commands can be
specified by configuring C<CommandExts> appropriately; the composite
command can be run by the user from a factory in the same way as
ordinary commands.

=item * options, parameters, switches and filespecs

An option is any command-line option; i.e., a specification set off by
a command-line by a specifier (like C<-v> or C<--outfile>). Parameters
are command-line options that accept a value (C<-title mydb>);
switches are boolean flags (C<--no-filter>). Filespecs are barewords
at the end of the command line that usually indicate input or output
files. In this module, this includes files that capture STDIN, STDOUT,
or STDERR via redirection.

=item * pseudo-program

A "pseudo-program" is a way to refer to a collection of related
applications that are run independently from the command line, rather
than via a frontend program. The C<blast+> suite of programs is an
example: C<blastn>, C<makeblastdb>, etc. C<CommandExts> can be
configured to create a single factory for a suite of related,
independent programs that treats each independent program as a
"pseudo-program" command.

=back

This module essentially adds the non-assembler-specific wrapper
machinery of fangly's L<Bio::Tools::Run::AssemblerBase> to the
L<Bio::Tools::Run::WrapperBase> namespace, adding the general
command-handling capability of L<Bio::Tools::Run::BWA>. It creates run
factories that are automatically Bio::ParameterBaseI compliant,
meaning that C<available_parameters()>, C<set_parameters()>,
C<get_parameters>, C<reset_parameters()>, and C<parameters_changed()>
are available.

=head1 DEVELOPER INTERFACE

C<CommandExts> is currently set up to read particular package globals
which define the program, the commands available, command-line options
for those commands, and human-readable aliases for those options.

The easiest way to use C<CommandExts> is probably to create two modules:

 Bio::Tools::Run::YourRunPkg
 Bio::Tools::Run::YourRunPkg::Config

The package globals should be defined in the C<Config> module, and the
run package itself should begin with the following mantra:

 use YourRunPkg::Config;
 use Bio::Tools::Run::WrapperBase;
 use Bio::Tools::Run::WrapperBase::CommandExts;
 sub new {
     my $class = shift;
     my @args = @_;
     my $self = $class->SUPER::new(@args);
     ...
     return $self;
 }

The following globals can/should be defined in the C<Config> module:

  $program_name
  $program_dir
  $use_dash
  $join
  @program_commands 
  %command_prefixes
  @program_params
  @program_switches 
  %param_translation
  %composite_commands
  %command_files

See L</Config Globals> for detailed descriptions.

The work of creating a run wrapper with C<CommandExts> lies mainly in
setting up the globals. The key methods for the developer interface are:

=over 

=item * program_dir($path_to_programs)

Set this to point the factory to the executables.

=item * _run(@file_args)

Runs an instantiated factory with the given file args. Use in the
 C<run()> method override.

=item *  _create_factory_set()

Returns a hash of instantiated factories for each true command from a
composite command factory. The hash keys are the true command names, so
you could do

 $cmds = $composite_fac->_create_factory_set;
 for (@true_commands) {
    $cmds->{$_}->_run(@file_args);
 }

=item * executables($cmd,[$fullpath])

For pseudo-programs, this gets/sets the full path to the executable of
the true program corresponding to the command C<$cmd>.

=back

=head2 Implementing Composite Commands

=head2 Implementing Pseudo-programs

To indicate that a package wraps disparate programs under a single pseudo program, use an asterisk before the program name:

 package Bio::Tools::Run::YourPkg::Config;
 ...
 our $program_name = '*blast+';

and C<_run> will know what to do. Specify the rest of the globals as
if the desired programs were commands. Use the basename of the
programs for the command names.

If all the programs can be found in a single directory, just specify
that directory in C<program_dir()>. If not, use C<executables()> to set the paths to each program explicitly:

 foreach (keys %cmdpaths) {
    $self->executables($_, $cmdpaths{$_});
 }

=head2 Config Globals

Here is an example config file. Further details in prose are below.

 package Dummy::Config;
 use strict;
 use warnings;
 no warnings qw(qw);
 use Exporter;
 our (@ISA, @EXPORT, @EXPORT_OK);
 push @ISA, 'Exporter';
 @EXPORT = qw(
              $program_name
              $program_dir
              $use_dash
              $join
              @program_commands
              %command_prefixes
              @program_params
              @program_switches
              %param_translation
              %command_files
              %composite_commands
             );

 our $program_name = '*flurb';
 our $program_dir = 'C:\cygwin\usr\local\bin';
 our $use_dash = 'mixed';
 our $join = ' ';
 
 our @program_commands = qw(
  rpsblast
  find
  goob
  blorb
  multiglob
   );

 our %command_prefixes = (
     blastp => 'blp',
     tblastn => 'tbn',
     goob => 'g',
     blorb => 'b',
     multiglob => 'm'
     );

 our @program_params = qw(
     command
     g|narf
     g|schlurb
     b|scroob
     b|frelb
     m|trud
 );
 
 our @program_switches = qw(
     g|freen
     b|klep
 );
 
 our %param_translation = (
     'g|narf'     => 'n',
     'g|schlurb'  => 'schlurb',
     'g|freen'    => 'f',
     'b|scroob'   => 's',
     'b|frelb'    => 'frelb'
     );
 
 our %command_files = (
     'goob'       => [qw( fas faq )],
     );
 
 our %composite_commands = (
     'multiglob' => [qw( blorb goob )]
     );
 1;

C<$use_dash> can be one of C<single>, C<double>, or C<mixed>. See L<Bio::Tools::Run::WrapperBase>.

There is a syntax for the C<%command_files> specification. The token
matching C<[a-zA-Z0-9_]+> in each element of each arrayref becomes the
named filespec parameter for the C<_run()> method in the wrapper
class. Additional symbols surrounding this token indicate how this
argument should be handled. Some examples:

 >out  : stdout is redirected into the file 
         specified by (..., -out => $file,... )
 <in   : stdin is accepted from the file 
         specified by (..., -in => $file,... )
 2>log : stderr is redirected into the file
         specified by (..., -log => $file,... )
 #opt  : this filespec argument is optional
         (no throw if -opt => $option is missing)
 2>#log: if -log is not specified in the arguments, the stderr() 
         method will capture stderr
 *lst  : this filespec can take multiple arguments,
         specify using an arrayref (..., -lst => [$file1, $file2], ...)
 *#lst : an optional list

The tokens above are examples; they can be anything matching the above regexp.

=head1 USER INTERFACE

Using a wrapper created with C<Bio::Tools::Run::WrapperBase::CommandExts>:

=over 

=item * Getting a list of available commands, parameters, and filespecs:

To get a list of commands, simply:

 @commands = Bio::Tools::Run::ThePkg->available_commands;

The wrapper will generally have human-readable aliases for each of the
command-line options for the wrapped program and commands. To obtain a
list of the parameters and switches available for a particular
command, do

 $factory = Bio::Tools::Run::ThePkg->new( -command => 'glurb' );
 @params = $factory->available_parameters('params');
 @switches = $factory->available_parameters('switches');
 @filespec = $factory->available_parameters('filespec');
 @filespec = $factory->filespec; # alias

=item * Create factories

The factory is a handle on the program and command you wish to
run. Create a factory using C<new> to set command-line parameters:

 $factory = Bio::Tools::Run::ThePkg->new( -command => 'glurb', 
                                          -freen => 1,
                                          -furschlugginer => 'vreeble' );

A shorthand for this is:
 
 $factory = Bio::Tools::Run::ThePkg->new_glurb( 
                                       -freen => 1, 
                                       -furschlugginer => 'vreeble' );

=item * Running programs

To run the program, use the C<run> method, providing filespecs as arguments

 $factory = Bio::Tools::Run::ThePkg->new_assemble( -min_qual => 63 );
 $factory->run( -faq1 => 'read1.fq', -faq2 => 'read2.fq', 
                -ref => 'refseq.fas', -out => 'new.sam' );
 # do another
 $factory->run( -faq1 => 'read-old1.fq', -faq2 => 'read-old2.fq', 
                -ref => 'refseq.fas', -out => 'old.sam' ); 

Messages on STDOUT and STDERR are dumped into their respective attributes:

 $stdout = $factory->stdout;
 $stderr = $factory->stderr;

unless STDOUT and/or STDERR are part of the named files in the filespec.

=item * Setting/getting/resetting/polling parameters.

A C<CommandExts>-based factory is always L<Bio::ParameterBaseI>
compliant. That means that you may set, get, and reset parameters
using C<set_parameters()>, C<get_parameters()>, and
C<reset_parameters>. You can ask whether parameters have changed since
they were last accessed by using the predicate
C<parameters_changed>. See L<Bio::ParameterBaseI> for more details.

Once set, parameters become attributes of the factory. Thus, you can get their values as follows:

 if ($factory->freen) { 
    $furs = $factory->furshlugginer;
    #...
 }

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support

Please direct usage questions or support issues to the mailing list:

L<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Mark A. Jensen

Email maj -at- fortinbras -dot- us

Describe contact details here

=head1 CONTRIBUTORS

Dan Kortschak ( dan -dot- kortschak -at- adelaide -dot- edu -dot- au )

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::Run::WrapperBase; # need these methods in WrapperBase/maj
use strict;
use warnings;
no warnings qw(redefine);

use Bio::Root::Root;
use File::Spec;
use IPC::Run;
use base qw(Bio::Root::Root Bio::ParameterBaseI);

our $AUTOLOAD;

=head2 new()

 Title   : new
 Usage   : 
 Function: constructor for WrapperBase::CommandExts ; 
           correctly binds configuration variables
           to the WrapperBase object
 Returns : Bio::Tools::Run::WrapperBase object with command extensions
 Args    : 
 Note    : this method subsumes the old _register_program_commands and
           _set_program_options, leaving out the assembler-specific
           parms ($qual_param and out_type())

=cut

sub new {
    my ($class, @args) = @_;
    my $self = bless ({}, $class);
    # pull in *copies* of the Config variables from the caller namespace:
    my ($pkg, @goob) = caller();
    my ($commands,
	$prefixes,
	$params,
	$switches,
	$translation,
	$use_dash,
	$join,
	$name,
	$dir,
	$composite_commands,
	$files);
    for (qw( @program_commands 
             %command_prefixes
             @program_params
             @program_switches 
             %param_translation
             $use_dash
             $join
             $program_name
             $program_dir
             %composite_commands
             %command_files ) ) {
	my ($sigil, $var) = m/(.)(.*)/;
	my $qualvar = "${sigil}${pkg}::${var}";
	for ($sigil) {
	    /\@/ && do { $qualvar = "\[$qualvar\]" };
	    /\%/ && do { $qualvar = "\{$qualvar\}" };
	}
	my $locvar = "\$${var}";
	$locvar =~ s/program_|command_|param_//g;
	eval "$locvar = $qualvar";
    }
    # set up the info registry hash
    my %registry;
    if ($composite_commands) {
	$self->_register_composite_commands($composite_commands,
					    $params,
					    $switches,
					    $prefixes);
    }
    @registry{qw( _commands _prefixes _files 
                  _params _switches _translation
                  _composite_commands )} =
	($commands, $prefixes, $files,
	 $params, $switches, $translation, 
	 $composite_commands);
    $self->{_options} = \%registry;
    if (not defined $use_dash) {
	$self->{'_options'}->{'_dash'}      = 1;
    } else {
	$self->{'_options'}->{'_dash'}      = $use_dash;
    }
    if (not defined $join) {
	$self->{'_options'}->{'_join'}      = ' ';
    } else {
	$self->{'_options'}->{'_join'}      = $join;
    }
    if ($name =~ /^\*/) {
	$self->is_pseudo(1);
	$name =~ s/^\*//;
    }
    $self->program_name($name) if not defined $self->program_name();
    $self->program_dir($dir) if not defined $self->program_dir();
    $self->set_parameters(@args);
    $self->parameters_changed(1); # set on instantiation, per Bio::ParameterBaseI
    return $self;
}

=head2 program_name

 Title   : program_name
 Usage   : $factory->program_name($name)
 Function: get/set the executable name
 Returns:  string
 Args    : string

=cut

sub program_name {
    my ($self, $val) = @_;
    $self->{'_program_name'} = $val if $val;
    return $self->{'_program_name'};
}

=head2 program_dir

 Title   : program_dir
 Usage   : $factory->program_dir($dir)
 Function: get/set the program dir
 Returns:  string
 Args    : string

=cut

sub program_dir {
    my ($self, $val) = @_;
    $self->{'_program_dir'} = $val if $val;
    return $self->{'_program_dir'};
}

=head2 _register_program_commands()

 Title   : _register_program_commands
 Usage   : $factory->_register_program_commands( \@commands, \%prefixes )
 Function: Register the commands a program accepts (for programs that act
           as frontends for a set of commands, each command having its own
           set of params/switches)
 Returns : true on success
 Args    : arrayref to a list of commands (scalar strings),
           hashref to a translation table of the form
           { $prefix1 => $command1, ... } [optional]
 Note    : To implement a program with this kind of calling structure, 
           include a parameter called 'command' in the 
           @program_params global
 Note    : The translation table is used to associate parameters and 
           switches specified in _set_program_options with the correct
           program command. In the globals @program_params and
           @program_switches, specify elements as 'prefix1|param' and 
           'prefix1|switch', etc.

=cut

=head2 _set_program_options

 Title   : _set_program_options
 Usage   : $factory->_set_program_options( \@ args );
 Function: Register the parameters and flags that an assembler takes.
 Returns : 1 for success
 Args    : - arguments passed by the user
           - parameters that the program accepts, optional (default: none)
           - switches that the program accepts, optional (default: none)
           - parameter translation, optional (default: no translation occurs)
           - dash option for the program parameters, [1|single|double|mixed],
             optional (default: yes, use single dashes only)
           - join, optional (default: ' ')

=cut

=head2 _translate_params

 Title   : _translate_params
 Usage   : @options = $assembler->_translate_params( );
 Function: Translate the Bioperl arguments into the arguments to pass to the
           program on the command line
 Returns : Arrayref of arguments
 Args    : none

=cut

sub _translate_params {
  my ($self)   = @_;
  # Get option string
  my ($params, $switches, $join, $dash, $translat) =
      @{$self->{_options}}{qw(_params _switches _join _dash _translation)};

  # access the multiple dash choices of _setparams...
  my @dash_args;
  $dash ||= 1; # default as advertised
  for ($dash) {
      $_ eq '1' && do {
	  @dash_args = ( -dash => 1 );
	  last;
      };
      /^s/ && do { #single dash only
	  @dash_args = ( -dash => 1);
	  last;
      };
      /^d/ && do { # double dash only
	  @dash_args = ( -double_dash => 1);
	  last;
      };
      /^m/ && do { # mixed dash: one-letter opts get -,
                  # long opts get --
	  @dash_args = ( -mixed_dash => 1);
	  last;
      };
      do { 
	  $self->warn( "Dash spec '$dash' not recognized; using 'single'" );
	  @dash_args = ( -dash => 1 );
      };
  }
  my $options  = $self->_setparams(
    -params    => $params,
    -switches  => $switches,
    -join      => $join,
    @dash_args
  );

  # Translate options
  my @options  = split(/(\s|$join)/, $options);
  for (my $i = 0; $i < scalar @options; $i++) {
    my ($prefix, $name) = ( $options[$i] =~ m/^(-{0,2})(.+)$/ );
    if (defined $name) {
	if ($name =~ /command/i) {
	    $name = $options[$i+2]; # get the command
	    splice @options, $i, 4;
	    $i--;
	    # don't add the command if this is a pseudo-program
	    unshift @options, $name unless ($self->is_pseudo); # put command first
	}
	elsif (defined $$translat{$name}) {
	    $options[$i] = $prefix.$$translat{$name};
	}
    } 
    else {
	splice @options, $i, 1;
	$i--;
    }
  }
  $options = join('', @options);

  # this is a kludge for mixed options: the reason mixed doesn't 
  # work right on the pass through _setparams is that the 
  # *aliases* and not the actual params are passed to it. 
  # here we just rejigger the dashes
  if ($dash =~ /^m/) {
      $options =~ s/--([a-z0-9](?:\s|$))/-$1/gi;
  }

  # Now arrayify the options
  @options = split(' ', $options);

  return \@options;
}

=head2 executable()

 Title   : executable
 Usage   : 
 Function: find the full path to the main executable,
           or to the command executable for pseudo-programs
 Returns : full path, if found
 Args    : [optional] explicit path to the executable
           (will set the appropriate command exec if
            applicable)
           [optional] boolean flag whether or not to warn when exe no found
 Note    : overrides WrapperBase.pm
            
=cut

sub executable {
    my $self = shift;
    my ($exe, $warn) = @_;
    if ($self->is_pseudo) {
	return $self->{_pathtoexe} = $self->executables($self->command,$exe);
    }

    # otherwise
    # setter
    if (defined $exe) {
	$self->throw("binary '$exe' does not exist") unless -e $exe;
	$self->throw("'$exe' is not executable") unless -x $exe;
	return $self->{_pathtoexe} = $exe;
    }

    # getter
    return $self->{_pathtoexe} if defined $self->{_pathstoexe};

    # finder
    return $self->{_pathtoexe} = $self->_find_executable($exe, $warn);
}

=head2 executables()

 Title   : executables
 Usage   : 
 Function: find the full path to a command's executable
 Returns : full path (scalar string)
 Args    : command (scalar string), 
           [optional] explicit path to this command exe
           [optional] boolean flag whether or not to warn when exe no found

=cut

sub executables {
    my $self = shift;
    my ($cmd, $exe, $warn) = @_;
    # for now, barf if this is not a pseudo program
    $self->throw("This wrapper represents a single program with commands, not multiple programs; can't use executables()") unless $self->is_pseudo;
    $self->throw("Command name required at arg 1") unless defined $cmd;
    $self->throw("The desired executable '$cmd' is not registered as a command") unless grep /^$cmd$/, @{$self->{_options}->{_commands}};

    # setter
    if (defined $exe) {
	$self->throw("binary '$exe' does not exist") unless -e $exe;
	$self->throw("'$exe' is not executable") unless -x $exe;
	$self->{_pathstoexe} = {} unless defined $self->{_pathstoexe};
	return $self->{_pathstoexe}->{$cmd} = $exe;
    }

    # getter
    return $self->{_pathstoexe}->{$cmd} if defined $self->{_pathstoexe}->{$cmd};
    
    $exe ||= $cmd;
    # finder
    return $self->{_pathstoexe}->{$cmd} = $self->_find_executable($exe, $warn);
}

=head2 _find_executable()

 Title   : _find_executable
 Usage   : my $exe_path = $fac->_find_executable($exe, $warn);
 Function: find the full path to a named executable,
 Returns : full path, if found
 Args    : name of executable to find
           [optional] boolean flag whether or not to warn when exe no found
 Note    : differs from executable and executables in not
           setting any object attributes

=cut

sub _find_executable {
    my $self = shift;
    my ($exe, $warn) = @_;

    if ($self->is_pseudo && !$exe) {
	if (!$self->command) {
	    # this throw probably appropriate
	    # the rest are now warns if $warn.../maj
	    $self->throw( 
		"The ".__PACKAGE__." wrapper represents several different programs;".
		"arg1 to _find_executable must be specified explicitly,".
		"or the command() attribute set");
	}
	else {
	    $exe = $self->command;
	}
    }
    $exe ||= $self->program_path;

    my $path;
    if ($self->program_dir) {
	$path = File::Spec->catfile($self->program_dir, $exe);
    } else {
	$path = $exe;
	$self->warn('Program directory not specified; use program_dir($path).') if $warn;
    }

    # use provided info - we are allowed to follow symlinks, but refuse directories
    map { return $path.$_ if ( -x $path.$_ && !(-d $path.$_) ) } ('', '.exe') if defined $path;

    # couldn't get path to executable from provided info, so use system path
    $path = $path ? " in $path" : undef;
    $self->warn("Executable $exe not found$path, trying system path...") if $warn;
    if ($path = $self->io->exists_exe($exe)) {
	return $path;
    } else {
	$self->warn("Cannot find executable for program '".($self->is_pseudo ? $self->command : $self->program_name)."'") if $warn;
	return;
    }
}

=head2 _register_composite_commands()

 Title   : _register_composite_commands
 Usage   : 
 Function: adds subcomand params and switches for composite commands
 Returns : true on success
 Args    : \%composite_commands,
           \@program_params,
           \@program_switches

=cut

sub _register_composite_commands {
    my $self = shift;
    my ($composite_commands, $program_params, 
	$program_switches, $command_prefixes) = @_;
    my @sub_params;
    my @sub_switches;
    foreach my $cmd (keys %$composite_commands) {
	my $pfx = $command_prefixes->{$cmd} || $cmd;
	foreach my $subcmd ( @{$$composite_commands{$cmd}} ) {
	    my $spfx = $command_prefixes->{$subcmd} || $subcmd;
	    my @sub_program_params = grep /^$spfx\|/, @$program_params;
	    my @sub_program_switches = grep /^$spfx\|/, @$program_switches;
	    for (@sub_program_params) {
		m/^$spfx\|(.*)/;
		push @sub_params, "$pfx\|${spfx}_".$1;
	    }
	    for (@sub_program_switches) {
		m/^$spfx\|(.*)/;
		push @sub_switches, "$pfx\|${spfx}_".$1;
	    }
	}
    }
    push @$program_params, @sub_params;
    push @$program_switches, @sub_switches;
    # translations for subcmd params/switches not necessary
    return 1;
}

=head2 _create_factory_set()

 Title   : _create_factory_set
 Usage   : @facs = $self->_create_factory_set
 Function: instantiate a set of individual command factories for
           a given composite command
           Factories will have the correct parameter fields set for
           their own subcommand
 Returns : hash of factories: ( $subcmd_prefix => $subcmd_factory, ... )
 Args    : none

=cut

sub _create_factory_set {
    my $self = shift;
    $self->throw('command not set') unless $self->command;
    my $cmd = $self->command;
    $self->throw('_create_factory_set only works on composite commands') 
	unless grep /^$cmd$/, keys %{$self->{_options}->{_composite_commands}};
    my %ret;
    my $class = ref $self;
    my $subargs_hash = $self->_collate_subcmd_args($cmd);
    for (keys %$subargs_hash) {
	$ret{$_} = $class->new( -command => $_,  @{$$subargs_hash{$_}} );
    }
    return %ret;
}

=head2 _collate_subcmd_args()

 Title   : _collate_subcmd_args
 Usage   : $args_hash = $self->_collate_subcmd_args
 Function: collate parameters and switches into command-specific
           arg lists for passing to new()
 Returns : hash of named argument lists
 Args    : [optional] composite cmd prefix (scalar string) 
           [default is 'run']

=cut

sub _collate_subcmd_args {
    my $self = shift;
    my $cmd = shift;
    my %ret;
    # default command is 'run'
    $cmd ||= 'run';
    return unless $self->{'_options'}->{'_composite_commands'};
    return unless $self->{'_options'}->{'_composite_commands'}->{$cmd};
    my @subcmds = @{$self->{'_options'}->{'_composite_commands'}->{$cmd}};

    my $cur_options = $self->{'_options'};
    # collate
    foreach my $subcmd (@subcmds) {
	# find the composite cmd form of the argument in 
	# the current params and switches
	# e.g., map_max_mismatches
	my $pfx = $self->{_options}->{_prefixes}->{$subcmd} || $subcmd;
	my @params = grep /^${pfx}_/, @{$$cur_options{'_params'}};
	my @switches = grep /^${pfx}_/, @{$$cur_options{'_switches'}};
	$ret{$subcmd} = [];
	# create an argument list suitable for passing to new() of
	# the subcommand factory...
	foreach my $opt (@params, @switches) {
	    my $subopt = $opt; 
	    $subopt =~ s/^${pfx}_//; 
	    push(@{$ret{$subcmd}}, '-'.$subopt => $self->$opt) if defined $self->$opt;
	}
    }
    return \%ret;
}

=head2 _run

 Title   : _run
 Usage   : $fac->_run( @file_args )
 Function: Run a command as specified during object contruction
 Returns : true on success
 Args    : a specification of the files to operate on according
           to the filespec

=cut

sub _run {
    my ($self, @args) = @_;
    # _translate_params will provide an array of command/parameters/switches
    # -- these are set at object construction
    # to set up the run, need to add the files to the call
    # -- provide these as arguments to this function
    my $cmd = $self->command if $self->can('command');
    my $opts = $self->{_options};
    my %args; 
    $self->throw("No command specified for the object") unless $cmd;
    # setup files necessary for this command
    my $filespec = $opts->{'_files'}->{$cmd};
    my @switches;
    my ($in, $out, $err);
    # some applications rely completely on switches
    if (defined $filespec && @$filespec) {
	# parse args based on filespec
	# require named args
	$self->throw("Named args are required") unless !(@args % 2);
	s/^-// for @args;
	%args = @args;
	# validate
	my @req = map { 
	    my $s = $_;
	    $s =~ s/^-.*\|//;
	    $s =~ s/^[012]?[<>]//;
	    $s =~ s/[^a-zA-Z0-9_]//g; 
	    $s
	} grep !/[#]/, @$filespec;
	!defined($args{$_}) && $self->throw("Required filearg '$_' not specified") for @req;
	# set up redirects and file switches
	for (@$filespec) {
	    m/^1?>#?(.*)/ && do {
		defined($args{$1}) && ( open $out, '>', $args{$1} or $self->throw("Could not write file '$args{$1}': $!") );
		next;
	    };
	    m/^2>#?(.*)/ && do {
		defined($args{$1}) && ( open $err, '>', $args{$1} or $self->throw("Could not write file '$args{$1}': $!") );
		next;
	    };
	    m/^<#?(.*)/ && do {
		defined($args{$1}) && ( open $in, '<', $args{$1} or $self->throw("Could not read file '$args{$1}': $!") );
		next;
	    };
	    if (m/^-(.*)\|/) {
		push @switches, $self->_dash_switch($1);
	    } else {
		push @switches, undef;
            }
	}
    }
    my $dum;
    $in || ($in = \$dum);
    $out || ($out = \$self->{'stdout'});
    $err || ($err = \$self->{'stderr'});
    
    # Get program executable
    my $exe = $self->executable;
    $self->throw("Can't find executable for '".($self->is_pseudo ? $self->command : $self->program_name)."'; can't continue") unless $exe;

    # Get command-line options
    my $options = $self->_translate_params();
    # Get file specs sans redirects in correct order
    my @specs = map { 
	my $s = $_; 
	$s =~ s/^-.*\|//;
	$s =~ s/[^a-zA-Z0-9_]//g; 
	$s
    } grep !/[<>]/, @$filespec;
    my @files = @args{@specs};
    # expand arrayrefs
    my $l = $#files;
    
    # Note: below code block may be brittle, see link on this:
    # http://lists.open-bio.org/pipermail/bioperl-l/2010-June/033439.html
    
    for (0..$l) {
	if (ref($files[$_]) eq 'ARRAY') {
	    splice(@switches, $_, 1, ($switches[$_]) x @{$files[$_]});
	    splice(@files, $_, 1, @{$files[$_]});
	}
    }
    
    
    @files = map {
        my $s = shift @switches;
        defined $_ ? ($s, $_): ()
    } @files;
    @files = map { defined $_ ? $_ : () } @files; # squish undefs
    my @ipc_args = ( $exe, @$options, @files );
    $self->{_last_execution} = join( $self->{'_options'}->{'_join'}, @ipc_args );
    eval {
	IPC::Run::run(\@ipc_args, $in, $out, $err) or
	    die ("There was a problem running $exe : ".$$err);
    };

    if ($@) {
	$self->throw("$exe call crashed: $@") unless $self->no_throw_on_crash;
	return 0;
    }

     return 1;
}



=head2 no_throw_on_crash()

 Title   : no_throw_on_crash
 Usage   : 
 Function: prevent throw on execution error
 Returns : 
 Args    : [optional] boolean

=cut

sub no_throw_on_crash {
    my $self = shift;
    return $self->{'_no_throw'} = shift if @_;
    return $self->{'_no_throw'};
}

=head2 last_execution()

 Title   : last_execution
 Usage   : 
 Function: return the last executed command with options
 Returns : string of command line sent to IPC::Run
 Args    : 

=cut

sub last_execution {
    my $self = shift;
    return $self->{'_last_execution'};
}

=head2 _dash_switch()

 Title   : _dash_switch
 Usage   : $version = $fac->_dash_switch( $switch )
 Function: Returns an appropriately dashed switch for the executable
 Args    : A string containing a switch without dashes
 Returns : string containing an appropriately dashed switch for the current executable

=cut

sub _dash_switch {
	my ($self, $switch) = @_;

	my $dash = $self->{'_options'}->{'_dash'};
	for ($dash) {
		$_ eq '1' && do {
			$switch = '-'.$switch;
			last;
		};
		/^s/ && do { #single dash only
			$switch = '-'.$switch;
			last;
		};
		/^d/ && do { # double dash only
			$switch = '--'.$switch;
			last;
		};
		/^m/ && do { # mixed dash: one-letter opts get -,
			$switch = '-'.$switch;
			$switch =~ s/^(-[a-z0-9](?:\w+))$/-$1/i;
			last;
		};
		do { 
			$self->warn( "Dash spec '$dash' not recognized; using 'single'" );
			$switch = '-'.$switch;
		};
	}

	return $switch;
}

=head2 stdout()

 Title   : stdout
 Usage   : $fac->stdout()
 Function: store the output from STDOUT for the run, 
           if no file specified in _run arguments
 Example : 
 Returns : scalar string
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub stdout {
    my $self = shift;
    return $self->{'stdout'} = shift if @_;
    return $self->{'stdout'};
}

=head2 stderr()

 Title   : stderr
 Usage   : $fac->stderr()
 Function: store the output from STDERR for the run, 
           if no file is specified in _run arguments
 Example : 
 Returns : scalar string
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub stderr {
    my $self = shift;
    return $self->{'stderr'} = shift if @_;
    return $self->{'stderr'};
}

=head2 is_pseudo()

 Title   : is_pseudo
 Usage   : $obj->is_pseudo($newval)
 Function: returns true if this factory represents
           a pseudo-program
 Example : 
 Returns : value of is_pseudo (boolean)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub is_pseudo {
    my $self = shift;
    
    return $self->{'is_pseudo'} = shift if @_;
    return $self->{'is_pseudo'};
}

=head2 AUTOLOAD

AUTOLOAD permits 

 $class->new_yourcommand(@args);

as an alias for

 $class->new( -command => 'yourcommand', @args );

=cut

sub AUTOLOAD {
    my $class = shift;
    my $tok = $AUTOLOAD;
    my @args = @_;
    $tok =~ s/.*:://;
    unless ($tok =~ /^new_/) {
	$class->throw("Can't locate object method '$tok' via package '".ref($class)?ref($class):$class); 
    }
    my ($cmd) = $tok =~ m/new_(.*)/;
    return $class->new( -command => $cmd, @args );
}

=head1 Bio:ParameterBaseI compliance

=head2 set_parameters()

 Title   : set_parameters
 Usage   : $pobj->set_parameters(%params);
 Function: sets the parameters listed in the hash or array
 Returns : true on success
 Args    : [optional] hash or array of parameter/values.  

=cut

sub set_parameters {
    my ($self, @args) = @_;

    # currently stored stuff
    my $opts = $self->{'_options'};
    my $params = $opts->{'_params'};
    my $switches = $opts->{'_switches'};
    my $translation = $opts->{'_translation'};
    my $use_dash = $opts->{'_dash'};
    my $join = $opts->{'_join'};
    unless (($self->can('command') && $self->command) 
	    || (grep /command/, @args)) {
	push @args, '-command', 'run';
    }
    my %args = @args;
    my $cmd = $args{'-command'} || $args{'command'} || ($self->can('command') && $self->command);
    if ($cmd) {
	my (@p,@s, %x);
	$self->warn('Command present, but no commands registered') unless $self->{'_options'}->{'_commands'};
	$self->throw("Command '$cmd' not registered") unless grep /^$cmd$/, @{$self->{'_options'}->{'_commands'}};
	$cmd = $self->{_options}->{_prefixes}->{$cmd} || $cmd;
	
	@p = (grep(!/^.*?\|/, @$params), grep(/^${cmd}\|/, @$params));
	@s = (grep(!/^.*?\|/, @$switches), grep(/^${cmd}\|/, @$switches));
	s/.*?\|// for @p;
	s/.*?\|// for @s;
	@x{@p, @s} = @{$translation}{
	    grep( !/^.*?\|/, @$params, @$switches),
	    grep(/^${cmd}\|/, @$params, @$switches) };
	$opts->{_translation} = $translation = \%x;
	$opts->{_params} = $params = \@p;
	$opts->{_switches} = $switches = \@s;
    }
    $self->_set_from_args(
	\@args,
	-methods => [ @$params, @$switches, 'program_name', 'program_dir', 'out_type' ],
	-create =>  1,
	# when our parms are accessed, signal parameters are unchanged for
	# future reads (until set_parameters is called)
	-code => 
	' my $self = shift; 
          $self->parameters_changed(0);
          return $self->{\'_\'.$method} = shift if @_;
          return $self->{\'_\'.$method};'
	);
    # the question is, are previously-set parameters left alone when
    # not specified in @args?
    $self->parameters_changed(1);
    return 1;
}

=head2 reset_parameters()

 Title   : reset_parameters
 Usage   : resets values
 Function: resets parameters to either undef or value in passed hash
 Returns : none
 Args    : [optional] hash of parameter-value pairs

=cut

sub reset_parameters {
    my ($self, @args) = @_;

    my @reset_args;
    # currently stored stuff
    my $opts = $self->{'_options'};
    my $params = $opts->{'_params'};
    my $switches = $opts->{'_switches'};
    my $translation = $opts->{'_translation'};
    my $qual_param = $opts->{'_qual_param'};
    my $use_dash = $opts->{'_dash'};
    my $join = $opts->{'_join'};

    # handle command name
    my %args = @args;
    my $cmd = $args{'-command'} || $args{'command'} || $self->command;
    $args{'command'} = $cmd;
    delete $args{'-command'};
    @args = %args;
    # don't like this, b/c _set_program_args will create a bunch of
    # accessors with undef values, but oh well for now /maj

    for my $p (@$params) {
	push(@reset_args, $p => undef) unless grep /^[-]?$p$/, @args;
    }
    for my $s (@$switches) {
	push(@reset_args, $s => undef) unless grep /^[-]?$s$/, @args;
    }
    push @args, @reset_args;
    $self->set_parameters(@args);
    $self->parameters_changed(1);
}

=head2 parameters_changed()

 Title   : parameters_changed
 Usage   : if ($pobj->parameters_changed) {...}
 Function: Returns boolean true (1) if parameters have changed
 Returns : Boolean (0 or 1)
 Args    : [optional] Boolean

=cut

sub parameters_changed {
    my $self = shift;
    return $self->{'_parameters_changed'} = shift if @_;
    return $self->{'_parameters_changed'};
}

=head2 available_parameters()

 Title   : available_parameters
 Usage   : @params = $pobj->available_parameters()
 Function: Returns a list of the available parameters
 Returns : Array of parameters
 Args    : 'params' for settable program parameters
           'switches' for boolean program switches
           default: all 

=cut

sub available_parameters {
    my $self = shift;
    my $subset = shift;
    my $opts = $self->{'_options'};
    my @ret;
    for ($subset) {
	(!defined || /^a/) && do {
	    @ret = (@{$opts->{'_params'}}, @{$opts->{'_switches'}});
	    last;
	};
	m/^p/i && do {
	    @ret = @{$opts->{'_params'}};
	    last;
	};
	m/^s/i && do {
	    @ret = @{$opts->{'_switches'}};
	    last;
	};
	m/^c/i && do {
	    @ret = @{$opts->{'_commands'}};
	    last;
	};
	m/^f/i && do { # get file spec
	    return @{$opts->{'_files'}->{$self->command}};
	};
	do { #fail
	    $self->throw("available_parameters: unrecognized subset");
	};
    }
    return @ret;
}

sub available_commands { shift->available_parameters('commands') }
sub filespec { shift->available_parameters('filespec') }

=head2 get_parameters()

 Title   : get_parameters
 Usage   : %params = $pobj->get_parameters;
 Function: Returns list of key-value pairs of parameter => value
 Returns : List of key-value pairs
 Args    : [optional] A string is allowed if subsets are wanted or (if a
           parameter subset is default) 'all' to return all parameters

=cut

sub get_parameters {
    my $self = shift;
    my $subset = shift;
    $subset ||= 'all';
    my @ret;
    my $opts = $self->{'_options'};
    for ($subset) {
	m/^p/i && do { #params only
	    for (@{$opts->{'_params'}}) {
		push(@ret, $_, $self->$_) if $self->can($_) && defined $self->$_;
	    }
	    last;
	};
	m/^s/i && do { #switches only
	    for (@{$opts->{'_switches'}}) {
		push(@ret, $_, $self->$_) if $self->can($_) && defined $self->$_;
	    }
	    last;
	};
	m/^a/i && do { # all
	    for ((@{$opts->{'_params'}},@{$opts->{'_switches'}})) {
		push(@ret, $_, $self->$_) if $self->can($_) && defined $self->$_;
	    }
	    last;
	};
	do {
	    $self->throw("get_parameters: unrecognized subset");
	};
    }
    return @ret;
}

1;
