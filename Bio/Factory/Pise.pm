=head1 NAME

Bio::Factory::Pise

=head1 SYNOPSIS

=head1 DESCRIPTION

    A class to create Pise application objects.

       my $factory = new Bio::Factory::Pise(-email => 'me@myhome');

    Then you can create an application object (Pise::Run::Tools::PiseApplication):
       my $program = $factory->program('genscan');

    The email is mandatory, as in the Web interface. This is due to
    the fact that your program might enter infinite loops, or just run
    many jobs: the Pise server maintainer needs a contact (s/he
    could of course cancel any requests from your address...).
    If you plan to run a lot of heavy jobs, or to do a course with many 
    students, please ask the maintainer before.

    The remote parameter stands for the actual CGI location, except when 
    set at the factory creation step, where it is rather the root of all CGI.
    There are default values for most of Pise programs. 

    You can either set remote:
       1) at factory creation
           my $factory = Bio::Factory::Pise->new(-remote => 'http://somewhere/Pise/cgi-bin',
						 -email => 'me@myhome');
       2) at program creation:
           my $program = $factory->program('water', 
					   -remote => 'http://somewhere/Pise/cgi-bin/water.pl'
					   )
       3) at any time before running:
           $program->remote('http://somewhere/Pise/cgi-bin/water.pl');
           $job = $program->run();

       3) when running:
           $job = $program->run(-remote => 'http://somewhere/Pise/cgi-bin/water.pl');


    You can also retrieve a previous job results by providing its url:
           $job = $factory->job($url);
    You get the url of a job by:
           $job->jobid;


=cut

# Let the code begin...

package Bio::Factory::Pise;

use vars qw($AUTOLOAD @ISA %REMOTE);
use strict;

use Bio::Root::Root;
use Bio::Tools::Run::PiseApplication;
use Bio::Tools::Run::PiseJob;
use Bio::Factory::ApplicationFactoryI;

@ISA = qw(Bio::Root::Root Bio::Factory::ApplicationFactoryI );

%REMOTE = (
    'default' => 'http://bioweb.pasteur.fr/cgi-bin/seqanal',
    'clustalw' => 'http://bioweb.pasteur.fr/cgi-bin/seqanal/clustalw.pl'
);

=head2 new

 Title   : new()
 Usage   : my $program = Bio::Factory::Pise->new(-remote => 'http://somewhere/cgi-bin/Pise', -email -> $email);
 Function: Creates a Bio::Factory::Pise object, which function is to 
           create interface object (Bio::Tools::Run::program) for programs.
           Email is mandatory.
 Example :
 Returns : An instance of Bio::Factory::Pise.

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($remote) =
	$self->_rearrange([qw(REMOTE )],
			  @args);
    my ($email) =
	$self->_rearrange([qw(EMAIL )],
			@args);
    
    my ($verbose) =
	$self->_rearrange([qw(VERBOSE )],
			  @args);
    
    if (defined $remote) {
	$self->{REMOTE} = $remote;
    }

    if (defined $email) {
	$self->{EMAIL} = $email;
    } else {
	$self->throw("Email is mandatory.")
    }
    if (defined $verbose) {
	$self->{VERBOSE} = $verbose;
    } else {
	$self->{VERBOSE} = 0;
    }
    return $self;
}

=head2 program

 Title   : program()
 Usage   : my $program = Bio::Factory::Pise->program($program, -remote => 'http://somewhere/cgi-bin/Pise', -email -> $email, @params);
 Function: Creates a representation of a single Pise program.
 Example :
 Returns : An instance of Bio::Tools::Run::PiseApplication::$program.

=cut

sub program {
    my ($self, $program, @args) = @_;

    my ($remote) =
      $self->_rearrange([qw(REMOTE )],
			@args);
    my ($email) =
	$self->_rearrange([qw(EMAIL )],
			  @args);
    
    my ($verbose) =
	$self->_rearrange([qw(VERBOSE )],
			  @args);
    if (! $remote) {
	if (defined $self->{REMOTE}) {
	    if ($self->{REMOTE} =~ /$program/) {
		$remote = $self->{REMOTE};
	    } else {
		$remote = $self->{REMOTE} . "/$program.pl";
	    }
	} else {
	    if (defined $REMOTE{$program}) {
		$remote = $REMOTE{$program};
	    } else {
		$remote = $REMOTE{'default'} . "/$program.pl";
	    }
	}
    }
    if (! $email) {
	$email = $self->{EMAIL};
    }
    if (! $verbose) {
	$verbose= $self->{VERBOSE};
    }

    no strict "subs";

    my $package = "Bio::Tools::Run::PiseApplication::$program";  
    my $pise_program;

    eval ("use $package");
    $self->throw("Problem to load Bio::Tools::Run::PiseApplication::${program}\n\n$@")
	if $@;

    eval($pise_program = $package->new($remote, $email) );
    use strict "subs";

    foreach my $param ($pise_program->parameters_order) {
	my $param_name = $param;
	$param_name =~ tr/a-z/A-Z/;
	#print STDERR "setting $param_name ...?\n";
	my ($value) =
	    $self->_rearrange([$param_name],
			      @args);
	if ($value) {
	    #print STDERR "setting $param to $value\n";
	    $pise_program->$param($value);
	}
    }

    #print STDERR "creation of $pise_program ", ref($pise_program), "\n"; 
    return $pise_program;

}

=head2 job

 Title   : job(url)
 Usage   : my $job = Bio::Factory::Pise->job('http://somewhere/cgi-bin/Pise/tmp/dnapars/A3459687595869098');
 Function: Creates a previously run job by providing its jobid (url of results).
 Example :
 Returns : An instance of Bio::Tools::Run::PiseJob.

=cut

sub job {
    my ($self, $jobid) = @_;
    my $job = Bio::Tools::Run::PiseJob->job($jobid);
    return $job;
}

1;
