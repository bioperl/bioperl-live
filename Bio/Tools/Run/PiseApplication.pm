# $Id$
#

=head1 NAME

Bio::Tools::Run::PiseApplication

=head1 SYNOPSIS

  #

=head1 DESCRIPTION

A class to manage Pise programs information, configuring parameters
and submit jobs. It is the super-class of all the 
Bio::Tools::Run::PiseApplication::program classes.

This class is preferably created through the Bio::Factory::Pise factory:

  my $factory = Bio::Factory::Pise->new(-email => 'me@myhome');
  my $program = $factory->program('mfold');

By submitting a job, you create a Bio::Tools::Run::PiseJob instance with 
the parameters you have just set. Bio::Tools::Run::PiseJob class handles 
a specific job state and results.

  my $factory = Bio::Factory::Pise->new(-email => 'me@myhome');
  my $program = $factory->program('water', 
				    -sequencea => $seqa,
				    -seqall => $seqb);

  # run: submit and waits for completion
  my $job = $program->run();

  # for long jobs
  my $job = $program->submit(); # only submit the request
  my $jobid = $job->jobid;
  # later, from another script
  my $job = Bio::Tools::Run::PiseJob->fromUrl($jobid);
  if ($job->terminated) {
      print $job->stdout;
  }


=head2 Pise parameters setting.

The @params list should contain a list of -parameter =E<gt> value pairs.

  my @params = (-query => $file, -protein_db => "genpept");
  my $program = $factory->program('blast2', @params);

or directly :

  my $program = $factory->program('blast2', query => $file, 
     protein_db => "genpept");

Each program parameter is described in the documentation of the
corresponding Bio::Tools::Run::PiseApplication::program documentation.

You can change parameters at any time by calling the corresponding
method, e.g, changing the parameter E in blast2:

  my $program = $factory->program('blast2', -protein_db => "genpept");
  $program->query($seq);
  $program->Expect($value);

Parameter of Pise type "Sequence" and "InFile" may be given as string,
filename, or filehandle.  Parameter of type "Sequence" may also be
given as Bio::Seq or Bio::SimpleAlign objects.

=head2 Job output

See Bio::Tools::Run::PiseJob for how to fetch results and chain programs.

=head2 Remote and email parameters:

Email must be set at factory creation.

The remote parameter stands for the actual CGI location.  There are
default values for most of Pise programs.

You can either set remote at:

=over 3

=item 1 factory creation

  my $factory = Bio::Factory::Pise->new(
         -remote = 'http://somewhere/Pise/cgi-bin',
	 -email => 'me@myhome');

=item 2 program creation

  my $program = $factory->program('water',
	  -remote = 'http://somewhere/Pise/cgi-bin/water.pl');

=item 3 any time before running:

  $program->remote('http://somewhere/Pise/cgi-bin/water.pl');
  $program->run();

=back

=cut

#'

package Bio::Tools::Run::PiseApplication;

use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::Root::Root;
use Bio::Tools::Run::PiseJob;

@ISA = qw(Bio::Root::Root);

=head2 new

 Title   : new()
 Usage   : my $program = Bio::Tools::Run::PiseApplication->new($remote, $email);
 Function: Creates a Bio::Tools::Run::PiseApplication::program object, 
           where program stands for any 
           of the Pise programs.
           This method should not be used directly, but rather by 
           a Bio::Factory::Pise instance.
 Example :
 Returns : An instance of Bio::Tools::Run::PiseApplication::program.

=cut

sub new {
    my ($class, $remote, $email, $verbose) = @_;

    my $self = $class->SUPER::new;
    if (! defined $remote) {
	$self->throw(ref($self) . ": You must provide a Pise CGI url (-remote).");
    }
    $self->{REMOTE} = $remote;
    if (! defined $email) {
	$self->throw(ref($self) . ": You must provide a valid email.");
    }
    $self->{EMAIL} = $email;
    if (defined $verbose) {
	$self->{VERBOSE} = $verbose;
    } else {
	$self->{VERBOSE} = 0;
    }
    $self->{RESULTS_TYPE} = "url";
    return $self;
}

=head2 remote

 Title   : remote
 Usage   : my $remote = $program->remote;
 Function: Called from Bio::Tools::Run::PiseJob to get/set program/Pise 
           configuration informations from the 
           Bio::Tools::Run::PiseApplication::program class.
 Example :
 Returns : A string containing the url of the Pise server CGI.

=cut

sub remote {
    my $self = shift;
    if (@_) { 	$self->{REMOTE} = shift ;    }
    return $self->{REMOTE} ;
}

=head2 email

 Title   : email()
 Usage   : my $email = $program->email();
 Function: Called from Bio::Tools::Run::PiseJob to get/set program/Pise 
           configuration informations from the 
           Bio::Tools::Run::PiseApplication::program class.
 Example :
 Returns : A string containing the user email for submissions.

=cut

sub email {
    my $self = shift;
    if (@_) { 	$self->{EMAIL} = shift ;    }
    return $self->{EMAIL} ;
}

=head2 verbose

 Title   : verbose()
 Usage   : $program->verbose(1);
 Function: Ask the object to tells more.
 Example :
 Returns : Actual value.

=cut

sub verbose {
    my $self = shift;
    if (@_) { 	$self->{VERBOSE} = shift ;    }
    return $self->{VERBOSE} ;
}

=head2 param_type

 Title   : param_type()
 Usage   : my $type = $program->param_type($param);
 Function: Tells the Pise parameter type of $param (e.g: Sequence, 
           String, Excl, ...).
 Example :
 Returns : A string containing the Pise parameter type.

=cut

sub param_type {
    my $self = shift;
    my $param = shift;
    my $type;
    if ($param =~ /_data$/) {
	my $p;
	($p = $param) =~ s/_data//;
	$type = $self->type($p) ;
#	print STDERR "param: $param => $p type=$type\n" if ($self->{DEBUG});
    } else {
	$type = $self->type($param) ;
    }
    return $type;
}

=head2 run

 Title   : run()
 Usage   : $program->run();
           $program->run(10);
 Function: Submit the job and waits for completion. You may provide an
           interval for completion checking.
 Example :
 Returns : The instance of Bio::Tools::Run::Pisejob that has been run.

=cut

sub run {
    my ($self, @args) = @_;

    my ($remote) =
	$self->_rearrange([qw(REMOTE )],
			  @args);
    if (defined $remote) {
	$self->{REMOTE} = $remote;
    }

    my ($interval) =
	$self->_rearrange([qw(INTERVAL )],
			  @args);
    if (! defined $interval) {
	$interval = 10;
    }

    foreach my $param ($self->parameters_order) {
	my $param_name = $param;
	$param_name =~ tr/a-z/A-Z/;
	my ($value) =
	    $self->_rearrange([$param_name],
			      @args);
	if ($value) {
	    print STDERR "run: setting $param to $value\n" if $self->{VERBOSE};
	    $self->$param($value);
	}
    }

    my $pisejob = $self->submit;
    if (! defined $pisejob) {
	$self->throw(ref($self) . "::run: no job created");
    }
    if ($pisejob->error) {
	print STDERR ref($self) . "::run: error while submitting:", $pisejob->error_message,"\n" if $self->{VERBOSE};
	return $pisejob;
    }
	
    my $jobid = $pisejob->jobid();

    if ( ! ($pisejob->terminated) ) {
	$pisejob->results_type($self->{RESULTS_TYPE});
    }
    while ( ! ($pisejob->terminated) ) {
	print STDERR ref($self), "::run: waiting for completion...($jobid)\n" if $self->{VERBOSE};
	sleep $interval;
	last if ($pisejob->error);
    }
    
    $self->{_LASTJOBID} = $jobid;
    return $pisejob;
}

=head2 submit

 Title   : submit()
 Usage   : $program->submit();
 Function: Submit the job.
 Example :
 Returns : The instance of Bio::Tools::Run::Pisejob that has been run.

=cut

sub submit {
    my ($self, @args) = @_;

    my ($remote) =
	$self->_rearrange([qw(REMOTE )],
			  @args);
    if (defined $remote) {
	$self->{REMOTE} = $remote;
    }

    foreach my $param ($self->parameters_order) {
	my $param_name = $param;
	$param_name =~ tr/a-z/A-Z/;
	my ($value) =
	    $self->_rearrange([$param_name],
			      @args);
	if ($value) {
	    print STDERR "submit: setting $param to $value\n" if $self->{VERBOSE};
	    $self->$param($value);
	}
    }

    my $pisejob = Bio::Tools::Run::PiseJob->new($self, $self->{VERBOSE});

    if (! defined $pisejob) {
	$self->throw(ref($self) . "::submit: no job created");
    }
    if ($pisejob->error) {
	return $pisejob;
    }
    if ( ! ($pisejob->terminated) ) {
	$pisejob->results_type($self->{RESULTS_TYPE});
    }
    my $jobid = $pisejob->jobid();
    $self->{_LASTJOBID} = $jobid;

    print STDERR ref($self), "::submit: job running, url: $jobid\n" if $self->{VERBOSE};

    return $pisejob;
}

=head2 results_type

 Title   : results_type()
 Usage   : $program->results_type($type);
 Function: Enables to change result delivery from one email per file
           to url notification or attached files. $type is either: url, 
           attachment, email. This information will be provided to the job
           when detached and submitted through the run method.
 Example :
 Returns : 

=cut

sub results_type {
    my ($self, $type) = @_;
    $self->{RESULTS_TYPE} = $type;
    print STDERR ref($self), "::results_type: results type changed to: $type\n" if $self->{VERBOSE};
}


=head2 paraminfo

 Title   : paraminfo()
 Usage   : $program->paraminfo();
 Function: Displays parameters and prompts.
 Example :
 Returns : 

=cut

sub paraminfo {
    my $self = shift;
    my $prompt;
    foreach my $param ($self->parameters) {
	$prompt = $self->prompt($param);
	if ($prompt) {
	    print "$param\t\t$prompt\n";
	}
    }
    return;
}

=head2 _init_params

 Title   : _init_params
 Usage   : $self->_init_params(@params);
 Function: Internal. To be called from Pise::program::new method after
           all the data structures have been initialized.
 Example :
 Returns : 

=cut

sub _init_params {
    my $self = shift;
    my @params = @_;

    my ($param, $value);
    while (@params)  {
	$param =   shift @params;
	$value =  shift @params;
	next if( $param =~ /^-/ ); # don't want named parameters
	print STDERR "init_params $param to $value\n" if $self->{VERBOSE};
	$self->$param($value);	
    }
}

=head2 _OK_FIELD

 Title   : _OK_FIELD()
 Usage   : if ($self->_OK_FIELD($param)) ...
 Function: Checks if $param is a known parameter for the specific program.
 Example :
 Returns : TRUE/FALSE

=cut

sub _OK_FIELD {
    my ($self, $param) = @_;
    print STDERR "_OK_FIELD: $param\n" if $self->{DEBUG};
    if (grep /^$param$/, $self->parameters) {
	return 1;
    }
    if ($param =~ /(\w+)_data$/) {
	$param = $1;
	my $type = $self->param_type($param);
	if ($type eq "InFile" || $type eq "Sequence") {
	    return 1;
	}
    }

}

sub AUTOLOAD {
    my $self = shift;
    my $param = $AUTOLOAD;
    $param =~ s/.*:://;

    #print STDERR "AUTOLOAD: $param\n";
    $self->throw("Unallowed parameter or unknown procedure: $param !") unless $self->_OK_FIELD($param);
    $self->{$param} = shift if @_;
    print STDERR ref($self), "::$param: $param set to: ",$self->{$param}, "\n" if $self->{DEBUG};
    return $self->{$param};
}

#
# from here on, all the methods are accessors to Pise data structures.
#

sub parameters {
    my $self = shift;
    return ($self->parameters_order);
}

sub command {
    my $self = shift;
    return $self->{COMMAND} ;
}

sub program {
    my $self = shift;
    return ( $self->command ) ;
}

sub version {
    my $self = shift;
    return $self->{VERSION} ;
}

sub title {
    my $self = shift;
    if (@_) { 	$self->{TITLE} = shift ;    }
    return $self->{TITLE} ;
}

sub description {
    my $self = shift;
    if (@_) { 	$self->{DESCRIPTION} = shift ;    }
    return $self->{DESCRIPTION} ;
}

sub authors {
    my $self = shift;
    if (@_) { 	$self->{AUTHORS} = shift ;    }
    return $self->{AUTHORS} ;
}

sub doclink {
    my $self = shift;
    if (@_) { 	$self->{DOCLINK} = shift ;    }
    return $self->{DOCLINK} ;
}

sub reference {
    my $self = shift;
    return @{ $self->{REFERENCE} };
}

sub seqinput {
    my $self = shift;
    if (@_) { 	$self->{SEQINPUT} = shift ;    }
    return $self->{SEQINPUT} ;
}

sub seqtype {
    my $self = shift;
    if (@_) { 	$self->{SEQTYPE} = shift ;    }
    return $self->{SEQTYPE} ;
}

sub top_parameters {
    my $self = shift;
    return @{ $self->{TOP_PARAMETERS} };
}

sub parameters_order {
    my $self = shift;
    return @{ $self->{PARAMETERS_ORDER} };
}

sub by_group_parameters {
    my $self = shift;
    return @{ $self->{BY_GROUP_PARAMETERS} };
}

sub type {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{TYPE}{$param} = shift ;    }
        return  $self->{TYPE}{$param};
    } else {
        return  %{ $self->{TYPE} };
    }
}

sub format {
    my $self = shift;
    if (@_) { 
        my $param = shift;
        if (@_) {
            my $language = shift;
	    if (@_) { 	$self->{FORMAT}{$param}{$language} = shift ;    }
            return  $self->{FORMAT}{$param}{$language};
        } else {
	    if (defined $self->{FORMAT}{$param}) {
		return  %{ $self->{FORMAT}{$param} };
	    } else {
		return $self->{FORMAT}{$param};
	    }
        }
    } else {
        return  %{ $self->{FORMAT} };
    }
}

sub filenames {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{FILENAMES}{$param} = shift ;    }
        return  $self->{FILENAMES}{$param};
    } else {
        return  %{ $self->{FILENAMES} };
    }
}

sub seqfmt {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{SEQFMT}{$param} = shift ;    }
        return  @{ $self->{SEQFMT}{$param} };
    } else {
        return  %{ $self->{SEQFMT} };
    }
}

sub size {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{SIZE}{$param} = shift ;    }
        return  $self->{SIZE}{$param};
    } else {
        return  %{ $self->{SIZE} };
    }
}

sub group {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{GROUP}{$param} = shift ;    }
        return  $self->{GROUP}{$param};
    } else {
        return  %{ $self->{GROUP} };
    }
}

sub ishidden {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{ISHIDDEN}{$param} = shift ;    }
        return  $self->{ISHIDDEN}{$param};
    } else {
        return  %{ $self->{ISHIDDEN} };
    }
}

sub iscommand {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{ISCOMMAND}{$param} = shift ;    }
        return  $self->{ISCOMMAND}{$param};
    } else {
        return  %{ $self->{ISCOMMAND} };
    }
}

sub ismandatory {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{ISMANDATORY}{$param} = shift ;    }
        return  $self->{ISMANDATORY}{$param};
    } else {
        return  %{ $self->{ISMANDATORY} };
    }
}

sub isstandout {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{ISSTANDOUT}{$param} = shift ;    }
        return  $self->{ISSTANDOUT}{$param};
    } else {
        return  %{ $self->{ISSTANDOUT} };
    }
}

sub _interface_standout {
    my $self = shift;
    if (@_) { 	$self->{_INTERFACE_STANDOUT} = shift ;    }
    return  $self->{_INTERFACE_STANDOUT};
}

sub _standout_file {
    my $self = shift;
    if (@_) { 	$self->{_STANDOUT_FILE} = shift ;    }
    return  $self->{_STANDOUT_FILE};
}

sub prompt {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{PROMPT}{$param} = shift ;    }
        return  $self->{PROMPT}{$param};
    } else {
        return  %{ $self->{PROMPT} };
    }
}

sub vlist {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{VLIST}{$param} = shift ;    }
        return  @{ $self->{VLIST}{$param} };
    } else {
        return  %{ $self->{VLIST} };
    }
}

sub flist {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 
	    my $value = shift;
	    if (@_) { 	$self->{FLIST}{$param}{$value} = shift ;    }
	    return $self->{FLIST}{$param}{$value};
	} else {
	    if (defined $self->{FLIST}{$param}) {
		return  %{ $self->{FLIST}{$param} };
	    } else {
		return  $self->{FLIST}{$param} ;
	    }
	}
    } else {
        return  %{ $self->{FLIST} };
    }
}

sub separator {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{SEPARATOR}{$param} = shift ;    }
        return  $self->{SEPARATOR}{$param};
    } else {
        return  %{ $self->{SEPARATOR} };
    }
}

sub vdef {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{VDEF}{$param} = shift ;    }
	# may be a list (List type parameters)
	# must be casted by user
	return $self->{VDEF}{$param};
    } else {
        return  %{ $self->{VDEF} };
    }
}

sub precond {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 
            my $language = shift;
	    if (@_) { 	$self->{PRECOND}{$param}{$language} = shift ;    }
            return  $self->{PRECOND}{$param}{$language};
	} else {
	    if (defined $self->{PRECOND}{$param}) {
		return  %{ $self->{PRECOND}{$param} };
	    } else {
		return  $self->{PRECOND}{$param} ;
	    }
	}
    } else {
        return  %{ $self->{PRECOND} };
    }
}

sub ctrl {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 
	    my $language = shift;
	    if (@_) { 
		my $test = shift;
		if (@_) { 	
		    $self->{CTRL}{$param}{$language}{$test} = shift ; 
		}
		return $self->{CTRL}{$param}{$language}{$test};
	    } else {
		if (defined $self->{CTRL}{$param}{$language}) {
		    return %{ $self->{CTRL}{$param}{$language} };
		} else {
		    return $self->{CTRL}{$param}{$language};
		}
	    }
	} else {
	    if (defined $self->{CTRL}{$param}) {
		return %{ $self->{CTRL}{$param} };
	    } else {
		return $self->{CTRL}{$param};
	    }
	}
    } else {
        return  %{ $self->{CTRL} };
    }
}

sub pipeout {
    my $self = shift;
    if (@_) { 
        my $param = shift;
        if (@_) {
            my $test = shift;
	    if (@_) { 	$self->{PIPEOUT}{$param}{$test} = shift ;    }
            return  $self->{PIPEOUT}{$param}{$test} ;
	} else {
	    if (defined $self->{PIPEOUT}{$param}) {
		return  %{ $self->{PIPEOUT}{$param} };
	    } else {
		return  $self->{PIPEOUT}{$param};
	    }
	}
    } else {
        return  %{ $self->{PIPEOUT} };
    }
}

sub withpipeout {
    my $self = shift;
    if (@_) { 
        my $param = shift;
        if (@_) {
            my $test = shift;
	    if (@_) { 	$self->{WITHPIPEOUT}{$param}{$test} = shift ;    }
            return  @{ $self->{WITHPIPEOUT}{$param}{$test} };
	} else {
	    if (defined $self->{WITHPIPEOUT}{$param}) {
		return  %{ $self->{WITHPIPEOUT}{$param} };
	    } else {
		return  $self->{WITHPIPEOUT}{$param} ;
	    }
	}
    } else {
        return  %{ $self->{WITHPIPEOUT} };
    }
}

sub pipein {
    my $self = shift;
    if (@_) { 
        my $param = shift;
        if (@_) {
            my $type = shift;
	    if (@_) { 	$self->{PIPEIN}{$param}{$type} = shift ;    }
            return   $self->{PIPEIN}{$param}{$type} ;
	} else {
	    if (defined $self->{PIPEIN}{$param}) {
		return  %{ $self->{PIPEIN}{$param} };
	    } else {
		return  $self->{PIPEIN}{$param};
	    }
	}
    } else {
        return  %{ $self->{PIPEIN} };
    }
}

sub withpipein {
    my $self = shift;
    if (@_) { 
        my $param = shift;
        if (@_) {
            my $type = shift;
	    if (@_) { 	$self->{WITHPIPEIN}{$param}{$type} = shift ;    }
            return  @{ $self->{WITHPIPEIN}{$param}{$type} };
	} else {
	    if (defined $self->{WITHPIPEIN}{$param}) {
		return  %{ $self->{WITHPIPEIN}{$param} };
	    } else {
		return  $self->{WITHPIPEIN}{$param} ;
	    }
	}
    } else {
        return  %{ $self->{WITHPIPEIN} };
    }
}

sub isclean {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{ISCLEAN}{$param} = shift ;    }
        return  $self->{ISCLEAN}{$param};
    } else {
        return  %{ $self->{ISCLEAN} };
    }
}

sub issimple {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{ISSIMPLE}{$param} = shift ;    }
        return  $self->{ISSIMPLE}{$param};
    } else {
        return  %{ $self->{ISSIMPLE} };
    }
}

sub paramfile {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{PARAMFILE}{$param} = shift ;    }
        return  $self->{PARAMFILE}{$param};
    } else {
        return  %{ $self->{PARAMFILE} };
    }
}

sub comment {
    my $self = shift;
    if (@_) { 
        my $param = shift;
        if (@_) {       
	    $self->{COMMENT}{$param} = [ @{$_[0]} ] ;    
	}
	if (defined $self->{COMMENT}{$param} ) {
            return  @{ $self->{COMMENT}{$param} };
	} else {
	    return $self->{COMMENT}{$param}
	}
    } else {
        return  %{ $self->{COMMENT} };
    }
}


sub scalemin {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{SCALEMIN}{$param} = shift ;    }
        return  $self->{SCALEMIN}{$param};
    } else {
        return  %{ $self->{SCALEMIN} };
    }
}

sub scalemax {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{SCALEMAX}{$param} = shift ;    }
        return  $self->{SCALEMAX}{$param};
    } else {
        return  %{ $self->{SCALEMAX} };
    }
}

sub scaleinc {
    my $self = shift;
    if (@_) { 
        my $param = shift;
	if (@_) { 	$self->{SCALEINC}{$param} = shift ;    }
        return  $self->{SCALEINC}{$param};
    } else {
        return  %{ $self->{SCALEINC} };
    }
}


1;


