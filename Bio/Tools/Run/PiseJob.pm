# $Id$
#

=head1 NAME

Bio::Tools::Run::PiseJob

=head1 SYNOPSIS

=head1 DESCRIPTION

  Bio::Tools::Run::PiseJob class handles a specific job state and results.
  A Bio::Tools::Run::PiseJob instance should be created by a subclass of 
  Bio::Tools::Run::PiseApplication  class, e.g
  Bio::Tools::Run::PiseApplication::genscan or 
  Bio::Tools::Run::PiseApplication::dnapars, ... (see 
  Bio::Tools::Run::PiseApplication class) :

    my $job = Bio::Tools::Run::PiseJob->new($self, $self->{VERBOSE});

  This class may also be used as a mean to get informations
  about a running job, or to get results after a long computation:

    my $job = Bio::Factory::Pise->job($url);
    print $job->content('infile.aln');

  Once the job is created, you can get results:

    foreach my $result ($job->get_results) {
	print $job->content($result);
	$job->save($result, "myfile"); # $job->save($result) keeps the name
	print $job->stdout;            # print job standard output
	print $job->stderr;            # print job standard error
    }

  You can feed a result file as a filehandle to a bioperl parser :

     my $parser = Bio::Tools:Genscan->new (-fh => $job->fh('genscan.out'));
     my $parser = Bio::Tools:BPlite->new (-fh => $job->fh('blast2.txt'));

   ... or to another pise job:

     my $neighbor = $factory->program ('neighbor',
					infile => $job->fh('outfile'));

  You can lookup up for a type of result that could be piped to another
  Pise program:

     my $matrix = $job->lookup_piped_file('phylip_dist');

  returns the url of the just calculated Phylip distances matrix file, 
  produced by e.g DNADIST or PROTDIST.

  All the available pipe types may be obtained by:
    $job->lookup_piped_files;

=cut

#'

package Bio::Tools::Run::PiseJob;

use vars qw(@ISA);
use strict;
use Bio::Root::Root;
use Bio::AlignIO;
use Bio::Tools::Run::PiseJobParser;
use XML::Parser::PerlSAX;
use LWP::UserAgent;
use HTTP::Request::Common;
use POSIX;

@ISA = qw(Bio::Root::Root);


=head2 new

 Title   : new()
 Usage   : $job = Bio::Tools::Run::PiseJob->new($application, $verbose);
 Function: Creates a Bio::Tools::Run::PiseJob object. 
           This is normally called by an application object
           - i.e a subclass of the Bio::Tools::Run::PiseApplication class, 
           for submitting a job. 
           This method actually submit the job and parse results.
 Example :
 Returns : An instance of Bio::Tools::Run::PiseJob.

=cut

sub new {
    my ($class, $application, $verbose) = @_;
    my $self = $class->SUPER::new();

    $self->{APPLICATION} = $application;
    $self->{VERBOSE} = $verbose;
    $self->{DEBUG} = 0;

    $self->_init;
    $self->_submit;
    return $self;
}

=head2 verbose

 Title   : verbose()
 Usage   : $program->verbose(1);
 Function: Ask the object to tells more.
 Example :
 Returns : 

=cut

sub verbose {
    my $self = shift;
    if (@_) { 	$self->{VERBOSE} = shift ;    }
    return $self->{VERBOSE} ;
}

=head2 job

 Title   : job()
 Usage   : $job = Bio::Tools::Run::PiseJob->job(url);
 Function: Creates a Bio::Tools::Run::PiseJob object from an already 
           run job by giving the url of the job result page.
           May also be called through Bio::Factory::Pise->job(url);
 Example :
 Returns : An instance of Bio::Tools::Run::PiseJob.

=cut

sub job {
    my ($class, $url, $verbose) = @_; 
    my $self = Bio::Tools::Run::PiseJob->SUPER::new();

    $self->{JOBID} = $url;
    $self->{VERBOSE} = $verbose;
    $self->{ERROR} = undef;
    $self->{ERROR_MESSAGE} = undef;
    $self->{TERMINATED} = 0;
    $self->{RESULT_FILES} = undef;
    $self->{RESULTS} = undef;
    $self->{SCRATCH_DIR} = undef;
    $self->{DEBUG} = 0;
    $self->{PIPES} = {};
    $self->{TMPFILES} = [];
    $self->{PIPED_FILE_TYPE} = {};

    my $ua = $self->_get_ua;
    my $res = $ua->request(GET $self->{JOBID});

    if ($res->is_success) {
	$self->{RESULTS} = $res->content;
	if ($self->_parse($res->content) < 0) {
	    $self->{ERROR} = 1;
	    $self->{ERROR_MESSAGE} = ref($self) . " _fromUrl: parsing error";
	}
    } else {
	$self->{ERROR} = 1;
	$self->{ERROR_MESSAGE} = ref($self) . " _fromUrl: " . $res->message;
	$self->throw(ref($self) . " _fromUrl: " . $res->message);
    }

    return $self;
}

=head2 jobid

 Title   : jobid()
 Usage   : $job->jobid();
 Function: Returns the url of the job result page.
 Example :
 Returns : 

=cut

sub jobid {
    my $self = shift;
    return $self->{JOBID};
}

=head2 error

 Title   : error()
 Usage   : $job->error();
 Function: Tells if the job has been successfully run. This is the case 
           when the job has been submitted, but the Pise server has 
           detected user errors (missing mandatory parameter, unallowed 
           value,...). This also happen when the user provided an
           invalid url, or the http request could not be submitted.
           See method error_message().

 Example :
 Returns : TRUE/FALSE

=cut

sub error {
    my $self = shift;
    return $self->{ERROR};
}

=head2 error_message

 Title   : error_message()
 Usage   : $job->error_message();
 Function: Returns the error message.
 Example :
 Returns : A string.

=cut

sub error_message {
    my $self = shift;
    return $self->{ERROR_MESSAGE};
}

=head2 get_results

 Title   : get_results()
 Usage   : $job->get_results();
 Function: Provides the urls of the result files.
 Example :
 Returns : A list of urls.

=cut

sub get_results {
    my $self = shift;

    if (! $self->{JOBID}) {
	$self->throw("Bio::Tools::Run::PiseJob::get_results: your job has no jobid");
    }
    if (! $self->{TERMINATED}) {
	$self->throw("Bio::Tools::Run::PiseJob::get_results: your job is not terminated");
    }

    return @{ $self->{RESULT_FILES} };
}

=head2 get_pipes

 Title   : get_pipes()
 Usage   : $job->get_pipes($result);
 Function: Provides the names of the programs that can use this type of
           result. $result is an url, that can be provided through the
           get_results method.
 Example :
 Returns : A list of program names.

=cut

sub get_pipes {
    my $self = shift;
    my $result_file = shift;

    if (! $self->{JOBID}) {
	$self->throw("Bio::Tools::Run::PiseJob::get_pipes: your job has no jobid");
    }
    if (! $self->{TERMINATED}) {
	$self->throw("Bio::Tools::Run::PiseJob::get_pipes: your job is not terminated");
    }

    my %pipes = %{ $self->{PIPES}};
    if (defined $pipes{$result_file}) {
	my @pipes = @{ $pipes{$result_file} };
	return @pipes;
    } else {
	return undef;
    }
}

=head2 piped_file_type

 Title   : piped_file_type()
 Usage   : $job->piped_file_type($result);
 Function: Provides the Pise type of $result. $result is an url, 
           that can be provided through the get_results method.
 Example :
 Returns : A Pise pipetype name.

=cut

sub piped_file_type {
    my $self = shift;
    my $result_file = shift;
    if (! $self->{JOBID}) {
	$self->throw("Bio::Tools::Run::PiseJob::piped_file_type: your job has no jobid");
    }
    if (! $self->{TERMINATED}) {
	$self->throw("Bio::Tools::Run::PiseJob::piped_file_type: your job is not terminated");
    }
    return $self->{PIPED_FILE_TYPE}{$result_file};
}

=head2 lookup_piped_files

 Title   : lookup_piped_files()
 Usage   : $pipe_types = $job->lookup_piped_files();
 Function: Returns the pipe types produced by the job
           (e.g:  phylip_tree, seqsfile, readseq_ok_alig, ...). 
           You have to call lookup_piped_file($type) to get the actual
           correponding result file.
 Example :
 Returns : A string.

=cut

sub lookup_piped_files {
    my $self = shift;
    my $pipe_type = shift;
    if (! $self->{JOBID}) {
	$self->throw(ref($self) . " lookup_piped_files: your job has no jobid");
    }
    if (! $self->{TERMINATED}) {
	$self->throw(ref($self) . " lookup_piped_files: your job is not terminated");
    }
    return (values %{$self->{PIPED_FILE_TYPE}});

}

=head2 lookup_piped_file

 Title   : lookup_piped_file(type)
 Usage   : $result = $job->lookup_piped_file($type);
 Function: Returns the name of the result file of pipe type $type 
           (e.g:  phylip_tree, seqsfile, readseq_ok_alig, ...). $result 
           is an url.
 Example :
 Returns : A string (an url).

=cut

sub lookup_piped_file {
    my $self = shift;
    my $pipe_type = shift;
    if (! $self->{JOBID}) {
	$self->throw("Bio::Tools::Run::PiseJob::lookup_piped_file: your job has no jobid");
    }
    if (! $self->{TERMINATED}) {
	$self->throw("Bio::Tools::Run::PiseJob::lookup_piped_file: your job is not terminated");
    }
    foreach my $result (@{ $self->{RESULT_FILES}} ) {
	if ($self->{PIPED_FILE_TYPE}{$result} eq $pipe_type) {
	    return $result;
	}
    }
}

=head2 terminated

 Title   : terminated()
 Usage   : $job->terminated();
 Function: Tells whether the job has terminated.
 Example :
 Returns : TRUE/FALSE.

=cut

sub terminated {
    my $self = shift;
    my $jobid = shift;
    if (! defined $jobid) {
	$jobid = $self->{JOBID};
    }
    if  (! defined $jobid) {
	$self->{ERROR} = 1;
	$self->{ERROR_MESSAGE} = ref($self) . " terminated: no jobid?";
    }
    my $ua = $self->_get_ua;

    my $res = $ua->request(GET $jobid);

    if ($res->is_success) {
	$self->{RESULTS} = $res->content;
	if ($self->_parse($res->content) < 0) {
	    $self->{ERROR} = 1;
	    $self->{ERROR_MESSAGE} = ref($self) . " terminated: parsing error";
	}
    } else {
	$self->{ERROR} = 1;
	$self->{ERROR_MESSAGE} = ref($self) . " terminated: " . $res->message;
	$self->throw(ref($self) .  " terminated: " . $res->message);
    }

    if ($self->{TERMINATED}) {
	return 1;
    }

    return 0;
}

=head2 save

 Title   : save()
 Usage   : $filename = $job->save($result);
           $filename = $job->save($result, $name);
 Function: Save the result in a file. $result is an url, 
           that can be provided through the get_results method. You can
           provide your own filename. By default, the file name will be
           the same as the result name.
 Example :
 Returns : A file name.

=cut

sub save {
    my $self = shift;
    my $jobid;
    my $url;
    my $file;
    my $result;

    if (! $self->{JOBID}) {
	$self->throw("Bio::Tools::Run::PiseJob::save: your job has no jobid");
    }
    if (! $self->{TERMINATED}) {
	$self->throw("Bio::Tools::Run::PiseJob::save: your job is not finished");
    }
    if (@_) {
	$result = shift;
    } else {
	$self->throw("Bio::Tools::Run::PiseJob::save: you must provide the result url");
    }

    my $tmp_url = $self->{JOBID};
    if (@_) {
	$file = shift;
    } else {
	$file = $result;
	if ($file =~ /http/) {
	    $file =~ s/$tmp_url//;
	    if (defined $self->{PROGRAM}) {
		my $cmd = $self->{PROGRAM};
		$file =~ s/$cmd//;
		$file =~ s/\w?\d+\///;
		$file =~ s/\///g;
	    } else {
		$file =~ s/\w+\/\w?\d+\///;
		$file =~ s/\///g;
	    }
	}
    }

    my $ua = $self->_get_ua;

    foreach $url (@{ $self->{RESULT_FILES}}) {
	if ($self->{DEBUG}) {
	    print STDERR ref($self), "::save: $url (",$self->{PROGRAM},")\n";
	}
	if ($url =~ /$result/) {
	    my $res = $ua->request(GET $url);
	    
	    if ($res->is_success) {
		open(FILE,"> $file") || die "cannot open $file: $!";
		print FILE $res->content;
		close FILE;
		return $file;
	    } else {
		$self->{ERROR} = 1;
		$self->{ERROR_MESSAGE} = ref($self) . " save: " . $res->message;
		$self->throw(ref($self) . " save: " . $res->message);
	    }
	}
    }
}

=head2 content

 Title   : content()
 Usage   : $s = $job->content($result);
 Function: Provides the content of $result. $result is an url, 
           that can be provided through the get_results method. 
           By default, $result is the standard output.
 Example :
 Returns : A string.

=cut

sub content {
    my $self = shift;
    my $jobid;
    my $url;
    my $file;

    if (@_) {
	$file = shift;
    } else {
	$file = $self->{PROGRAM} . ".out";
    }

    if (! $self->{JOBID}) {
	$self->throw("Bio::Tools::Run::PiseJob::content: your job has no jobid");
    }
    if (! $self->{TERMINATED}) {
	$self->throw("Bio::Tools::Run::PiseJob::content: your job is not terminated");
    }

    my $ua = $self->_get_ua;

    foreach $url (@{ $self->{RESULT_FILES}}) {
	if ($self->{DEBUG}) {
	    print STDERR ref($self) . " content: $url (",$self->{PROGRAM},")\n";
	}
	if ($url =~ /$file/) {
	    if ($self->{DEBUG}) {
		print STDERR ref($self) . " content: this one!\n";
	    }
	    my $res = $ua->request(GET $url);
	    
	    if ($res->is_success) {
		return $res->content;
	    } else {
		$self->{ERROR} = 1;
		$self->{ERROR_MESSAGE} = ref($self) . " content: " . $res->message;
		$self->throw(ref($self) . " content: " . $res->message);
	    }
	}
    }
}

=head2 stdout

 Title   : stdout()
 Usage   : print $job->stdout();
 Function: Provides the content of the job standard output. 
 Example :
 Returns : A string.

=cut

sub stdout {
    my $self = shift;

    if (! $self->{JOBID}) {
	$self->throw("Bio::Tools::Run::PiseJob::stdout: your job has no jobid");
    }
    if (! $self->{TERMINATED}) {
	$self->throw("Bio::Tools::Run::PiseJob::stdout: your job has not terminated");
    }

    return $self->content($self->{PROGRAM} . ".out");

}

sub output {
    my $self = shift;
    return($self->stdout);
}

=head2 stderr

 Title   : stderr()
 Usage   : print $job->stderr();
 Function: Provides the content of the job standard error. 
 Example :
 Returns : A string.

=cut

sub stderr {
    my $self = shift;

    if (! $self->{JOBID}) {
	$self->throw("Bio::Tools::Run::PiseJob::stderr: your job has no jobid");
    }
    if (! $self->{TERMINATED}) {
	$self->throw("Bio::Tools::Run::PiseJob::stderr: your job has not terminated");
    }

    return $self->content($self->{PROGRAM} . ".err");

}

=head2 fh

 Title   : fh()
 Usage   : $fh = $job->fh($result);
 Function: Provides a filhandle for a result.
           $result is an url, that can be provided through the 
           get_results method. 

           Be aware that you must re-ask for it for a subsequent use. For
           instance, if you first use it for an input parameter:
             my $program = Pise::program->new ( ...,
		  			      file => $previous_job->fh('..'),
					      );
             my $job = $program->run;

           A subsequent run of the same object: will need a re-initialization:
             $program->file($previous_job->fh('..'));
             my $job2 = $program->run;

 Example :
 Returns : A filehandle.

=cut

sub fh {
    my $self = shift;
    my $jobid;
    my $url;
    my $file;

    if (! $self->{JOBID}) {
	$self->throw("Bio::Tools::Run::PiseJob::fh: your job has no jobid");
    }
    if (! $self->{TERMINATED}) {
	$self->throw("Bio::Tools::Run::PiseJob::fh: your job has not terminated");
    }
    if (@_) {
	$file = shift;
    } else {
	$file = $self->{PROGRAM} . ".out";
    }

    my $ua = $self->_get_ua;

    foreach $url (@{ $self->{RESULT_FILES}}) {
	if ($self->{DEBUG}) {
	    print STDERR "DEBUG> Bio::Tools::Run::PiseJob fh: $url (",$self->{PROGRAM},")\n";
	}
	if ($url =~ /$file/) {
	    if ($self->{DEBUG}) {
		print STDERR "Bio::Tools::Run::PiseJob::fh: this one ($file)!\n";
	    }
	    my $res = $ua->request(GET $url);
	    
	    if ($res->is_success) {
		@{ $self->{FH_DATA} } = split( "\n", $res->content);
		my $class = ref($self) || $self;
		my $s = Symbol::gensym;
		tie $$s,$class,$self;
		return $s;
	    } else {
		$self->{ERROR} = 1;
		$self->{ERROR_MESSAGE} = ref($self) . " fh: " . $res->message;
		$self->throw(ref($self) . " fh: " . $res->message);
	    }
	}
    }


}

=head2 results_type

 Title   : results_type()
 Usage   : $job->results_type($type);
 Function: Enables to change result delivery from one email per file
           to url notification or attached files. $type is either: url, 
           attachment, email. 
 Example :
 Returns : 1 if success, 0 if job already terminated. 

=cut

sub results_type {
    my $self = shift;
    my $results_type;
    if (! $self->{JOBID}) {
	$self->throw("Bio::Tools::Run::PiseJob::results_type: your job has no jobid");
    }
    if ($self->{TERMINATED}) {
	print STDERR "Bio::Tools::Run::PiseJob::results_type: job already terminated\n" if $self->{VERBOSE};
	return 0;
    }
    if (@_) {
	$results_type = shift;
    } else {
	$results_type = "url";
    }

    my $jobid = $self->{JOBID};
    my $application = $self->{APPLICATION};
    
    my $scratch_dir = (defined $self->{SCRATCH_DIR}) ? $self->{SCRATCH_DIR} : "" ;
    my $command = $application->program;
    if ($scratch_dir eq "") {
	($scratch_dir = $jobid) =~ s/http.+\/(\w?\d+)/$1/;
	$scratch_dir =~ s/index.html//;
	$scratch_dir = "$command/$scratch_dir";
    }

    my $ua = $self->_get_ua;
    
    my $remote = $self->{REMOTE};
    $remote =~ s/$command\.pl//;
    $remote .= "lib/results.pl";
    print STDERR "Bio::Tools::Run::PiseJob::results_type: running $remote to change results type ($results_type scratch_dir: $scratch_dir)\n" if $self->{VERBOSE};

    my $res = $ua->request(POST $remote, [command => $command, email => $ENV{'USER'}, results_type => $results_type, scratch_dir => $scratch_dir]);

    if ($res->is_success) {
	return 1;
    } else {
	$self->throw("Bio::Tools::Run::PiseJob::results_type: " . $res->message);
    }
}

=head2 value

 Title   : value(param)
 Usage   : $job->value(param);
 Function: 
 Example :
 Returns : value of parameter param, if available.

=cut

sub value {
    my $self = shift;
    my $param;
    if (! $self->{JOBID}) {
	$self->throw("Bio::Tools::Run::PiseJob::value: your job has no jobid");
    }
    if (! $self->{TERMINATED}) {
	$self->throw("Bio::Tools::Run::PiseJob::value: the job has not terminated");
    }
    if (@_) {
	$param = shift;
    } else {
	return;
    }

    if (exists $self->{VALUE}{$param}) {
	return $self->{VALUE}{$param};
    }
}

=head2 _init

 Title   : _init()
 Usage   : $self->_init;
 Function: Internal. Initializes parameters. Called by new.
 Example :
 Returns : 

=cut

sub _init {
    my $self = shift;
    my $application = $self->{APPLICATION};

    $self->{PROGRAM} = $application->program;
    $self->{REMOTE} = $application->remote;
    $self->{EMAIL} = $application->email;
    $self->{JOBID} = undef;
    $self->{ERROR} = undef;
    $self->{ERROR_MESSAGE} = undef;
    $self->{TERMINATED} = 0;
    $self->{ARGS} = undef;
    $self->{RESULT_FILES} = undef;
    $self->{RESULTS} = undef;
    $self->{SCRATCH_DIR} = undef;
    $self->{PIPES} = {};
    $self->{PIPED_FILE_TYPE} = {};
    $self->{UA} = undef;

    foreach my $param ($application->parameters) { 
	my $value;
	print STDERR "Bio::Tools::Run::PiseJob::_init param type: ", $application->param_type($param), "\n" if $self->{DEBUG};
	$value = $application->$param();
	if (defined $value) {
	    print STDERR "Bio::Tools::Run::PiseJob::_init param value: $value, ref: ",ref($value),"\n" if $self->{DEBUG};

	    if ($application->param_type($param) eq "Sequence" || $application->param_type($param) eq "InFile") {
		if (ref($value)) {
		    print STDERR ref($self), "::_init: ",ref($value), "\n" if $self->{DEBUG};
		    if (ref($value) eq "GLOB" || $value->isa('IO::Handle')) {
			print STDERR "Bio::Tools::Run::PiseJob::_init got filehandle ",ref($value),"\n" if $self->{DEBUG};
			while (<$value>) {
			    $self->{ARGS}{$param . "_data"} .= $_;
			}
		    } elsif ($application->param_type($param) eq "Sequence" && 
			$value->isa("Bio::PrimarySeqI")) {
			$self->{ARGS}{$param . "_data"} = $value->seq;
		    } elsif ($application->param_type($param) eq "Sequence" && 
			     $value->isa("Bio::SimpleAlign")) {
			#my $tmpfile = POSIX::tmpnam;
			my $tmpfile = $param . ".fasta";
			# bioperl 1.0
			my $out = Bio::AlignIO->new(-file => ">$tmpfile", '-format' => 'fasta');
			$out->write_aln($value);
			close(TMP);
			push (@{$self->{TMPFILES}}, $tmpfile);
			print STDERR "Bio::Tools::Run::PiseJob::_init written alignment to $tmpfile\n" if $self->{VERBOSE};
			$self->{ARGS}{$param} = $tmpfile;
		    }
		} else {
		    if (ref(\$value) eq "SCALAR" && -f $value) {
			$self->{ARGS}{$param} = $value;
			print STDERR "Bio::Tools::Run::PiseJob::_init got file ($value)\n" if $self->{DEBUG};
		    } else {
			$self->{ARGS}{$param . "_data"} = $value;
		    }
		}
	    } else {
		$self->{ARGS}{$param} = $value;
	    }
	}
    }

    $self->{ARGS}{'email'} = $self->{EMAIL};
}

=head2 _submit

 Title   : _submit()
 Usage   : $self->_submit();
 Function: Internal. Sends the http request on a Pise server. Called by new.
 Example :
 Returns : -1 if an error has occured
           jobid else
 Exceptions: when the job has already been submitted.

=cut

sub _submit {

    my $self = shift;

    if (defined $self->{JOBID}) {
	print STDERR ref($self) . " submit: this job has been already setup and launched\n";
	$self->{ERROR} = 1;
	$self->{ERROR_MESSAGE} = ref($self) . " _submit: this job has been already setup and launched";
	$self->throw(ref($self) . " _submit: this job has been already setup and launched");

    }

    my $remote = $self->{REMOTE};

    my $application = $self->{APPLICATION};
    my $type;
    my $value;
    my $vdef;
    my @content;

    foreach my $param (keys %{ $self->{ARGS} }) {
	$type = $application->param_type($param);
        $value = $self->{ARGS}{$param};
	if ($type eq "InFile" || $type eq "Sequence") {
	    if ($param !~ /_data$/) {
		stat($value);
		if (-e _) {
		    push (@content, $param => [$value]);
		    print STDERR "_submit(1): $param: file $value\n" if ($self->{DEBUG});
		} else {
		    push (@content, $param => $value);
		    print STDERR "_submit(1): $param: not file (1)\n" if ($self->{DEBUG});
		}
	    } else {
		push (@content, $param => $value);
		print STDERR "_submit(1): $param: not file ($value)(2)\n" if ($self->{DEBUG});
	    }
	} elsif ($type eq "Switch") {
	    if ($value) {
		push (@content, $param => "on");
	    }
	} elsif ($type eq "List") {
	    foreach my $v (@{ $value }) {
		push (@content, $param => $v);
	    }
	} else {
	    push (@content, $param => $value);
	}
	
#	print STDERR "$param ($type): $content{$param}\n" if ($self->{DEBUG});;
    }
    
    # dealing with default values
    # they are more or less assumed by the Pise system, so it's better to
    # fill them
    foreach my $param ($application->parameters_order) {
	$type = $application->param_type($param) ;
	if (! defined $self->{ARGS}{$param}) {
	    $vdef = $application->vdef($param) ;
	    if ($vdef && $vdef ne "\"\"") {
		if ($type eq "Switch") {
		    push (@content, $param => "on");
		} elsif ($type eq "List") {
		    foreach my $v (@{ $vdef }) {
			push (@content, $param => $v);
		    }
		} else {
		    print STDERR "_submit(2): setting $param to vdef $vdef\n" if $self->{DEBUG};
		    push (@content, $param => $vdef);
		}
	    }
	}
    }

    if ($self->{DEBUG}) {
	my $i;
	for ($i=0; $i <= scalar(@content); $i++) {
	    print STDERR "PiseJob _submit(3): $content[$i]\n";
	}
    }

    print STDERR ref($self), "::_submit: submitting request ($remote)...\n"  if $self->{VERBOSE};

    my $ua = $self->_get_ua;

    my $res = $ua->request(POST $remote,
			   Content_Type => 'form-data', 
			   Content      => \@content);

    foreach my $tmpfile (@{ $self->{TMPFILES}} ) {
	print STDERR "removing $tmpfile\n" if $self->{VERBOSE};
	unlink $tmpfile;
    }

    if ($res->is_success) {
#	if ($self->{DEBUG}) {
#	    print STDERR "submit:\n", $res->content;
#	}
	$self->{RESULTS} = $res->content;
	if ($self->_parse($res->content) >= 0) {
	    return $self->jobid;
	} else {
	    print STDERR ref($self) . " _submit: parse error, result content: " . $res->content if $self->{VERBOSE};
 	    return $self->jobid;
	}
    } else {
	$self->{ERROR} = 1;
	$self->{ERROR_MESSAGE} = ref($self) . " _submit: " . $res->message;
	$self->{TERMINATED} = 1;
	return -1;
    }
}

=head2 _parse

 Title   : _parse()
 Usage   : $self->_parse();
 Function: Internal. Parses Pise XHTML results page and fills data structures.
           Called by frmoUrl or by _submit.
 Example :
 Returns : 

=cut

sub _parse {
    my $self = shift;
    my $content;
    if (@_) {
	$content = shift;
    } elsif (defined $self->{RESULTS}) {
	$content = $self->{RESULTS};
    } else {
	print STDERR "parse: you must provide the REMOTE results page\n";
	return -1;
    }
    my $handler;
    if ($self->{VERBOSE}) {
 	$handler = Bio::Tools::Run::PiseJobParser->new(1);
   } else {
	$handler = Bio::Tools::Run::PiseJobParser->new;
    }
    my $parser = XML::Parser::PerlSAX->new  (Handler => $handler);
    $self->{PARSER} = $parser;
    my $content = $self->_clean_content($content);

    eval {$parser->parse($content)};

    if ($@) {
	print STDERR "parse: cannot parse this job:\n$@\n";
	print STDERR $content;
	return -1;
    } else {
	if (! $self->{JOBID}) {
	    $self->{JOBID} = $handler->bioweb_result;
	}
	$self->{SCRATCH_DIR} = $handler->scratch_dir;
	my @results_files = $handler->hrefs;
	$self->{RESULT_FILES} = [@results_files];
	foreach my $result (@results_files) {
	    $self->{PIPED_FILE_TYPE}{$result} = $handler->piped_file_type($result);
	}
	my %pipes = $handler->pipes;
        if (defined %pipes) {
	    foreach my $f (keys %pipes) {
		if (defined $pipes{$f}) {
		    my @p = @{ $pipes{$f} };
		    foreach my $p (@p) {
			push (@{$self->{PIPES}{$f}}, $p);
		    }
		}
	    }
	}

	# parameters hidden values
	;
	foreach my $param (keys %{ $handler->{value}}) {
	    $self->{VALUE}{$param} = $handler->{value}{$param};
	    #print STDERR "DEBUG> Bio::Tools::Run::PiseJob _parse: $param => ", $self->{VALUE}{$param}, "\n" if $self->{VERBOSE};
	}

	$self->{TERMINATED} = $handler->terminated;
	if ($handler->error) {
	    $self->{ERROR} = $handler->error;
	    $self->{ERROR_MESSAGE} = $handler->error_message;
	    print STDERR ref($self) . " _parse: an error has occured (", $self->{PROGRAM}, ") : ",$handler->error_message, "\n" if $self->{VERBOSE};
	    return -1;
	}
    }

}

sub _get_ua {
    my $self = shift;
    my $ua;
    if (defined $self->{UA}) {
	$ua = $self->{UA};
    } else {
	$ua = LWP::UserAgent->new;
	$ua->agent("Pise/" . $self->{VERSION} . $ua->agent);
	$self->{UA} = $ua;
    }
    return $ua;
}

=head2 READLINE

 Title   : READLINE()
 Usage   : 
 Function: Internal - see perltie.
 Example :
 Returns : A string.

=cut

sub READLINE {
  my $self = shift;
  if (scalar(@{ $self->{pisejob}->{FH_DATA} }) > 0) {
      my $line = shift @{ $self->{pisejob}->{FH_DATA} };
      return "$line\n";
  } else {
      return undef;
  }
}

=head2 TIEHANDLE

 Title   : TIEHANDLE()
 Usage   : 
 Function: Internal - see perltie.
 Example :
 Returns : 

=cut

sub TIEHANDLE {
  my $class = shift;
  return bless {pisejob => shift}, $class;
}

=head2 _clean_content

 Title   : _clean_content()
 Usage   : my $content = $self->_clean_content($content);
 Function: Internal. Useful to call before XML parsing.
 Example :
 Returns :

=cut

sub _clean_content {
    my $self = shift;
    my $content = shift;

    $content =~ s/\&/\&amp;/g;
#    $content =~ s/</&lt;/g;
#    $content =~ s/>/&gt;/g;
    my $title;
    my $head;
    my $foot;
#    if ($content !~ /<\?xml/) {
#	$head = "<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n";
#    }
    if ($content !~ /DOCTYPE/) {
	$head .= "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\"
    \"http://www.w3.org/TR/xhtml1/DTD/strict.dtd\">\n";
    }
    if ($content !~ /<html>/i) { 
	if (defined $self->{APPLICATION}) {
	    my $application = $self->{APPLICATION};
	    $title = $application->title;
	} else {
	    $title = "unknown title";
	}
	$head .= "<HTML>
<HEAD><TITLE>$title</TITLE><h1>$title</h1>
</HEAD>
<BODY>
";
    }
    $content = $head . $content;

    if ($content !~ /<\/html>/i) {
	$foot = "
</BODY></HTML>
";
	$content = $content . $foot;
    }
#    print STDERR "clean_content:\n",$content;

    return $content;

}


1;
