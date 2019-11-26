#
# BioPerl module for Bio::Tools::Run::AnalysisFactory
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Martin Senger <martin.senger@gmail.com>
# For copyright and disclaimer see below.
#

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::AnalysisFactory - A directory of analysis tools

=head1 SYNOPSIS

  # list all available analyses from the default location,
  # using a default (SOAP) access method
  use Bio::Tools::Run::AnalysisFactory;
  my $list = Bio::Tools::Run::AnalysisFactory->new();
                ->available_analyses;
  use Data::Dumper; print Dumper ($list);

  # ditto, but from a different location
  use Bio::Tools::Run::AnalysisFactory;
  my $list =
     Bio::Tools::Run::AnalysisFactory->new(-location => 'http://somewhere/something')
                ->available_analyses;

  # ...and using a different access method
  # (this example is not yet impelmented)
  use Bio::Tools::Run::AnalysisFactory;
  my $list =
     Bio::Tools::Run::AnalysisFactory->new(-location => 'http://somewhere/something',
                                           -access => 'novella')
                ->available_analyses;

  # list available categories of analyses
  use Bio::Tools::Run::AnalysisFactory;
  my $categories =
     Bio::Tools::Run::AnalysisFactory->new();
                ->available_categories;
  use Data::Dumper; print Dumper ($categories);

  # show all analyses group by categories
  use Bio::Tools::Run::AnalysisFactory;
  my $factory = Bio::Tools::Run::AnalysisFactory->new();
  foreach $cat ( @{ $factory->available_categories } ) {
    my @sublist = @{ $factory->available_analyses ($cat) };
    print "$cat:\n\t",
          join ("\n\t", @{ $factory->available_analyses ($cat) }),
          "\n";
  }

  # create an analysis object
  use Bio::Tools::Run::AnalysisFactory;
  $service = Bio::Tools::Run::AnalysisFactory->new();
                 ->create_analysis ('edit.seqret');
  $service->run (
                #...
                )->results;

=head1 DESCRIPTION

The module represents a list of available analysis tools from a given
location using a given access method. Additionally, for any of the
available analyses, it can create an object of type C<Bio::Tools::Run::Analysis>.

The module is a higher-level abstraction whose main job is to load a
'real-work-doing' implementation. Which one is used, it depends on the
C<-access> parameter. The same design is used here as for
C<Bio::Tools::Run::Analysis> module.

There is available a I<SOAP> access to almost all EMBOSS applications,
running at European Bioinformatics Institute.

The documentation of all C<public> methods are to be found
in C<Bio::Factory::AnalysisI>.

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

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR

Martin Senger (martin.senger@gmail.com)

=head1 COPYRIGHT

Copyright (c) 2003, Martin Senger and EMBL-EBI.
All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 SEE ALSO

=over 4

=item *

http://www.ebi.ac.uk/soaplab/Perl_Client.html

=back

=head1 APPENDIX

Here is the rest of the object methods.  Internal methods are preceded
with an underscore _.

=cut


# Let the code begin...

package Bio::Tools::Run::AnalysisFactory;

use vars qw(@ISA $Revision);
use strict;

use Bio::Root::Root;
use Bio::Factory::AnalysisI;
@ISA = qw(Bio::Root::Root Bio::Factory::AnalysisI);


BEGIN {
    $Revision = q$Id$;
}

# -----------------------------------------------------------------------------

# Available (understood) parameters:
# -access
#  (+ parameters used in guessing an access)

# -----------------------------------------------------------------------------

=head2 new

 Usage   : my $factory =
             Bio::Tools::Run::AnalysisFactory->new(-access => 'soap',
                                                   -location => 'http://...');
 Returns : a new Bio::Tools::Run::AnalysisFactory object representing a list
           of available analyses
 Args    : There may be additional arguments which are specific
           to the access method (see methods 'new' or '_initialize'
           of the access-specific implementations (such as module
	   Bio::Tools::Run::AnalysisFactory::soap for a SOAP-based access).

           The recognised and used arguments are:
             -access
             -location
             -httpproxy
             -timeout

It builds, populates and returns a new C<Bio::Tools::Run::AnalysisFactory> object. This
is how it is seen from the outside. But in fact, it builds, populates
and returns a more specific lower-level object, for example
C<Bio::Tools::Run::AnalysisFactory::soap> object - which one it is it depends on the C<-access>
parameter.

=over 4

=item -access

It indicates what lower-level module to load.  Default is 'soap'.
Other (but future) possibilities are:

   -access => 'novella'
   -access => 'local'

=item -location

A location of the service. The contents is access-specific (see
details in the lower-level implementation modules).

Default is C<http://www.ebi.ac.uk/soaplab/services> (there are
services running at European Bioinformatics Institute on top of most
of EMBOSS analyses, and on some others).

=item -httpproxy

In addition to the I<location> parameter, you may need to specify also
a location/URL of an HTTP proxy server (if your site requires
one). The expected format is C<http://server:port>.  There is no
default value. It is also an access-specific parameter which may not
be used by all access methods.

=item -timeout

For long(er) running jobs the HTTP connection may be time-outed. In
order to avoid it (or, vice-versa, to call timeout sooner) you may
specify C<timeout> with the number of seconds the connection will be
kept alive. Zero means to keep it alive forever. The default value is
two minutes.

=back

=cut

sub new {
    my ($caller,@args) = @_;
    my $class = ref($caller) || $caller;
  
    if ($class eq 'Bio::Tools::Run::AnalysisFactory') {

	# this is called only the first time when somebody calls: 'new
	# Bio::Tools::Run::AnalysisFactory (...)', and it actually loads a
	# 'real-work-doing' module and call this new() method again
	# (unless the loaded module has its own new() method)

	my %param = @args;
	@param { map { lc $_ } keys %param } = values %param; # lowercase keys
	my $access =
	    $param {'-access'} ||                  # use -access parameter
	    $class->_guess_access ( \%param ) ||   # or guess from other parameters
	    'soap';                                # or use a default access method
	$access = "\L$access";	# normalize capitalization to lower case

	# remember the access method (putting it into @args means that the
	# object - when created - will remember it)
	push (@args, (-access => $access)) unless $param {'-access'};

	# load module with the real implementation - as defined in $access
	return undef unless (&_load_access_module ($access));

	# this calls this same method new() - but now its object part
	# (see the upper branche above) is called
	return "Bio::Tools::Run::AnalysisFactory::$access"->new (@args);

    } else {

	# if $caller is an object, or if it is an underlying
	# 'real-work-doing' class (e.g. Bio::Tools::Run::AnalysisFactory::soap)
	# then we want to call SUPER to create and bless a new object

	my ($self) = $class->SUPER::new (@args);

	# now the $self is an empty object - we will populate it from
	# the $caller - if $caller is an object (so we do cloning here)

	if (ref ($caller)) {
	    %{ $self } = %{ $caller };
	}

	# and finally add values from '@args' into the newly created
	# object (the values will overwrite the values copied above);
	# this is done by calling '_initialize' of the 'real-work-doing'
	# class (if there is no one there, there is always an empty one
	# in Bio::Root::Root)

	$self->_initialize (@args);
	return $self;
    }

}

# -----------------------------------------------------------------------------

=head2 _load_access_module

 Usage   : $class->_load_access_module ($access)
 Returns : 1 on success, undef on failure
 Args    : 'access' should contain the last part of the
           name of a module who does the real implementation

It does (in the run-time) a similar thing as

   require Bio::Tools::Run::AnalysisFactory::$access

It prints an error on STDERR if it fails to find and load the module
(for example, because of the compilation errors in the module).

=cut

sub _load_access_module {
  my ($access) = @_;

  my $load = "Bio/Tools/Run/AnalysisFactory/$access.pm";
  eval {
    require $load;
  };

  if ( $@ ) {
    Bio::Root::Root->throw (<<END);
$load: $access cannot be found or loaded
Exception $@
For more information about the Analysis system please see the Bio::Tools::Run::AnalysisFactory docs.
END
  ;
    return;
  }
  return 1;
}

# -----------------------------------------------------------------------------

=head2 _guess_access

 Usage   : $class->_guess_access ($rh_params)
 Returns : string with a guessed access protocol (e.g. 'soap'),
           or undef if the guessing failed
 Args    : 'rh_params' is a hash reference containing parameters given
           to the 'new' method.

It makes an expert guess what kind of access/transport protocol should
be used to access the underlying analysis. The guess is based on the
parameters in I<rh_params>. Remember that this method is called only
if there was no I<-access> parameter which could tell directly what
access method to use.

=cut

sub _guess_access {
   my ($class, $rh_params) = @_;
   return undef;
}

# -----------------------------------------------------------------------------

=head2 VERSION and Revision

 Usage   : print $Bio::Tools::Run::AnalysisFactory::VERSION;
           print $Bio::Tools::Run::AnalysisFactory::Revision;

=cut

1;
__END__
