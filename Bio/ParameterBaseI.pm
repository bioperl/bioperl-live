#
# BioPerl module for Bio::ParameterBaseI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chris Fields <cjfields at bioperl dot org>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::ParameterBaseI - Simple interface class for any parameter-related data such
as IDs, database name, program arguments, and other odds and ends.

=head1 SYNOPSIS

  # Bio::DB::MyParams implements Bio::ParameterBaseI

  @params = (-db   => 'protein',
             -id   => \@ids,
             -retmax => 10);

  $pobj->Bio::DB::MyDBParams->new();

  # sets only parameters passed; results in a state change if any parameter
  # passed is new or differs from previously set value

  $pobj->set_params(@params);

  # reset all parameters (sets to undef); results in a state change

  $pobj->reset_params();

  # resets parameters to those in %param (sets all others to undef); resets the
  # object state to indicate change.

  $pobj->reset_params(@params);

  # direct get/set; results in a state change if any parameter passed is new or
  # differs from previously set value

  $pobj->db('nucleotide');
  @ids = $pobj->id();

  # retrieve list containing set defined parameters

  %myparams = $pobj->get_parameters();

  # checks whether the state of the object has changed (i.e. parameter has
  # changed, so on)

  if ($pobj->parameters_changed) {
     # run new search
  } else {
     # return cached search
  }

  # available parameters

  @params = $pobj->available_parameters();

  # retrieve string (URI, query, etc); calling to* methods changes object state
  # to indicate data hasn't changed (so future calls to parameters_changed()
  # will return FALSE)

  $query = $pobj->to_string(); # returns raw string
  $uri = $pobj->to_uri(); #  returns URI-based object
  $uri = $pobj->to_my_data_struct(); #  returns implemenation-specific data structure
  ...

=head1 DESCRIPTION

This is a class interface which focuses on common parameter-related tasks such
as building simple database queries, URI-related requests, program arguments,
etc.

Implementing classes use the following ways to set parameters:

1) Create a new instance of a ParameterBaseI-implementing object.

  $pobj->Bio::DB::MyParamClass->new(-db => 'local', -id => \@ids);

2) Pass the parameters as a hash or array to set_parameters(), which sets the
parameters listed in the hash but leaves all others as is.

  $pobj->set_parameters(-retmax => 100, -retstart => 20); 

3) Pass the parameters as a hash or array to reset_parameters(), which sets the
parameters listed in the hash and resets everything else.

  $pobj->reset_parameters(-term => 'pyrimidine'); # sets db and id to undef

4) Pass values using specific getter/setters.

  $pobj->id(\@ids); # sets IDs

There is no restriction on what one uses to set up individual parameter
getter/setters, though there are some other options implemented in BioPerl (for
instance, Bio::Root::RootI::_set_from_args()).

A key requirement is there be a way to detect changes in the state of the
ParameterBaseI object so that any object with a Bio::ParameterBaseI can decide
whether to submit a new request or return cached data. State changes are
revealed by the returned values of the parameters_changed() method, which is a
simple boolean set to TRUE when the object is first instantiated or parameters
have changed. 

When retrieving anything using the implementation-specific to_* methods (such as
to_query, to_string, to_uri, to_request, etc), the ParameterBaseI object state
is set to FALSE to indicate the data has been accessed and indicate reaccessing
will retrieve the same value. The observing object can then independently decide
whether to rerun the cached query or return a previously cached result. 

One can also use indiviual getter/setters to retrieve single parameter values as
well as use parameter_hash() to retrieve all of the parameters in one go as a
hash. To check which parameters are available use available_parameters().  Args
passed to 

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Email cjfields at bioperl dot org

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::ParameterBaseI;
use strict;
use warnings;

use base qw(Bio::Root::RootI);

=head2 set_parameters

 Title   : set_parameters
 Usage   : $pobj->set_parameters(%params);
 Function: sets the parameters listed in the hash or array
 Returns : None
 Args    : [optional] hash or array of parameter/values.  

=cut

sub set_parameters {
    shift->throw_not_implemented;
}

=head2 reset_parameters

 Title   : reset_parameters
 Usage   : resets values
 Function: resets parameters to either undef or value in passed hash
 Returns : none
 Args    : [optional] hash of parameter-value pairs

=cut

sub reset_parameters {
    shift->throw_not_implemented;
}

=head2 parameters_changed

 Title   : parameters_changed
 Usage   : if ($pobj->parameters_changed) {...}
 Function: Returns boolean true (1) if parameters have changed
 Returns : Boolean (0 or 1)
 Args    : [optional] Boolean

=cut

sub parameters_changed {
    shift->throw_not_implemented;
}

=head2 available_parameters

 Title   : available_parameters
 Usage   : @params = $pobj->available_parameters()
 Function: Returns a list of the available parameters
 Returns : Array of parameters
 Args    : [optional, implementation-dependent] string for returning subset of
           parameters

=cut

sub available_parameters {
    shift->throw_not_implemented;
}

=head2 get_parameters

 Title   : get_parameters
 Usage   : %params = $pobj->get_parameters;
 Function: Returns list of key-value pairs of parameter => value
 Returns : List of key-value pairs
 Args    : [optional] A string is allowed if subsets are wanted or (if a
           parameter subset is default) 'all' to return all parameters

=cut

sub get_parameters {
    shift->throw_not_implemented;
}

=head1 to* methods

All to_* methods are implementation-specific

=cut

1;

