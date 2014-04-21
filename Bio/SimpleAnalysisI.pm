#
# BioPerl module for Bio::SimpleAnalysisI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Martin Senger <martin.senger@gmail.com>
# For copyright and disclaimer see below.
#

# POD documentation - main docs before the code

=head1 NAME

Bio::SimpleAnalysisI - A simple interface to any (local or remote) analysis tool

=head1 SYNOPSIS

This is an interface module - you do not instantiate it.
Use other modules instead (those that implement this interface).

=head1 DESCRIPTION

This interface contains public methods for accessing and controlling
local and remote analysis tools. It is meant to be used on the client
side. The interface consists only of a necessary set of methods for
synchronous invocation of analysis tools. For more complex set,
including an asynchronous access, see interface C<Bio::AnalysisI>
(which inherits from this one, by the way).

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

  https://github.com/bioperl/bioperl-live/issues

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

=over

=item *

http://www.ebi.ac.uk/Tools/webservices/soaplab/guide

=back

=head1 APPENDIX

This is actually the main documentation...

If you try to call any of these methods directly on this
C<Bio::SimpleAnalysisI> object you will get a I<not implemented> error
message.

=cut


# Let the code begin...

package Bio::SimpleAnalysisI;
use strict;

use base qw(Bio::Root::RootI);

# -----------------------------------------------------------------------------

=head2 analysis_name

 Usage   : $tool->analysis_name;
 Returns : a name of this analysis
 Args    : none

=cut

sub analysis_name { shift->throw_not_implemented(); }

# -----------------------------------------------------------------------------

=head2 analysis_spec

 Usage   : $tool->analysis_spec;
 Returns : a hash reference describing this analysis
 Args    : none

The returned hash reference uses the following keys (not all of them always
present, perhaps others present as well): C<name>, C<type>, C<version>,
C<supplier>, C<installation>, C<description>.

=cut

sub analysis_spec { shift->throw_not_implemented(); }

# -----------------------------------------------------------------------------

=head2 input_spec

 Usage   : $tool->input_spec;
 Returns : an array reference with hashes as elements
 Args    : none

The analysis input data are named, and can be also associated with a
default value, with allowed values and with few other attributes. The
names are important for feeding the analysis with the input data (the
inputs are given to methods C<run> and C<wait_for> as name/value
pairs).

=cut

sub input_spec { shift->throw_not_implemented(); }

# -----------------------------------------------------------------------------

=head2 result_spec

 Usage   : $tool->result_spec;
 Returns : a hash reference with result names as keys
           and result types as values
 Args    : none

An analysis can produce several results, or the same result in several
different formats. All such results are named and can be retrieved
using their names by metod C<result>.

Here is an example of the result specification:

  $result_spec = {
          'outseq' => 'String',
          'report' => 'String',
          'detailed_status' => 'String'
        };

=cut

sub result_spec { shift->throw_not_implemented(); }

# -----------------------------------------------------------------------------

=head2 run

 Usage   : $tool->run ( ['sequence=@my.seq', 'osformat=embl'] )
 Returns : $self
 Args    : data and parameters for this execution
           (in various formats)

Create a job, start it, and wait for its completion. The method is
identical to the method C<wait_for>. Why there are two methods doing
the same? Because it is expected that the sub-classes may implement
them differently (an example is an interface C<Bio::AnalysisI> which
uses method C<run> for an asynchronous execution and method
C<wait_for> for a synchronous one.

Usually, after this call, you ask for results of the finished job:

    $analysis->run (...)->result;

The input data and prameters for this execution can be specified in
various ways:

=over

=item array reference

The array has scalar elements of the form

   name = [[@]value]

where C<name> is the name of an input data or input parameter (see
method C<input_spec> for finding what names are recognized by this
analysis) and C<value> is a value for this data/parameter. If C<value>
is missing a 1 is assumed (which is convenient for the boolean
options). If C<value> starts with C<@> it is treated as a local
filename, and its contents is used as the data/parameter value.

=item hash reference

The same as with the array reference but now there is no need to use
an equal sign. The hash keys are input names and hash values their
data. The values can again start with a C<@> sign indicating a local
filename.

=back

=cut

sub run { shift->throw_not_implemented(); }

# -----------------------------------------------------------------------------

=head2 wait_for

 Usage   : $tool->wait_for ( { 'sequence' => '@my,file' } )
 Returns : $self
 Args    : the same as for method 'run'

Create a job, start it and wait for its completion. The method is
identical to the method C<run>. See details in the C<run> method.

=cut

sub wait_for { shift->throw_not_implemented(); }

# -----------------------------------------------------------------------------

=head2 status

 Usage   : $tool->status
 Returns : string describing a status of the execution
 Args    : none

It returns one of the following strings (and perhaps more if a server
implementation extended possible job states):

   CREATED              (not run yet)
   COMPLETED            (run and finished normally)
   TERMINATED_BY_ERROR  (run and finished with an error or a signal)

=cut

sub status { shift->throw_not_implemented(); }

# -----------------------------------------------------------------------------

=head2 result

 Usage   : $job->result (...)
 Returns : a result created by running an analysis
 Args    : none (but an implementation may choose
           to add arguments for instructions how to process
           the raw result)

The method returns a scalar representing a result of an executed
job. If the job was terminated by an error the result may contain an
error message instead of the real data (or both, depending on the
implementation).

=cut

sub result { shift->throw_not_implemented(); }

# -----------------------------------------------------------------------------

1;
__END__

