#
# BioPerl module for wrapping runtime parameters
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chad Matsalla (bioinformatics1 at dieselwurks dot com)
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::ParametersI - A Base object for the parameters used to run programs

=head1 SYNOPSIS

  # do not use this object directly, it provides the following methods
  # for its subclasses

  my $void   = $obj->set_parameter("parameter_name","parameter_value"); 
  my $value  = $obj->get_parameter("parameter_name");

=head1 DESCRIPTION

This is a basic container to hold the parameters used to run a program.

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

=head1 AUTHOR - Chad Matsalla

Email bioinformatics1 at dieselwurks dot com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Run::ParametersI;
use strict;

# Object preamble - inherits from Bio::Root::Root


use base qw(Bio::Root::RootI);

=head2 get_parameter

 Title   : get_parameter
 Usage   : $parameter_object->get_parameter($param_name);
 Function: Get the value of a parameter named $param_name
 Returns : A scalar that should be a string
 Args    : A scalar that should be a string

=cut

sub get_parameter {
   my ($self,$arg) = @_;
     $self->throw_not_implemented;
}


=head2 set_parameter

 Title   : set_parameter
 Usage   : $parameter_object->set_parameter($param_name => $param_value);
 Function: Set the value of a parameter named $param_name to $param_value
 Returns : Void
 Args    : A hash containing name=>value pairs

=cut

sub set_parameter {
   my ($self,$name,$value) = @_;
     $self->throw_not_implemented;
}



1;
