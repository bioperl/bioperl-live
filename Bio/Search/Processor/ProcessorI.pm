
#
# BioPerl module for Bio::Search::Processor::ProcessorI
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Processor::ProcessorI - Abstract Interface for Processor Objects

=head1 SYNOPSIS

use Bio::Search::Processor

my $processor = new Bio::Search::Processor -file      => 'mysearchrun',
                                           -algorithm => 'Blast';

while ($result = $processor->next_result()) {

    $id = $result->get_query_id();
    $lib = $result->get_library_name();
    $size = $result->get_library_size();

    foreach $hit ( $result->get_hits() ) {
        $matchid = $hit->get_id();
        $matchdesc = $hit->get_desc();
        # etc, etc, do stuff.
    }
}

=head1 DESCRIPTION

A Processor object is used to generate Bio::Search::Result::* objects, given
a source of Search data (a file or filehandle).  The Processor object works
very much like the SeqIO system: once initialized with the new() method, the
Processor object will continue to yield as many Result objects as are available
from the data source (for single "runs" this is often only one Result object).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::Processor::ProcessorI;

use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Exporter);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
    my($self,@args) = @_;

    my $make = $self->SUPER::_initialize;

# set stuff in self from @args

    my($file, $fh, $algorithm) =
	$self->_rearrange( [ qw(
                                FILE
                                FH
                                ALGORITHM
                               )
                           ], @args);

    defined $file || defined $fh or $self->throw("Must supply either a file or fh!");
    defined $algorithm or $self->throw("Must specify algorithm!");

    $self->_algorithm($algorithm);

    $fh = new IO::File $file if defined $file;
    $self->_filehandle($fh);

    return $make; # success - we hope!
}

=head2 next_result

 Title   : next_result
 Usage   : $result = $processor->next_result()
 Function: Returns the next Bio::Search::Result::* object available from
           the provided data stream
 Returns : Bio::Search::Result object
 Args    : <none>


=cut

sub next_result{
   my ($self,@args) = @_;

   $self->throw("Abstract processor call of next_result() - your processor type has not implemented this method!");

}

=head2 _filehandle

 Title   : _filehandle
 Usage   : $obj->_filehandle($newval)
 Function: Internal function to handle handles
 Returns : value of _filehandle
 Args    : newvalue (optional)


=cut

sub _filehandle{
    my ($obj,$value) = @_;

    if( defined $value) {
	$obj->{'_filehandle'} = $value;
    }
    return $obj->{'_filehandle'};
}

=head2 _algorithm

 Title   : _algorithm
 Usage   : $obj->_algorithm($newval)
 Function: Internal function to specify algorithm
 Returns : value of _algorithm
 Args    : newvalue (optional)


=cut

sub _algorithm{
    my ($obj,$value) = @_;

    if( defined $value) {
	$obj->{'_algorithm'} = $value;
    }
    return $obj->{'_algorithm'};
}


1;

__END__


