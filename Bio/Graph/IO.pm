# $Id$
#
# BioPerl module for Bio::Graph::IO
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1  NAME

Bio::Graph::IO - Class for reading /writing biological graph data.

=head1  SYNOPSIS

  # This is a class for reading /writing biological data that can
  # be represented by graphs e.g., protein interaction data.

  # e.g., a graph reformatter..
  my $graph_in = Bio::Graph::IO->new(-file =>'myfile.dat',
                                     -format=>'dip' );
  my $graph = $graph_in->next_graph();
  my $graph_out = Bio::Graph::IO->new(-file =>'outfile.dat',
                                      -format=>'psixml') ;
  $graph_out->write_graph($graph);

=head1  DESCRIPTION

This class is analagous to the SeqIO and AlignIO classes. To read in a
file of a particular format, file and format are given as key/value
pairs as arguments.  The Bio::Graph::IO checks that the appropriate
module is available and loads it.

At present only the DIP tab delimited format and PSI XML format are supported

=head1 METHODS

The main methods are:

=head2  $graph = $io-E<gt>next_graph()

The next_graph method does not imply that multiple graphs are
contained in a file, more to maintain the consistency of nomenclature
with the $seqio-E<gt>next_seq() and $alnio-E<gt>next_aln() methods.

=head2  $io-E<gt>write_graph($graph) (not implemented yet).

Writes the graph data to file in requested format.


=head1  REQUIREMENTS

To read or write from XML you will need the XML::Twig module available
from CPAN.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.

Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Richard Adams

Email richard.adams@ed.ac.uk

=cut

use strict;
package Bio::Graph::IO;
use base qw(Bio::Root::IO);

=head2  new

 Name       : new
 Usage      : $io = Bio::Graph::IO->new(-file => 'myfile.dat', 
                                        -format => 'dip');
 Returns    : A Bio::Graph::IO stream initialised to the appropriate format.
 Args       : Named parameters: 
              -file      => $filename
              -format    => format
	      -threshold => a confidence score for the interaction, optional

=cut

sub new {
	my ($caller, @args) = @_;
	my $class           = ref($caller) || $caller;
	if ($class =~ /Bio::Graph::IO::(\S+)/){
		my $self = $class->SUPER::new(@args);
		$self->_initialize_io(@args);
		return $self;
	} else {
		my %param = @args;
		@param{ map { lc $_ } keys %param } = values %param;
		if (!exists($param{'-format'})) {
			Bio::Root::Root->throw("Must specify a valid format!");
		} 
		my $format = $param{'-format'};
		$format    = "\L$format";	
		return unless ($class->_load_format_module($format)); 
		return "Bio::Graph::IO::$format"->new(@args);
	}
}

=head2    next_graph

 Name       : next_graph
 Usage      : $gr = $io->next_graph().
 Returns    : A Bio::Graph::ProteinGraph object.
 Args       : None

=cut

sub next_graph {
   my ($self, $gr) = @_;
   $self->throw("Sorry, you cannot read from a generic Bio::Graph::IO object.");
}

=head2    write_graph

 Name       : write_graph
 Usage      : $gr = $io->write_graph($graph).
 Args       : A Bio::Graph object.
 Returns    : None

=cut

sub write_graph {
   my ($self, $gr) = @_;
   $self->throw("Sorry, you can't write from a generic Bio::GraphIO object.");
}


=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL BioGraphIO stuff*
 Function: Loads up (like use) a module at run time on demand
 Returns :
 Args    :

=cut

sub _load_format_module {

my ($self, $format) = @_;
    my $module = "Bio::Graph::IO::" . $format;
    my $ok;

    eval {
	$ok = $self->_load_module($module);
    };
    if ( $@ ) {
    print STDERR <<END
$self: $format cannot be found
Exception $@
For more information about the Bio::Graph::IO system please see the Bio:Graph::IO docs.
END
  ;
  }
  return $ok;

}

sub _initialize_io {

	my ($self, @args) = @_;
	$self->SUPER::_initialize_io(@args);
	my ($th) = $self->_rearrange( [qw(THRESHOLD)], @args);
	$self->{'_th'} = $th;
	return $self;

}

1;
