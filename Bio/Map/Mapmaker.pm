# BioPerl module for Bio::Mapmaker
#
# Cared for by Jason Stajich <bioinformatics1@dieselwurks.com>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Mapmaker - Parser for Mapmaker Linkage Map output files

=head1 SYNOPSIS

=head1 DESCRIPTION

Parse Mapmaker files.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Chad Matsalla

Email bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::Mapmaker;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Root::IO;
use Bio::Map::Marker;
use Bio::Map::OrderedPositionWithDistance;

@ISA = qw(Bio::Root::Root Bio::Root::IO);

=head2 new(-file => $filename)

 Title   : new(-file => $filename)
 Usage   : my $obj = new Bio::Map::Mapmaker(-file => $filename);
 Function: Parses a Mapmaker object
 Returns : A Bio::Map::Mapmaker object
 Args    : -file => the name of a file

=cut

sub new {
	my($class,@args) = @_;
	my $self = {};
	my %param = @args;
	my $filename = $param{'-file'};
	my $type = $param{'-type'};
	$self->{'filename'} = $filename;
	$self->{'type'} = $type;
	$self->{fh} = Bio::Root::IO->new(-file=>"$filename");
	return bless $self;
}


=head2 next_marker()

 Title   : next_marker()
 Usage   : 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub next_marker {
	my $self = shift;
	if ($self->{done_markers}) { return; }
	my $line;
	if (!$self->{processed_top}) {
		until (($line = $self->{fh}->_readline()) =~ /^  Markers          Distance/) { }
		$self->{processed_top} = 1;
	}
	$line = $self->{fh}->_readline();
	chomp $line;
	if ($line =~ /----------/) {
		$self->{done_markers} = 1;
	}
	$line =~ /\s+(\S+)\s+(\S+)\s+(\S+)/;
	my ($marker_number,$marker_name,$marker_distance);
	$marker_number = $1;
        $marker_name = $2;
        $marker_distance = $3;
	my $o_marker = new Bio::Map::Marker(-name=> $marker_name,
		-position => new Bio::Map::OrderedPositionWithDistance(
			-positions => $marker_number,
                	-distance => $marker_distance
                	)
		);
}

=head2 summary_info()

 Title   : summary_info()
 Usage   : $o_mapmaker->summary_info();
 Function:
 Example :
 Returns :
 Args    :


=cut

sub summary_info {
	my $self = shift;
	if (!$self->{done_markers}) {
		$self->warn("You can't get the summary info until all markers are read from the mapmaker file. Sorry.\n");
		return;
	}
	my $line = $self->{fh}->_readline();
	chomp $line;
	return $line;
}






=head2 attach_EventHandler

 Title   : attach_EventHandler
 Usage   : $parser->attatch_EventHandler($handler)
 Function: Adds an event handler to listen for events
 Returns : none
 Args    : Bio::Event::EventHandlerI

=cut

sub attach_EventHandler{
    my ($self,$handler) = @_;
    return if( ! $handler );
    if( ! $handler->isa('Bio::Event::EventHandlerI') ) {
	$self->warn("Ignoring request to attatch handler ".ref($handler). ' because it is not a Bio::Event::EventHandlerI');
    }
    $self->{'_handler'} = $handler;
    return;
}

=head2 _eventHandler

 Title   : _eventHandler
 Usage   : private
 Function: Get the EventHandler
 Returns : Bio::Event::EventHandlerI
 Args    : none


=cut

sub _eventHandler{
   my ($self) = @_;
   return $self->{'_handler'};
}

sub _initialize {
    my($self, @args) = @_;
    $self->{'_handler'} = undef;
    # not really necessary unless we put more in RootI
    $self->SUPER::_initialize(@args);
    
    # initialize the IO part
    $self->_initialize_io(@args);
    $self->attach_EventHandler(new Bio::TreeIO::TreeEventBuilder());
}

=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL TreeIO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_format_module {
  my ($format) = @_;
  my ($module, $load, $m);

  $module = "_<Bio/TreeIO/$format.pm";
  $load = "Bio/TreeIO/$format.pm";

  return 1 if $main::{$module};
  eval {
    require $load;
  };
  if ( $@ ) {
    print STDERR <<END;
$load: $format cannot be found
Exception $@
For more information about the TreeIO system please see the TreeIO docs.
This includes ways of checking for formats at compile time, not run time
END
  ;
    return;
  }
  return 1;
}


=head2 _guess_format

 Title   : _guess_format
 Usage   : $obj->_guess_format($filename)
 Function:
 Example :
 Returns : guessed format of filename (lower case)
 Args    :

=cut

sub _guess_format {
   my $class = shift;
   return unless $_ = shift;
   return 'newick'   if /\.(dnd|newick|nh)$/i;
   return 'phyloxml' if /\.(xml)$/i;
}

sub DESTROY {
    my $self = shift;

    $self->close();
}

1;
