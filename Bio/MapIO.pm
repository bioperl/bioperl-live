#
# BioPerl module for Bio::MapIO
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::MapIO - A Map Factory object

=head1 SYNOPSIS

    use Bio::MapIO;
    my $mapio = Bio::MapIO->new(-format => "mapmaker",
			       -file   => "mapfile.map");

    while( my $map = $mapio->next_map ) { 
	# get each map
	foreach my $marker ( $map->each_element ) {
	    # loop through the markers associated with the map
	}
    }

=head1 DESCRIPTION

This is the Factory object for reading Maps from a data stream or file.

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
of the bugs and their resolution. Bug reports can be submitted the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::MapIO;
use strict;


use base qw(Bio::Root::Root Bio::Root::IO Bio::Factory::MapFactoryI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::MapIO->new();
 Function: Builds a new Bio::MapIO object 
 Returns : Bio::MapIO
 Args    :


=cut

sub new {
  my($caller,@args) = @_;

  my $class = ref($caller) || $caller;
  
  # or do we want to call SUPER on an object if $caller is an
  # object?
  if( $class =~ /Bio::MapIO::(\S+)/ ) {
	my ($self) = $class->SUPER::new(@args);	
	$self->_initialize(@args);
	return $self;
    } else { 
	
	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	my $format = $param{'-format'} || 
	    $class->_guess_format( $param{'-file'} || $ARGV[0] ) ||
		'mapmaker';
	$format = "\L$format";	# normalize capitalization to lower case

	# normalize capitalization
	return unless( $class->_load_format_module($format) );
	return "Bio::MapIO::$format"->new(@args);
    }

}


=head2 format

 Title   : format
 Usage   : $format = $stream->format()
 Function: Get the map format
 Returns : map format
 Args    : none

=cut

# format() method inherited from Bio::Root::IO


=head2 Bio::Factory::MapFactoryI methods

=cut

=head2 next_map

 Title   : next_tree
 Usage   : my $map = $factory->next_map;
 Function: Get a map from the factory
 Returns : L<Bio::Map::MapI>
 Args    : none


=head2 write_map

 Title   : write_tree
 Usage   : $factory->write_map($map);
 Function: Write a map out through the factory
 Returns : none
 Args    : L<Bio::Map::MapI>

=cut


=head2 attach_EventHandler

 Title   : attach_EventHandler
 Usage   : $parser->attatch_EventHandler($handler)
 Function: Adds an event handler to listen for events
 Returns : none
 Args    : L<Bio::Event::EventHandlerI>

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
 Returns : L<Bio::Event::EventHandlerI>
 Args    : none


=cut

sub _eventHandler{
   my ($self) = @_;
   return $self->{'_handler'};
}

sub _initialize {
    my($self, @args) = @_;
    $self->{'_handler'} = undef;
    
    # initialize the IO part
    $self->_initialize_io(@args);
#    $self->attach_EventHandler(Bio::MapIO::MapEventBuilder->new());
}

=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL MapIO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_format_module {
  my ($self,$format) = @_;
  my $module = "Bio::MapIO::" . $format;
  my $ok;  
  eval {
      $ok = $self->_load_module($module);
  };
  if ( $@ ) {
    print STDERR <<END;
$self: $format cannot be found
Exception $@
For more information about the MapIO system please see the MapIO docs.
This includes ways of checking for formats at compile time, not run time
END
  ;
  }
  return $ok;
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
   return 'mapmaker'   if /\.(map)$/i;
   return 'mapxml' if /\.(xml)$/i;
}

sub DESTROY {
    my $self = shift;

    $self->close();
}

1;
