# $Id$
#
# BioPerl module for Bio::SearchIO
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO - Driver for parsing Sequence Database Searches (Blast,FASTA,...)

=head1 SYNOPSIS

    use Bio::SearchIO;
    my $searchio = new Bio::SearchIO( -format => 'blastxml',
                                      -file   => 'blastout.xml' );
    while ( my $result = $in->next_result() ) {
       while( my $hit = $result->next_hit ) {
	# process the Bio::Search::HitI object
           while( my $align = $sbjct->next_align ) { 
	    # process the Bio::Search::HSPI object
	}
    }
=head1 DESCRIPTION

This is a driver for instantiating a parser for report files from
sequence database searches.

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

=head1 AUTHOR - Jason Stajich & Steve Chervitz

Email jason@bioperl.org
Email sac@bioperl.org

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SearchIO;
use strict;
use vars qw(@ISA);

# Object preamble - inherits from Bio::Root::IO

use Bio::Root::IO;
use Bio::Event::EventGeneratorI;
use Bio::SearchIO::SearchResultEventBuilder;
use Bio::AnalysisParserI;
use Symbol();

@ISA = qw( Bio::Root::IO Bio::Event::EventGeneratorI Bio::AnalysisParserI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO();
 Function: Builds a new Bio::SearchIO object 
 Returns : Bio::SearchIO initialized with the correct format
 Args    : Args    : -file => $filename
           -format => format
           -fh => filehandle to attach to
           -result_factory => Object implementing Bio::Factory::ResultFactoryI
           -hit_factory    => Object implementing Bio::Factory::HitFactoryI
           -writer         => Object implementing Bio::SearchIO::SearchWriterI

=cut

sub new {
  my($caller,@args) = @_;
  my $class = ref($caller) || $caller;
    
  # or do we want to call SUPER on an object if $caller is an
  # object?
  if( $class =~ /Bio::SearchIO::(\S+)/ ) {
    my ($self) = $class->SUPER::new(@args);	
    $self->_initialize(@args);
    return $self;
  } else { 
    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
    my $format = $param{'-format'} ||
      $class->_guess_format( $param{'-file'} || $ARGV[0] ) || 'blast';
    
    # normalize capitalization to lower case
    $format = "\L$format";
    
    return undef unless( &_load_format_module($format) );
    return "Bio::SearchIO::${format}"->new(@args);
  }
}

=head2 newFh

 Title   : newFh
 Usage   : $fh = Bio::SearchIO->newFh(-file=>$filename,
                                      -format=>'Format')
 Function: does a new() followed by an fh()
 Example : $fh = Bio::SearchIO->newFh(-file=>$filename,
                                      -format=>'Format')
           $result = <$fh>;   # read a ResultI object
           print $fh $result; # write a ResultI object
 Returns : filehandle tied to the Bio::SearchIO::Fh class
 Args    :

=cut

sub newFh {
  my $class = shift;
  return unless my $self = $class->new(@_);
  return $self->fh;
}

=head2 fh

 Title   : fh
 Usage   : $obj->fh
 Function:
 Example : $fh = $obj->fh;      # make a tied filehandle
           $result = <$fh>;     # read a ResultI object
           print $fh $result;   # write a ResultI object
 Returns : filehandle tied to the Bio::SearchIO::Fh class
 Args    :

=cut


sub fh {
  my $self = shift;
  my $class = ref($self) || $self;
  my $s = Symbol::gensym;
  tie $$s,$class,$self;
  return $s;
}

=head2 attach_EventHandler

 Title   : attach_EventHandler
 Usage   : $parser->attatch_EventHandler($handler)
 Function: Adds an event handler to listen for events
 Returns : none
 Args    : Bio::SearchIO::EventHandlerI

=cut

sub attach_EventHandler{
    my ($self,$handler) = @_;
    return if( ! $handler );
    if( ! $handler->isa('Bio::SearchIO::EventHandlerI') ) {
	$self->warn("Ignoring request to attatch handler ".ref($handler). ' because it is not a Bio::SearchIO::EventHandlerI');
    }
    $self->{'_handler'} = $handler;
    return;
}

=head2 _eventHandler

 Title   : _eventHandler
 Usage   : private
 Function: Get the EventHandler
 Returns : Bio::SearchIO::EventHandlerI
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
    #$self->SUPER::_initialize(@args);
    
    # initialize the IO part
    $self->_initialize_io(@args);
    $self->attach_EventHandler(new Bio::SearchIO::SearchResultEventBuilder());
    $self->{'_reporttype'} = '';

    my ( $writer, $rfactory, $hfactory, $use_factories ) =
      $self->_rearrange([qw(WRITER 
			    RESULT_FACTORY 
			    HIT_FACTORY
			    USE_FACTORIES)], @args);

    $self->writer( $writer ) if $writer;

    # TODO: Resolve this issue:
    # The $use_factories flag is a temporary hack to allow factory-based and 
    # non-factory based SearchIO objects to co-exist. 
    # steve --- Sat Dec 22 04:41:20 2001
    if( $use_factories) {
      if( not defined $self->{'_result_factory'}) {
	$self->result_factory( $rfactory || $self->default_result_factory_class->new );
      }
      if( not defined $self->{'_hit_factory'}) {
	$self->hit_factory( $hfactory || $self->default_hit_factory_class->new );
      }
    }
}

=head2 next_result

 Title   : next_result
 Usage   : $result = stream->next_result
 Function: Reads the next ResultI object from the stream and returns it.

           Certain driver modules may encounter entries in the stream that
           are either misformatted or that use syntax not yet understood
           by the driver. If such an incident is recoverable, e.g., by
           dismissing a feature of a feature table or some other non-mandatory
           part of an entry, the driver will issue a warning. In the case
           of a non-recoverable situation an exception will be thrown.
           Do not assume that you can resume parsing the same stream after
           catching the exception. Note that you can always turn recoverable
           errors into exceptions by calling $stream->verbose(2) (see
           Bio::Root::RootI POD page).
 Returns : A Bio::Search::Result::ResultI object
 Args    : n/a

=cut

sub next_result {
   my ($self) = @_;
   $self->throw_not_implemented;
}

=head2 write_result

 Title   : write_result
 Usage   : $stream->write_result($result_result, @other_args)
 Function: Writes data from the $result_result object into the stream.
         : Delegates to the to_string() method of the associated 
         : WriterI object.
 Returns : 1 for success and 0 for error
 Args    : Bio::Search:Result::ResultI object,
         : plus any other arguments for the Writer
 Throws  : Bio::Root::Exception if a Writer has not been set.

=cut

sub write_result {
   my ($self, $result, @args) = @_;

   if( not ref($self->{'_result_writer'}) ) {
       $self->throw("ResultWriter not defined.");
   }
   my $str = $self->writer->to_string( $result, @args );
   #print "Got string: \n$str\n";
   $self->_print( "$str" );
   return 1;
}


=head2 writer

 Title   : writer
 Usage   : $writer = $stream->writer;
 Function: Sets/Gets a SearchWriterI object to be used for this searchIO.
 Returns : 1 for success and 0 for error
 Args    : Bio::SearchIO::SearchWriterI object (when setting)
 Throws  : Bio::Root::Exception if a non-Bio::SearchIO::SearchWriterI object
           is passed in.

=cut

sub writer {
    my ($self, $writer) = @_;
    if( ref($writer) and $writer->isa( 'Bio::SearchIO::SearchWriterI' )) {
        $self->{'_result_writer'} = $writer;
    }
    elsif( defined $writer ) {
        $self->throw("Can't set ResultWriter. Not a Bio::SearchIO::SearchWriterI: $writer");
    }
    return $self->{'_result_writer'};
}


=head2 hit_factory

 Title   : hit_factory
 Usage   : $hit_factory = $stream->hit_factory;  (get)
         : $stream->hit_factory( $factory );     (set)
 Function: Sets/Gets a factory object to create hit objects for this SearchIO
 Returns : Bio::Factory::HitFactoryI object 
 Args    : Bio::Factory::HitFactoryI object (when setting)
 Throws  : Bio::Root::Exception if a non-Bio::Factory::HitFactoryI object  
           is passed in.
 Comments: A SearchIO implementation should provide a default hit factory.

=cut

sub hit_factory {
    my ($self, $fact) = @_;
    if( ref $fact and $fact->isa( 'Bio::Factory::HitFactoryI' )) {
    	   $self->{'_hit_factory'} = $fact;
    }
    elsif( defined $fact ) {
        $self->throw("Can't set HitFactory. Not a Bio::Factory::HitFactoryI: $fact");
    }
    return $self->{'_hit_factory'};
}

=head2 result_factory

 Title   : result_factory
 Usage   : $result_factory = $stream->result_factory;  (get)
         : $stream->result_factory( $factory );        (set)
 Function: Sets/Gets a factory object to create result objects for this SearchIO.
 Returns : Bio::Factory::ResultFactoryI object 
 Args    : Bio::Factory::ResultFactoryI object (when setting)
 Throws  : Bio::Root::Exception if a non-Bio::Factory::ResultFactoryI object
           is passed in.
 Comments: A SearchIO implementation should provide a default result factory.

=cut

sub result_factory {
    my ($self, $fact) = @_;
    if( ref $fact and $fact->isa( 'Bio::Factory::ResultFactoryI' )) {
    	   $self->{'_result_factory'} = $fact;
    }
    elsif( defined $fact ) {
        $self->throw("Can't set ResultFactory. Not a Bio::Factory::ResultFactoryI: $fact");
    }
    return $self->{'_result_factory'};
}

=head2 default_hit_factory_class

 Title   : default_hit_factory_class
 Usage   : $res_factory = $obj->default_hit_factory_class()->new( @args )
 Function: Returns the name of the default class to be used for creating
           Bio::Search::Hit::HitI objects.
 Example :
 Returns : A string containing the name of a class that implements 
           the Bio::Search::Hit::HitI interface.
 Args    : none
 Comments: Bio::SearchIO does not implement this method. It throws a NotImplemented
           exception

=cut

sub default_hit_factory_class {
    my $self = shift;
# TODO: Uncomment this when Jason's SearchIO code conforms
#    $self->throw_not_implemented;
}


=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL SearchIO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_format_module {
  my ($format) = @_;
  my ($module, $load, $m);

  $module = "_<Bio/SearchIO/$format.pm";
  $load = "Bio/SearchIO/$format.pm";

  return 1 if $main::{$module};
  eval {
    require $load;
  };
  if ( $@ ) {
    print STDERR <<END;
$load: $format cannot be found
Exception $@
For more information about the SearchIO system please see the SearchIO docs.
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
   return 'blast'   if (/blast/i or /\.bl\w$/i);
   return 'fasta' if (/fasta/i or /\.fas$/i);
   return 'blastxml' if (/blast/i and /\.xml$/i);
}


sub DESTROY {
    my $self = shift;

    $self->close();
}

sub TIEHANDLE {
  my $class = shift;
  return bless {processor => shift}, $class;
}

sub READLINE {
  my $self = shift;
  return $self->{'processor'}->next_result() unless wantarray;
  my (@list, $obj);
  push @list, $obj while $obj = $self->{'processor'}->next_result();
  return @list;
}

sub PRINT {
  my $self = shift;
  $self->{'processor'}->write_result(@_);
}

1;
__END__
