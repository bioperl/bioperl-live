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
    my $searchio = new Bio::SearchIO(-format => 'blast', -file => 'file.bls' );
    while( my $sbjct = $searchio->next_subject ) {
	# process the Bio::Search::SubjectI object
	while( my $hsp = $sbjct->next_hsp ) { 
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

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SearchIO;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Root::IO;
use Bio::SearchIO::EventGeneratorI;
use Bio::SearchIO::SearchResultEventBuilder;

@ISA = qw(Bio::Root::Root Bio::Root::IO Bio::SearchIO::EventGeneratorI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SearchIO();
 Function: Builds a new Bio::SearchIO object 
 Returns : Bio::SearchIO initialized with the correct format
 Args    : Args    : -file => $filename
           -format => format
           -fh => filehandle to attach to

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
	    $class->_guess_format( $param{'-file'} || $ARGV[0] ) ||
		'blast';
	$format = "\L$format";	# normalize capitalization to lower case

	# normalize capitalization
	return undef unless( &_load_format_module($format) );
	return "Bio::SearchIO::$format"->new(@args);
    }
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
    $self->SUPER::_initialize(@args);
    
    # initialize the IO part
    $self->_initialize_io(@args);
    $self->attach_EventHandler(new Bio::SearchIO::SearchResultEventBuilder());
    $self->{'_reporttype'} = '';
    
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
   return 'blast'   if /\.(bls|blast)$/i;
   return 'fasta' if /\.(fas|fasta)$/i;
   return 'blastxml'     if /\.xml$/i;
}


sub DESTROY {
    my $self = shift;

    $self->close();
}

1;
