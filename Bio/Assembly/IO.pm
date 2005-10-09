# $Id$
#
# BioPerl module for Bio::Assembly::IO
#
#   based on the Bio::SeqIO module
#       by Ewan Birney <birney@sanger.ac.uk>
#       and Lincoln Stein  <lstein@cshl.org>
#
# Copyright Robson Francisco de Souza
#
# You may distribute this module under the same terms as perl itself
#
# _history

# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::IO - Handler for Assembly::IO Formats

=head1 SYNOPSIS

    use Bio::Assembly::IO;

    $in  = Bio::Assembly::IO->new(-file=>"<inputfilename",
                                -format=>'phrap');
    $out = Bio::Assembly::IO->new(-file=>">outputfilename",
                                -format=>'phrap');

    while ( my $seq = $in->next_seq() ) {
       $out->write_seq($seq);
    }

=head1 DESCRIPTION

Bio::Assembly::IO is a handler module for formats in the Assembly::IO set
(e.g. Bio::Assembly::IO::phrap).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bio.perl.org/MailList.html     - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR

Robson Francisco de Souza

E-mail: rfsouza@citri.iq.usp.br

=head1 CONTRIBUTORS

#

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Assembly::IO;

use Bio::Root::Root;
use Bio::Root::IO;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : Bio::Assembly::IO->new(-file =>$filename,-format=>'format')
 Function: Returns a new assembly stream
 Returns : A Bio::Assembly::IO::Handler initialised 
           with the appropriate format
 Args    : -file => $filename
           -format => format

=cut

sub new {
    my ($caller,@args) = @_;
    my $class = ref($caller) || $caller;
    
    # or do we want to call SUPER on an object if $caller is an
    # object?
    if( $class =~ /Bio::Assembly::IO::(\S+)/ ) {
	my ($self) = $class->SUPER::new(@args);	
	$self->_initialize(@args);
	return $self;
    } else { 

	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys

	$class->throw("Need at least a file name to proceed!")
	    unless (defined $param{'-file'} || defined $ARGV[0]);

	my $format = $param{'-format'} || 
	    $class->_guess_format( $param{-file} || $ARGV[0] );
	$format = "\L$format";	# normalize capitalization to lower case

	# normalize capitalization
	return undef unless( $class->_load_format_module($format) );
	return "Bio::Assembly::IO::$format"->new(@args);
    }
}

# _initialize is chained for all SeqIO classes

sub _initialize {
    my($self, @args) = @_;
    # initialize the IO part
    $self->_initialize_io(@args);
}

=head2 next_assembly

 Title   : next_assembly
 Usage   : $cluster = $stream->next_assembly()
 Function: Reads the next assembly object from the stream and returns it.
 Returns : a Bio::Assembly::ScaffoldI compliant object
 Args    : none

=cut

sub next_assembly {
   my ($self, $seq) = @_;
   $self->throw("Sorry, you cannot read from a generic Bio::Assembly::IO object.");
}


=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL Assembly::IO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_format_module {
  my ($self,$format) = @_;
  my $module = "Bio::Assembly::IO::" . $format;
  my $ok;

  eval {
      $ok = $self->_load_module($module);
  };
  if ( $@ ) {
    print STDERR <<END;
$self: could not load $format - for more details on supported formats please see the Assembly::IO docs
Exception $@
END
  ;
  }
  return $ok;
}

=head2 _guess_format

 Title   : _guess_format
 Usage   : $obj->_guess_format($filename)
 Function: guess format based on file suffix
 Example :
 Returns : guessed format of filename (lower case)
 Args    :
 Notes   : formats that _filehandle() will guess includes
           only phrap, by now.

=cut

sub _guess_format {
   my $class = shift;
   my $arg   = shift;

   return unless defined($arg);
   return 'ace' if ($arg =~ /\.ace\.\d+$/i);
   return 'phrap' if ($arg =~ /\.phrap\.out$/i);
}

sub DESTROY {
    my $self = shift;

    $self->close();
}

# I need some direction on these!! The module works so I haven't fiddled with them!
# Me neither! (rfsouza)

sub TIEHANDLE {
    my ($class,$val) = @_;
    return bless {'seqio' => $val}, $class;
}

sub READLINE {
  my $self = shift;
  return $self->{'seqio'}->next_seq() unless wantarray;
  my (@list, $obj);
  push @list, $obj while $obj = $self->{'seqio'}->next_seq();
  return @list;
}

sub PRINT {
  my $self = shift;
  $self->{'seqio'}->write_seq(@_);
}

1;

