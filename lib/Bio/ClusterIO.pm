#
# BioPerl module for Bio::ClusterIO.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Andrew Macgregor <andrew@anatomy.otago.ac.nz>
#
# Copyright Andrew Macgregor, Jo-Ann Stanton, David Green
# Molecular Embryology Group, Anatomy & Structural Biology, University of Otago
# http://anatomy.otago.ac.nz/meg
#
# You may distribute this module under the same terms as perl itself
#
# _history
#
# May 7, 2002 - changed from UniGene.pm to more generic ClusterIO.pm
# by Andrew Macgregor
#
# April 17, 2002 - Initial implementation by Andrew Macgregor
# POD documentation - main docs before the code

=head1 NAME

Bio::ClusterIO - Handler for Cluster Formats

=head1 SYNOPSIS

  #NB: This example is unigene specific

  use Bio::ClusterIO;

  $stream  = Bio::ClusterIO->new('-file' => "Hs.data", 
                                 '-format' => "unigene");
  # note: we quote -format to keep older perl's from complaining.

  while ( my $in = $stream->next_cluster() ) {
      print $in->unigene_id() . "\n";
      while ( my $sequence = $in->next_seq() ) {
          print $sequence->accession_number() . "\n";
      }
  }
  # Parsing errors are printed to STDERR.

=head1 DESCRIPTION

The ClusterIO module works with the ClusterIO format module to read
various cluster formats such as NCBI UniGene.


=head1 CONSTRUCTORS

=head2 Bio::ClusterIO-E<gt>new()

   $str = Bio::ClusterIO->new(-file => 'filename',
                              -format=>$format);

The new() class method constructs a new Bio::ClusterIO object.  The
returned object can be used to retrieve or print cluster
objects. new() accepts the following parameters:

=over 4

=item -file

A file path to be opened for reading.

=item -format

Specify the format of the file.  Supported formats include:

   unigene		*.data	UniGene build files.
   dbsnp		*.xml	dbSNP XML files

If no format is specified and a filename is given, then the module
will attempt to deduce it from the filename.  If this is unsuccessful,
the main UniGene build format is assumed.

The format name is case insensitive.  'UNIGENE', 'UniGene' and
'unigene' are all supported, as are dbSNP, dbsnp, and DBSNP

=back

=head1 OBJECT METHODS

See below for more detailed summaries.  The main methods are:

=head2 $cluster = $str-E<gt>next_cluster()

Fetch the next cluster from the stream.


=head2 TIEHANDLE(), READLINE(), PRINT()

These I've left in here because they were in the SeqIO
module. Feedback appreciated. There they provide the tie interface.
See L<perltie> for more details.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Andrew Macgregor

Email andrew@anatomy.otago.ac.nz

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...

package Bio::ClusterIO;

use strict;


use base qw(Bio::Root::Root Bio::Root::IO);



=head2 new

 Title   : new
 Usage   : Bio::ClusterIO->new(-file => $filename, -format => 'format')
 Function: Returns a new cluster stream
 Returns : A Bio::ClusterIO::Handler initialised with the appropriate format
 Args    : -file => $filename
           -format => format

=cut


my $entry = 0;

sub new {
    my ($caller,@args) = @_;
    my $class = ref($caller) || $caller;
    
    # or do we want to call SUPER on an object if $caller is an
    # object?
    if( $class =~ /Bio::ClusterIO::(\S+)/ ) {
	my ($self) = $class->SUPER::new(@args);	
	$self->_initialize(@args);
	return $self;
    } else { 

	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	my $format = $param{'-format'} || 
	    $class->_guess_format( $param{-file} || $ARGV[0] );
	$format = "\L$format";	# normalize capitalization to lower case

	return unless( $class->_load_format_module($format) );
	return "Bio::ClusterIO::$format"->new(@args);
    }
}

=head2 format

 Title   : format
 Usage   : $format = $stream->format()
 Function: Get the cluster format
 Returns : cluster format
 Args    : none

=cut

# format() method inherited from Bio::Root::IO


# _initialize is chained for all ClusterIO classes

sub _initialize {
    my($self, @args) = @_;
    # initialize the IO part
    $self->_initialize_io(@args);
}

=head2 next_cluster

 Title   : next_cluster
 Usage   : $cluster = $stream->next_cluster()
 Function: Reads the next cluster object from the stream and returns it.
 Returns : a L<Bio::ClusterI> compliant object
 Args    : none


=cut

sub next_cluster {
   my ($self, $seq) = @_;
   $self->throw("Sorry, you cannot read from a generic Bio::ClusterIO object.");
}

=head2 cluster_factory

 Title   : cluster_factory
 Usage   : $obj->cluster_factory($newval)
 Function: Get/set the object factory to use for creating the cluster
           objects.
 Example : 
 Returns : a L<Bio::Factory::ObjectFactoryI> compliant object
 Args    : on set, new value (a L<Bio::Factory::ObjectFactoryI> 
           compliant object or undef, optional)


=cut

sub cluster_factory{
    my $self = shift;

    return $self->{'cluster_factory'} = shift if @_;
    return $self->{'cluster_factory'};
}

=head2 object_factory

 Title   : object_factory
 Usage   : $obj->object_factory($newval)
 Function: This is an alias to cluster_factory with a more generic name.
 Example : 
 Returns : a L<Bio::Factory::ObjectFactoryI> compliant object
 Args    : on set, new value (a L<Bio::Factory::ObjectFactoryI> 
           compliant object or undef, optional)


=cut

sub object_factory{
    return shift->cluster_factory(@_);
}

=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL ClusterIO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_format_module {
  my ($self,$format) = @_;
  my $module = "Bio::ClusterIO::" . $format;
  my $ok;
  
  eval {
      $ok = $self->_load_module($module);
  };
  if ( $@ ) {
    print STDERR <<END;
$self: could not load $format - for more details on supported formats please see the ClusterIO docs
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
 Notes   : formats that _filehandle() will guess include unigene and dbsnp

=cut

sub _guess_format {
   my $class = shift;
   return unless $_ = shift;
   return 'unigene'   if /\.(data)$/i;
   return 'dbsnp'     if /\.(xml)$/i;
}

sub DESTROY {
    my $self = shift;

    $self->close();
}

# I need some direction on these!! The module works so I haven't fiddled with them!

sub TIEHANDLE {
    my ($class,$val) = @_;
    return bless {'seqio' => $val}, $class;
}

sub READLINE {
  my $self = shift;
  return $self->{'seqio'}->next_seq() || undef unless wantarray;
  my (@list, $obj);
  push @list, $obj while $obj = $self->{'seqio'}->next_seq();
  return @list;
}

sub PRINT {
  my $self = shift;
  $self->{'seqio'}->write_seq(@_);
}

1;

