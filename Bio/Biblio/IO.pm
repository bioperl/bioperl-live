# $Id$
#
# BioPerl module for Bio::Biblio::IO
#
# Cared for by Martin Senger <senger@ebi.ac.uk>
# For copyright and disclaimer see below.

# POD documentation - main docs before the code

=head1 NAME

Bio::Biblio::IO - Handling the bibliographic references

=head1 SYNOPSIS

    use Bio::Biblio::IO;

    $in = Bio::Biblio::IO->new(-fh => $io_handle , '-format' => 'medlinexml');

    while ( my $citation = $in->next() ) {
	print $citation;
    }

Now, to actually get an I<$io_handle> representing a citation in an XML format,
use I<Bio::Biblio> module which returns an XML string and convert this string
into an IO handle:

    use Bio::Biblio;
    my $xml = new Bio::Biblio->get_by_id ('94033980');
    my $reader = Bio::Biblio::IO->new ('-fh' => IO::String->new ($xml),
				       '-format' => 'medlinexml',
				       '-result' => 'bibref');
    while (my $citation = $reader->next()) {
       ... do something here with $citation
       }

=head1 DESCRIPTION

 # to be written

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/MailList.shtml      - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR

Martin Senger (senger@ebi.ac.uk)

=head1 COPYRIGHT

Copyright (c) 2002 European Bioinformatics Institute. All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are preceded with a _

=cut


# Let the code begin...

package Bio::Biblio::IO;

use strict;
use vars qw(@ISA);

use Bio::Root::Root;
use Bio::Root::IO;
use Symbol();

@ISA = qw(Bio::Root::Root Bio::Root::IO);

my $entry = 0;

sub new {
    my ($caller, @args) = @_;
    my $class = ref ($caller) || $caller;
    
    # if $caller is an object, or if it is an underlying
    # 'real-work-doing' class (e.g. Bio::Biblio::IO::medlinexml) then
    # we want to call SUPER to create and bless an object
    if( $class =~ /Bio::Biblio::IO::(\S+)/ ) {
	my ($self) = $class->SUPER::new (@args);	
	$self->_initialize (@args);
	return $self;

    # this is called only the first time when somebody calls: 'new
    # Bio::Biblio::IO (...)', and it actually loads a 'real-work-doing'
    # module and call this new() method again (unless the loaded
    # module has its own new() method)
    } else { 

	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	my $format = $param{'-format'} || 
	    $class->_guess_format( $param{-file} || $ARGV[0] ) ||
		'medlinexml';
	$format = "\L$format";	# normalize capitalization to lower case

	# load module with the real implementation - as defined in $format
	return undef unless (&_load_format_module ($format));

	# this will call this same method new() - but rather its
	# upper (object) branche
	return "Bio::Biblio::IO::$format"->new(@args);
    }
}

sub newFh {
  my $class = shift;
  return unless my $self = $class->new(@_);
  return $self->fh;
}


sub fh {
  my $self = shift;
  my $class = ref($self) || $self;
  my $s = Symbol::gensym;
  tie $$s,$class,$self;
  return $s;
}

# _initialize is chained for all Bio::Biblio::IO classes

sub _initialize {
    my ($self, @args) = @_;
    # initialize the IO part
    $self->_initialize_io (@args);
}

=head2 next_bibref

 Usage   : $citation = stream->next_bibref
 Function: Reads the next citation object from the stream and returns it.

 Returns : a Bio::Biblio::IO sequence object
 Args    : none

=cut

sub next_bibref {
   my ($self) = shift;
   $self->throw ("Sorry, you cannot read from a generic Bio::Biblio::IO object.");
}

# -----------------------------------------------------------------------------

=head2 _load_format_module

 Usage   : $class->_load_format_module ($format)
 Returns : 1 on success, undef on failure
 Args    : 'format' should contain the last part of the
           name of a module who does the real implementation

It does (in run-time) a similar thing as

   require Bio::Biblio::IO::$format

It prints an error on STDERR if it fails to find and load the module
(for example, because of the compilation errors in the module).

=cut

sub _load_format_module {
  my ($format) = @_;
  my ($module, $load, $m);

  $module = "_<Bio/Biblio/IO/$format.pm";
  $load = "Bio/Biblio/IO/$format.pm";

  return 1 if $main::{$module};
  eval {
    require $load;
  };
  if ( $@ ) {
    print STDERR <<END;
$load: $format cannot be found or loaded
Exception $@
For more information about the Biblio system please see the Bio::Biblio::IO docs.
This includes ways of checking for formats at compile time, not run time
END
  ;
    return;
  }
  return 1;
}

=head2 _guess_format

 Usage   : $class->_guess_format ($filename)
 Returns : string with a guessed format of the input data (e.g. 'medlinexml')
 Args    : a file name whose extension can help to guess its format

It makes an expert guess what kind of data are in the given file
(but be prepare that $filename may be empty).

=cut

sub _guess_format {
   my $class = shift;
   return unless $_ = shift;
   return 'medlinexml'   if /\.(xml|medlinexml)$/i;
   return;
}

sub DESTROY {
    my $self = shift;

    $self->close();
}

#sub TIEHANDLE {
#    my ($class,$val) = @_;
#    return bless {'seqio' => $val}, $class;
#}

#sub READLINE {
#  my $self = shift;
#  return $self->{'seqio'}->next_seq() unless wantarray;
#  my (@list, $obj);
#  push @list, $obj while $obj = $self->{'seqio'}->next_seq();
#  return @list;
#}

1;
