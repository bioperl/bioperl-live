# $Id: IO.pm 16690 2010-01-14 07:27:29Z kortsch $
#
# BioPerl module for Bio::Assembly::IO
#
#   based on the Bio::SeqIO module
#       by Ewan Birney <birney@ebi.ac.uk>
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

    while ( my $scaffold = $in->next_assembly() ) {
       # do something with Bio::Assembly::Scaffold instance
       # ...
       $out->write_assembly(-scaffold => $scaffold);
    }

    $in->close;
    $out->close;

=head1 DESCRIPTION

Bio::Assembly::IO is a handler module for formats in the Assembly::IO set
(e.g. Bio::Assembly::IO::phrap).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

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


use strict;

use base qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : $stream = Bio::Assembly::IO->new( -file   => $filename,
                                             -format =>'format'    )
 Function: Returns a new assembly stream
 Returns : A Bio::Assembly::IO::Handler initialised
           with the appropriate format
 Args    : -file   => $filename
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

        if ($format =~ /-/) {
            ($format, my $variant) = split('-', $format, 2);
            push @args, (-variant => $variant);
        }

        return unless( $class->_load_format_module($format) );
        return "Bio::Assembly::IO::$format"->new(@args);
    }
}


=head2 format

 Title   : format
 Usage   : $format = $stream->format()
 Function: Get the assembly format
 Returns : assembly format
 Args    : none

=cut

# format() method inherited from Bio::Root::IO


# _initialize is chained for all SeqIO classes
sub _initialize {
    my($self, @args) = @_;
    # initialize the IO part
    $self->_initialize_io(@args);
}


=head2 next_assembly

 Title   : next_assembly
 Usage   : $scaffold = $stream->next_assembly()
 Function: Reads the next assembly object from the stream and returns it.
 Returns : a Bio::Assembly::ScaffoldI compliant object
 Args    : none

=cut

sub next_assembly {
   my ($self) = @_;
   $self->throw("Cannot read from a generic Bio::Assembly::IO object.");
}


=head2 next_contig

 Title   : next_contig
 Usage   : $contig = $stream->next_contig()
 Function: Reads the next contig or singlet from the stream and returns it.
 Returns : a Bio::Assembly::Contig or Bio::Contig::Assembly::Singlet
 Args    : none

=cut

sub next_contig {
   my ($self) = @_;
   $self->throw("Cannot read from a generic Bio::Assembly::IO object.");
}


=head2 write_assembly

  Title   : write_assembly
  Usage   : $stream->write_assembly($assembly)
  Function: Write the assembly object in desired format. This method calls
            write_header(), write_contigs() and write_footer() internally.
  Returns : 1 on success, 0 for error
  Args    : A Bio::Assembly::Scaffold object

=cut

sub write_assembly {
    my ($self, @args) = @_;
    my ($scaf, $write_singlets) = $self->_rearrange([qw(SCAFFOLD SINGLETS)], @args);

    # Sanity check
    if ( !$scaf || !$scaf->isa('Bio::Assembly::ScaffoldI') ) {
        $self->throw("Must provide a Bio::Assembly::Scaffold object when calling write_assembly");
    }

    # Write header
    $self->write_header($scaf);

    # ID-sorted contig and read entries
    my @contig_ids  = $scaf->get_contig_ids;
    if ($write_singlets) {
        push @contig_ids, $scaf->get_singlet_ids;
    }
    @contig_ids = _sort(@contig_ids);

    # Write contigs
    for my $contig_id ( @contig_ids ) {
        my $contig = $scaf->get_contig_by_id($contig_id) ||
                     $scaf->get_singlet_by_id($contig_id);
        $self->write_contig($contig);
    }

    # Write footer
    $self->write_footer($scaf);

    return 1;      
}


=head2 write_header

  Title   : write_header
  Usage   : $stream->write_header($assembly)
  Function: Write the start of the assembly file. It can be called at any time,
            not when starting to write the assembly file. 
  Returns : 1 on success, 0 for error
  Args    : A Bio::Assembly::Scaffold object or ... (check the specific format
            driver for more details)

=cut

sub write_header {
   my ($self) = @_;
   $self->throw("Cannot write from a generic Bio::Assembly::IO object.");
}


=head2 write_contig

  Title   : write_contig
  Usage   : $stream->write_contig($contig)
  Function: Write a contig object in the desired format.
  Returns : 1 on success, 0 for error
  Args    : A Bio::Assembly::Contig object

=cut

sub write_contig {
   my ($self) = @_;
   $self->throw("Cannot write from a generic Bio::Assembly::IO object.");
}


=head2 write_footer

  Title   : write_footer
  Usage   : $stream->write_footer($assembly)
  Function: Write the start of the assembly file.
  Returns : 1 on success, 0 for error
  Args    : A Bio::Assembly::Scaffold object or ... (check the specific format
            driver for more details)

=cut

sub write_footer {
   my ($self) = @_;
   $self->throw("Cannot write from a generic Bio::Assembly::IO object.");
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
           ace, phrap and tigr at the moment

=cut

sub _guess_format {
   my $class = shift;
   my $arg   = shift;

   return unless defined($arg);
   return 'ace'    if ($arg =~ /\.ace/i);
   return 'phrap'  if ($arg =~ /\.phrap/i);
   return 'tigr'   if ($arg =~ /\.tigr/i);
   return 'maq'    if ($arg =~ /\.maq/i);
   return 'sam'    if ($arg =~ /\.[bs]am/i);
   return 'bowtie' if ($arg =~ /\.bowtie/i);

}


=head2 _sort

    Title   : _sort
    Usage   : @sorted_values = $ass_io->_sort(@values)
    Function: Sort a list of values naturally if Sort::Naturally is installed
              (nicer), lexically otherwise (not as nice, but safe)
    Returns : array of sorted values
    Args    : array of values to sort

=cut

sub _sort {
    my @arr = @_;
    my @sorted_arr;
    if (eval { require Sort::Naturally }) {
        @sorted_arr = Sort::Naturally::nsort( @arr ); # natural sort (better)
    } else {
        @sorted_arr = sort @arr; # lexical sort (safe)
    }
    return @sorted_arr;
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
