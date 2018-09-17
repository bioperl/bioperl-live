#
# BioPerl module for Bio::TreeIO
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

Bio::TreeIO - Parser for Tree files

=head1 SYNOPSIS

  {
      use Bio::TreeIO;
      my $treeio = Bio::TreeIO->new(-format => 'newick',
                                    -file   => 'globin.dnd');
      while( my $tree = $treeio->next_tree ) {
		print "Tree is ", $tree->number_nodes, "\n";
      }
  }

=head1 DESCRIPTION

This is the driver module for Tree reading from data streams and
flatfiles.  This is intended to be able to create Bio::Tree::TreeI
objects.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 CONTRIBUTORS

Allen Day E<lt>allenday@ucla.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO;
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::TreeIO::TreeEventBuilder;

use base qw(Bio::Root::Root Bio::Root::IO Bio::Event::EventGeneratorI Bio::Factory::TreeFactoryI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::TreeIO->new();
 Function: Builds a new Bio::TreeIO object 
 Returns : Bio::TreeIO
 Args    : a hash.  useful keys:
   -format : Specify the format of the file.  Supported formats:

     newick             Newick tree format
     nexus              Nexus tree format
     nhx                NHX tree format
     svggraph           SVG graphical representation of tree
     tabtree            ASCII text representation of tree
     lintree            lintree output format

=cut

sub new {
  my($caller,@args) = @_;
  my $class = ref($caller) || $caller;

    # or do we want to call SUPER on an object if $caller is an
    # object?n

  my $obj;
  if( $class =~ /Bio::TreeIO::(\S+)/ ) {
    $obj = $class->SUPER::new(@args);
    $obj->_initialize(@args);
  } else {     
    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
    my $format = $param{'-format'} || 
      $class->_guess_format( $param{'-file'} || $ARGV[0] ) ||
      'newick';
    $format = "\L$format";  # normalize capitalization to lower case
      
    # normalize capitalization
    return undef unless( $class->_load_format_module($format) );
    $obj = "Bio::TreeIO::$format"->new(@args);
  }
  return $obj;
}


=head2 format

 Title   : format
 Usage   : $format = $obj->format()
 Function: Get the tree format
 Returns : tree format
 Args    : none

=cut

# format() method inherited from Bio::Root::IO


=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $treeio->next_tree;
 Function: Gets the next tree off the stream
 Returns : Bio::Tree::TreeI or undef if no more trees
 Args    : none

=cut

sub next_tree{
   my ($self) = @_;
   $self->throw("Cannot call method next_tree on Bio::TreeIO object must use a subclass");
}

=head2 write_tree

 Title   : write_tree
 Usage   : $treeio->write_tree($tree);
 Function: Writes a tree onto the stream
 Returns : none
 Args    : Bio::Tree::TreeI


=cut

sub write_tree{
   my ($self,$tree) = @_;
   $self->throw("Cannot call method write_tree on Bio::TreeIO object must use a subclass");
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
    $self->warn("Ignoring request to attach handler ".ref($handler). ' because it is not a Bio::Event::EventHandlerI');
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

    $self->get_params; # Initialize the default parameters.

    my ($nen,$ini) = $self->_rearrange
    ([qw(NEWLINE_EACH_NODE INTERNAL_NODE_ID)],@args);
    $self->set_param('newline_each_node',$nen);
    $self->set_param('internal_node_id',$ini);

    $self->attach_EventHandler(Bio::TreeIO::TreeEventBuilder->new
                   (-verbose => $self->verbose(), @args));
    $self->_initialize_io(@args);
    #$self->debug_params;
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
  my ($self,$format) = @_;
  my $module = "Bio::TreeIO::" . $format;
  my $ok;
  
  eval {
      $ok = $self->_load_module($module);
  };

  if ( $@ ) {
    print STDERR <<END;
$self: $format cannot be found
Exception $@
For more information about the TreeIO system please see the TreeIO docs.
This includes ways of checking for formats at compile time, not run time
END
  ;
  }
  return $ok;
}

sub param {
  my $self = shift;
  my $param = shift;
  my $value = shift;

  if (defined $value) {
    $self->get_params->{$param} = $value;
  }
  return $self->get_params->{$param};
}

sub set_param {
  my $self = shift;
  my $param = shift;
  my $value = shift;

  #print STDERR "[$param] -> [undef]\n" if (!defined $value);
  return unless (defined $value);
  #print STDERR "[$param] -> [$value]\n";

  $self->get_params->{$param} = $value;
  return $self->param($param);
}

sub params {
  my $self = shift;
  return $self->get_params;
}
sub get_params {
  my $self = shift;

  if (!defined $self->{_params}) {
    $self->{_params} = $self->get_default_params;
  }

  return $self->{_params};
}

sub set_params {
  my $self = shift;
  my $params = shift;

  # Apply all the passed parameters to our internal parm hashref.
  my $cur_params = $self->get_params;
  $self->{_params} = { %$cur_params, %$params };

  return $self->get_params;
}

sub get_default_params {
  my $self = shift;
  
  return {};
}

sub debug_params {
  my $self = shift;

  my $params = $self->get_params;

  print STDERR "{\n";
  foreach my $param (keys %$params) {
    my $value = $params->{$param};
    print STDERR "  [$param] -> [$value]\n";
  }
  print STDERR "}\n";
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
   return 'nhx'   if /\.(nhx)$/i;
   return 'phyloxml' if /\.(xml)$/i;
   return 'svggraph' if /\.svg$/i;
   return 'lintree'  if( /\.(lin|lintree)$/i );
}

sub DESTROY {
    my $self = shift;

    $self->close();
}

sub TIEHANDLE {
  my $class = shift;
  return bless {'treeio' => shift},$class;
}

sub READLINE {
  my $self = shift;
  return $self->{'treeio'}->next_tree() || undef unless wantarray;
  my (@list,$obj);
  push @list,$obj  while $obj = $self->{'treeio'}->next_tree();
  return @list;
}

sub PRINT {
  my $self = shift;
  $self->{'treeio'}->write_tree(@_);
}

1;
