package Bio::Graphics::SimpleConfigurator;

# $Id$
#
# BioPerl module for Bio::Graphics::SimpleConfigurator
#
# Cared for by Paul Edlefsen <paul@systemsbiology.org>
#
# Copyright Paul Edlefsen
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Graphics::SimpleConfigurator - Interface for describing
Configurator objects in bioperl

=head1 SYNOPSIS

    # get a SimpleConfigurator somehow
    my $fg_color = $configurator->get('fgcolor');

=head1 DESCRIPTION

This object contains the various configuration parameters.  It
is devided up into sections and tags.  There is also the concept
of a default section which is referenced when no section is
passed to the object's methods.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org            - General discussion
http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Robert Hubley

Email rhubley@systemsbiology.org

=head1 CONTRIBUTORS

Paul Edlefsen, pedlefsen@systemsbiology.org
Lincoln Stein, lstein@cshl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

use strict;
use vars qw( @ISA );
use Bio::Root::Root;
use Bio::Graphics::ConfiguratorI;

@ISA = qw( Bio::Root::Root Bio::Graphics::ConfiguratorI );

use Carp;

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Graphics::SimpleConfigurator();
 Function: Builds a new Bio::Graphics::SimpleConfigurator object
 Returns : A Bio::Graphics::SimpleConfigurator
 Args    : None

=cut

sub new {
  my( $caller, @args ) = @_;

  my $self = $caller->SUPER::new( @args );
  $self->_initialize_simple_configurator( @args );
  return $self;
} # new(..)

sub _initialize_simple_configurator {
  my $self = shift;
  my @args = @_;

  return if( $self->{ '_simple_configurator_initialized' } );
  $self->{ '_sections' } = { '_general' => {} };
  $self->{ '_section_order' } = [];

  $self->{ '_simple_configurator_initialized' }++;
  return $self;
} # _initialize_simple_configurator(..)

=head2 get_sections

 Title   : get_sections
 Usage   : my @values = $configurator->get_sections();
 Function: Returns a list of the valid sections except
           the default or undef.
 Returns : A list of the sections which can be queried.
 Args    : (optional section as string, tag as string)

  The sections are guaranteed to be returned in the order in which
  they were added.

=cut

sub get_sections {
  return @{ shift->{ '_section_order' } };
} # get_sections(..)

=head2 get_tags

 Title   : get_tags
 Usage   : my @values = $configurator->get_tags();
           or
           my @values = $configurator->get_tags('dna');
 Function: Returns a list of tags for a given section
           or only the default tags section if no section
           is given.
 Returns : A scalar list of tags
 Args :    A section name or undef, 'general', or 'default' for default
           settings.

  Tags are NOT guaranteed to be returned in the order in which they were added.

=cut

sub get_tags {
  my ( $self, $section ) = @_;
  if ( defined $section ) {
    if( ( $section =~ /^default$/i ) || ( $section =~ /^general$/i ) ) {
      $section = '_general';
    }
    return ( keys %{ $self->{ '_sections' }->{ $section } } );
  } else {
    return ( keys %{ $self->{ '_sections' }->{ '_general' } } );
  }
} # get_tags(..)

=head2 get

 Title   : get
 Usage   : my $value = $configurator->get('height');
           or
           my $value = $configurator->get('dna','height');
 Function: Returns a tag value from a configurator from the
           either the default "_general" section or from
           a specified section or undef.
 Returns : A scalar value for the tag
 Args    : The tag name is required.  If there are two arguments then the
           first will be interpreted as the section name.  If it is
           undef or 'general' or 'default' then the default section
           will be used.

=cut

sub get {
  my ( $self, @params ) = @_;
  
  return unless ( defined @params );
  
  if ( $#params == 0 ) {
    return $self->{ '_sections' }->{ '_general' }->{ $params[ 0 ] };
  } elsif ( $#params == 1 ) {
    my $section = $params[ 0 ];
    if( ( $section =~ /^default$/i ) || ( $section =~ /^general$/i ) ) {
      $section = '_general';
    }
    return $self->{ '_sections' }->{ $section }->{ $params[ 1 ] };
  } else {
    return;
  }
} # get(..)

=head2 set

 Title   : set
 Usage   : $configurator->set('fgcolor','chartreuse');
           or
           $configurator->set('EST','fgcolor','chartreuse');
 Function: Set a value for a tag
 Returns : The old value of the tag
 Args :    The tag name and new value are required.  If there are two
           arguments then the first will be interpreted as the section
           name.  If it is undef or 'general' or 'default' then the
           default section will be used.

  If the given value is the string 'undef' then the tag will be
  undefined.  If a section has no defined tags then it will be
  removed.  If the named section is new then it will be added at the
  end of the list of sections.

=cut

sub set {
  my ( $self, @params ) = @_;

  return unless ( defined @params && ( $#params > 0 ) );

  my $old_value = "";
  if ( $#params == 1 ) {
    $old_value = $self->{ '_sections' }->{ '_general' }->{ $params[ 0 ] };
    if( $params[ 1 ] eq 'undef' ) {
      delete $self->{ '_sections' }->{ '_general' }->{ $params[ 0 ] };
    } else {
      $self->{ '_sections' }->{ '_general' }->{ $params[ 0 ] } = $params[ 1 ];
    }
  } elsif ( $#params == 2 ) {
    my $section = $params[ 0 ];
    if( ( $section =~ /^default$/i ) || ( $section =~ /^general$/i ) ) {
      $section = '_general';
    }
    if( $params[ 2 ] eq 'undef' ) {
      $old_value = $self->{ '_sections' }->{ $section }->{ $params[ 1 ] };
      delete $self->{ '_sections' }->{ $section }->{ $params[ 1 ] };
      if( ( $section ne '_general' ) &&
          !scalar( keys( %{ $self->{ '_sections' }->{ $section } } ) ) ) {
        # Erp, the section is no more.
        delete $self->{ '_sections' }->{ $section };
        @{ $self->{ '_section_order' } } =
          grep { $_ ne $section } @{ $self->{ '_section_order' } };
      }
    } else {
      if( ( $section ne '_general' ) &&
          !defined( $self->{ '_sections' }->{ $section } ) ) {
        # Keep track of the order.
        push( @{ $self->{ '_section_order' } }, $section );
        undef $old_value;
      } else {
        $old_value = $self->{ '_sections' }->{ $section }->{ $params[ 1 ] };
      }
      $self->{ '_sections' }->{ $section }->{ $params[ 1 ] } = $params[ 2 ];
    }
  } else {
    return;
  }
  return $old_value;
} # set(..)

=head2 get_and_eval

 Title   : get_and_eval
 Usage   : my $value = $configurator->get_and_eval('height');
           or
           my $value = $configurator->get_and_eval('dna','height');
 Function: This works like get() except that it is
           also able to evaluate code references.  These are
           options whose values begin with the characters
           "sub {".  In this case the value will be passed to
           an eval() and the resulting codereference returned.
 Returns : A value of the tag or undef.
 Args    : The tag name is required.  If there are two arguments then the
           first will be interpreted as the section name.  If it is
           undef or 'general' or 'default' then the default section
           will be used.

  THIS COULD BE DANGEROUS!  Arbitrarily eval'ing user code is unwise.
  You have been warned.

=cut

sub get_and_eval {
  my $val = shift->get( @_ );
  return $val unless $val =~ /^sub\s*\{/;
  my $coderef = eval $val;
  warn $@ if $@;
  return $coderef;
} # get_and_eval(..)

1;

__END__
