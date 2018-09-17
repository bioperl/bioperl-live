#
# BioPerl module for Bio::Phenotype::Measure
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Christian M. Zmasek <czmasek-at-burnham.org> or <cmzmasek@yahoo.com>
#
# (c) Christian M. Zmasek, czmasek-at-burnham.org, 2002.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2002.
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
#
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Phenotype::Measure - Representation of context/value(-range)/unit triplets

=head1 SYNOPSIS

  use Bio::Phenotype::Measure;

  my $measure = Bio::Phenotype::Measure->new( -context     => "length",
                                              -description => "reduced length in 4(Tas1r3)",
                                              -start       => 0,
                                              -end         => 15,
                                              -unit        => "mm",
                                              -comment     => "see also Miller et al" );

  print $measure->context();
  print $measure->description();
  print $measure->start();
  print $measure->end();
  print $measure->unit();
  print $measure->comment();

  print $measure->to_string();

=head1 DESCRIPTION

Measure is for biochemically defined phenotypes or any other types of measures.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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

Christian M. Zmasek

Email: czmasek-at-burnham.org  or  cmzmasek@yahoo.com

WWW:   http://monochrome-effect.net/

Address: 

  Genomics Institute of the Novartis Research Foundation
  10675 John Jay Hopkins Drive
  San Diego, CA 92121

=head1 APPENDIX

The rest of the documentation details each of the object
methods.

=cut


# Let the code begin...

package Bio::Phenotype::Measure;
use strict;

use base qw(Bio::Root::Root);


=head2 new

 Title   : new
 Usage   : my $me = Bio::Phenotype::Measure->new( -context     => "length",
                                                  -description => "reduced length in 4(Tas1r3)",
                                                  -start       => 0,
                                                  -end         => 15,
                                                  -unit        => "mm",
                                                  -comment     => "see Miller also et al" );                      
 Function: Creates a new Measure object.
 Returns : A new Measure object.
 Args    : -context     => the context
           -description => a description
           -start       => the start value
           -end         => the end value
           -unit        => the unit
           -comment     => a comment

=cut

sub new {
    my( $class, @args ) = @_;
    
    my $self = $class->SUPER::new( @args );

    my ( $con, $desc, $start, $end, $unit, $comment )
    = $self->_rearrange( [ qw( CONTEXT
                               DESCRIPTION
                               START
                               END
                               UNIT
                               COMMENT ) ], @args );

    $self->init(); 
 
    $con     && $self->context( $con );
    $desc    && $self->description( $desc );
    $start   && $self->start( $start );
    $end     && $self->end( $end );
    $unit    && $self->unit( $unit );
    $comment && $self->comment( $comment );
                           
    return $self;
    
} # new




=head2 init

 Title   : init()
 Usage   : $measure->init();   
 Function: Initializes this Measure to all "".
 Returns : 
 Args    :

=cut

sub init {
    my( $self ) = @_;

    $self->context( "" );
    $self->description( "" );
    $self->start( "" );
    $self->end( "" );
    $self->unit( "" );
    $self->comment( "" );
  
} # init




=head2 context

 Title   : context
 Usage   : $measure->context( "Ca-conc" );
           or 
           print $measure->context(); 
 Function: Set/get for the context of this Measure.
 Returns : The context.
 Args    : The context (optional).

=cut

sub context {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_context" } = $value;
    }
   
    return $self->{ "_context" };
    
} # context




=head2 description

 Title   : description
 Usage   : $measure->description( "reduced in 4(Tas1r3)" );
           or 
           print $measure->description(); 
 Function: Set/get for the description of this Measure.
 Returns : A description.
 Args    : A description (optional).

=cut

sub description {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_description" } = $value;
    }
   
    return $self->{ "_description" };
    
} # description




=head2 start 

 Title   : start
 Usage   : $measure->start( 330 );
           or 
           print $measure->start(); 
 Function: Set/get for the start value of this Measure.
 Returns : The start value.
 Args    : The start value (optional).

=cut

sub start {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_start" } = $value;
    }
   
    return $self->{ "_start" };
    
} #  start




=head2 end 

 Title   : end 
 Usage   : $measure->end( 459 );
           or 
           print $measure->end(); 
 Function: Set/get for the end value of this Measure.
 Returns : The end value.
 Args    : The end value (optional).

=cut

sub end {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_end" } = $value;
    }
   
    return $self->{ "_end" };
    
} # end




=head2 unit

 Title   : unit
 Usage   : $measure->unit( "mM" );
           or 
           print $measure->unit(); 
 Function: Set/get for the unit of this Measure.
 Returns : The unit.
 Args    : The unit (optional).

=cut

sub unit {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_unit" } = $value;
    }
   
    return $self->{ "_unit" };
    
} # unit




=head2 comment

 Title   : comment
 Usage   : $measure->comment( "see also Miller et al" );
           or 
           print $measure->comment();
 Function: Set/get for an arbitrary comment about this Measure.
 Returns : A comment.
 Args    : A comment (optional).

=cut

sub comment {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_comment" } = $value;
    }
   
    return $self->{ "_comment" };
    
} # comment




=head2 to_string

 Title   : to_string()
 Usage   : print $measure->to_string();
 Function: To string method for Measure objects.
 Returns : A string representations of this Measure.
 Args    :

=cut

sub to_string {
    my ( $self ) = @_;

    my $s = "";
    
    $s .= "-- Context:\n";
    $s .= $self->context()."\n";
    $s .= "-- Description:\n";
    $s .= $self->description()."\n";
    $s .= "-- Start:\n";
    $s .= $self->start()."\n";
    $s .= "-- End:\n";
    $s .= $self->end()."\n";
    $s .= "-- Unit:\n";
    $s .= $self->unit()."\n";
    $s .= "-- Comment:\n";
    $s .= $self->comment();
    
    return $s;
    
} # to_string



1;
