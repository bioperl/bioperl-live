#
# BioPerl module for Bio::Phenotype::OMIM::MiniMIMentry
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

Bio::Phenotype::OMIM::MiniMIMentry - Representation of a Mini MIM entry

=head1 SYNOPSIS

  use Bio::Phenotype::OMIM::MiniMIMentry;

  $mm = Bio::Phenotype::OMIM::MiniMIMentry->new( -description  => "The central form of ...",
                                                 -created      => "Victor A. McKusick: 6/4/1986",
                                                 -contributors => "Kelly A. Przylepa - revised: 03/18/2002",
                                                 -edited       => "alopez: 06/03/1997" );


=head1 DESCRIPTION

This class representats of Mini MIM entries.
This class is intended to be used together with a OMIM entry class. 


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

package Bio::Phenotype::OMIM::MiniMIMentry;
use strict;

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : $mm = Bio::Phenotype::OMIM::MiniMIMentry->new( -description  => "The central form of ...",
                                                          -created      => "Victor A. McKusick: 6/4/1986",
                                                          -contributors => "Kelly A. Przylepa - revised: 03/18/2002",
                                                          -edited       => "alopez: 06/03/1997" );

 Function: Creates a new MiniMIMentry object.
 Returns : A new MiniMIMentry object.
 Args    : -description  => a description
           -created      => name(s) and date(s) (free form)
           -contributors => name(s) and date(s) (free form)
           -edited       => name(s) and date(s) (free form)

=cut

sub new {

    my( $class, @args ) = @_;
    
    my $self = $class->SUPER::new( @args );

    my ( $desc, $created, $contributors, $edited )
    = $self->_rearrange( [ qw( DESCRIPTION
                               CREATED
                               CONTRIBUTORS
                               EDITED ) ], @args );

    $self->init(); 

    $desc         && $self->description( $desc );
    $created      && $self->created( $created );
    $contributors && $self->contributors( $contributors );
    $edited       && $self->edited( $edited );

    return $self;

} # new




=head2 init

 Title   : init()
 Usage   : $mm->init();   
 Function: Initializes this MiniMIMentry to all "".
 Returns : 
 Args    :

=cut

sub init {

    my( $self ) = @_;
   
    $self->description( "" );
    $self->created( "" );
    $self->contributors( "" );
    $self->edited( "" );
    
  
} # init




=head2 description

 Title   : description
 Usage   : $mm->description( "The central form of ..." );
           or
           print $mm->description();
 Function: Set/get for the description field of the Mini MIM database.
 Returns : The description.
 Args    : The description (optional).

=cut

sub description {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_description" } = $value;
    }

    return $self->{ "_description" };

} # description




=head2 created

 Title   : created
 Usage   : $mm->created( "Victor A. McKusick: 6/4/1986" );
           or
           print $mm->created();
 Function: Set/get for the created field of the Mini MIM database.
 Returns : Name(s) and date(s) [scalar - free form].
 Args    : Name(s) and date(s) [scalar - free form] (optional).

=cut

sub created {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_created" } = $value;
    }

    return $self->{ "_created" };

} # created




=head2 contributors

 Title   : contributors
 Usage   : $mm->contributors( "Kelly A. Przylepa - revised: 03/18/2002" );
           or
           print $mm->contributors();
 Function: Set/get for the contributors field of the Mini MIM database.
 Returns : Name(s) and date(s) [scalar - free form].
 Args    : Name(s) and date(s) [scalar - free form] (optional).

=cut

sub contributors {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_contributors" } = $value;
    }

    return $self->{ "_contributors" };

} # contributors




=head2 edited

 Title   : edited
 Usage   : $mm->edited( "alopez: 06/03/1997" );
           or
           print $mm->edited();
 Function: Set/get for the edited field of the Mini MIM database.
 Returns : Name(s) and date(s) [scalar - free form].
 Args    : Name(s) and date(s) [scalar - free form] (optional).

=cut

sub edited {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_edited" } = $value;
    }

    return $self->{ "_edited" };

} # edited




=head2 to_string

 Title   : to_string()
 Usage   : print $mm->to_string();
 Function: To string method for MiniMIMentry objects.
 Returns : A string representations of this MiniMIMentry.
 Args    :

=cut

sub to_string {
    my ( $self ) = @_;

    my $s = "";
    
    $s .= "-- Description:\n";
    $s .= $self->description()."\n";
    $s .= "-- Created:\n";
    $s .= $self->created()."\n";
    $s .= "-- Contributors:\n";
    $s .= $self->contributors()."\n";
    $s .= "-- Edited:\n";
    $s .= $self->edited();
  
    return $s;
    
} # to_string 


1;
