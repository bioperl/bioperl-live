# $Id: TranscriptionFactor.pm,v 1.6 2006/07/17 14:16:53 sendu Exp $
#
# BioPerl module for Bio::Map::TranscriptionFactor
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
# 
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::TranscriptionFactor - A transcription factor modelled as a mappable
element

=head1 SYNOPSIS

  use Bio::Map::TranscriptionFactor;
  use Bio::Map::GeneMap;
  use Bio::Map::Position;

  # model a TF that binds 500bp upstream of the BRCA2 gene in humans and
  # 250bp upstream of BRCA2 in mice
  my $tf = Bio::Map::TranscriptionFactor->get(-universal_name => 'tf1');
  my $map1 = Bio::Map::GeneMap->get(-universal_name => "BRCA2",
                                    -species => "human");
  my $map2 = Bio::Map::GeneMap->get(-universal_name => "BRCA2",
                                    -species => "mouse");
  Bio::Map::Position->new(-map => $map1,
                          -element => $tf,
                          -start => -500,
                          -length => 10);
  Bio::Map::Position->new(-map => $map2,
                          -element => $tf,
                          -start => -250,
                          -length => 10);

  # Find out where the transcription factor binds
  foreach $pos ($tf->get_positions) {
    print $tf->universal_name, " binds at position " $pos->value, " relative to ",
          $pos->relative->description, " of gene ",
          $pos->map->universal_name, " in species ", $pos->map->species, "\n";
  }

=head1 DESCRIPTION

A transcription factor modelled as a mappable element. It can have mulitple
binding sites (positions) near multiple genes (maps).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

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

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::TranscriptionFactor;
use strict;

use base qw(Bio::Map::Mappable);

our $TFS = {};

=head2 new

 Title   : new
 Usage   : my $tf = Bio::Map::TranscriptionFactor->new();
 Function: Builds a new Bio::Map::TranscriptionFactor object
 Returns : Bio::Map::TranscriptionFactor
 Args    : -universal_name => string name of the TF (in a form common to all
                              species that have the TF, but unique amongst
                              non-orthologous TFs), REQUIRED
           -description    => string, free text description of the TF

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($u_name, $desc) = $self->_rearrange([qw(UNIVERSAL_NAME DESCRIPTION)], @args);
    $u_name || $self->throw("You must supply a -universal_name");
    $self->universal_name($u_name);
    
    defined $desc && $self->description($desc);
    
    return $self;
}

=head2 get

 Title   : get
 Usage   : my $obj = Bio::Map::TranscriptionFactor->get();
 Function: Builds a new Bio::Map::TranscriptionFactor object (like new()), or
           gets a pre-existing one that shares the same universal_name.
 Returns : Bio::Map::TranscriptionFactor
 Args    : -universal_name => string name of the TF (in a form common to all
                              species that have the TF, but unique amongst
                              non-orthologous TFs), REQUIRED
           -description    => string, free text description of the TF

=cut

sub get {
    my ($class, @args) = @_;
    my ($u_name) = Bio::Root::Root->_rearrange([qw(UNIVERSAL_NAME)], @args);
    
    if ($u_name && defined $TFS->{$u_name}) {
        return $TFS->{$u_name};
    }
    
    return $class->new(@args);
}

=head2 universal_name

 Title   : universal_name
 Usage   : my $name = $obj->universal_name
 Function: Get/Set TF name, corresponding to the name of the TF in a form shared
           by orthologous versions of the TF in different species, but otherwise
           unique.
 Returns : string
 Args    : none to get, OR string to set

=cut

sub universal_name {
    my ($self, $value) = @_;
    if (defined $value) {
        delete $TFS->{$self->{'_uname'}} if $self->{'_uname'};
        $self->{'_uname'} = $value;
        $TFS->{$value} = $self;
    }
    return $self->{'_uname'};
}

=head2 description

 Title   : description
 Usage   : my $desc = $obj->description
 Function: Get/Set a description of the TF.
 Returns : string
 Args    : none to get, OR string to set

=cut

sub description {
    my $self = shift;
    if (@_) { $self->{desc} = shift }
    return $self->{desc} || '';
}

1;
