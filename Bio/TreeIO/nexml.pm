#
# BioPerl module for Bio::TreeIO::nexml
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chase Miller <chmille4@gmail.com>
#
# Copyright Chase Miller
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::nexml - A TreeIO driver module for parsing NeXML tree files

=head1 SYNOPSIS

  use Bio::TreeIO;
  my $in = Bio::TreeIO->new(-file => 'data.nexml' -format => 'Nexml');
  while( my $tree = $in->next_tree ) {
  }

=head1 DESCRIPTION

This is a driver module for parsing tree data in a NeXML format. For
more information on NeXML, visit L<http://www.nexml.org>.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Chase Miller

Email chmille4@gmail.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::TreeIO::nexml;
use strict;

use lib '../..';
use Bio::Event::EventGeneratorI;
use IO::String;
use Bio::Nexml::Factory;
use Bio::Phylo::IO qw (parse unparse);


use base qw(Bio::TreeIO);


sub _initialize {
    my $self = shift;
    $self->SUPER::_initialize(@_);
    $self->{_doc} = undef;
}

=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $treeio->next_tree
 Function: Gets the next tree in the stream
 Returns : L<Bio::Tree::TreeI>
 Args    : none


=cut

sub next_tree {
    my ($self) = @_;
    unless ( $self->{'_parsed'} ) {
        $self->_parse;
    }
    return $self->{'_trees'}->[ $self->{'_treeiter'}++ ];
}

=head2 doc

 Title   : doc
 Usage   : $treeio->doc
 Function: Returns the biophylo nexml document object
 Returns : Bio::Phylo::Project
 Args    : none or Bio::Phylo::Project object

=cut

sub doc {
	my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_doc'} = $value;
	}
	return $obj->{'_doc'};
}


=head2 rewind

 Title   : rewind
 Usage   : $treeio->rewind
 Function: Resets the stream
 Returns : none
 Args    : none

=cut

sub rewind {
    my $self = shift;
    $self->{'_treeiter'} = 0;
}

sub _parse {
    my ($self) = @_;
    
    $self->{'_parsed'}   = 1;
    $self->{'_treeiter'} = 0;
    my $fac = Bio::Nexml::Factory->new();
    
    
    $self->doc(parse(
 	'-file'       => $self->{'_file'},
 	'-format'     => 'nexml',
 	'-as_project' => '1'
 	));
 	
 	$self->{'_trees'} = $fac->create_bperl_tree($self);
}

=head2 write_tree

 Title   : write_tree
 Usage   : $treeio->write_tree($tree);
 Function: Writes a tree onto the stream
 Returns : none
 Args    : L<Bio::Tree::TreeI>


=cut

sub write_tree {
	my ($self, $bp_tree) = @_;
	
	my $fac = Bio::Nexml::Factory->new();
	my $taxa = $fac->create_bphylo_taxa($bp_tree);
	my ($tree) = $fac->create_bphylo_tree($bp_tree, $taxa);
	
	my $forest = Bio::Phylo::Factory->create_forest();
	$self->doc(Bio::Phylo::Factory->create_project());
	
	$forest->set_taxa($taxa);
	$forest->insert($tree);
	
	$self->doc->insert($forest);
	
	my $ret = $self->_print($self->doc->to_xml());
	$self->flush;
	return $ret;
}


1;
