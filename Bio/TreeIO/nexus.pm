# $Id$
#
# BioPerl module for Bio::TreeIO::nexus
#
# Cared for by Jason Stajich <jason-at-open-bio-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::nexus - A TreeIO driver module for parsing Nexus tree output from PAUP

=head1 SYNOPSIS

  use Bio::TreeIO;
  my $in = new Bio::TreeIO(-file => 't/data/cat_tre.tre');
  while( my $tree = $in->next_tree ) {
  }

=head1 DESCRIPTION

This is a driver module for parsing PAUP Nexus tree format which
basically is just a remapping of trees.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-open-bio-dot-org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::nexus;
use vars qw(@ISA);
use strict;

use Bio::TreeIO;
use Bio::Event::EventGeneratorI;
use IO::String;

@ISA = qw(Bio::TreeIO );


=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $treeio->next_tree
 Function: Gets the next tree in the stream
 Returns : Bio::Tree::TreeI
 Args    : none


=cut

sub next_tree {
    my ($self) = @_;
    unless ( $self->{'_parsed'} ) { 
	$self->_parse;
    }
    return $self->{'_trees'}->[$self->{'_treeiter'}++];
}

sub rewind { 
    shift->{'_treeiter'} = 0;
}

sub _parse {
   my ($self) = @_;

   $self->{'_parsed'} = 1;
   $self->{'_treeiter'} = 0;

   while( defined ( $_ = $self->_readline ) ) {
       next if /^\s+$/;
       last;
   }
   return unless( defined $_ );
   
   unless( /^\#NEXUS/i ) {
       $self->warn("File does not start with #NEXUS"); #'
	   return;
   }
   my $state = 0;
   my %translate;
   while( defined ( $_ = $self->_readline ) ) {
       next if ( /^\[/);
       if( /^\[/ ) {
	   $state = 1 if $state != 0;
       } elsif( /^\]/ ) {
	   $state = 0 if $state == 1;
       } elsif( /^\s*Translate/i ) { 
	   $state = 3;
       } elsif( $state == 3 ) {
	   if( /^\s+(\S+)\s+([^\s\,\;]+)\s*([\,\;])?\s*$/ ) {
	       $translate{$1} = $2;
	       $state = 1 if( defined $3 && $3 eq ';' );
	   } elsif( /^\s+;/) {
	       $state = 1;
	   }
       } elsif( /^\s*tree\s+(\S+)\s+\=\s+(?:\[\S+\])?\s*([^\;]+;)\s*$/ ) {
	   my ($tree_name,$tree_str) = ($1,$2);

	   # MrBayes does not print colons for node label
	   # $tree_str =~ s/\)(\d*\.\d+)\)/:$1/g;
	   my $buf = new IO::String($tree_str);	   
	   my $treeio = new Bio::TreeIO(-format => 'newick',
					-fh     => $buf);
	   my $tree = $treeio->next_tree;
	   foreach my $node ( grep { $_->is_Leaf } $tree->get_nodes ) {
	       my $id = $node->id;
	       my $lookup = $translate{$id};
	       $node->id($lookup || $id);
	   }
	   $tree->id($tree_name) if defined $tree_name;
	   push @{$self->{'_trees'}},$tree;
       } elsif( /^\s*Begin(\s+trees)?;/io ) {
	   $state = 1;
       } elsif( /^\s*End(\s+trees)?;/io ) {
	   $state = 0;
	   return;
       } else { 
	   $self->debug("$state, $_");
       }
   }
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


1;
