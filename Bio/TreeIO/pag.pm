# $Id$
#
# BioPerl module for Bio::TreeIO::pag
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::pag - Bio::TreeIO driver for Pagel format

=head1 SYNOPSIS

  use Bio::TreeIO;
  my $in = Bio::TreeIO->new(-format => 'nexus',
                            -file   => 't/data/adh.mb_tree.nexus');

  my $out = Bio::TreeIO->new(-format => 'pag');
  while( my $tree = $in->next_tree ) {
    $out->write_tree($tree);
  }
  
=head1 DESCRIPTION

Convert a Bio::TreeIO to Pagel format. Currenty  
More information here http://sapc34.rdg.ac.uk/meade/Mark/

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

  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TreeIO::pag;
use vars qw(@ISA);
use strict;

use Bio::TreeIO;

@ISA = qw(Bio::TreeIO );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::TreeIO::pag();
 Function: Builds a new Bio::TreeIO::pag object 
 Returns : an instance of Bio::TreeIO::pag
 Args    : -file/-fh for filename or filehandles
 

=cut


=head2 write_tree

 Title   : write_tree
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub write_tree{
   my ($self,@args) = @_;
   for my $t ( @args ) {
       for my $node ( $t->get_nodes ) {
	   next unless defined $node->ancestor; #skip root
	   # for now assumes that characters have been stored
	   # as tag-values
	   my @tags = sort $node->get_all_tags;
	   my @charstates = map { ($node->get_tag_values($_))[0] } @tags;
	   $self->_print(join(", ", ($node->id || 
				     sprintf("node%d",$node->internal_id)),
			      sprintf("node%d",$node->ancestor->internal_id),
			      sprintf("%.6f",$node->branch_length),
			      @charstates),"\n");
       }
   }
}

=head2 next_tree

 Title   : next_tree
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub next_tree{
   my ($self,@args) = @_;
   $self->throw("parsing not implemented yet");
}


1;
