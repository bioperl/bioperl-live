# $Id$
#
# BioPerl module for Bio::Tools::Phylo::PAML::Codeml
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich, Aaron J Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Phylo::PAML::Codeml - Parses output from the PAML program codeml. 

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

This object is used to parse the output from the PAML program codeml.
You can use the Bio::Tools::Run::Phylo::Phylo::PAML::Codeml module to RUN
codeml.

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
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich, Aaron Mackey

Email jason@bioperl.org
Email amackey@virginia.edu

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Phylo::PAML::Codeml;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Root::IO;
use Bio::TreeIO;
use IO::String;

@ISA = qw(Bio::Root::Root Bio::Root::IO );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Phylo::PAML::Codeml();
 Function: Builds a new Bio::Tools::Phylo::PAML::Codeml object 
 Returns : Bio::Tools::Phylo::PAML::Codeml
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->_initialize_io(@args);
  $self->_parse_mlc();
  return $self;
}

=head2 get_trees

 Title   : get_trees
 Usage   : my @trees = $codemlparser->get_trees();
 Function: Returns a list of trees (if any) are in the output file
 Returns : List of L<Bio::Tree::TreeI> objects
 Args    : none


=cut

sub get_trees{
   my ($self) = @_;


}

=head2 get_statistics

 Title   : get_statistics
 Usage   : my $data = $codemlparser->get_statistics
 Function: Retrieves the set of pairwise comparisons 
 Returns : Hash Reference keyed as 'seqname' -> 'seqname' -> 'datatype'
 Args    : none


=cut

sub get_statistics {
   my ($self) = @_;
   

}


# parse the mlc file

sub _parse_mlc {
    my ($self) = @_;
    my %data;
    while( defined( $_ = $self->_readline) ) {
	print;
	# Aaron this is where the parsing should begin

	# I'll do the Tree objects if you like - 
	# I'd do it by building an IO::String for the
	# the tree data 
	# or does it make more sense to parse this out of a collection of 
	# files?
	if( /^TREE/ ) {
	    # ...
	    while( defined($_ = $self->_readline) ) {
		if( /^\(/) {
		    my $treestr = new IO::String($_);
		    my $treeio = new Bio::TreeIO(-fh => $treestr,
						 -format => 'newick');
		    # this is very tenative here!!
		    push @{$self->{'_trees'}}, $treeio->next_tree;
		}
	    }
	}
    }
}

1;
