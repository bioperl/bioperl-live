#
#  BioPerl module for Bio::Assembly::ScaffoldI
#
# Copyright by Robson F. de Souza
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::ScaffoldI - Abstract Inteface of Sequence Assemblies

=head1 SYNOPSIS

    # get a Bio::Assembly::ScaffoldI object somehow

    foreach my $contig ($assembly->all_contigs) {
       # do something (see Bio::Assembly::Contig)
    }

=head1 DESCRIPTION

This interface defines the basic set of methods an object should have
to manipulate assembly data.

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

=head1 AUTHOR - Robson Francisco de Souza

Email: rfsouza@citri.iq.usp.br

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#
# Now, let's code!


package Bio::Assembly::ScaffoldI;

use strict;
use Carp;

# Inheritance

use base qw(Bio::Root::RootI);

#
# Methods

=head1 Accessing general assembly data

=cut

=head2 get_nof_contigs

	Title   : get_nof_contigs
	Usage   : $assembly->get_nof_contigs()
	Function: Get the number of contigs included in the assembly
	Returns : integer
	Args    : none

=cut

sub get_nof_contigs {
    my $self = shift;

    $self->throw_not_implemented();
}

=head2 get_nof_singlets

	Title   : get_nof_singlets
	Usage   : $assembly->get_nof_singlets()
	Function: Get the number of singlets included in the assembly
	Returns : integer
	Args    : none

=cut

sub get_nof_singlets {
    my $self = shift;

    $self->throw_not_implemented();
}

=head2 get_contig_ids

	Title   : get_contig_ids
	Usage   : $assembly->get_contig_ids()
	Function: Access list of contig IDs from assembly
	Returns : an array if there are any contigs in the assembly.
                  undef otherwise
	Args    : an array of contig IDs

=cut

sub get_contig_ids {
    my $self = shift;

    $self->throw_not_implemented();
}

=head2 get_singlet_ids

	Title   : get_singlet_ids
	Usage   : $assembly->get_singlet_ids()
	Function: Access list of singlet IDs from assembly
	Returns : an array if there are any singlets in the assembly.
                  undef otherwise
	Args    : an array of singlet IDs

=cut

sub get_singlet_ids {
    my $self = shift;

    $self->throw_not_implemented();
}

=head2 get_contig_by_id

    Title   : get_contig_by_id
    Usage   : $assembly->get_contig_by_id($id)
    Function: Get a reference for a contig from the assembly
    Returns : a Bio::Assembly::Contig object or undef
    Args    : [string] contig unique identifier (ID)

=cut

sub get_contig_by_id {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 get_singlet_by_id

    Title   : get_singlet_by_id
    Usage   : $assembly->get_singlet_by_id()
    Function: Get a reference for a singlet from the assembly
    Returns : Bio::Assembly::Singlet object or undef
    Args    : [string] a singlet ID

=cut

sub get_singlet_by_id {
    my $self = shift;
    $self->throw_not_implemented();
}

=head1 Modifier methods

Implementation of these methods is optional in the sense that
read-only implementations may not have these. If an object implements
one of them, it should however implement all.

=cut

=head2 add_contig

	Title   : add_contig
	Usage   : $assembly->add_contig($contig)
	Function: Add another contig to the Bio::Assembly::ScaffoldI object
	Returns : 1 on success, 0 otherwise
	Args    : a Bio::Assembly:Contig object

    See Bio::Assembly::Contig for more information

=cut

#---------------------
sub add_contig {
#---------------------
	my ($self) = @_;
	$self->throw_not_implemented();
}

=head2 add_singlet

	Title   : add_singlet
	Usage   : $assembly->add_singlet($seq)
	Function: Add another singlet to the Bio::Assembly::ScaffoldI object
	Returns : 1 on success, 0 otherwise
	Args    : a Bio::Assembly::Singlet object

=cut

#---------------------
sub add_singlet {
#---------------------
	my ($self) = @_;
	$self->throw_not_implemented();
}

=head2 remove_contigs

        Title   : remove_contigs
	Usage   : $assembly->remove_contigs(1..4)
	Function: Remove contig from assembly object
	Returns : a Bio::Assembly::Contig object
	Args    : a list of contig IDs

    See function get_contig_ids() above

=cut

#---------------------
sub remove_contigs {
#---------------------
	my ($self) = @_;
	$self->throw_not_implemented();
}

=head2 remove_singlets

        Title   : remove_singlets
	Usage   : $assembly->remove_singlets(1..4)
	Function: Remove singlets from assembly object
	Returns : an array of Bio::Assembly::Singlet objects
	Args    : an array of singlet IDs 

    See function get_singlet_ids() above

=cut

#---------------------
sub remove_singlets {
#---------------------
	my ($self) = @_;
	$self->throw_not_implemented();
}

=head1 Contig and singlet selection methos

=cut

=head2 select_contigs

	Title   : select_contig
	Usage   : $assembly->select_contig
	Function: Selects an array of contigs from the assembly
	Returns : an array of Bio::Assembly::Contig objects
	Args    : an array of contig ids

    See function get_contig_ids() above

=cut

#---------------------
sub select_contigs {
#---------------------
	my ($self) = @_;
	$self->throw_not_implemented();
}

=head2 select_singlets

	Title   : select_singlets
	Usage   : $assembly->select_singlets(@list)
	Function: Selects an array of singlets from the assembly
	Returns : an array of Bio::Assembly::Singlet objects
	Args    : an array of singlet ids

    See function get_singlet_ids() above

=cut

#---------------------
sub select_singlets {
#---------------------
	my ($self) = @_;
	$self->throw_not_implemented();
}

=head2 all_contigs

	Title   : all_contigs
	Usage   : my @contigs = $assembly->all_contigs
	Function: Returns a list of all contigs in this assembly.
		  Contigs are both clusters and alignments of one
		  or more reads, with an associated consensus
		  sequence.
	Returns : array of Bio::Assembly::Contig
	Args    : none

=cut

#---------------------
sub all_contigs {
#---------------------
	my ($self) = @_;
	$self->throw_not_implemented();
}

=head2 all_singlets

    Title   : all_singlets
    Usage   : my @singlets = $assembly->all_singlets
    Function: Returns a list of all singlets in this assembly.
	      Singlets are isolated reads, without non-vector
	      matches to any other read in the assembly.
    Returns : array of Bio::Assembly::Singlet objects
    Args    : none

=cut

#---------------------
sub all_singlets {
#---------------------
	my ($self) = @_;
	$self->throw_not_implemented();
}

1;
