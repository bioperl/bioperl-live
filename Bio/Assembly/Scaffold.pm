# $Id$
#
#  BioPerl module for Bio::Assembly::Scaffold
#
# Copyright by Robson F. de Souza
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::Scaffold - Perl module to hold and manipulate sequence assembly data.

=head1 SYNOPSYS

    # Module loading
    use Bio::Assembly::IO;

    # Assembly loading methods
    my $aio = new Bio::Assembly::IO(-file=>"test.ace.1", -format=>'phrap');
    my $assembly = $aio->next_assembly;

    foreach my $contig ($assembly->all_contigs) {
        # do something... (see Bio::Assembly::Contig)
    }

=head1 DESCRIPTION

Bio::Assembly::Scaffold was developed to store and manipulate data
from sequence assembly programs like Phrap. It implements the
ScaffoldI interface and intends to be generic enough to be used by
Bio::Assembly::IO drivers written to programs other than Phrap.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bio.perl.org/MailList.html     - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Robson Francisco de Souza

rfsouza@citri.iq.usp.br

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Assembly::Scaffold;

use strict;
use vars qw(@ISA $VERSION);

use Bio::Root::Root;
use Bio::Assembly::ScaffoldI;
use Bio::Annotation::Collection;

$VERSION = '0.0.1';
@ISA = qw(Bio::Root::Root Bio::Assembly::ScaffoldI);

=head2 new ()

    Title   : new
    Usage   : $assembly = new (-source=>'program_name',
			       -contigs=>\@contigs,
			       -id=>"assembly 1");
    Function: creates a new assembly object
    Returns : 
    Args    : 
              -source  : [string] sequence assembly program
              -contigs : reference to array of 
                         Bio::Assembly::Contig objects
              -id      : [string] assembly name

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($src,$contigs,$id) = $self->_rearrange([qw(SOURCE CONTIGS ID)], @args);

  $self->{'_contigs'} = {};
  $self->{'_singlets'} = {};
  $self->{'_seqs'} = {};
  $self->{'_annotation'} = Bio::Annotation::Collection->new();
  $self->{'_id'} = 'NoName';

  if (defined $contigs && ref($contigs = 'ARRAY')) {
      foreach my $contig (@{$contigs}) {
	  $self->add_contig($contig);
      }
  }

  $self->{'_id'} = $id if (defined $id);

  return $self;
}

=head1 Accessing general assembly data

=head2 id

    Title   : id
    Usage   : $assembly->id()
    Function: Get/Set assembly ID
    Returns : string or undef
    Args    : string

=cut

sub id {
    my $self = shift;
    my $id   = shift;

    $self->{'_id'} = $id if (defined $id);

    return $self->{'_id'};
}

=head2 annotation

    Title   : annotation
    Usage   : $assembly->annotation()
    Function: Get/Set assembly annotation object
    Returns : Bio::Annotation::Collection
    Args    : none

=cut

sub annotation {
    my ($self,$ref) = shift;

    $self->{'_annotation'} = $ref if (defined $ref);
    return $self->{'_annotation'};
}

=head2 get_nof_contigs

    Title   : get_nof_contigs
    Usage   : $assembly->get_nof_contigs()
    Function: Get the number of contigs included in the assembly
    Returns : integer
    Args    : none 

=cut

sub get_nof_contigs {
    my $self = shift;

    return scalar( $self->get_contig_ids() );
}

=head2 get_nof_sequences_in_contigs

    Title   : get_nof_sequences_in_contigs
    Usage   : $assembly->get_nof_sequences_in_contigs()
    Function: 

              Get the number of sequences included in the
              assembly. This number refers only to the sequences used
              to build contigs in this assembly. It does not includes
              contig consensus sequences or singlets.

    Returns : integer
    Args    : none

=cut

sub get_nof_sequences_in_contigs {
    my $self = shift;

    my $nof_seqs = 0;
    foreach my $contig ($self->all_contigs) {
	$nof_seqs += scalar( $contig->get_seq_ids() );
    }

    return $nof_seqs;
}

=head2 get_nof_singlets

    Title   : nof_singlets
    Usage   : $assembly->nof_singlets()
    Function: Get the number of singlets included in the assembly
    Returns : integer
    Args    : none

=cut

sub get_nof_singlets {
    my $self = shift;

    return scalar( $self->get_singlet_ids() );
}

=head2 get_seq_ids

    Title   : get_seq_ids
    Usage   : $assembly->get_seq_ids()
    Function: 

              Get the ID of sequences from all contigs.  This list
              refers only to the aligned sequences in contigs. It does
              not includes contig consensus sequences or singlets.

    Returns : array of strings
    Args    : none

=cut

sub get_seq_ids {
    my $self = shift;

    return keys %{ $self->{'_seqs'} };
}

=head2 get_contig_ids

    Title   : get_contig_ids
    Usage   : $assembly->get_contig_ids()
    Function: Access list of contig IDs from assembly
    Returns : an array, if there are any contigs in the
              assembly. An empty array otherwise
    Args    : none

=cut

sub get_contig_ids {
    my $self = shift;

    return sort keys %{$self->{'_contigs'}};
}

=head2 get_singlet_ids

    Title   : get_singlet_ids
    Usage   : $assembly->get_singlet_ids()
    Function: Access list of singlet IDs from assembly
    Returns : array of strings if there are any singlets
              otherwise an empty array
    Args    : none

=cut

sub get_singlet_ids {
    my $self = shift;

    return sort keys %{$self->{'_singlets'}};
}

=head2 get_seq_by_id

    Title   : get_seq_by_id
    Usage   : $assembly->get_seq_by_id($id)
    Function: 

              Get a reference for an aligned sequence
              This sequence must be part of a contig
              in the assembly.

    Returns : a Bio::LocatableSeq object
              undef if sequence $id is not found
              in any contig
    Args    : [string] sequence identifier (id)

=cut

sub get_seq_by_id {
    my $self = shift;
    my $seqID = shift;

    return undef unless (exists $self->{'_seqs'}{$seqID});

    return $self->{'_seqs'}{$seqID}->get_seq_by_name($seqID);
}

=head2 get_contig_by_id

    Title   : get_contig_by_id
    Usage   : $assembly->get_contig_by_id($id)
    Function: Get a reference for a contig
    Returns : a Bio::Assembly::Contig object or undef
    Args    : [string] contig unique identifier (ID)

=cut

sub get_contig_by_id {
    my $self = shift;
    my $contigID = shift;

    return undef unless (exists $self->{'_contigs'}{$contigID});

    return $self->{'_contigs'}{$contigID};
}

=head2 get_singlet_by_id

    Title   : get_singlet_by_id
    Usage   : $assembly->get_singlet_by_id()
    Function: Get a reference for a singlet
    Returns : Bio::PrimarySeqI object or undef
    Args    : [string] a singlet ID

=cut

sub get_singlet_by_id {
    my $self = shift;

    my $singletID = shift;

    return undef unless (exists $self->{'_singlets'}{$singletID});

    return $self->{'_singlets'}{$singletID};
}

=head1 Modifier methods

=head2 add_contig

    Title   : add_contig
    Usage   : $assembly->add_contig($contig)
    Function: Add a contig to the assembly
    Returns : 1 on success
    Args    : a Bio::Assembly::Contig object
	      order (optional)

=cut

sub add_contig {
    my $self = shift;
    my $contig = shift;

    if( !ref $contig || ! $contig->isa('Bio::Assembly::Contig') ) {
	$self->throw("Unable to process non Bio::Assembly::Contig object [", ref($contig), "]");
    }
    my $contigID  = $contig->id();
    if( !defined $contigID ) {
	$contigID = 'Unknown_' . ($self->get_nof_contigs() + 1);
	$contig->id($contigID);
	$self->warn("Attributing ID $contigID to unidentified Bio::Assembly::Contig object.");
    }

    $self->warn("Replacing contig $contigID with a new contig object")
	if (exists $self->{'_contigs'}{$contigID});
    $self->{'_contigs'}{$contigID} = $contig;

    foreach my $seqID ($contig->get_seq_ids()) {
	if (exists $self->{'_seqs'}{$seqID}) {
	    $self->warn( "Sequence $seqID already assigned to contig ".
			 $self->{'_seqs'}{$seqID}->id().". Moving to contig $contigID")
		unless ($self->{'_seqs'}{$seqID} eq $contig);
	}
	$self->{'_seqs'}{$seqID} = $contig;
    }

    return 1;
}

=head2 add_singlet

    Title   : add_singlet
    Usage   : $assembly->add_singlet($seq)
    Function: Add a singlet to the assembly
    Returns : 1 on success, 0 otherwise
    Args    : a Bio::PrimarySeqI object
		  order (optional)

=cut

sub add_singlet {
    my $self = shift;
    my $singlet = shift;

    if( !ref $singlet || ! $singlet->isa('Bio::PrimarySeqI') ) {
	$self->warn("Unable to process non Bio::SeqI object [", ref($singlet), "]");
	return 0;
    }

    my $singletID = $singlet->id();
    $self->warn("Replacing singlet $singletID wih a new sequence object")
	if (exists $self->{'_contigs'}{$singletID});
    $self->{'_singlets'}{$singletID} = $singlet;

    return 1;
}

=head2 update_seq_list

    Title   : update_seq_list
    Usage   : $assembly->update_seq_list()
    Function: 

              Synchronizes the assembly registry for sequences in
              contigs and contig actual aligned sequences content. You
              probably want to run this after you remove/add a
              sequence from/to a contig in the assembly.

    Returns : nothing
    Args    : none 

=cut

sub update_seq_list {
    my $self = shift;

    $self->{'_seqs'} = {};
    foreach my $contig ($self->all_contigs) {
	foreach my $seqID ($contig->get_seq_ids) {
	    $self->{'_seqs'}{$seqID} = $contig;
	}
    }

    return 1;
}

=head2 remove_contigs

    Title   : remove_contigs
    Usage   : $assembly->remove_contigs(1..4)
    Function: Remove contig from assembly object
    Returns : an array of removed Bio::Assembly::Contig
              objects
    Args    : an array of contig IDs 

    See function get_contig_ids() above

=cut

#---------------------
sub remove_contigs {
#---------------------
    my ($self,@args) = @_;

    my @ret = ();
    foreach my $contigID (@args) {
	foreach my $seqID ($self->get_contig_by_id($contigID)->get_seq_ids()) {
	    delete $self->{'_seqs'}{$seqID};
	}
	push(@ret,$self->{'_contigs'}{$contigID});
	delete $self->{'_contigs'}{$contigID};
    }

    return @ret;
}

=head2 remove_singlets

    Title   : remove_singlets
    Usage   : $assembly->remove_singlets(@singlet_ids)
    Function: Remove singlet from assembly object
    Returns : the Bio::SeqI objects removed
    Args    : a list of singlet IDs

    See function get_singlet_ids() above

=cut

#---------------------
sub remove_singlets {
#---------------------
    my ($self,@args) = @_;

    my @ret = ();
    foreach my $singletID (@args) {
	push(@ret,$self->{'_singlets'}{$singletID});
	delete $self->{'_singlets'}{$singletID};
    }

    return @ret;
}

=head1 Contig and singlet selection methos

=head2 select_contigs

    Title   : select_contigs
    Usage   : $assembly->select_contigs(@list)
    Function: Select an array of contigs from the assembly
    Returns : an array of Bio::Assembly::Contig objects
    Args    : an array of contig ids

    See function get_contig_ids() above

=cut

#---------------------
sub select_contigs {
#---------------------
    my ($self,@args) = @_;

    my @contigs = ();
    foreach my $contig (@args) {
	unless (exists $self->{'_contigs'}{$contig}) {
	    $self->warn("$contig contig not found. Ignoring...");
	    next;
	}
	push(@contigs, $self->{'_contigs'}{$contig});
    }

    return @contigs;
}

=head2 select_singlets

    Title   : select_singlets
    Usage   : $assembly->select_singlets(@list)
    Function: Selects an array of singlets from the assembly
    Returns : an array of Bio::SeqI objects
    Args    : an array of singlet ids

    See function get_singlet_ids() above

=cut

#---------------------
sub select_singlets {
#---------------------
    my ($self,@args) = @_;

    my @singlets = ();
    foreach my $singlet (@args) {
	unless (exists $self->{'_singlets'}{$singlet}) {
	    $self->warn("$singlet singlet not found. Ignoring...");
	    next;
	}
	push(@singlets, $self->{'_singlets'}{$singlet});
    }

    return @singlets;
}

=head2 all_contigs

    Title   : all_contigs
    Usage   : my @contigs = $assembly->all_contigs
    Function: 

              Returns a list of all contigs in this assembly.  Contigs
              are both clusters and alignments of one or more reads,
              with an associated consensus sequence.

    Returns : array of Bio::Assembly::Contig (in lexical id order)
    Args    : none

=cut

#---------------------
sub all_contigs {
#---------------------
    my ($self) = @_;

    my @contigs = ();
    foreach my $contig (sort { $a cmp $b } keys %{ $self->{'_contigs'} }) {
	push(@contigs, $self->{'_contigs'}{$contig});
    }

    return @contigs;
}

=head2 all_singlets

    Title   : all_singlets
    Usage   : my @singlets = $assembly->all_singlets
    Function: 

              Returns a list of all singlets in this assembly.
	      Singlets are isolated reads, without non-vector
	      matches to any other read in the assembly.

    Returns : array of Bio::SeqI (in lexical order by id)
    Args    : none

=cut

#---------------------
sub all_singlets {
#---------------------
    my ($self) = @_;

    my @singlets = ();
    foreach my $singlet (sort { $a cmp $b } keys %{ $self->{'_singlets'} }) {
	push(@singlets, $self->{'_singlets'}{$singlet});
    }

    return @singlets;
}

# =head1 Internal Methods

1;
