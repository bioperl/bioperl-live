#
# BioPerl module for Bio::Assembly::Scaffold
#
# Copyright by Robson F. de Souza
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::Scaffold - Perl module to hold and manipulate sequence assembly
data.

=head1 SYNOPSIS
# 
    # Module loading
    use Bio::Assembly::IO;

    # Assembly loading methods
    my $aio = Bio::Assembly::IO->new(-file=>"test.ace.1", -format=>'phrap');
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

rfsouza@citri.iq.usp.br

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Assembly::Scaffold;

use strict;
use Bio::Annotation::Collection;

use base qw(Bio::Root::Root Bio::Assembly::ScaffoldI);

=head2 new ()

    Title   : new
    Usage   : $scaffold = new ( -id       => "assembly 1",
                                -source   => 'program_name',
                                -contigs  => \@contigs,
                                -singlets => \@singlets );
    Function: creates a new scaffold object
    Returns : Bio::Assembly::Scaffold
    Args    : -id       : [string] scaffold name
              -source   : [string] sequence assembly program
              -contigs  : reference to array of Bio::Assembly::Contig objects
              -singlets : reference to array of Bio::Assembly::Singlet objects


=cut

sub new {
  my($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($id, $src, $contigs, $singlets) = $self->_rearrange(
      [qw(ID SOURCE CONTIGS SINGLETS)], @args);

  # Scaffold defaults
  $self->{'_id'} = 'NoName';
  $self->{'_source'} = 'Unknown';
  $self->{'_contigs'} = {};
  $self->{'_singlets'} = {};
  $self->{'_seqs'} = {};
  $self->{'_annotation'} = Bio::Annotation::Collection->new();

  # Import manual info
  $self->{'_id'} = $id if (defined $id);
  $self->{'_source'} = $src if (defined $src);
  
  # Add contigs and singlets to scaffold
  if (defined $contigs && ref($contigs = 'ARRAY')) {
    for my $contig (@{$contigs}) {
      if( ! ref $contig || ! $contig->isa('Bio::Assembly::Contig') ) {
        $self->throw("Bio::Assembly::Scaffold::new is unable to process non".
          "Bio::Assembly::Contig object [", ref($contig), "]");
      }
      $self->add_contig($contig);
    }
  }
  if (defined $singlets && ref($singlets = 'ARRAY')) {
    for my $singlet (@{$singlets}) {
      if( ! ref $singlet || ! $singlet->isa('Bio::Assembly::Singlet') ) {
        $self->throw("Bio::Assembly::Scaffold::new is unable to process non".
          "Bio::Assembly::Singlet object [", ref($singlet), "]");
      }
      $self->add_singlet($singlet);
    }
  }
  
  return $self;
}

=head1 Accessing general assembly data

=cut

=head2 id

    Title   : id
    Usage   : $assembly->id()
    Function: Get/Set assembly ID
    Returns : string or undef
    Args    : string

=cut

sub id {
    my ($self, $id) = @_;
    return $self->{'_id'} = $id if (defined $id);
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
    my ($self, $ref) = shift;
    $self->{'_annotation'} = $ref if (defined $ref);
    return $self->{'_annotation'};
}

=head2 get_nof_contigs

    Title   : get_nof_contigs
    Usage   : $assembly->get_nof_contigs()
    Function: Get the number of contigs included in the scaffold
    Returns : integer
    Args    : none

=cut

sub get_nof_contigs {
    my $self = shift;
    return scalar( $self->get_contig_ids() );
}

=head2 get_nof_contig_seqs

    Title   : get_nof_contig_seqs
    Usage   : $assembly->get_nof_contig_seqs()
    Function: Get the number of sequences included in contigs of the 
              scaffold (no consensus sequences or singlets)
    Returns : integer
    Args    : none

=cut

sub get_nof_contig_seqs {
    my $self = shift;

    my $nof_seqs = 0;
    foreach my $contig ($self->all_contigs) {
        $nof_seqs += scalar( $contig->get_seq_ids() );
    }

    return $nof_seqs;
}
# function alias for backward compatibility
*get_nof_sequences_in_contigs = \&get_nof_contig_seqs;


=head2 get_nof_singlets (get_nof_singlet_seqs)

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
*get_nof_singlet_seqs = \&get_nof_singlets;

=head2 get_all_seq_ids

    Title   : get_all_seq_ids
    Usage   : $assembly->get_all_seq_ids()
    Function: Get the ID of all sequences making up the scaffold
              (sequences from contigs and singlets, not consensus).
    Returns : array of strings
    Args    : none

=cut

sub get_all_seq_ids {
    my $self = shift;
    return keys %{ $self->{'_seqs'} };
}

=head2 get_nof_seqs

    Title   : get_nof_seqs
    Usage   : $assembly->get_nof_seqs()
    Function: Get total number of sequences making up the scaffold
              (sequences from contigs and singlets, not consensus).
    Returns : integer
    Args    : none

=cut

sub get_nof_seqs {
    my $self = shift;
    return scalar $self->get_all_seq_ids;
}

=head2 get_contig_seq_ids

    Title   : get_contig_seq_ids
    Usage   : $assembly->get_contig_seq_ids()
    Function: Get the ID of all sequences in contigs
    Returns : array of strings
    Args    : none

=cut

sub get_contig_seq_ids {
    my $self = shift;
    my @ids;
    for my $contig ( $self->all_contigs ) {
        push @ids, $contig->get_seq_ids;
    }
    return @ids;
}
# function alias for backward compatibility
*get_seq_ids = \&get_contig_seq_ids; 

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

    return wantarray
        ? sort keys %{$self->{'_contigs'}}
        : scalar keys %{$self->{'_contigs'}};
}

=head2 get_singlet_ids (get_singlet_seq_ids)

    Title   : get_singlet_ids
    Usage   : $assembly->get_singlet_ids()
    Function: Access list of singlet IDs from assembly
    Returns : array of strings if there are any singlets
              otherwise an empty array
    Args    : none

=cut

sub get_singlet_ids {
    my $self = shift;

    return wantarray
        ? sort keys %{$self->{'_singlets'}}
        : scalar keys %{$self->{'_singlets'}};
}
*get_singlet_seq_ids = \&get_singlet_ids;

=head2 get_seq_by_id

    Title   : get_seq_by_id
    Usage   : $assembly->get_seq_by_id($id)
    Function: Get a reference for an sequence making up the scaffold 
              (from a contig or singlet, not consensus)
    Returns : a Bio::LocatableSeq object
              undef if sequence $id is not found in the scaffold
    Args    : [string] sequence identifier (id)

=cut

sub get_seq_by_id {
    my $self = shift;
    my $seqID = shift;

    return unless (exists $self->{'_seqs'}{$seqID});

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

    return unless (exists $self->{'_contigs'}{$contigID});

    return $self->{'_contigs'}{$contigID};
}

=head2 get_singlet_by_id

    Title   : get_singlet_by_id
    Usage   : $assembly->get_singlet_by_id()
    Function: Get a reference for a singlet
    Returns : Bio::Assembly::Singlet object or undef
    Args    : [string] a singlet ID

=cut

sub get_singlet_by_id {
    my $self = shift;

    my $singletID = shift;

    return unless (exists $self->{'_singlets'}{$singletID});

    return $self->{'_singlets'}{$singletID};
}

=head1 Modifier methods

=cut

=head2 add_contig

    Title   : add_contig
    Usage   : $assembly->add_contig($contig)
    Function: Add a contig to the assembly
    Returns : 1 on success
    Args    : a Bio::Assembly::Contig object
          order (optional)

=cut

sub add_contig {
    my ($self, $contig) = @_;

    # Input check
    if( !ref $contig || ! $contig->isa('Bio::Assembly::Contig') ) {
        $self->throw("Bio::Assembly::Scaffold::add_contig is unable to process".
            " non Bio::Assembly::Contig object [", ref($contig), "]");
    }
    
    # Create and attribute contig ID
    my $contigID  = $contig->id();
    if( !defined $contigID ) {
        $contigID = 'Unknown_' . ($self->get_nof_contigs() + 1);
        $contig->id($contigID);
        $self->warn("Attributing ID $contigID to unnamed Bio::Assembly::Contig".
            " object.");
    }

    # Adding contig to scaffold
    $self->warn("Replacing contig $contigID with a new contig object")
        if (exists $self->{'_contigs'}{$contigID});
    $self->{'_contigs'}{$contigID} = $contig;
    $contig->assembly($self); # weak circular reference

    # Put contig sequences in the list of sequences belonging to the scaffold
    foreach my $seqID ($contig->get_seq_ids()) {
        if (exists $self->{'_seqs'}{$seqID} &&
            not($self->{'_seqs'}{$seqID} eq $contig) ) {
            $self->warn( "Sequence $seqID already assigned to object ".
                $self->{'_seqs'}{$seqID}->id().". Moving to contig $contigID");
        }
        $self->{'_seqs'}{$seqID} = $contig;
    }
    
    return 1;
}

=head2 add_singlet

    Title   : add_singlet
    Usage   : $assembly->add_singlet($seq)
    Function: Add a singlet to the assembly
    Returns : 1 on success
    Args    : a Bio::Assembly::Singlet object
              order (optional)

=cut

sub add_singlet {
    my ($self, $singlet) = @_;

    # Input check
    if ( !ref $singlet || ! $singlet->isa('Bio::Assembly::Singlet') ) {
        $self->throw("Bio::Assembly::Scaffold::add_singlet is unable to process".
            " non Bio::Assembly::Singlet object [", ref($singlet), "]");
    }
    
    # Create and attribute singlet ID
    my $singletID = $singlet->id();
    if( !defined $singletID ) {
        $singletID = 'Unknown_' . ($self->get_nof_singlets() + 1);
        $singlet->id($singletID);
        $self->warn("Attributing ID $singletID to unnamed Bio::Assembly::".
            "Singlet object.");
    }
    
    # Adding singlet to scaffold
    $self->warn("Replacing singlet $singletID with a new singlet object")
        if (exists $self->{'_singlets'}{$singletID});
    $self->{'_singlets'}{$singletID} = $singlet;
    $singlet->assembly($self); # weak circular reference

    # Put singlet sequence in the list of sequences belonging to the scaffold
    my $seqID = $singlet->seqref->id();
    if (exists $self->{'_seqs'}{$seqID} &&
        not($self->{'_seqs'}{$seqID} eq $singlet) ) {
        $self->warn( "Sequence $seqID already assigned to object ".
            $self->{'_seqs'}{$seqID}->id().". Moving to singlet $singletID");
    }
    $self->{'_seqs'}{$seqID} = $singlet;

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

    Returns : 1 for success
    Args    : none 

=cut

sub update_seq_list {
    my $self = shift;
    
    $self->{'_seqs'} = {};

    # Put sequences in contigs in list of sequences belonging to the scaffold
    foreach my $contig ($self->all_contigs) {
        my $contigID = $contig->id();
        foreach my $seqID ($contig->get_seq_ids) {
            if (exists $self->{'_seqs'}{$seqID} &&
                not($self->{'_seqs'}{$seqID} eq $contig) ) {
                $self->warn( "Sequence $seqID already assigned to object ".
                    $self->{'_seqs'}{$seqID}->id().". Moving to contig $contigID");
            }
            $self->{'_seqs'}{$seqID} = $contig;
        }
    }
    
    # Put singlet sequences in the list of sequences belonging to the scaffold
    foreach my $singlet ($self->all_singlets) {
        my $singletID = $singlet->id();
        my $seqID     = $singlet->seqref->id();
        if (exists $self->{'_seqs'}{$seqID} &&
            not($self->{'_seqs'}{$seqID} eq $singlet) ) {
            $self->warn( "Sequence $seqID already assigned to object ".
                $self->{'_seqs'}{$seqID}->id().". Moving to singlet $singletID");
        }
        $self->{'_seqs'}{$seqID} = $singlet;
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

sub remove_contigs {
    my ($self, @args) = @_;

    my @ret = ();
    foreach my $contigID (@args) {
        foreach my $seqID ($self->get_contig_by_id($contigID)->get_seq_ids()) {
            delete $self->{'_seqs'}{$seqID};
        }
        push(@ret, $self->{'_contigs'}{$contigID});
        delete $self->{'_contigs'}{$contigID};
    }

    return @ret;
}

=head2 remove_singlets

    Title   : remove_singlets
    Usage   : $assembly->remove_singlets(@singlet_ids)
    Function: Remove singlet from assembly object
    Returns : the Bio::Assembly::Singlet objects removed
    Args    : a list of singlet IDs

    See function get_singlet_ids() above

=cut

sub remove_singlets {
    my ($self,@args) = @_;

    my @ret = ();
    foreach my $singletID (@args) {
        push(@ret,$self->{'_singlets'}{$singletID});
        delete $self->{'_singlets'}{$singletID};
    }

    return @ret;
}

=head2 remove_features_collection

    Title   : remove_features_collection
    Usage   : $assembly->remove_features_collection()
    Function: Removes the collection of features associated to every
              contig and singlet of the scaffold. This can be useful
              to save some memory (when contig and singlet features are
              not needed).
    Returns   : none
    Argument  : none

=cut

sub remove_features_collection {
    my ($self) = @_;
    for my $obj ( $self->all_contigs, $self->all_singlets ) {
        $obj->remove_features_collection;
    }
    return;
}


=head1 Contig and singlet selection methods

=cut

=head2 select_contigs

    Title   : select_contigs
    Usage   : $assembly->select_contigs(@list)
    Function: Select an array of contigs from the assembly
    Returns : an array of Bio::Assembly::Contig objects
    Args    : an array of contig ids

    See function get_contig_ids() above

=cut

sub select_contigs {
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
    Returns : an array of Bio::Assembly::Singlet objects
    Args    : an array of singlet ids

    See function get_singlet_ids() above

=cut

sub select_singlets {
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

sub all_contigs {
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

    Returns : array of Bio::Assembly::Singlet objects (in lexical order by id)
    Args    : none

=cut

sub all_singlets {
    my ($self) = @_;

    my @singlets = ();
    foreach my $singlet (sort { $a cmp $b } keys %{ $self->{'_singlets'} }) {
    push(@singlets, $self->{'_singlets'}{$singlet});
    }

    return @singlets;
}


1;
