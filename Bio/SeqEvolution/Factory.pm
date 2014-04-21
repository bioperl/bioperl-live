#
# BioPerl module for Bio::SeqEvolution::Factory
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki at bioperl dot org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqEvolution::Factory - Factory object to instantiate sequence evolving classes

=head1 SYNOPSIS

    # not an instantiable class

=head1 DESCRIPTION

This is the factory class that can be used to call for a specific
model to mutate a sequence.

Bio::SeqEvolution::DNAPoint is the default for nucleotide sequences
and the only implementation at this point.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

  Heikki Lehvaslaiho E<lt>heikki at bioperl dot orgE<gt>

=head1 CONTRIBUTORS

Additional contributor's names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqEvolution::Factory;
use strict;
use Bio::Root::Root;
use Bio::SeqEvolution::EvolutionI;
use base qw(Bio::Root::Root Bio::SeqEvolution::EvolutionI);

=head2 new

  Title   : new
  Usage   : my $obj = Bio::SeqEvolution::Factory->new();
  Function: Builds a new Bio:SeqEvolution::EvolutionI object
  Returns : Bio:SeqEvolution::EvolutionI object
  Args    : -type           => class name

See L<Bio:SeqEvolution::EvolutionI>

=cut

sub new {
    my($caller,@args) = @_;
    my $class = ref($caller) || $caller;

    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param; # lowercase keys

    if ( $class eq 'Bio::SeqEvolution::Factory') {
        #my %param = @args;
        #@param{ map { lc $_ } keys %param } = values %param; # lowercase keys

        if (exists $param{'-type'}) {
#            $self->type($param{'-type'});
        } else {
            $param{'-type'} = 'Bio::SeqEvolution::DNAPoint';
            #$self->type('Bio::SeqEvolution::DNAPoint'} unless $seq->alphabet == 'protein'
        }
        my $type = $param{'-type'};
        return unless( $class->_load_format_module($param{'-type'}) );
        return $type->new(%param);
    } else {
        my ($self) = $class->SUPER::new(%param);
        $self->_initialize(%param);
        return $self;
    }
}

sub _initialize {
    my($self, @args) = @_;

    $self->SUPER::_initialize(@args);
    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param; # lowercase keys

    exists $param{'-seq'} && $self->seq($param{'-seq'});
    exists $param{'-set_mutated_seq'} && $self->set_mutated_seq($param{'-set_mutated_seq'});
    exists $param{'-identity'} && $self->identity($param{'-identity'});
    exists $param{'-pam'} && $self->pam($param{'-pam'});
    exists $param{'-mutation_count'} && $self->mutation_count($param{'-mutation_count'});

}


=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL SeqIO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_format_module {
	my ($self, $format) = @_;
	my $module =  $format;
	my $ok;

	eval {
		$ok = $self->_load_module($module);
	};
	if ( $@ ) {
		print STDERR <<END;
$self: $format cannot be found
Exception $@
END
		;
	}
	return $ok;
}


=head2 type

 Title   : type
 Usage   : $obj->type($newval)
 Function: Set used evolution model. It is set by giving a
           valid Bio::SeqEvolution::* class name 
 Returns : value of type
 Args    : newvalue (optional)

Defaults to Bio::SeqEvolution::DNAPoint.

=cut

sub type{
   my $self = shift;
   if (@_) {
       $self->{'_type'} = shift @_;
       $self->_load_module($self->{'_type'});
   }
   return $self->{'_type'} || 'Bio::SeqEvolution::DNAPoint';
}

=head1 mutation counters

The next three methods set a value to limit the number of mutations
introduced the the input sequence.

=cut

=head2 identity

 Title   : identity
 Usage   : $obj->identity($newval)
 Function: Set the desired identity between original and mutated sequence
 Returns : value of identity
 Args    : newvalue (optional)

=cut

sub identity{
   my $self = shift;
   $self->{'_identity'} = shift @_ if @_;
   return $self->{'_identity'};
}


=head2 pam

 Title   : pam
 Usage   : $obj->pam($newval)
 Function: Set the wanted Percentage of Accepted Mutations, PAM
 Returns : value of PAM
 Args    : newvalue (optional)

When you are measuring sequence divergence, PAM needs to be
estimated. When you are generating sequences, PAM is simply the count
of mutations introduced to the reference sequence normalised to the
original sequence length.

=cut

sub pam{
   my $self = shift;
   $self->{'_pam'} = shift @_ if @_;
   return $self->{'_pam'};
}

=head2 mutation_count

 Title   : mutation_count
 Usage   : $obj->mutation_count($newval)
 Function: Set the number of wanted mutations to the sequence 
 Returns : value of mutation_count
 Args    : newvalue (optional)

=cut

sub mutation_count{
   my $self = shift;
   $self->{'_mutation_count'} = shift @_ if @_;
   return $self->{'_mutation_count'};
}



=head2 seq

 Title   : seq
 Usage   : $obj->seq($newval)
 Function: Set the sequence object for the original sequence
 Returns : The sequence object
 Args    : newvalue (optional)

Setting this will reset mutation and generated mutation counters.

=cut

sub seq {
   my $self = shift;
   if (@_) {
       $self->{'_seq'} = shift @_ ;
       return $self->{'_seq'};
       $self->reset_mutation_counter;
       $self->reset_sequence_counter;
   }
   return $self->{'_seq'};
}

=head2 seq_type

 Title   : seq_type
 Usage   : $obj->seq_type($newval)
 Function: Set the returned seq_type to one needed
 Returns : value of seq_type
 Args    : newvalue (optional)

Defaults to Bio::PrimarySeq.

=cut

sub seq_type{
   my $self = shift;
   if (@_) {
       $self->{'_seq_type'} = shift @_;
       $self->_load_module($self->{'_seq_type'});
   }
   return $self->{'_seq_type'} || 'Bio::PrimarySeq';
}


=head2 get_mutation_counter

 Title   : get_mutation_counter
 Usage   : $obj->get_mutation_counter()
 Function: Get the count of sequences created
 Returns : value of counter
 Args    : -

=cut

sub get_mutation_counter{
   return shift->{'_mutation_counter'};
}


=head2 reset_mutation_counter

 Title   : reset_mutation_counter
 Usage   : $obj->reset_mutation_counter()
 Function: Resert the counter of mutations
 Returns : value of counter
 Args    : -

=cut

sub reset_mutation_counter{
   shift->{'_mutation_counter'} = 0;
   return 1;
}


=head2 get_sequence_counter

 Title   : get_sequence_counter
 Usage   : $obj->get_sequence_counter()
 Function: Get the count of sequences created
 Returns : value of counter
 Args    : -

=cut

sub get_sequence_counter{
   return shift->{'_sequence_counter'};
}

=head2 reset_sequence_counter

 Title   : reset_sequence_counter
 Usage   : $obj->reset_sequence_counter()
 Function: Resert the counter of sequences created
 Returns : value of counter
 Args    : -

This is called when ever mutated sequences are reassigned new values
using methods seq() and mutated_seq().  As a side affect, this method
also recreates the intermal alignment that is used to calculate the
sequence identity.

=cut

sub reset_sequence_counter{
   my $self = shift;
   $self->{'_sequence_counter'} = 0;
   $self->_init_alignment;
   return 1;
}



=head2 each_seq

 Title   : each_seq
 Usage   : $obj->each_seq($int)
 Function:
 Returns : an array of sequences mutated from the reference sequence
           according to evolutionary parameters given
 Args    : -

=cut

sub each_seq{
   my $self = shift;
   my $number = shift;

   $self->throw("[$number] ". ' should be a positive integer')
       unless $number =~ /^[+\d]+$/;

   my @array;
   for (my $count=1; $count<$number; $count++) {
       push @array, $self->next_seq();

   }
   return @array;
}



=head2 each_mutation

  Title   : each_mutation
  Usage   : $obj->each_mutation
  Function: return the mutations leading to the last generated 
            sequence in objects 
  Returns : an array of Bio::Variation::DNAMutation objects
  Args    : optional argument to return an array of  stringified names

=cut

sub each_mutation {
    my $self = shift;
    my $string = shift;

    return @{$self->{'_mutations'}} if $string;

    return map {
        /(\d+)(\w*)>(\w*)/;
#        print;
        my $dnamut = Bio::Variation::DNAMutation->new
            ('-start'         => $1,
             '-end'           => $1,
             '-length'        => 1,
             '-isMutation'    => 1
            );
        $dnamut->allele_ori( Bio::Variation::Allele->new(-seq => $2,
                                                         -alphabet => 'dna') );
        $dnamut->add_Allele( Bio::Variation::Allele->new(-seq => $3,
                                                         -alphabet => 'dna') );
        $dnamut;
    } @{$self->{'_mutations'}}
}


sub get_alignment_identity  {
    my $self = shift;
    return $self->{'_align'}->overall_percentage_identity;
}


sub get_alignmet {
   my $self = shift;
   return $self->{'_align'}->remove_gaps('-', 'all-gaps');
}


=head1 Internal methods

=cut


=head2 _increase_mutation_counter

 Title   : _increase_mutation_counter
 Usage   : $obj->_increase_mutation_counter()
 Function: Internal method to increase the counter of mutations performed
 Returns : value of counter
 Args    : -

=cut

sub _increase_mutation_counter{
   return shift->{'_mutation_counter'}++;
}



=head2 _increase_sequence_counter

 Title   : _increase_sequence_counter
 Usage   : $obj->_increase_sequence_counter()
 Function: Internal method to increase the counter of sequences created
 Returns : value of counter
 Args    : -

=cut

sub _increase_sequence_counter{
   return shift->{'_sequence_counter'}++;
}


1;

