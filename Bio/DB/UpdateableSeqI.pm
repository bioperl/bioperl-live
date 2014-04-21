#
#
# BioPerl module for Bio::DB::UpdateableSeqI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
#
# _history
# June 18, 2000 - module begun
#
# POD Doc - main docs before code

=head1 NAME

Bio::DB::UpdateableSeqI - An interface for writing to a database of sequences.

=head1 SYNOPSIS 

    # get a Bio::DB::UpdateableSeqI somehow
    eval {
	my ( @updatedseqs, @newseqs, @deadseqs);
	my $seq = $db->get_Seq_by_id('ROA1_HUMAN');
	$seq->desc('a new description');

	push @updatedseqs, $seq;

	$db->write_seq(\@updatedseqs, \@newseqs, \@deadseqs);
    };
    if( $@ ) {
	print STDERR "an error when trying to write seq : $@\n";
    }

=head1 DESCRIPTION

This module seeks to provide a simple method for pushing sequence changes 
back to a Sequence Database - which can be an SQL compliant database, a file 
based database, AceDB, etc.

=head1 AUTHOR

Jason Stajich E<lt>jason@bioperl.orgE<gt>

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues           

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

#Lets start some code

package Bio::DB::UpdateableSeqI;

use strict;



use base qw(Bio::DB::SeqI);

=head2 write_seq

  Title   : write_seq
  Usage   : write_seq(\@updatedseqs, \@addedseqs, \@deadseqs)
  Function: updates sequences in first array,
            adds sequences in the second array,
            and removes sequences in the third array.
  Example :
  Returns :
  Args    : arrays of sequence objects that must be obtained from
            Bio::DB::UpdateableSeqI.

=cut

sub write_seq {
    my ($self) = @_;
    
    $self->throw("Abstract database call of write_seq. Your database has not implemented this method!");

}

=head2 _add_seq

 Title   : _add_seq
 Usage   : _add_seq($seq)
 Function: Adds a new sequence
 Example : 
 Returns : will throw an exception if
           sequences accession number already exists
 Args    : a new seq object - should have an accession number

=cut

sub _add_seq {
    my ($self ) = @_;
    
    $self->throw("Abstract database call of _add_seq. Your database has not implemented this method!");

}

=head2 _remove_seq

 Title   : _remove_seq
 Usage   : _remove_seq($seq)
 Function: Removes an existing sequence
 Example : 
 Returns : will throw an exception if
           sequence does not exists for the primary_id
 Args    : a seq object that was retrieved from Bio::DB::UpdateableSeqI

=cut

sub _remove_seq {
    my ($self) = @_;
    
    $self->throw("Abstract database call of _remove_seq. Your database has not implemented this method!");

}

=head2 _update_seq

 Title   : _update_seq
 Usage   : _update_seq($seq)
 Function: Updates a sequence
 Example : 
 Returns : will throw an exception if
           sequence is out of sync from expected val.
 Args    : a seq object that was retrieved from Bio::DB::UpdateableSeqI

=cut

sub _update_seq {
    my ($self) = @_;
    
    $self->throw("Abstract database call of _update_seq. Your database has not implemented this method!");

}


=head1 Methods inherieted from Bio::DB::RandomAccessI

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id('ROA1_HUMAN')
 Function: Gets a Bio::Seq object by its name
 Returns : a Bio::Seq object
 Args    : the id (as a string) of a sequence
 Throws  : "id does not exist" exception


=cut

=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   : $seq = $db->get_Seq_by_acc('X77802');
 Function: Gets a Bio::Seq object by accession number
 Returns : A Bio::Seq object
 Args    : accession number (as a string)
 Throws  : "acc does not exist" exception


=cut

=head1 Methods inheirited from Bio::DB::SeqI

=head2 get_PrimarySeq_stream

 Title   : get_PrimarySeq_stream
 Usage   : $stream = get_PrimarySeq_stream
 Function: Makes a Bio::DB::SeqStreamI compliant object
           which provides a single method, next_primary_seq
 Returns : Bio::DB::SeqStreamI
 Args    : none


=cut

=head2 get_all_primary_ids

 Title   : get_all_ids
 Usage   : @ids = $seqdb->get_all_primary_ids()
 Function: gives an array of all the primary_ids of the 
           sequence objects in the database. These
           maybe ids (display style) or accession numbers
           or something else completely different - they
           *are not* meaningful outside of this database
           implementation.
 Example :
 Returns : an array of strings
 Args    : none


=cut

=head2 get_Seq_by_primary_id

 Title   : get_Seq_by_primary_id
 Usage   : $seq = $db->get_Seq_by_primary_id($primary_id_string);
 Function: Gets a Bio::Seq object by the primary id. The primary
           id in these cases has to come from $db->get_all_primary_ids.
           There is no other way to get (or guess) the primary_ids
           in a database.

           The other possibility is to get Bio::PrimarySeqI objects
           via the get_PrimarySeq_stream and the primary_id field
           on these objects are specified as the ids to use here.
 Returns : A Bio::Seq object
 Args    : accession number (as a string)
 Throws  : "acc does not exist" exception


=cut

1;



