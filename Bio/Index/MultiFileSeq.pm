
#
# BioPerl module for Bio::DB::MultiFileSeq
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::MultiFileSeq - Base class for MultiFileSeq's #'

=head1 SYNOPSIS

  # not really appropiate - just providing a common base
  # class for multiple sequence files built using the Bio::Index::Abstract
  # system

=head1 DESCRIPTION


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Index::MultiFileSeq;
use vars qw(@ISA);
use strict;

use Bio::SeqIO;
use Bio::SeqIO::MultiFile;
use Bio::Index::Abstract;
use Bio::DB::SeqI;


@ISA = qw(Bio::Index::Abstract Bio::DB::SeqI);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 _file_format

 Title   : _file_format
 Usage   : $self->_file_format
 Function: Derived classes should override this
           method (it throws an exception here)
           to give the file format of the files used
 Example :
 Returns : 
 Args    :


=cut

sub _file_format{
   my ($self,@args) = @_;

   $self->throw("Class must provide a file format method correctly");
}

=head2 fetch

  Title   : fetch
  Usage   : $index->fetch( $id )
  Function: Returns a Bio::Seq object from the index
  Example : $seq = $index->fetch( 'dJ67B12' )
  Returns : Bio::Seq object
  Args    : ID

=cut

sub fetch {
    my( $self, $id ) = @_;
    my ($desc,$acc,$out);
    my $db = $self->db();
    my $seq;

    if (my $rec = $db->{ $id }) {
	my ($file, $begin, $end) = $self->unpack_record( $rec );
        
        # Get the (possibly cached) filehandle
        my $fh = $self->_file_handle( $file );

        # move to start
        seek($fh, $begin, 0);

	# make a new SeqIO object

	my $seqio = Bio::SeqIO->new( '-format' => $self->_file_format,
				  -fh => $fh);
	
	$seq = $seqio->next_seq();
	$seqio->close();
    }

    # we essentially assumme that the primary_id for the database
    # is the display_id
    $seq->primary_id($seq->display_id());

    return $seq;
}

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id()
 Function: retrieves a sequence object, identically to
           ->fetch, but here behaving as a Bio::DB::BioSeqI
 Returns : new Bio::Seq object
 Args    : string represents the id


=cut

sub get_Seq_by_id{
   my ($self,$id) = @_;

   return $self->fetch($id);
}

=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   : $seq = $db->get_Seq_by_acc()
 Function: retrieves a sequence object, identically to
           ->fetch, but here behaving as a Bio::DB::BioSeqI
 Returns : new Bio::Seq object
 Args    : string represents the accession number


=cut

sub get_Seq_by_acc {
   my ($self,$id) = @_;

   return $self->fetch($id);
}

=head2 get_PrimarySeq_stream

 Title   : get_PrimarySeq_stream
 Usage   : $stream = get_PrimarySeq_stream
 Function: Makes a Bio::DB::SeqStreamI compliant object
           which provides a single method, next_primary_seq
 Returns : Bio::DB::SeqStreamI
 Args    : none


=cut

sub get_PrimarySeq_stream{
    my $self = shift;
    my $num  = $self->_file_count() || 0;
    my @file;
    
    for (my $i = 0; $i < $num; $i++) {
        my( $file, $stored_size ) = $self->unpack_record( $self->db->{"__FILE_$i"} );
	push(@file,$file);
    }
   
    my $out = Bio::SeqIO::MultiFile->new( '-format' => $self->_file_format , -files => \@file);
    return $out;
}

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

sub get_all_primary_ids{
   my ($self,@args) = @_;
    my $db = $self->db;
   
   # the problem is here that we have indexed things both on
   # accession number and byte position. 

   # We could take two options
   # here - loop over the database, returnin g only one copy of each
   # id that points to the same byte position, or we rely on semantics
   # of accession numbers.

   # someone is going to index a database with no accession numbers.
   # doh!. We have to uniquify the index...

    my( %bytepos );
   while (my($id, $rec) = each %$db) {
       my ($file, $begin, $end) = $self->unpack_record( $rec );
       $bytepos{$begin} = $id;
   }

   return values %bytepos;
}


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
 Args    : primary id (as a string)
 Throws  : "acc does not exist" exception


=cut

sub get_Seq_by_primary_id {
   my ($self,$id) = @_;
   return $self->fetch($id);
}

1;
