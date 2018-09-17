#
# BioPerl module for Bio::DB::InMemoryCache
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::DB::InMemoryCache - Abstract interface for a sequence database

=head1 SYNOPSIS

  $cachedb = Bio::DB::InMemoryCache->new( -seqdb => $real_db,
                                          -number => 1000);
  #
  # get a database object somehow using a concrete class
  #

  $seq = $cachedb->get_Seq_by_id('ROA1_HUMAN');

  #
  # $seq is a Bio::Seq object
  #

=head1 DESCRIPTION

This is a memory cache system which saves the objects returned by
Bio::DB::RandomAccessI in memory to a hard limit of sequences.

=head1 CONTACT

Ewan Birney E<lt>birney@ebi.ac.ukE<gt>

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

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::DB::InMemoryCache;


use strict;

use Bio::Seq;

use base qw(Bio::Root::Root Bio::DB::SeqI);

sub new {
    my ($class,@args) = @_;

    my $self = Bio::Root::Root->new();
    bless $self,$class;

    my ($seqdb,$number,$agr) =
		$self->_rearrange([qw(SEQDB NUMBER AGRESSION)],@args);

    if( !defined $seqdb || !ref $seqdb ||
		  !$seqdb->isa('Bio::DB::RandomAccessI') ) {
       $self->throw("Must be a RandomAccess database not a [$seqdb]");
    }

    if( !defined $number ) {
        $number = 1000;
    }

    $self->seqdb($seqdb);
    $self->number($number);
    $self->agr($agr);

    # we consider acc as the primary id here
    $self->{'_cache_number_hash'} = {};
    $self->{'_cache_id_hash'}     = {};
    $self->{'_cache_acc_hash'}    = {};
    $self->{'_cache_number'}      = 1;

    return $self;
}

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id('ROA1_HUMAN')
 Function: Gets a Bio::Seq object by its name
 Returns : a Bio::Seq object
 Args    : the id (as a string) of a sequence
 Throws  : "id does not exist" exception

=cut

sub get_Seq_by_id{
   my ($self,$id) = @_;

   if( defined $self->{'_cache_id_hash'}->{$id} ) {
		my $acc = $self->{'_cache_id_hash'}->{$id};
		my $seq = $self->{'_cache_acc_hash'}->{$acc};
		$self->{'_cache_number_hash'}->{$seq->accession} =
		  $self->{'_cache_number'}++;
		return $seq;
   } else {
		return $self->_load_Seq('id',$id);
   }
}

=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   : $seq = $db->get_Seq_by_acc('X77802');
 Function: Gets a Bio::Seq object by accession number
 Returns : A Bio::Seq object
 Args    : accession number (as a string)
 Throws  : "acc does not exist" exception

=cut

sub get_Seq_by_acc{
   my ($self,$acc) = @_;

   #print STDERR "In cache get for $acc\n";
   if( defined $self->{'_cache_acc_hash'}->{$acc} ) {
       #print STDERR "Returning cached $acc\n";
       my $seq = $self->{'_cache_acc_hash'}->{$acc};
       $self->{'_cache_number_hash'}->{$seq->accession} =
			$self->{'_cache_number'}++;
       return $seq;
   } else {
     return $self->_load_Seq('acc',$acc);
   }
}



sub number {
    my ($self, $number) = @_;
    if ($number) {
        $self->{'number'} = $number;
    } else {
        return $self->{'number'};
    }
}

sub seqdb {
    my ($self, $seqdb) = @_;
    if ($seqdb) {
        $self->{'seqdb'} = $seqdb;
    } else {
        return $self->{'seqdb'};
    }
}

sub agr {
    my ($self, $agr) = @_;
    if ($agr) {
        $self->{'agr'} = $agr;
    } else {
        return $self->{'agr'};
    }
}


sub _load_Seq {
  my ($self,$type,$id) = @_;

  my $seq;

  if( $type eq 'id') {
    $seq = $self->seqdb->get_Seq_by_id($id);
  }elsif ( $type eq 'acc' ) {
    $seq = $self->seqdb->get_Seq_by_acc($id);
  } else {
    $self->throw("Bad internal error. Don't understand $type");
  }
  if( ! $seq ) {
      # warding off bug #1628
      $self->debug("could not find seq $id in seqdb\n");
      return;
  }

  if( $self->agr() ) {
      #print STDERR "Pulling out into memory\n";
      my $newseq = Bio::Seq->new( -display_id => $seq->display_id,
				  -accession_number  => $seq->accession,
				  -seq        => $seq->seq,
				  -desc       => $seq->desc,
				  );
      if( $self->agr() == 1 ) {
	  foreach my $sf ( $seq->top_SeqFeatures() ) {
	      $newseq->add_SeqFeature($sf);
	  }

	  $newseq->annotation($seq->annotation);
      }
      $seq = $newseq;
  }

  if( $self->_number_free < 1 ) {
    # remove the latest thing from the hash
    my @accs = sort { $self->{'_cache_number_hash'}->{$a} <=>
								$self->{'_cache_number_hash'}->{$b} }
                          keys %{$self->{'_cache_number_hash'}};

    my $acc = shift @accs;
    # remove this guy
    my $seq = $self->{'_cache_acc_hash'}->{$acc};

    delete $self->{'_cache_number_hash'}->{$acc};
    delete $self->{'_cache_id_hash'}->{$seq->id};
    delete $self->{'_cache_acc_hash'}->{$acc};
  }

  # up the number, register this sequence into the hash.
  $self->{'_cache_id_hash'}->{$seq->id} = $seq->accession;
  $self->{'_cache_acc_hash'}->{$seq->accession} = $seq;
  $self->{'_cache_number_hash'}->{$seq->accession} = $self->{'_cache_number'}++;

  return $seq;
}


sub _number_free {
  my $self = shift;

  return $self->number - scalar(keys %{$self->{'_cache_number_hash'}});
}

=head2 get_Seq_by_version

 Title   : get_Seq_by_version
 Usage   : $seq = $db->get_Seq_by_version('X77802.1');
 Function: Gets a Bio::Seq object by sequence version
 Returns : A Bio::Seq object
 Args    : accession.version (as a string)
 Throws  : "acc.version does not exist" exception

=cut

sub get_Seq_by_version{
   my ($self,@args) = @_;
   $self->throw("Not implemented it");
}

## End of Package

1;
