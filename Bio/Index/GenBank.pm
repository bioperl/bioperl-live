#
# $Id$
#
# BioPerl module for Bio::Index::Abstract
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::GenBank - Interface for indexing one or more GenBank
files (i.e. flat file GenBank format).

=head1 SYNOPSIS

    # Complete code for making an index for several
    # GenBank files
    use Bio::Index::GenBank;
    use strict;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::GenBank->new(-filename => $Index_File_Name, 
				       -write_flag => 'WRITE');
    $inx->make_index(@ARGV);

    # Print out several sequences present in the index
    # in gcg format
    use Bio::Index::GenBank;
    use Bio::SeqIO;
    use strict;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::GenBank->new(-filename => $Index_File_Name);
    my $seqio = new Bio::SeqIO(-format => 'gcg');
    foreach my $id (@ARGV) {
        my $seq = $inx->fetch($id); # Returns Bio::Seq object
        $seqio->write_seq($seq);
    }

    # alternatively
    my ($id, $acc);
    my $seq1 = $inx->get_Seq_by_id($locus);
    my $seq2 = $inx->get_Seq_by_acc($acc);

=head1 DESCRIPTION

By default the index that's created uses the LOCUS, ACCESSION, and
VERSION identifiers as keys. Inherits functions for managing dbm 
files from Bio::Index::Abstract.pm, and provides the basic 
functionality for indexing GenBank files, and retrieving the 
sequence from them. For best results 'use strict'.

You can also set or customize the unique key used to retrieve by 
writing your own function and calling the id_parser() method.
For example:

   $inx->id_parser(\&get_id);
   # make the index
   $inx->make_index($file_name);

   # here is where the retrieval key is specified
   sub get_id {
      my $line = shift;
      $line =~ /clone="(\S+)"/;
      $1;
   }

Details on configuration and additional example code are available in the
biodatabases.pod file.


=head1 FEED_BACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Ewan Birney

Email - birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let's begin the code...

package Bio::Index::GenBank;

use vars qw(@ISA);
use strict;

use Bio::Index::AbstractSeq;
use Bio::Seq;

@ISA = qw(Bio::Index::AbstractSeq);

sub _type_stamp {
    return '__GenBank_FLAT__'; # What kind of index are we?
}

sub _version {
    return 0.1;
}

=head2 _index_file

  Title   : _index_file
  Usage   : $index->_index_file($file_name, $i)
  Function: Specialized function to index GenBank format files.
            Is provided with a filename and an integer
            by make_index in its SUPER class.
  Example : 
  Returns : 
  Args    : 

=cut

sub _index_file {
	my( $self,
		 $file, # File name
		 $i     # Index-number of file being indexed
	  ) = @_;

	my $begin = 0;

	my $id_parser = $self->id_parser;

	open GENBANK,$file or 
	  $self->throw("Can't open file for read : $file");

	while (<GENBANK>) {
		if (/^LOCUS/) {
			$begin = tell(GENBANK) - length($_);
		}
		for my $id (&$id_parser($_)) {
			$self->add_record($id, $i, $begin) if $id;
		}
	}
	close GENBANK;
	1;
}

=head2 id_parser

  Title   : id_parser
  Usage   : $index->id_parser( CODE )
  Function: Stores or returns the code used by record_id to
            parse the ID for record from a string.
            Returns \&default_id_parser (see below) if not
            set. An entry will be added to
            the index for each string in the list returned.
  Example : $index->id_parser( \&my_id_parser )
  Returns : reference to CODE if called without arguments
  Args    : CODE

=cut

sub id_parser {
	my ($self,$code) = @_;

	if ($code) {
		$self->{'_id_parser'} = $code;
	}
	return $self->{'_id_parser'} || \&default_id_parser;
}

=head2 default_id_parser

  Title   : default_id_parser
  Usage   : $id = default_id_parser($line)
  Function: The default parser for GenBank.pm
  Returns : Array of specified id's
  Args    : a line string

=cut

sub default_id_parser {
	my $line = shift;
	my @accs = ();
	if ( $line =~ /^LOCUS\s+(\S+)/ ) {
		push @accs,$1;
	} elsif ( $line =~ /^ACCESSION(.*)/ ) {
		(@accs) = $1 =~ /\s*(\S+)/g;
	} elsif ( /^VERSION\s+(.*)/) {
		my $a = $1;
		$a =~ s/\s+$//;
		$a =~ s/GI\://;
		push @accs, split(/\s+/,$a);
	}
	@accs;
}

=head2 _file_format

 Title   : _file_format
 Usage   : Internal function for indexing system
 Function: Provides file format for this database
 Example :
 Returns : 
 Args    :

=cut

sub _file_format{
   my ($self,@args) = @_;
   return 'GenBank';
}

1;

__END__
