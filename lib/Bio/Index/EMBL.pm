#
# BioPerl module for Bio::Index::EMBL
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::EMBL - Interface for indexing (multiple) EMBL/Swissprot
.dat files (i.e. flat file EMBL/Swissprot format).

=head1 SYNOPSIS

    # Complete code for making an index for several
    # EMBL files
    use Bio::Index::EMBL;
    use strict;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::EMBL->new(-filename => $Index_File_Name,
				    -write_flag => 'WRITE');
    $inx->make_index(@ARGV);

    # Print out several sequences present in the index
    # in Fasta format
    use Bio::Index::EMBL;
    use strict;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::EMBL->new(-filename => $Index_File_Name);
    my $out = Bio::SeqIO->new(-format => 'Fasta',-fh => \*STDOUT);

    foreach my $id (@ARGV) {
        my $seq = $inx->fetch($id); # Returns Bio::Seq object
	$out->write_seq($seq);
    }

    # alternatively
    my ($id, $acc);
    my $seq1 = $inx->get_Seq_by_id($id);
    my $seq2 = $inx->get_Seq_by_acc($acc);

=head1 DESCRIPTION

Inherits functions for managing dbm files from Bio::Index::Abstract.pm,
and provides the basic funtionallity for indexing EMBL files, and
retrieving the sequence from them. Heavily snaffled from James Gilbert
and his Fasta system. Note: for best results 'use strict'.

The keys are the identifiers in the ID and AC lines.

=head1 FEED_BACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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

=head1 AUTHOR - Ewan Birney

Email - birney@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let's begin the code...


package Bio::Index::EMBL;

use strict;
use Bio::Seq;

use base qw(Bio::Index::AbstractSeq);

sub _type_stamp {
    return '__EMBL_FLAT__'; # What kind of index are we?
}


sub _version {
    return 0.1;
}

=head2 _index_file

  Title   : _index_file
  Usage   : $index->_index_file( $file_name, $i )
  Function: Specialist function to index EMBL format files.
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

    my( $begin, # Offset from start of file of the start
                # of the last found record.
        $id,    # ID of last found record.
	@accs,   # accession of last record. Also put into the index
        );

    $begin = 0;

    open my $EMBL, '<', $file or $self->throw("Could not read file '$file': $!");

    # In Windows, text files have '\r\n' as line separator, but when reading in
    # text mode Perl will only show the '\n'. This means that for a line "ABC\r\n",
    # "length $_" will report 4 although the line is 5 bytes in length.
    # We assume that all lines have the same line separator and only read current line.
    my $init_pos   = tell($EMBL);
    my $curr_line  = <$EMBL>;
    my $pos_diff   = tell($EMBL) - $init_pos;
    my $correction = $pos_diff - length $curr_line;
    seek $EMBL, $init_pos, 0; # Rewind position to proceed to read the file

    # Main indexing loop
    $id = undef;
    @accs = ();
    while (<$EMBL>) {
	if( m{^//} ) {
	    if( ! defined $id ) {
		$self->throw("Got to a end of entry line for an EMBL flat file with no parsed ID. Considering this a problem!");
		next;
	    }
	    if( ! @accs ) {
		$self->warn("For id [$id] in embl flat file, got no accession number. Storing id index anyway");
	    }

	    $self->add_record($id, $i, $begin);

	    foreach my $acc (@accs) {
		if( $acc ne $id ) {
		    $self->add_record($acc, $i, $begin);
		}
	    }
	} elsif (/^ID\s+(\S+)/) {
	    $id = $1;
	    # not sure if I like this. Assummes tell is in bytes.
	    # we could tell before each line and save it.
            $begin = tell($EMBL) - length( $_ ) - $correction;
	
	} elsif (/^AC\s+(.*)?/) {
            push @accs , split (/[; ]+/, $1);
	} else {
	    # do nothing
	}
    }

    close $EMBL;
    return 1;
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

   return 'EMBL';
}



1;
