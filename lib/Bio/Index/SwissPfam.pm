#
# BioPerl module for Bio::Index::SwissPfam
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#

# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::SwissPfam - Interface for indexing swisspfam files

=head1 SYNOPSIS

    use Bio::Index::SwissPfam;
    use strict;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::SwissPfam->new('-filename' => $Index_File_Name, 
                         		 '-write_flag' => 'WRITE');
    $inx->make_index(@ARGV);

    use Bio::Index::SwissPfam;
    use strict;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::SwissPfam->new('-filename' => $Index_File_Name);

    foreach my $id (@ARGV) {
        my $seq = $inx->fetch($id); # Returns stream
	     while( <$seq> ) {
	         if(/^>/) {
	    	       print;
		          last;
	         }
	     }
    }


=head1 DESCRIPTION

SwissPfam is one of the flat files released with Pfam. This modules
provides a way of indexing this module.

Inherits functions for managing dbm files from Bio::Index::Abstract.pm, and 
provides the basic funtionallity for indexing SwissPfam files.  Only 
retrieves FileStreams at the moment. Once we have something better 
(ie, an object!), will use that. Heavily snaffled from Index::Fasta system of 
James Gilbert. Note: for best results 'use strict'.

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

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let's begin the code...

package Bio::Index::SwissPfam;

use strict;
use Bio::Seq;

use base qw(Bio::Index::Abstract);

sub _type_stamp {
    return '__SWISSPFAM_FLAT__'; # What kind of index are we?
}

sub _version {
    return 0.1;
}

=head2 _index_file

  Title   : _index_file
  Usage   : $index->_index_file( $file_name, $i )
  Function: Specialist function to index swisspfam format files.
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
        $end,   # Offset from start of file of the end
                # of the last found record.
        $id,    # ID of last found record.
	$acc,   # accession of last record. Also put into the index
	$nid, $nacc, # new ids for the record just found
        );

    $begin = 0;
    $end   = 0;

    open my $SP, '<', $file or $self->throw("Could not read file '$file': $!");

    # In Windows, text files have '\r\n' as line separator, but when reading in
    # text mode Perl will only show the '\n'. This means that for a line "ABC\r\n",
    # "length $_" will report 4 although the line is 5 bytes in length.
    # We assume that all lines have the same line separator and only read current line.
    my $init_pos   = tell($SP);
    my $curr_line  = <$SP>;
    my $pos_diff   = tell($SP) - $init_pos;
    my $correction = $pos_diff - length $curr_line;
    seek $SP, $init_pos, 0; # Rewind position to proceed to read the file

    # Main indexing loop
    while (<$SP>) {
        if (/^>(\S+)\s+\|=*\|\s+(\S+)/) {
	    $nid = $1;
	    $nacc = $2;
            my $new_begin = tell($SP) - length( $_ ) - $correction;
            $end = $new_begin - 1;

	    if( $id ) {
		$self->add_record($id, $i, $begin, $end);
		if( $acc ne $id ) {
		    $self->add_record($acc, $i, $begin, $end);
		}
	    }
            $begin = $new_begin;
	    $id = $nid;
	    $acc = $nacc;
        }
    }
    # Don't forget to add the last record
    $end = tell($SP);
    $self->add_record($id, $i, $begin, $end) if $id;

    close $SP;
    return 1;
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
    my $desc;
    my $db = $self->db();
    if (my $rec = $db->{ $id }) {
        my( @record );
        my ($file, $begin, $end) = $self->unpack_record( $rec );
        # Get the (possibly cached) filehandle
        my $fh = $self->_file_handle( $file );

        # move to start
        seek($fh, $begin, 0);

        return $fh;
    } else {
	$self->throw("Unable to find a record for $id in SwissPfam flat file index");
    }
}

1;
