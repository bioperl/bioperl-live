# BioPerl module for Bio::Index::Fastq
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Tony Cox <avc@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::Fastq - Interface for indexing (multiple) fastq files

=head1 SYNOPSIS

    # Complete code for making an index for several
    # fastq files
    use Bio::Index::Fastq;
    use strict;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::Fastq->new(
        '-filename' => $Index_File_Name,
        '-write_flag' => 1);
    $inx->make_index(@ARGV);

    # Print out several sequences present in the index
    # in Fastq format
    use Bio::Index::Fastq;
    use strict;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::Fastq->new('-filename' => $Index_File_Name);
    my $out = Bio::SeqIO->new('-format' => 'Fastq','-fh' => \*STDOUT);

    foreach my $id (@ARGV) {
        my $seq = $inx->fetch($id); # Returns Bio::Seq::Quality object
	$out->write_seq($seq);
    }

    # or, alternatively
    my $id;
    my $seq = $inx->get_Seq_by_id($id); #identical to fetch

=head1 DESCRIPTION

Inherits functions for managing dbm files from Bio::Index::Abstract.pm,
and provides the basic funtionallity for indexing fastq files, and
retrieving the sequence from them. Note: for best results 'use strict'.

Bio::Index::Fastq supports the Bio::DB::BioSeqI interface, meaning
it can be used as a Sequence database for other parts of bioperl

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

=head1 AUTHOR - Tony Cox

Email - avc@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Index::Fastq;

use strict;

use Bio::Seq;

use base qw(Bio::Index::AbstractSeq);

#
# Suggested fix by Michael G Schwern <schwern@pobox.com> to
# get around a clash with CPAN shell...
#

sub _version {
    return 0.2;
}

=head2 _file_format

 Title   : _file_format
 Function: The file format for this package, which is needed
           by the SeqIO system when reading the sequence.
 Returns : 'Fastq'

=cut

sub _file_format {
    return 'Fastq';
}



=head2 _index_file

  Title   : _index_file
  Usage   : $index->_index_file( $file_name, $i )
  Function: Specialist function to index FASTQ format files.
            Is provided with a filename and an integer
            by make_index in its SUPER class.
  Example : 
  Returns : 
  Args    : 

=cut

sub _index_file {
    my( $self,
        $file, # File name
        $i,    # Index-number of file being indexed
        ) = @_;
    
    my( $begin,     # Offset from start of file of the start
                    # of the last found record.
        );

    $begin = 0;

    my $id_parser = $self->id_parser;
    my $c = 0;
    open my $FASTQ, '<', $file or $self->throw("Could not read file '$file': $!");

    # In Windows, text files have '\r\n' as line separator, but when reading in
    # text mode Perl will only show the '\n'. This means that for a line "ABC\r\n",
    # "length $_" will report 4 although the line is 5 bytes in length.
    # We assume that all lines have the same line separator and only read current line.
    my $init_pos   = tell($FASTQ);
    my $curr_line  = <$FASTQ>;
    my $pos_diff   = tell($FASTQ) - $init_pos;
    my $correction = $pos_diff - length $curr_line;
    seek $FASTQ, $init_pos, 0; # Rewind position to proceed to read the file

    # Main indexing loop
    while (<$FASTQ>) {
        if (/^@/) {
            my $begin = tell($FASTQ) - length( $_ ) - $correction;
            foreach my $id (&$id_parser($_)) {
                $self->add_record($id, $i, $begin);
                $c++;
            }
        }
    }

    close $FASTQ;
    return ($c);
}

=head2 id_parser

  Title   : id_parser
  Usage   : $index->id_parser( CODE )
  Function: Stores or returns the code used by record_id to
            parse the ID for record from a string.  Useful
            for (for instance) specifying a different
            parser for different flavours of FASTQ file. 
            Returns \&default_id_parser (see below) if not
            set. If you supply your own id_parser
            subroutine, then it should expect a fastq
            description line.  An entry will be added to
            the index for each string in the list returned.
  Example : $index->id_parser( \&my_id_parser )
  Returns : ref to CODE if called without arguments
  Args    : CODE

=cut

sub id_parser {
    my( $self, $code ) = @_;
    
    if ($code) {
        $self->{'_id_parser'} = $code;
    }
    return $self->{'_id_parser'} || \&default_id_parser;
}



=head2 default_id_parser

  Title   : default_id_parser
  Usage   : $id = default_id_parser( $header )
  Function: The default Fastq ID parser for Fastq.pm
            Returns $1 from applying the regexp /^>\s*(\S+)/
            to $header.
  Returns : ID string
  Args    : a fastq header line string

=cut

sub default_id_parser {    
    if ($_[0] =~ /^@\s*(\S+)/) {
        return $1;
    } else {
        return;
    }
}

1;
