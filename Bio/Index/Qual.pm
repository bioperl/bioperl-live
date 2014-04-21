#
#
# BioPerl module for Bio::Index::Qual
#
# Copied almost verbatim from James Gilbert's Bio::Index::Fasta 
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::Qual - Interface for indexing (multiple) fasta qual files

=head1 SYNOPSIS

    # Complete code for making an index for several
    # qual files
    use Bio::Index::Qual;
    use strict;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::Qual->new(
        '-filename' => $Index_File_Name,
        '-write_flag' => 1);
    $inx->make_index(@ARGV);

    # Print out several sequences present in the index
    # in Qual format
    use Bio::Index::Qual;
    use strict;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::Qual->new('-filename' => $Index_File_Name);
    my $out = Bio::SeqIO->new('-format' => 'qual','-fh' => \*STDOUT);

    foreach my $id (@ARGV) {
        my $seq = $inx->fetch($id); # Returns Bio::Seq object
	$out->write_seq($seq);
    }

    # or, alternatively
    my $id;
    my $seq = $inx->get_Seq_by_id($id); #identical to fetch

=head1 DESCRIPTION

Inherits functions for managing dbm files from Bio::Index::Abstract.pm,
and provides the basic funtionallity for indexing qual files, and
retrieving the sequence from them. For best results 'use strict'.

Bio::Index::Qual supports the Bio::DB::BioSeqI interface, meaning
it can be used as a Sequence database for other parts of bioperl

Additional example code is available in scripts/index/*PLS and in 
the Bioperl Tutorial (L<http://www.bioperl.org/wiki/Bptutorial.pl>).

Note that by default the key for the sequence will be the first continuous
string after the 'E<gt>' in the qual header. If you want to use a specific
substring of the qual header you must use the id_parser() method.

You can also set or customize the unique key used to retrieve by 
writing your own function and calling the id_parser() method.
For example:

   $inx->id_parser(\&get_id);
   # make the index
   $inx->make_index($file_name);

   # here is where the retrieval key is specified
   sub get_id {
      my $line = shift;
      $line =~ /^(\d+)/;
      $1;
   }

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

=head1 AUTHOR - James Gilbert, Mark Johnson

Email - jgrg@sanger.ac.uk, johnsonm-at-gmail-dot-com 

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Index::Qual;

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
 Returns : 'qual'

=cut

sub _file_format {
	return 'qual';
}



=head2 _index_file

  Title   : _index_file
  Usage   : $index->_index_file( $file_name, $i )
  Function: Specialist function to index QUAL format files.
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

	open my $QUAL, '<', $file or $self->throw("Could not read file '$file': $!");

	# In Windows, text files have '\r\n' as line separator, but when reading in
	# text mode Perl will only show the '\n'. This means that for a line "ABC\r\n",
	# "length $_" will report 4 although the line is 5 bytes in length.
	# We assume that all lines have the same line separator and only read current line.
	my $init_pos   = tell($QUAL);
	my $curr_line  = <$QUAL>;
	my $pos_diff   = tell($QUAL) - $init_pos;
	my $correction = $pos_diff - length $curr_line;
	seek $QUAL, $init_pos, 0; # Rewind position to proceed to read the file

	# Main indexing loop
	while (<$QUAL>) {
		if (/^>/) {
			my $begin = tell($QUAL) - length( $_ ) + 1 - $correction;

			foreach my $id (&$id_parser($_)) {
				$self->add_record($id, $i, $begin);
			}
		}
	}

	close $QUAL;
	return 1;
}

=head2 id_parser

  Title   : id_parser
  Usage   : $index->id_parser( CODE )
  Function: Stores or returns the code used by record_id to
            parse the ID for record from a string.  Useful
            for (for instance) specifying a different
            parser for different flavours of Qual file. 
            Returns \&default_id_parser (see below) if not
            set. If you supply your own id_parser
            subroutine, then it should expect a qual
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
  Function: The default Qual ID parser for Qual.pm
            Returns $1 from applying the regexp /^>\s*(\S+)/
            to $header.
  Returns : ID string
  Args    : a qual header line string

=cut

sub default_id_parser {    
	if ($_[0] =~ /^>\s*(\S+)/) {
		return $1;
	} else {
		return;
	}
}

1;
