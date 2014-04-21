#
# BioPerl module for Bio::Index::Swissprot
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::Swissprot - Interface for indexing one or more
Swissprot files.

=head1 SYNOPSIS

  # Make an index for one or more Swissprot files:

    use Bio::Index::Swissprot;
    use strict;

    my $index_file_name = shift;
    my $inx = Bio::Index::Swissprot->new(
        -filename => $index_file_name,
        -write_flag => 1);
    $inx->make_index(@ARGV);

  # Print out several sequences present in the index in Genbank
  # format:

    use Bio::Index::Swissprot;
    use Bio::SeqIO;
    use strict;

    my $out = Bio::SeqIO->new( -format => 'genbank',
                               -fh => \*STDOUT );
    my $index_file_name = shift;
    my $inx = Bio::Index::Swissprot->new(-filename => $index_file_name);

    foreach my $id (@ARGV) {
        my $seq = $inx->fetch($id); # Returns a Bio::Seq object
        $out->write_seq($seq);
    }

    # alternatively
    my ($id, $acc);
    my $seq1 = $inx->get_Seq_by_id($id);
    my $seq2 = $inx->get_Seq_by_acc($acc);

=head1 DESCRIPTION

By default the index that is created uses the AC and ID identifiers
as keys. This module inherits functions for managing dbm files from 
Bio::Index::Abstract.pm, and provides the basic functionality 
for indexing Swissprot files and retrieving Sequence objects from 
them. For best results 'use strict'.

You can also set or customize the unique key used to retrieve by 
writing your own function and calling the id_parser() method.
For example:

   $inx->id_parser(\&get_id);
   # make the index
   $inx->make_index($index_file_name);

   # here is where the retrieval key is specified
   sub get_id {
      my $line = shift;
      $line =~ /^KW\s+([A-Z]+)/i;
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
the bugs and their resolution.  Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ewan Birney

Also lorenz@ist.org, bosborne at alum.mit.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let's begin the code...

package Bio::Index::Swissprot;

use strict;

use Bio::Seq;

use base qw(Bio::Index::AbstractSeq);

sub _type_stamp {
	return '__Swissprot_FLAT__'; # What kind of index are we?
}

sub _version {
	return 0.1;
}

=head2 _index_file

  Title   : _index_file
  Usage   : $index->_index_file( $file_name, $i )
  Function: Specialist function to index Swissprot format files.
            Is provided with a filename and an integer
            by make_index in its SUPER class.
  Example : 
  Returns : 
  Args    : 

=cut

sub _index_file {
	# $file is file name, $i is number of file being indexed
	my( $self, $file, $i ) = @_;

	# Offset from start of file
	my $begin = 0;

	my $id_parser = $self->id_parser;

	open my $SWISSPROT, '<', $file or $self->throw("Could not read file '$file': $!");
    
        my %done_ids;

	# In Windows, text files have '\r\n' as line separator, but when reading in
	# text mode Perl will only show the '\n'. This means that for a line "ABC\r\n",
	# "length $_" will report 4 although the line is 5 bytes in length.
	# We assume that all lines have the same line separator and only read current line.
	my $init_pos   = tell($SWISSPROT);
	my $curr_line  = <$SWISSPROT>;
	my $pos_diff   = tell($SWISSPROT) - $init_pos;
	my $correction = $pos_diff - length $curr_line;
	seek $SWISSPROT, $init_pos, 0; # Rewind position to proceed to read the file

	while (<$SWISSPROT>) {
		if (/^ID\s+\S+/) {
			$begin = tell($SWISSPROT) - length( $_ ) - $correction;
		}
		for my $id (&$id_parser($_)) {
                        next if exists $done_ids{$id};
  			$self->add_record($id, $i, $begin) if $id;
                        $done_ids{$id} = 1;
		}
        if (m{//}) {
            %done_ids = ();
        }
	}
	close $SWISSPROT;
	return 1;
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
  Usage   : $id = default_id_parser( $line )
  Function: The default parser for Swissprot.pm
            Returns $1 from applying the regexp /^ID\s*(\S+)/
            or /^AC\s+([A-Z0-9]+)/ to the current line.
  Returns : ID string
  Args    : a line string

=cut

sub default_id_parser {
	my $line = shift;
	if ($line =~ /^ID\s*(\S+)/) {
		return $1;
	} 
	elsif ($line =~ /^AC\s+([A-Z0-9]+)/) {
		return $1;
	}
}

=head2 _file_format

 Title   : _file_format
 Usage   : Internal function for indexing system
 Function: Provides file format for this database
 Example :
 Returns :
 Args    :


=cut

sub _file_format {
   return 'swiss';
}

1;

__END__
