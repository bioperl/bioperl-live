#
# $Id$
#
# BioPerl module for Bio::Index::Abstract
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::Swissprot - Interface for indexing one or more
Swissprot files.

=head1 SYNOPSIS

    # Complete code for making an index for several
    # Swissprot files
    use Bio::Index::Swissprot;
    use strict;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::Swissprot->new(
                           -filename => $Index_File_Name,
					            -write_flag => 'WRITE');
    $inx->make_index(@ARGV);

    # Print out several sequences present in the index
    # in gcg format
    use Bio::Index::Swissprot;
    use Bio::SeqIO;
    use strict;

    my $out = Bio::SeqIO->new( -format => 'gcg',
                               -fh => \*STDOUT );
    my $Index_File_Name = shift;
    my $inx = Bio::Index::Swissprot->new(-filename => $Index_File_Name);

    foreach my $id (@ARGV) {
        my $seq = $inx->fetch($ID); # Returns Bio::Seq object
        $out->write_seq($seq);
    }

    # alternatively
    my $seq1 = $inx->get_Seq_by_id($ID);
    my $seq2 = $inx->get_Seq_by_acc($AC);

=head1 DESCRIPTION

By default the index that's created uses the AC and ID identifiers
as keys. Inherits functions for managing dbm files from 
Bio::Index::Abstract.pm, and provides the basic functionality 
for indexing Swissprot files, and retrieving sequence objects from 
them. For best results 'use strict'.

You can also set or customize the unique key used to retrieve by 
writing your own function and calling the id_parser() method.
For example:

   $inx->id_parser(\&get_id);
   # make the index
   $inx->make_index($file_name);

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

  bioperl-l@bioperl.org                - General discussion
  http://bio.perl.org/MailList.html    - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Ewan Birney

Also lorenz@ist.org, bosborne@cognia.com

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let's begin the code...

package Bio::Index::Swissprot;

use vars qw( @ISA);
use strict qw(vars);

use Bio::Index::AbstractSeq;
use Bio::Seq;

@ISA = qw(Bio::Index::AbstractSeq);

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

	open SWISSPROT,$file or $self->throw("Can't read file: $file");

	while (<SWISSPROT>) {
		if (/^ID\s+\S+/) {
			$begin = tell(SWISSPROT) - length( $_ );
		}
		for my $id (&$id_parser($_)) {
			$self->add_record($id, $i, $begin) if $id;
		}
	}
	close SWISSPROT;
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
	} elsif ($line =~ /^AC\s+([A-Z0-9]+)/) {
		return $1;
	}
	#else {
		#return;
	#}
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
   return 'swiss';
}

1;

__END__
