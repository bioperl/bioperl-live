#
# BioPerl module for Bio::Index::BlastTable
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chris Fields <cjfields@uiuc.edu>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::BlastTable - Indexes tabular Blast reports (-m 8 or -m 9 format) and
supports retrieval based on query accession(s)

=head1 SYNOPSIS

    use strict;
    use Bio::Index::BlastTable;
    my ($indexfile,$file1,$file2,$query);
    my $index = Bio::Index::BlastTable->new(-filename => $indexfile,
				                          -write_flag => 1);
    $index->make_index($file1,$file2);

    my $data = $index->get_stream($query);

    my $blast_result = $index->fetch_report($query);
    print "query is ", $blast_result->query_name, "\n";
    while ( my $hit = $blast_result->next_hit ) {
        print "Name ", $hit->name,"\n";
        while ( my $hsp = $hit->next_hsp ) {
            print "Score ", $hsp->score;
        }
        print "\n";
    }

=head1 DESCRIPTION

This object allows one to build an index on a tabular BLAST file (or files)
and provide quick access to the blast report for that accession.  This also
allows for ID parsing using a callback:

   $inx->id_parser(\&get_id);
   # make the index
   $inx->make_index($file_name);

   # here is where the retrieval key is specified
   sub get_id {
      my $line = shift;
      $line =~ /^>.+gi\|(\d+)/;
      $1;
   }

The indexer is capable of indexing based on multiple IDs passed back from the
callback; this is assuming of course all IDs are unique.

Note: for best results 'use strict'.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Index::BlastTable;
use strict;

use IO::String;
use Bio::SearchIO;
use base qw(Bio::Index::Abstract Bio::Root::Root);

sub _version {
	return ${Bio::Root::Version::VERSION};
}

=head2 new

  Usage   : $index = Bio::Index::Abstract->new(
                -filename    => $dbm_file,
                -write_flag  => 0,
                -dbm_package => 'DB_File',
                -verbose     => 0);

  Function: Returns a new index object.  If filename is
            specified, then open_dbm() is immediately called. 
            Bio::Index::Abstract->new() will usually be called
            directly only when opening an existing index.
  Returns : A new index object
  Args    : -filename    The name of the dbm index file.
            -write_flag  TRUE if write access to the dbm file is
                         needed.
            -dbm_package The Perl dbm module to use for the
                         index.
            -verbose     Print debugging output to STDERR if
                         TRUE.

=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
}

=head2 Bio::Index::Blast implemented methods

=cut

=head2 fetch_report

 Title   : fetch_report
 Usage   : my $blastreport = $idx->fetch_report($id);
 Function: Returns a Bio::SearchIO report object 
           for a specific blast report
 Returns : Bio::SearchIO
 Args    : valid id

=cut

sub fetch_report{
	my ($self,$id) = @_;
	my $fh = $self->get_stream($id);
	my $report = Bio::SearchIO->new(-noclose => 1,
									-format => 'blasttable',
						            -fh => $fh);
	return $report->next_result;
}


=head2 Require methods from Bio::Index::Abstract

=cut

=head2 _index_file

  Title   : _index_file
  Usage   : $index->_index_file( $file_name, $i )
  Function: Specialist function to index BLAST report file(s).
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

	open my $BLAST, '<', $file or $self->throw("Could not read file '$file': $!");
	my $indexpoint = 0;
	my $lastline = 0;
    my $last_query = '';
    my $is_m9;

    # In Windows, text files have '\r\n' as line separator, but when reading in
    # text mode Perl will only show the '\n'. This means that for a line "ABC\r\n",
    # "length $_" will report 4 although the line is 5 bytes in length.
    # We assume that all lines have the same line separator and only read current line.
    my $init_pos   = tell($BLAST);
    my $curr_line  = <$BLAST>;
    my $pos_diff   = tell($BLAST) - $init_pos;
    my $correction = $pos_diff - length $curr_line;
    seek $BLAST, $init_pos, 0; # Rewind position to proceed to read the file

	while( <$BLAST> ) {
        if (m{^#}) {
            $is_m9 ||= 1;
            if(m{^#\s+T?BLAST[PNX]}i ) {
                $indexpoint = tell($BLAST) - length($_) - $correction;
            }
            next
        }
        
        if (/^(?:([^\t]+)\t)(?:[^\t]+\t){7,}/) {
            next if $last_query eq $1;
            $indexpoint = tell($BLAST) - length($_) - $correction unless $is_m9;
            foreach my $id ($self->id_parser()->($1)) {
				$self->debug("id is $id, begin is $indexpoint\n");
				$self->add_record($id, $i, $indexpoint);
			}
            $last_query = $1;
        }
	}
}

# shamelessly stolen from Bio::Index::Fasta

=head2 id_parser

  Title   : id_parser
  Usage   : $index->id_parser( CODE )
  Function: Stores or returns the code used by record_id to
            parse the ID for record from a string.  Useful
            for (for instance) specifying a different
            parser for different flavours of blast dbs. 
            Returns \&default_id_parser (see below) if not
            set. If you supply your own id_parser
            subroutine, then it should expect a fasta
            description line.  An entry will be added to
            the index for each string in the list returned.
  Example : $index->id_parser( \&my_id_parser )
  Returns : ref to CODE if called without arguments
  Args    : CODE

=cut

sub id_parser {
	my( $self, $code ) =@_;

	if ($code) {
		$self->{'_id_parser'} = $code;
	}
	return $self->{'_id_parser'} || \&default_id_parser;
}

=head2 default_id_parser

  Title   : default_id_parser
  Usage   : $id = default_id_parser( $header )
  Function: The default Blast Query ID parser for Bio::Index::Blast.pm
            Returns $1 from applying the regexp /^>\s*(\S+)/
            to $header.
  Returns : ID string
  Args    : a header line string

=cut

sub default_id_parser
{
	if ($_[0] =~ /^\s*(\S+)/) {
		return $1;
	} else {
		return;
	}
}

=head2 Bio::Index::Abstract methods

=cut

=head2 filename

 Title   : filename
 Usage   : $value = $self->filename();
           $self->filename($value);
 Function: Gets or sets the name of the dbm index file.
 Returns : The current value of filename
 Args    : Value of filename if setting, or none if
           getting the value.

=head2 write_flag

 Title   : write_flag
 Usage   : $value = $self->write_flag();
           $self->write_flag($value);
 Function: Gets or sets the value of write_flag, which
           is wether the dbm file should be opened with
           write access.
 Returns : The current value of write_flag (default 0)
 Args    : Value of write_flag if setting, or none if
           getting the value.

=head2 dbm_package

 Usage   : $value = $self->dbm_package();
           $self->dbm_package($value);

 Function: Gets or sets the name of the Perl dbm module used. 
           If the value is unset, then it returns the value of
           the package variable $USE_DBM_TYPE or if that is
           unset, then it chooses the best available dbm type,
           choosing 'DB_File' in preference to 'SDBM_File'. 
           Bio::Abstract::Index may work with other dbm file
           types.

 Returns : The current value of dbm_package
 Args    : Value of dbm_package if setting, or none if
           getting the value.


=head2 get_stream

 Title   : get_stream
 Usage   : $stream = $index->get_stream( $id );
 Function: Returns a file handle with the file pointer
           at the approprite place

           This provides for a way to get the actual
           file contents and not an object 

           WARNING: you must parse the record deliminter
           *yourself*. Abstract wont do this for you 
           So this code

           $fh = $index->get_stream($myid);
           while( <$fh> ) {
              # do something
           }
           will parse the entire file if you do not put in
           a last statement in, like

           while( <$fh> ) {
              /^\/\// && last; # end of record
              # do something
           }

 Returns : A filehandle object
 Args    : string represents the accession number
 Notes   : This method should not be used without forethought 


=head2 open_dbm

  Usage   : $index->open_dbm()
  Function: Opens the dbm file associated with the index
            object.  Write access is only given if explicitly
            asked for by calling new(-write => 1) or having set
            the write_flag(1) on the index object.  The type of
            dbm file opened is that returned by dbm_package(). 
            The name of the file to be is opened is obtained by
            calling the filename() method.

  Example : $index->_open_dbm()
  Returns : 1 on success


=head2 _version

  Title   : _version
  Usage   : $type = $index->_version()
  Function: Returns a string which identifes the version of an
            index module.  Used to permanently identify an index
            file as having been created by a particular version
            of the index module.  Must be provided by the sub class
  Example : 
  Returns : 
  Args    : none

=head2 _filename

  Title   : _filename
  Usage   : $index->_filename( FILE INT )
  Function: Indexes the file
  Example : 
  Returns : 
  Args    : 

=head2 _file_handle

  Title   : _file_handle
  Usage   : $fh = $index->_file_handle( INT )
  Function: Returns an open filehandle for the file
            index INT.  On opening a new filehandle it
            caches it in the @{$index->_filehandle} array.
            If the requested filehandle is already open,
            it simply returns it from the array.
  Example : $fist_file_indexed = $index->_file_handle( 0 );
  Returns : ref to a filehandle
  Args    : INT

=head2 _file_count

  Title   : _file_count
  Usage   : $index->_file_count( INT )
  Function: Used by the index building sub in a sub class to
            track the number of files indexed.  Sets or gets
            the number of files indexed when called with or
            without an argument.
  Example : 
  Returns : INT
  Args    : INT


=head2 add_record

  Title   : add_record
  Usage   : $index->add_record( $id, @stuff );
  Function: Calls pack_record on @stuff, and adds the result
            of pack_record to the index database under key $id.
            If $id is a reference to an array, then a new entry
            is added under a key corresponding to each element
            of the array.
  Example : $index->add_record( $id, $fileNumber, $begin, $end )
  Returns : TRUE on success or FALSE on failure
  Args    : ID LIST

=head2 pack_record

  Title   : pack_record
  Usage   : $packed_string = $index->pack_record( LIST )
  Function: Packs an array of scalars into a single string
            joined by ASCII 034 (which is unlikely to be used
            in any of the strings), and returns it. 
  Example : $packed_string = $index->pack_record( $fileNumber, $begin, $end )
  Returns : STRING or undef
  Args    : LIST

=head2 unpack_record

  Title   : unpack_record
  Usage   : $index->unpack_record( STRING )
  Function: Splits the sting provided into an array,
            splitting on ASCII 034.
  Example : ( $fileNumber, $begin, $end ) = $index->unpack_record( $self->db->{$id} )
  Returns : A 3 element ARRAY
  Args    : STRING containing ASCII 034

=head2 DESTROY

 Title   : DESTROY
 Usage   : Called automatically when index goes out of scope
 Function: Closes connection to database and handles to
           sequence files
 Returns : NEVER
 Args    : NONE


=cut

1;
