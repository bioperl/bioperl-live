
#
# BioPerl module for Bio::Index::Abstract
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#          and James Gilbert <jgrg@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::Abstract - Abstract interface for indexing a flat file

=head1 SYNOPSIS

You should not be using this module directly

=head1 USING DB_FILE

To use DB_FILE and not SDBM for this index, do the following

    $Bio::Index::Abstract::USE_DBM_TYPE = 'DB_File';

=head1 DESCRIPTION

This object provides the basic mechanism to associate positions
in files with names. The position and filenames are stored in DBM
which can then be accessed later on. It is the equivalent of flat
file indexing (eg, SRS or efetch).

This object is the guts to the mechanism, which will be used by the
specific objects inheriting from it.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via
email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney, James Gilbert

Email - birney@sanger.ac.uk, jgrg@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with an "_" (underscore).

=cut


# Let the code begin...

package Bio::Index::Abstract;

use strict;
use Fcntl qw( O_RDWR O_CREAT O_RDONLY );
use vars qw( $TYPE_AND_VERSION_KEY $AUTOLOAD
             @ISA @EXPORT_OK $USE_DBM_TYPE $DB_HASH );

$USE_DBM_TYPE = 'SDBM_File'; # Choose 'DB_File' or 'SDMB_File'

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Exporter);
@EXPORT_OK = qw();

# new() is inherited from Bio::Root::Object

=head2 _initialize

  Title   : _initialize
  Usage   : $index->_initialize
  Function: Initializes data structures in the index object,
            opens or creates the dbm file, and performs a number
            of pre-flight checks, which  depend on the existance
            of certain methods in the sub class.  In particular
            the sub class should provide the _version() and
            _type() methods (see below).
  Example : $index->SUPER::_initialize # add this call to the
                                       # sub class's _initialize()
                                       # function if the sub
                                       # class needs to do its
                                       # own initializing
  Returns : 
  Args    : 

=cut

sub _initialize {
    my($self, $index_file, $write_flag) = @_;
    
    $index_file           or $self->throw("Index file name not given");

    $self->{'_filename'}   = $index_file;
    $self->{'_filehandle'} = []; # Array in which to cache open filehandles
    $self->{'_DB'}         = {}; # Gets tied to the DBM file
    
    # Open database
    $self->_open_dbm($write_flag);

    # scary up-cast for abstract databases.
    if( ref $self eq "Bio::Index::Abstract" ) { 
	my $type = $self->_code_base();
	bless $self, $type;
    }

    # set verbose to 0 

    $self->verbose(0);

    # Check or set this is the right kind and version of index
    $self->_type_and_version();
    
    # Check files haven't changed size since they were indexed
    $self->_check_file_sizes();
}


=head2 filename

  Title   : filename
  Usage   : $name = $index->filename()
  Function: Returns the name of the index file
  Example : 
  Returns : STRING
  Args    : NONE

=cut

sub filename {
    return $_[0]->{'_filename'};
}


=head2 db

  Title   : db
  Usage   : $index->db
  Function: Returns a ref to the hash which is tied to the dbm
            file.  Used internally when adding and retrieving
            data from the database.
  Example : $db = $index->db();
            $db->{ $some_key } = $data
            $data = $index->db->{ $some_key };
  Returns : ref to HASH
  Args    : NONE

=cut

sub db {
    return $_[0]->{'_DB'};
}


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
           will parse the entire file if you don't put in
           a last statement in, like

           while( <$fh> ) {
              /^\/\// && last; # end of record
              # do something
           }

 Returns : A filehandle object
 Args    : string represents the accession number
 Notes   : 
    This method should not be used without forethought 
=cut

sub get_stream {
   my ($self,$id) = @_;

   my ($desc,$acc,$out);
   my $db = $self->db();

   if (my $rec = $db->{ $id }) {
       my( @record );
       
       my ($file, $begin, $end) = $self->unpack_record( $rec );
        
       # Get the (possibly cached) filehandle
       my $fh = $self->_file_handle( $file );
       
       # move to start
       seek($fh, $begin, 0);
       
       return $fh;
   }
   else {
       $self->throw("Unable to find a record for $id in the flat file index");
   }
}



=head2 _open_dbm

  Title   : _open_dbm
  Usage   : $index->_open_dbm()
  Function: Called by _initialize() to create or open an existing
            index dbm file for adding and retieving records.  Write
            access is only given if explicitly asked for by passing 
            the string 'WRITE' as the second argument to the new()
            function (which calls initialize).  The type of dbm file
            opened, DB_File or SDBM_File, depends upon wether the
            global variable $USE_DBM_TYPE is set to 'DB_File' or
            'SDBM_File' respectively.  The name of the file to be is
            opened is obtained by calling the filename() method.
  Example : $index->_open_dbm('WRITE')
  Returns : 1 on success
  Args    : 'WRITE' - optional

=cut

sub _open_dbm {
    my( $self, $write_flag ) = @_;
    my( $db_type );
    
    my $index_file = $self->filename();

    my $db = $self->db();
    
    # Do we use the Berkley DB for improved speed and
    # cross platform compatability?
    if ($USE_DBM_TYPE eq 'DB_File') {
    
        require DB_File;
        DB_File->import( qw($DB_HASH) );
    
        # Only give write access if specifically asked for
        if (defined($write_flag) and $write_flag eq 'WRITE') {
            tie %$db, 'DB_File', $index_file, O_RDWR|O_CREAT, 0644, $DB_HASH
                or $self->throw("Can't open DB_File [ $index_file ] : $!");
        } else {
            tie %$db, 'DB_File', $index_file, O_RDONLY,         0644, $DB_HASH
                or $self->throw("Can't open DB_File [ $index_file ] : $!");
        }
    }
    
    # No, we use SDBM which is guaranteed to come with Perl
    # and produces smaller files
    elsif ($USE_DBM_TYPE eq 'SDBM_File') {
    
        require SDBM_File;
    
        # Only give write access if specifically asked for
        if (defined($write_flag) and $write_flag eq 'WRITE') {
            tie %$db, 'SDBM_File', $index_file, O_RDWR|O_CREAT, 0644
                or $self->throw("Can't open SDBM_File [ $index_file ] : $!");
        } else {
            tie %$db, 'SDBM_File', $index_file, O_RDONLY,         0644
                or $self->throw("Can't open SDBM_File [ $index_file ] $!");
        }
    
    } else {
        $self->throw("Don't know how to open dbm file type [ $USE_DBM_TYPE ]");
    }
    return 1;
}


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

=cut

sub _version {
    my $self = shift;
    
    $self->throw("In Bio::Index::Abstract, no _version method in sub class");
}

=head2 _code_base

 Title   : _code_base
 Usage   : $code = $db->_code_base();
 Function:
 Example :
 Returns : Code package to be used with this 
 Args    :


=cut

sub _code_base{
   my ($self) = @_;
   my $code_key    = '__TYPE_AND_VERSION';
   my $record;

   $record = $self->db->{$code_key};

   my($code,$version) = $self->unpack_record($record);
   if( wantarray ) {
       return ($code,$version);
   } else {
       return $code;
   }
}


=head2 _type_and_version

  Title   : _type_and_version
  Usage   : Called by _initalize
  Function: Checks that the index opened is made by the same index
            module and version of that module that made it.  If the
            index is empty, then it adds the information to the
            database.
  Example : 
  Returns : 1 or exception
  Args    : none

=cut

sub _type_and_version {
    my $self    = shift;
    my $key     = '__TYPE_AND_VERSION';
    my $version = $self->_version();
    my $type    = ref $self;
    
    # Add type and version key if missing else run check
    if (my $rec = $self->db->{ $key }) {
        my( $db_type, $db_version ) = $self->unpack_record($rec);
        $self->throw("This index file is from version [$db_version] - You need to rebuild it to use module version [$version]")
            unless $db_version == $version;
        $self->throw("This index file is type [$db_type] - Can't access it with module for [$type]")
            unless $db_type eq $type;
    } else {
        $self->add_record( $key, $type, $version )
            or $self->throw("Can't add Type and Version record");
    }
    return 1;
}


=head2 

  Title   : _check_file_sizes
  Usage   : $index->_check_file_sizes()
  Function: Verifies that the files listed in the database
            are the same size as when the database was built,
            or throws an exception.  Called by the _initialize()
            function.
  Example : 
  Returns : 1 or exception
  Args    : 

=cut

sub _check_file_sizes {
    my $self = shift;
    my $num  = $self->_file_count() || 0;
    
    for (my $i = 0; $i < $num; $i++) {
        my( $file, $stored_size ) = $self->unpack_record( $self->db->{"__FILE_$i"} );
        my $size = -s $file;
        unless ($size = $stored_size) {
            $self->throw("file $i [ $file ] has changed size $stored_size -> $size");
        }
    }
    return 1;
}


=head2 make_index

  Title   : make_index
  Usage   : $index->make_index( FILE_LIST )
  Function: Takes a list of file names, checks that they are
            all fully qualified, and then calls _index_file() on
            each.  It supplies _index_file() with the name of the
            file, and an integer which is stored with each record
            created by _index_file().  Can be called multiple times,
            and can be used to add to an existing index file.
  Example : $index->make_index( '/home/seqs1', '/home/seqs2', '/nfs/pub/big_db' );
  Returns : Number of files indexed
  Args    : LIST OF FILES

=cut

sub make_index {
    my($self, @files) = @_;
    my $count = 0;

    # We're really fussy/lazy, expecting all file names to be fully qualified
    $self->throw("No files to index provided") unless @files;
    foreach my $file (@files) {
        $self->throw("File name not fully qualified : $file") unless $file =~ m|^/|;
        $self->throw("File does not exist: $file")            unless -e $file;
    }

    # Add each file to the index
    FILE :
    foreach my $file (@files) {

        my $i; # index for this file
    
        # Get new index for this file and increment file count
        if ( defined(my $count = $self->_file_count) ) {
            $i = $count;
        } else {
            $i = 0; $self->_file_count(0);
        }

	# see whether this file has been already indexed
	my ($record,$number,$size);

	if( ($record = $self->db->{"__FILENAME_$file"}) ) {
	    ($number,$size) = $self->unpack_record($record);

	    # if it is the same size - fine. Otherwise die 
	    if( -s $file == $size ) {
		print "File $file already indexed. Skipping...\n";
		next FILE;
	    } else {
		$self->throw("In index, $file has changed size ($size). Indicates that the index is out of date");
	    }
	}

	# index this file
	print STDOUT "Indexing file $file\n";

	# this is supplied by the subclass and does the serious work
        $self->_index_file( $file, $i ); # Specific method for each type of index


        # Save file name and size for this index
        $self->add_record("__FILE_$i", $file, -s $file)
            or $self->throw("Can't add data to file: $file");
        $self->add_record("__FILENAME_$file", $i, -s $file)
            or $self->throw("Can't add data to file: $file");

        # increment file lines
	$i++; $self->_file_count($i);
	my $temp;
	$temp = $self->_file_count();
	

    }
    return $count;
}

=head2 _index_file

  Title   : _index_file
  Usage   : $index->_index_file( FILE INT )
  Function: Indexes the file
  Example : 
  Returns : 
  Args    : 

=cut

sub _index_file {
    my $self = shift;
    
    $self->throw("In Bio::Index::Abstract, no _index_file method provided by sub class");
}



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

=cut

sub _file_handle {
    my( $self, $i ) = @_;
    
    if (my $fh = $self->{'_filehandle'}[$i]) {
        return $fh; # Return cached filehandle
    } else {
        local *FH;
        my @rec = $self->unpack_record($self->db->{"__FILE_$i"})
            or $self->throw("Can't get filename for index : $i");
        my $file = $rec[0];
        open FH, $file or $self->throw("Can't open file for read : $file");
        $self->{'_filehandle'}[$i] = *FH; # Cache filehandle
        return *FH;
    }
}


=head2 _file_count

  Title   : _file_count
  Usage   : $index->_file_count( INT )
  Function: Used by the index building sub in a sub class to
            track the number of files indexed.  Sets or gets
            the number of files indexed when called with or
            without an argument.
  Example : 
  Returns : 
  Args    : INT

=cut

sub _file_count {
    my $self = shift;
    if (@_) {
        $self->db->{'__FILE_COUNT'} = shift;
    } else {
        return $self->db->{'__FILE_COUNT'};
    }    
}


=head2 add_record

  Title   : add_record
  Usage   : $index->add_record( $id, @stuff );
  Function: Calls pack_record on @stuff, and adds the result
            of pack_record to the index database under key $id.
  Example : $index->add_record( $id, $fileNumber, $begin, $end )
  Returns : TRUE on success or FALSE on failure
  Args    : ID LIST

=cut

sub add_record {
    my( $self, $id, @rec ) = @_;
    if( $self->verbose != 0 ) {
	print STDERR "Adding key $id\n";
    }

    $self->db->{$id} = $self->pack_record( @rec );
}


=head2 pack_record

  Title   : pack_record
  Usage   : $packed_string = $index->pack_record( LIST )
  Function: Packs an array of scalars into a single string
            joined by ASCII 034 (which is unlikely to be used
            in any of the strings), and returns it. 
  Example : $packed_string = $index->pack_record( $fileNumber, $begin, $end )
  Returns : STRING or undef
  Args    : LIST

=cut

sub pack_record {
    my( $self, @args ) = @_;
    return join "\034", @args;
}

=head2 unpack_record

  Title   : unpack_record
  Usage   : $index->unpack_record( STRING )
  Function: Splits the sting provided into an array,
            splitting on ASCII 034.
  Example : ( $fileNumber, $begin, $end ) = $index->unpack_record( $self->db->{$id} )
  Returns : A 3 element ARRAY
  Args    : STRING containing ASCII 034

=cut

sub unpack_record {
    my( $self, @args ) = @_;
    return split /\034/, $args[0];
}

=head2 verbose

 Title   : verbose
 Usage   : $obj->verbose($newval)
 Function: sets whether a report to STDERR should be issued or not
           for each sequence indexed. Helps track errors
 Example : 
 Returns : value of verbose
 Args    : newvalue (optional)


=cut

sub verbose{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'verbose'} = $value;
    }
    return $obj->{'verbose'};

}


=head2 DESTROY

 Title   : DESTROY
 Usage   : Called automatically when index goes out of scope
 Function: Closes connection to database and handles to
           sequence files
 Returns : NEVER
 Args    : NONE


=cut

sub DESTROY {
    my $self = shift;
    
    foreach my $fh (@{$self->{'_filehandle'}}) {
        close $fh;
    }
    
    untie %{$self->db};
}

1;





