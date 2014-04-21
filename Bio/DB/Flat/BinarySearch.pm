#
# BioPerl module for Bio::DB::Flat::BinarySearch
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Michele Clamp <michele@sanger.ac.uk>>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Flat::BinarySearch - BinarySearch search indexing system for sequence files

=head1 SYNOPSIS

  TODO: SYNOPSIS NEEDED!

=head1 DESCRIPTION

This module can be used both to index sequence files and also to
retrieve sequences from existing sequence files.

This object allows indexing of sequence files both by a primary key
(say accession) and multiple secondary keys (say ids).  This is
different from the Bio::Index::Abstract (see L<Bio::Index::Abstract>)
which uses DBM files as storage.  This module uses a binary search to
retrieve sequences which is more efficient for large datasets.

=head2 Index creation

    my $sequencefile;  # Some fasta sequence file

Patterns have to be entered to define where the keys are to be indexed
and also where the start of each record.  E.g. for fasta

    my $start_pattern   = '^>';
    my $primary_pattern = '^>(\S+)';

So the start of a record is a line starting with a E<gt> and the
primary key is all characters up to the first space after the E<gt>

A string also has to be entered to defined what the primary key
(primary_namespace) is called.

The index can now be created using

    my $index = Bio::DB::Flat::BinarySearch->new(
             -directory         => "/home/max/",
             -dbname            => "mydb",
              -start_pattern     => $start_pattern,
              -primary_pattern   => $primary_pattern,
             -primary_namespace => "ID",
              -format            => "fasta" );

    my @files = ("file1","file2","file3");

    $index->build_index(@files);

The index is now ready to use.  For large sequence files the perl way
of indexing takes a *long* time and a *huge* amount of memory.  For
indexing things like dbEST I recommend using the DB_File indexer, BDB.

The formats currently supported by this module are fasta, Swissprot,
and EMBL.

=head2 Creating indices with secondary keys

Sometimes just indexing files with one id per entry is not enough.  For
instance you may want to retrieve sequences from swissprot using
their accessions as well as their ids.

To be able to do this when creating your index you need to pass in
a hash of secondary_patterns which have their namespaces as the keys
to the hash.

e.g. For Indexing something like

ID   1433_CAEEL     STANDARD;      PRT;   248 AA.
AC   P41932;
DT   01-NOV-1995 (Rel. 32, Created)
DT   01-NOV-1995 (Rel. 32, Last sequence update)
DT   15-DEC-1998 (Rel. 37, Last annotation update)
DE   14-3-3-LIKE PROTEIN 1.
GN   FTT-1 OR M117.2.
OS   Caenorhabditis elegans.
OC   Eukaryota; Metazoa; Nematoda; Chromadorea; Rhabditida; Rhabditoidea;
OC   Rhabditidae; Peloderinae; Caenorhabditis.
OX   NCBI_TaxID=6239;
RN   [1]

where we want to index the accession (P41932) as the primary key and the
id (1433_CAEEL) as the secondary id.  The index is created as follows

    my %secondary_patterns;

    my $start_pattern   = '^ID   (\S+)';
    my $primary_pattern = '^AC   (\S+)\;';

    $secondary_patterns{"ID"} = '^ID   (\S+)';

    my $index = Bio::DB::Flat::BinarySearch->new(
                -directory          => $index_directory,
                  -dbname             => "ppp",
                  -write_flag         => 1,
                -verbose            => 1,
                -start_pattern      => $start_pattern,
                -primary_pattern    => $primary_pattern,
                -primary_namespace  => 'AC',
                -secondary_patterns => \%secondary_patterns);

    $index->build_index($seqfile);

Of course having secondary indices makes indexing slower and use more
memory.

=head2 Index reading

To fetch sequences using an existing index first of all create your sequence
object

    my $index = Bio::DB::Flat::BinarySearch->new(
                  -directory => $index_directory);

Now you can happily fetch sequences either by the primary key or
by the secondary keys.

    my $entry = $index->get_entry_by_id('HBA_HUMAN');

This returns just a string containing the whole entry.  This is
useful is you just want to print the sequence to screen or write it to a file.

Other ways of getting sequences are

    my $fh = $index->get_stream_by_id('HBA_HUMAN');

This can then be passed to a seqio object for output or converting
into objects.

    my $seq = Bio::SeqIO->new(-fh     => $fh,
                                -format => 'fasta');

The last way is to retrieve a sequence directly.  This is the
slowest way of extracting as the sequence objects need to be made.

    my $seq = $index->get_Seq_by_id('HBA_HUMAN');

To access the secondary indices the secondary namespace needs to be known

    $index->secondary_namespaces("ID");

Then the following call can be used

    my $seq   = $index->get_Seq_by_secondary('ID','1433_CAEEL');

These calls are not yet implemented

    my $fh    = $index->get_stream_by_secondary('ID','1433_CAEEL');
    my $entry = $index->get_entry_by_secondary('ID','1433_CAEEL');

=head1 FEEDBACK

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

=head1 AUTHOR - Michele Clamp

Email - michele@sanger.ac.uk

=head1 CONTRIBUTORS

Jason Stajich, jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with an "_" (underscore).

=cut

package Bio::DB::Flat::BinarySearch;

use strict;

use Fcntl qw(SEEK_END SEEK_CUR);

# rather than using tell which might be buffered
sub systell { sysseek( $_[0], 0, SEEK_CUR ) }
sub syseof  { sysseek( $_[0], 0, SEEK_END ) }

use File::Spec;
use Bio::Root::RootI;
use Bio::SeqIO;
use Bio::Seq;

use base qw(Bio::DB::RandomAccessI);

use constant CONFIG_FILE_NAME => 'config.dat';
use constant HEADER_SIZE      => 4;
use constant DEFAULT_FORMAT   => 'fasta';

my @formats = [ 'FASTA', 'SWISSPROT', 'EMBL' ];

=head2 new

 Title   : new
 Usage   : For reading
             my $index = Bio::DB::Flat::BinarySearch->new(
                     -directory => '/Users/michele/indices/dbest',
             -dbname    => 'mydb',
                     -format    => 'fasta');

           For writing

             my %secondary_patterns = {"ACC" => "^>\\S+ +(\\S+)"}
             my $index = Bio::DB::Flat::BinarySearch->new(
             -directory          => '/Users/michele/indices',
                     -dbname             => 'mydb',
             -primary_pattern    => "^>(\\S+)",
                     -secondary_patterns => \%secondary_patterns,
             -primary_namespace  => "ID");

             my @files = ('file1','file2','file3');

             $index->build_index(@files);


 Function: create a new Bio::DB::Flat::BinarySearch object
 Returns : new Bio::DB::Flat::BinarySearch
 Args    : -directory          Root directory for index files
           -dbname             Name of subdirectory containing indices
                               for named database
           -write_flag         Allow building index
           -primary_pattern    Regexp defining the primary id
           -secondary_patterns A hash ref containing the secondary
                               patterns with the namespaces as keys
           -primary_namespace  A string defining what the primary key
                               is

 Status  : Public

=cut

sub new {
    my ( $class, @args ) = @_;

    my $self = $class->SUPER::new(@args);

    bless $self, $class;

    my ( $index_dir, $dbname, $format, $write_flag, $primary_pattern,
        $primary_namespace, $start_pattern, $secondary_patterns )
      = $self->_rearrange(
        [
            qw(DIRECTORY
              DBNAME
              FORMAT
              WRITE_FLAG
              PRIMARY_PATTERN
              PRIMARY_NAMESPACE
              START_PATTERN
              SECONDARY_PATTERNS)
        ],
        @args
      );

    $self->index_directory($index_dir);
    $self->dbname($dbname);

    if ( $self->index_directory && $self->read_config_file ) {

        my $fh           = $self->primary_index_filehandle;
        my $record_width = $self->read_header($fh);
        $self->record_size($record_width);
    }
    $format ||= DEFAULT_FORMAT;
    $self->format($format);
    $self->write_flag($write_flag);

    if ( $self->write_flag && !$primary_namespace ) {
        (
            $primary_namespace, $primary_pattern,
            $start_pattern,     $secondary_patterns
        ) = $self->_guess_patterns( $self->format );
    }

    $self->primary_pattern($primary_pattern);
    $self->primary_namespace($primary_namespace);
    $self->start_pattern($start_pattern);
    $self->secondary_patterns($secondary_patterns);

    return $self;
}

sub new_from_registry {
    my ( $self, %config ) = @_;

    my $dbname   = $config{'dbname'};
    my $location = $config{'location'};

    my $index = Bio::DB::Flat::BinarySearch->new(
        -dbname    => $dbname,
        -index_dir => $location,
    );
}

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $obj->get_Seq_by_id($newval)
 Function:
 Example :
 Returns : value of get_Seq_by_id
 Args    : newvalue (optional)

=cut

sub get_Seq_by_id {
    my ( $self, $id ) = @_;

    # too many uninit variables...
    local $^W = 0;

    my ( $fh, $length ) = $self->get_stream_by_id($id);

    unless ( defined( $self->format ) ) {
        $self->throw("Can't create sequence - format is not defined");
    }

    return unless $fh;

    unless ( defined( $self->{_seqio} ) ) {

        $self->{_seqio} = Bio::SeqIO->new(
            -fh     => $fh,
            -format => $self->format
        );
    }
    else {
        $self->{_seqio}->fh($fh);
    }

    return $self->{_seqio}->next_seq;
}

=head2 get_entry_by_id

 Title   : get_entry_by_id
 Usage   : $obj->get_entry_by_id($newval)
 Function: Get a Bio::SeqI object for a unique ID
 Returns : Bio::SeqI
 Args    : string


=cut

sub get_entry_by_id {
    my ( $self, $id ) = @_;

    my ( $fh, $length ) = $self->get_stream_by_id($id);

    my $entry;

    sysread( $fh, $entry, $length );

    return $entry;
}

=head2 get_stream_by_id

 Title   : get_stream_by_id
 Usage   : $obj->get_stream_by_id($id)
 Function: Gets a Sequence stream for an id
 Returns : Bio::SeqIO stream
 Args    : Id to lookup by


=cut

sub get_stream_by_id {
    my ( $self, $id ) = @_;

    unless ( $self->record_size ) {
        if ( $self->index_directory && $self->read_config_file ) {

            my $fh           = $self->primary_index_filehandle;
            my $record_width = $self->read_header($fh);
            $self->record_size($record_width);
        }
    }
    my $indexfh = $self->primary_index_filehandle;
    syseof($indexfh);

    my $filesize = systell($indexfh);

    $self->throw("file was not parsed properly, record size is empty")
      unless $self->record_size;

    my $end = ( $filesize - $self->{'_start_pos'} ) / $self->record_size;
    my ( $newid, $rest, $fhpos ) =
      $self->find_entry( $indexfh, 0, $end, $id, $self->record_size );

    my ( $fileid, $pos, $length ) = split( /\t/, $rest );

#print STDERR "BinarySearch Found id entry $newid $fileid $pos $length:$rest\n";

    if ( !$newid ) {
        return;
    }

    my $file = $self->{_file}{$fileid};

    open my $IN, '<', $file or $self->throw("Could not read file '$file': $!");

    my $entry;

    sysseek( $IN, $pos, 0 );

    return ( $IN, $length );
}

=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   : $obj->get_Seq_by_acc($acc)
 Function: Gets a Bio::SeqI object by accession number
 Returns : Bio::SeqI object
 Args    : string representing accession number


=cut

sub get_Seq_by_acc {
    my ( $self, $acc ) = @_;

    # too many uninit variables...
    local $^W = 0;

    if ( $self->primary_namespace eq "ACC" ) {
        return $self->get_Seq_by_id($acc);
    }
    else {
        return $self->get_Seq_by_secondary( "ACC", $acc );
    }
}

=head2 get_Seq_by_version

 Title   : get_Seq_by_version
 Usage   : $obj->get_Seq_by_version($version)
 Function: Gets a Bio::SeqI object by accession.version number
 Returns : Bio::SeqI object
 Args    : string representing accession.version number


=cut

sub get_Seq_by_version {
    my ( $self, $acc ) = @_;

    # too many uninit variables...
    local $^W = 0;

    if ( $self->primary_namespace eq "VERSION" ) {
        return $self->get_Seq_by_id($acc);
    }
    else {
        return $self->get_Seq_by_secondary( "VERSION", $acc );
    }
}

=head2 get_Seq_by_secondary

 Title   : get_Seq_by_secondary
 Usage   : $obj->get_Seq_by_secondary($namespace,$acc)
 Function: Gets a Bio::SeqI object looking up secondary accessions
 Returns : Bio::SeqI object
 Args    : namespace name to check secondary namespace and an id


=cut

sub get_Seq_by_secondary {
    my ( $self, $name, $id ) = @_;

    my @names = $self->secondary_namespaces;

    my $found = 0;
    foreach my $tmpname (@names) {
        if ( $name eq $tmpname ) {
            $found = 1;
        }
    }

    if ( $found == 0 ) {
        $self->throw("Secondary index for $name doesn't exist\n");
    }

    my $fh = $self->open_secondary_index($name);

    syseof($fh);

    my $filesize = systell($fh);

    my $recsize = $self->{'_secondary_record_size'}{$name};

    #    print "Name " . $recsize . "\n";

    my $end = ( $filesize - $self->{'_start_pos'} ) / $recsize;

    #    print "End $end $filesize\n";
    my ( $newid, $primary_id, $pos ) =
      $self->find_entry( $fh, 0, $end, $id, $recsize );

    sysseek( $fh, $pos, 0 );

    #    print "Found new id $newid $primary_id\n";
    # We now need to shuffle up the index file to find the top secondary entry

    my $record = $newid;

    while ( $record =~ /^$newid/ && $pos >= 0 ) {

        $record = $self->read_record( $fh, $pos, $recsize );
        $pos = $pos - $recsize;

        #	print "Up record = $record:$newid\n";
    }

    $pos += $recsize;

    #    print "Top position is $pos\n";

    # Now we have to shuffle back down again to read all the secondary entries

    my $current_id = $newid;
    my %primary_id;

    $primary_id{$primary_id} = 1;

    while ( $current_id eq $newid ) {
        $record = $self->read_record( $fh, $pos, $recsize );

        # print "Record is :$record:\n";
        my ( $secid, $primary_id ) = split( /\t/, $record, 2 );
        $current_id = $secid;

        if ( $current_id eq $newid ) {
            $primary_id =~ s/ //g;

            #    print "Primary $primary_id\n";
            $primary_id{$primary_id} = 1;

            $pos = $pos + $recsize;

            #   print "Down record = $record\n";
        }
    }

    if ( !defined($newid) ) {
        return;
    }

    my @entry;

    foreach my $id ( keys %primary_id ) {
        push @entry, $self->get_Seq_by_id($id);
    }
    return wantarray ? @entry : $entry[0];

}

=head2 read_header

 Title   : read_header
 Usage   : $obj->read_header($fhl)
 Function: Reads the header from the db file
 Returns : width of a record
 Args    : filehandle


=cut

sub read_header {
    my ( $self, $fh ) = @_;

    my $record_width;

    sysread( $fh, $record_width, HEADER_SIZE );

    $self->{'_start_pos'} = HEADER_SIZE;
    $record_width =~ s/ //g;
    $record_width = $record_width * 1;

    return $record_width;
}

=head2 read_record

 Title   : read_record
 Usage   : $obj->read_record($fh,$pos,$len)
 Function: Reads a record from a filehandle
 Returns : String
 Args    : filehandle, offset, and length


=cut

sub read_record {
    my ( $self, $fh, $pos, $len ) = @_;

    sysseek( $fh, $pos, 0 );

    my $record;

    sysread( $fh, $record, $len );

    return $record;

}

=head2 get_all_primary_ids

 Title   : get_all_primary_ids
 Usage   : @ids = $seqdb->get_all_primary_ids()
 Function: gives an array of all the primary_ids of the
           sequence objects in the database.
 Returns : an array of strings
 Args    : none

=cut

sub get_all_primary_ids {
    my $self = shift;

    my $fh = $self->primary_index_filehandle;
    syseof($fh);
    my $filesize = systell($fh);
    my $recsize  = $self->record_size;
    my $end      = $filesize;

    my @ids;
    for ( my $pos = $self->{'_start_pos'} ; $pos < $end ; $pos += $recsize ) {
        my $record = $self->read_record( $fh, $pos, $recsize );
        my ($entryid) = split( /\t/, $record );
        push @ids, $entryid;
    }
    @ids;
}

=head2 find_entry

 Title   : find_entry
 Usage   : $obj->find_entry($fh,$start,$end,$id,$recsize)
 Function: Extract an entry based on the start,end,id and record size
 Returns : string
 Args    : filehandle, start, end, id, recordsize


=cut

sub find_entry {
    my ( $self, $fh, $start, $end, $id, $recsize ) = @_;

    my $mid = int( ( $end + 1 + $start ) / 2 );
    my $pos = ( $mid - 1 ) * $recsize + $self->{'_start_pos'};

    my ($record) = $self->read_record( $fh, $pos, $recsize );
    my ( $entryid, $rest ) = split( /\t/, $record, 2 );
    $rest =~ s/\s+$//;

    #    print "Mid $recsize $mid $pos:$entryid:$rest:$record\n";
    #    print "Entry :$id:$entryid:$rest\n";

    my ( $first, $second ) =
      $id le $entryid ? ( $id, $entryid ) : ( $entryid, $id );

    if ( $id eq $entryid ) {

        return ( $id, $rest, $pos - $recsize );

    }
    elsif ( $first eq $id ) {

        if ( $end - $start <= 1 ) {
            return;
        }
        my $end = $mid;

        #      print "Moving up $entryid $id\n";
        $self->find_entry( $fh, $start, $end, $id, $recsize );

    }
    elsif ( $second eq $id ) {

        #	print "Moving down $entryid $id\n";
        if ( $end - $start <= 1 ) {
            return;
        }

        $start = $mid;

        $self->find_entry( $fh, $start, $end, $id, $recsize );
    }

}

=head2 build_index

 Title   : build_index
 Usage   : $obj->build_index(@files)
 Function: Build the index based on a set of files
 Returns : count of the number of entries
 Args    : List of filenames


=cut

sub build_index {
    my ( $self, @files ) = @_;
    $self->write_flag
      or $self->throw('Cannot build index unless -write_flag is true');

    my $rootdir = $self->index_directory;

    if ( !defined($rootdir) ) {
        $self->throw("No index directory set - can't build indices");
    }

    if ( !-d $rootdir ) {
        $self->throw(
            "Index directory [$rootdir] is not a directory. Cant' build indices"
        );
    }

    my $dbpath = File::Spec->catfile( $rootdir, $self->dbname );
    if ( !-d $dbpath ) {
        warn "Creating directory $dbpath\n";
        mkdir $dbpath, 0777 or $self->throw("Couldn't create $dbpath: $!");
    }

    unless (@files) {
        $self->throw("Must enter an array of filenames to index");
    }

    foreach my $file (@files) {
        $file = File::Spec->rel2abs($file)
          unless File::Spec->file_name_is_absolute($file);
        unless ( -e $file ) {
            $self->throw("Can't index file [$file] as it doesn't exist");
        }
    }

    if ( my $filehash = $self->{_dbfile} ) {
        push @files, keys %$filehash;
    }

    my %seen;
    @files = grep { !$seen{$_}++ } @files;

    # Lets index
    $self->make_config_file( \@files );
    my $entries = 0;
    foreach my $file (@files) {
        $entries += $self->_index_file($file);
    }

    # update alphabet if necessary
    $self->make_config_file( \@files );

    # And finally write out the indices
    $self->write_primary_index;
    $self->write_secondary_indices;

    $entries;
}

=head2 _index_file

 Title   : _index_file
 Usage   : $obj->_index_file($newval)
 Function:
 Example :
 Returns : value of _index_file
 Args    : newvalue (optional)

=cut

sub _index_file {
    my ( $self, $file ) = @_;
    my $v = $self->verbose;
    open my $FILE, '<', $file or $self->throw("Could not read file '$file': $!");

    my $recstart = 0;
    my $fileid   = $self->get_fileid_by_filename($file);
    my $found    = 0;
    my $id;
    my $count = 0;

    my $primary       = $self->primary_pattern;
    my $start_pattern = $self->start_pattern;

    my $pos = 0;

    my $new_primary_entry;

    my $length;

    my $fh = $FILE;

    my $done = -1;

    my @secondary_names = $self->secondary_namespaces;
    my %secondary_id;
    my $last_one;

    # In Windows, text files have '\r\n' as line separator, but when reading in
    # text mode Perl will only show the '\n'. This means that for a line "ABC\r\n",
    # "length $_" will report 4 although the line is 5 bytes in length.
    # We assume that all lines have the same line separator and only read current line.
    my $init_pos   = tell($fh);
    my $curr_line  = <$fh>;
    my $pos_diff   = tell($fh) - $init_pos;
    my $correction = $pos_diff - length $curr_line;
    seek $fh, $init_pos, 0; # Rewind position to proceed to read the file

    while (<$fh>) {
        $last_one = $_;
        $self->{alphabet} ||= $self->guess_alphabet($_);
        if ( $_ =~ /$start_pattern/ ) {
            if ( $done == 0 ) {
                $id = $new_primary_entry;
                $self->{alphabet} ||= $self->guess_alphabet($_);

                my $tmplen = ( tell $fh ) - length($_) - $correction;

                $length = $tmplen - $pos;

                unless ( defined($id) ) {
                    $self->throw("No id defined for sequence");
                }
                unless ( defined($fileid) ) {
                    $self->throw("No fileid defined for file $file");
                }
                unless ( defined($pos) ) {
                    $self->throw( "No position defined for " . $id . "\n" );
                }
                unless ( defined($length) ) {
                    $self->throw( "No length defined for " . $id . "\n" );
                }
                $self->_add_id_position( $id, $pos, $fileid, $length,
                    \%secondary_id );

                $pos = $tmplen;

                if ( $count > 0 && $count % 1000 == 0 ) {
                    $self->debug("Indexed $count ids\n") if $v > 0;
                }

                $count++;
            }
            else {
                $done = 0;
            }
        }

        if ( $_ =~ /$primary/ ) {
            $new_primary_entry = $1;
        }

        my $secondary_patterns = $self->secondary_patterns;

        foreach my $sec (@secondary_names) {
            my $pattern = $secondary_patterns->{$sec};

            if ( $_ =~ /$pattern/ ) {
                $secondary_id{$sec} = $1;
            }
        }

    }

    # Remember to add in the last one

    $id = $new_primary_entry;

    # my $tmplen = (tell $fh) - length($last_one);
    my $tmplen = ( tell $fh );

    $length = $tmplen - $pos;

    if ( !defined($id) ) {
        $self->throw("No id defined for sequence");
    }
    if ( !defined($fileid) ) {
        $self->throw("No fileid defined for file $file");
    }
    if ( !defined($pos) ) {
        $self->throw( "No position defined for " . $id . "\n" );
    }
    if ( !defined($length) ) {
        $self->throw( "No length defined for " . $id . "\n" );
    }

    $self->_add_id_position( $id, $pos, $fileid, $length, \%secondary_id );
    $count++;

    close $FILE;
    $count;
}

=head2 write_primary_index

 Title   : write_primary_index
 Usage   : $obj->write_primary_index($newval)
 Function:
 Example :
 Returns : value of write_primary_index
 Args    : newvalue (optional)


=cut

sub write_primary_index {
    my ($self) = @_;

    my @ids = keys %{ $self->{_id} };

    @ids = sort { $a cmp $b } @ids;

    open my $INDEX, '>', $self->primary_index_file
      or $self->throw(
        "Could not write primary index file '" . $self->primary_index_file . "': $!" );

    my $recordlength =
      $self->{_maxidlength} +
      $self->{_maxfileidlength} +
      $self->{_maxposlength} +
      $self->{_maxlengthlength} + 3;

    print $INDEX sprintf( "%04d", $recordlength );

    foreach my $id (@ids) {

        if ( !defined( $self->{_id}{$id}{_fileid} ) ) {
            $self->throw("No fileid for $id\n");
        }
        if ( !defined( $self->{_id}{$id}{_pos} ) ) {
            $self->throw("No position for $id\n");
        }
        if ( !defined( $self->{_id}{$id}{_length} ) ) {
            $self->throw("No length for $id");
        }

        my $record =
            $id . "\t"
          . $self->{_id}{$id}{_fileid} . "\t"
          . $self->{_id}{$id}{_pos} . "\t"
          . $self->{_id}{$id}{_length};

        print $INDEX sprintf( "%-${recordlength}s", $record );

    }
}

=head2 write_secondary_indices

 Title   : write_secondary_indices
 Usage   : $obj->write_secondary_indices($newval)
 Function:
 Example :
 Returns : value of write_secondary_indices
 Args    : newvalue (optional)


=cut

sub write_secondary_indices {
    my ($self) = @_;

    # These are the different
    my @names = keys( %{ $self->{_secondary_id} } );

    foreach my $name (@names) {

        my @seconds = keys %{ $self->{_secondary_id}{$name} };

        # First we need to loop over to get the longest record.
        my $length = 0;

        foreach my $second (@seconds) {
            my $tmplen = length($second) + 1;
            my @prims  = keys %{ $self->{_secondary_id}{$name}{$second} };

            foreach my $prim (@prims) {
                my $recordlen = $tmplen + length($prim);

                if ( $recordlen > $length ) {
                    $length = $recordlen;
                }
            }
        }

        # Now we can print the index

        my $fh = $self->new_secondary_filehandle($name);

        print $fh sprintf( "%04d", $length );
        @seconds = sort @seconds;

        foreach my $second (@seconds) {

            my @prims = keys %{ $self->{_secondary_id}{$name}{$second} };
            my $tmp   = $second;

            foreach my $prim (@prims) {
                my $record = $tmp . "\t" . $prim;
                if ( length($record) > $length ) {
                    $self->throw(
"Something has gone horribly wrong - length of record is more than we thought [$length]\n"
                    );
                }
                else {
                    print $fh sprintf( "%-${length}s", $record );
                }
            }
        }

        close($fh);
    }
}

=head2 new_secondary_filehandle

 Title   : new_secondary_filehandle
 Usage   : $obj->new_secondary_filehandle($newval)
 Function:
 Example :
 Returns : value of new_secondary_filehandle
 Args    : newvalue (optional)


=cut

sub new_secondary_filehandle {
    my ( $self, $name ) = @_;

    my $indexdir = $self->_config_path;

    my $secindex = File::Spec->catfile( $indexdir, "id_$name.index" );

    open my $fh, '>', $secindex or $self->throw("Could not write file '$secindex': $!");
    return $fh;
}

=head2 open_secondary_index

 Title   : open_secondary_index
 Usage   : $obj->open_secondary_index($newval)
 Function:
 Example :
 Returns : value of open_secondary_index
 Args    : newvalue (optional)


=cut

sub open_secondary_index {
    my ( $self, $name ) = @_;

    if ( !defined( $self->{_secondary_filehandle}{$name} ) ) {

        my $indexdir = $self->_config_path;
        my $secindex = $indexdir . "/id_$name.index";

        if ( !-e $secindex ) {
            $self->throw("Index is not present for namespace [$name]\n");
        }

        open my $newfh, '<', $secindex or $self->throw("Could not read file '$secindex': $!");
        my $reclen = $self->read_header($newfh);

        $self->{_secondary_filehandle}{$name}  = $newfh;
        $self->{_secondary_record_size}{$name} = $reclen;
    }

    return $self->{_secondary_filehandle}{$name};

}

=head2 _add_id_position

 Title   : _add_id_position
 Usage   : $obj->_add_id_position($newval)
 Function:
 Example :
 Returns : value of _add_id_position
 Args    : newvalue (optional)


=cut

sub _add_id_position {
    my ( $self, $id, $pos, $fileid, $length, $secondary_id ) = @_;

    if ( !defined($id) ) {
        $self->throw("No id defined. Can't add id position");
    }
    if ( !defined($pos) ) {
        $self->throw("No position defined. Can't add id position");
    }
    if ( !defined($fileid) ) {
        $self->throw("No fileid defined. Can't add id position");
    }
    if ( !defined($length) || $length <= 0 ) {
        $self->throw(
            "No length defined or <= 0 [$length]. Can't add id position");
    }

    $self->{_id}{$id}{_pos}    = $pos;
    $self->{_id}{$id}{_length} = $length;
    $self->{_id}{$id}{_fileid} = $fileid;

    # Now the secondary ids

    foreach my $sec ( keys(%$secondary_id) ) {
        my $value = $secondary_id->{$sec};
        $self->{_secondary_id}{$sec}{$value}{$id} = 1;
    }

    $self->{_maxidlength} = length($id)
      if !exists $self->{_maxidlength}
          or length($id) >= $self->{_maxidlength};

    $self->{_maxfileidlength} = length($fileid)
      if !exists $self->{_maxfileidlength}
          or length($fileid) >= $self->{_maxfileidlength};

    $self->{_maxposlength} = length($pos)
      if !exists $self->{_maxposlength}
          or length($pos) >= $self->{_maxposlength};

    $self->{_maxlengthlength} = length($length)
      if !exists $self->{_maxlengthlength}
          or length($length) >= $self->{_maxlengthlength};
}

=head2 make_config_file

 Title   : make_config_file
 Usage   : $obj->make_config_file($newval)
 Function:
 Example :
 Returns : value of make_config_file
 Args    : newvalue (optional)

=cut

sub make_config_file {
    my ( $self, $files ) = @_;

    my @files = @$files;

    my $configfile = $self->_config_file;

    open my $CON, '>', $configfile
      or $self->throw("Could not write config file '$configfile': $!");

    # First line must be the type of index - in this case flat
    print $CON "index\tflat/1\n";

    # Now the fileids

    my $count = 0;

    foreach my $file (@files) {

        my $size = -s $file;

        print $CON "fileid_$count\t$file\t$size\n";

        $self->{_file}{$count}  = $file;
        $self->{_dbfile}{$file} = $count;
        $self->{_size}{$count}  = $size;
        $count++;
    }

    # Now the namespaces

    print $CON "primary_namespace\t" . $self->primary_namespace . "\n";

    # Needs fixing for the secondary stuff

    my $second_patterns = $self->secondary_patterns;

    my @second = keys %$second_patterns;

    if ( (@second) ) {
        print $CON "secondary_namespaces";

        foreach my $second (@second) {
            print $CON "\t$second";
        }
        print $CON "\n";
    }

    # Now the config format

    unless ( defined( $self->format ) ) {
        $self->throw(
            "Format does not exist in module - can't write config file");
    }
    else {
        my $format   = $self->format;
        my $alphabet = $self->alphabet;
        my $alpha    = $alphabet ? "/$alphabet" : '';
        print $CON "format\t" . "$format\n";
    }
    close $CON;
}

=head2 read_config_file

 Title   : read_config_file
 Usage   : $obj->read_config_file($newval)
 Function:
 Example :
 Returns : value of read_config_file
 Args    : newvalue (optional)

=cut

sub read_config_file {
    my ($self) = @_;
    my $configfile = $self->_config_file;
    return unless -e $configfile;

    open my $CON, '<', $configfile
      or $self->throw("Could not read config file '$configfile': $!");

    # First line must be type
    my $line = <$CON>;
    chomp($line);
    my $version;

    # This is hard coded as we only index flatfiles here
    if ( $line =~ m{index\tflat/(\d+)} ) {
        $version = $1;
    }
    else {
        $self->throw(
"First line not compatible with flat file index.  Should be something like\n\nindex\tflat/1"
        );
    }

    $self->index_type("flat");
    $self->index_version($version);

    while (<$CON>) {
        chomp;

        # Look for fileid lines
        if ( $_ =~ /^fileid_(\d+)\t(.+)\t(\d+)/ ) {
            my $fileid   = $1;
            my $filename = $2;
            my $filesize = $3;

            if ( !-e $filename ) {
                $self->throw("File [$filename] does not exist!");
            }
            if ( -s $filename != $filesize ) {
                $self->throw(
"Flatfile size for $filename differs from what the index thinks it is. Real size ["
                      . ( -s $filename )
                      . "] Index thinks it is ["
                      . $filesize
                      . "]" );
            }

            $self->{_file}{$fileid}     = $filename;
            $self->{_dbfile}{$filename} = $fileid;
            $self->{_size}{$fileid}     = $filesize;
        }

        # Look for namespace lines
        if (/(.*)_namespaces?\t(.+)/) {
            if ( $1 eq "primary" ) {
                $self->primary_namespace($2);
            }
            elsif ( $1 eq "secondary" ) {
                $self->secondary_namespaces( split "\t", $2 );
            }
            else {
                $self->throw("Unknown namespace name in config file [$1");
            }
        }

        # Look for format lines
        if ( $_ =~ /format\t(\S+)/ ) {

            # Check the format here?
            my $format = $1;

            # handle LSID format
            if ( $format =~ /^URN:LSID:open-bio\.org:(\w+)(?:\/(\w+))?/ ) {
                $self->format($1);
                $self->alphabet($2);
            }
            else {    # compatibility with older versions
                $self->format($1);
            }
        }
    }

    close($CON);

    # Now check we have all that we need

    my @fileid_keys = keys( %{ $self->{_file} } );

    if ( !(@fileid_keys) ) {
        $self->throw(
"No flatfile fileid files in config - check the index has been made correctly"
        );
    }

    if ( !defined( $self->primary_namespace ) ) {
        $self->throw("No primary namespace exists");
    }

    if ( !-e $self->primary_index_file ) {
        $self->throw( "Primary index file ["
              . $self->primary_index_file
              . "] doesn't exist" );
    }

    1;
}

=head2 get_fileid_by_filename

 Title   : get_fileid_by_filename
 Usage   : $obj->get_fileid_by_filename($newval)
 Function:
 Example :
 Returns : value of get_fileid_by_filename
 Args    : newvalue (optional)

=cut

sub get_fileid_by_filename {
    my ( $self, $file ) = @_;

    if ( !defined( $self->{_dbfile} ) ) {
        $self->throw(
            "No file to fileid mapping present.  Has the fileid file been read?"
        );
    }

    return $self->{_dbfile}{$file};
}

=head2 get_filehandle_by_fileid

 Title   : get_filehandle_by_fileid
 Usage   : $obj->get_filehandle_by_fileid($newval)
 Function:
 Example :
 Returns : value of get_filehandle_by_fileid
 Args    : newvalue (optional)

=cut

sub get_filehandle_by_fileid {
    my ( $self, $fileid ) = @_;

    if ( !defined( $self->{_file}{$fileid} ) ) {
        $self->throw("ERROR: undefined fileid in index [$fileid]");
    }

    open my $fh, '<', $self->{_file}{$fileid} or $self->throw("Could not read file '$self->{_file}{$fileid}': $!");
    return $fh;
}

=head2 primary_index_file

 Title   : primary_index_file
 Usage   : $obj->primary_index_file($newval)
 Function:
 Example :
 Returns : value of primary_index_file
 Args    : newvalue (optional)


=cut

sub primary_index_file {
    my ($self) = @_;

    return File::Spec->catfile( $self->_config_path,
        "key_" . $self->primary_namespace . ".key" );
}

=head2 primary_index_filehandle

 Title   : primary_index_filehandle
 Usage   : $obj->primary_index_filehandle($newval)
 Function:
 Example :
 Returns : value of primary_index_filehandle
 Args    : newvalue (optional)


=cut

sub primary_index_filehandle {
    my ($self) = @_;

    unless ( defined( $self->{'_primary_index_handle'} ) ) {
        my $primary_file = $self->primary_index_file;
        open $self->{'_primary_index_handle'}, '<', $primary_file
          or self->throw("Could not read file '$primary_file': $!\n");
    }
    return $self->{'_primary_index_handle'};
}

=head2 format

 Title   : format
 Usage   : $obj->format($newval)
 Function:
 Example :
 Returns : value of format
 Args    : newvalue (optional)


=cut

sub format {
    my ( $obj, $value ) = @_;
    if ( defined $value ) {
        $obj->{'format'} = $value;
    }
    return $obj->{'format'};

}

sub alphabet {
    my ( $obj, $value ) = @_;
    if ( defined $value ) {
        $obj->{alphabet} = $value;
    }
    return $obj->{alphabet};
}

=head2 write_flag

 Title   : write_flag
 Usage   : $obj->write_flag($newval)
 Function:
 Example :
 Returns : value of write_flag
 Args    : newvalue (optional)


=cut

sub write_flag {
    my ( $obj, $value ) = @_;
    if ( defined $value ) {
        $obj->{'write_flag'} = $value;
    }
    return $obj->{'write_flag'};

}

=head2 dbname

 Title   : dbname
 Usage   : $obj->dbname($newval)
 Function: get/set database name
 Example :
 Returns : value of dbname
 Args    : newvalue (optional)

=cut

sub dbname {
    my $self = shift;
    my $d    = $self->{flat_dbname};
    $self->{flat_dbname} = shift if @_;
    $d;
}

=head2 index_directory

 Title   : index_directory
 Usage   : $obj->index_directory($newval)
 Function:
 Example :
 Returns : value of index_directory
 Args    : newvalue (optional)


=cut

sub index_directory {
    my ( $self, $arg ) = @_;

    if ( defined($arg) ) {
        if ( $arg !~ m{/$} ) {
            $arg .= "/";
        }
        $self->{_index_directory} = $arg;
    }
    return $self->{_index_directory};

}

sub _config_path {
    my $self   = shift;
    my $root   = $self->index_directory;
    my $dbname = $self->dbname;
    File::Spec->catfile( $root, $dbname );
}

sub _config_file {
    my $self = shift;
    my $path = $self->_config_path;
    File::Spec->catfile( $path, CONFIG_FILE_NAME );
}

=head2 record_size

 Title   : record_size
 Usage   : $obj->record_size($newval)
 Function:
 Example :
 Returns : value of record_size
 Args    : newvalue (optional)


=cut

sub record_size {
    my $self = shift;
    $self->{_record_size} = shift if @_;
    return $self->{_record_size};
}

=head2 primary_namespace

 Title   : primary_namespace
 Usage   : $obj->primary_namespace($newval)
 Function:
 Example :
 Returns : value of primary_namespace
 Args    : newvalue (optional)

=cut

sub primary_namespace {
    my $self = shift;
    $self->{_primary_namespace} = shift if @_;
    return $self->{_primary_namespace};
}

=head2 index_type

 Title   : index_type
 Usage   : $obj->index_type($newval)
 Function:
 Example :
 Returns : value of index_type
 Args    : newvalue (optional)


=cut

sub index_type {
    my $self = shift;
    $self->{_index_type} = shift if @_;
    return $self->{_index_type};
}

=head2 index_version

 Title   : index_version
 Usage   : $obj->index_version($newval)
 Function:
 Example :
 Returns : value of index_version
 Args    : newvalue (optional)


=cut

sub index_version {
    my $self = shift;
    $self->{_index_version} = shift if @_;
    return $self->{_index_version};
}

=head2 primary_pattern

 Title   : primary_pattern
 Usage   : $obj->primary_pattern($newval)
 Function:
 Example :
 Returns : value of primary_pattern
 Args    : newvalue (optional)


=cut

sub primary_pattern {
    my $obj = shift;
    $obj->{'primary_pattern'} = shift if @_;
    return $obj->{'primary_pattern'};
}

=head2 start_pattern

 Title   : start_pattern
 Usage   : $obj->start_pattern($newval)
 Function:
 Example :
 Returns : value of start_pattern
 Args    : newvalue (optional)


=cut

sub start_pattern {
    my $obj = shift;
    $obj->{'start_pattern'} = shift if @_;
    return $obj->{'start_pattern'};
}

=head2 secondary_patterns

 Title   : secondary_patterns
 Usage   : $obj->secondary_patterns($newval)
 Function:
 Example :
 Returns : value of secondary_patterns
 Args    : newvalue (optional)


=cut

sub secondary_patterns {
    my ( $obj, $value ) = @_;
    if ( defined $value ) {
        $obj->{'secondary_patterns'} = $value;

        my @names = keys %$value;

        foreach my $name (@names) {
            $obj->secondary_namespaces($name);
        }
    }
    return $obj->{'secondary_patterns'};

}

=head2 secondary_namespaces

 Title   : secondary_namespaces
 Usage   : $obj->secondary_namespaces($newval)
 Function:
 Example :
 Returns : value of secondary_namespaces
 Args    : newvalue (optional)


=cut

sub secondary_namespaces {
    my ( $obj, @values ) = @_;

    if (@values) {
        push( @{ $obj->{'secondary_namespaces'} }, @values );
    }
    return @{ $obj->{'secondary_namespaces'} || [] };
}

## These are indexing routines to index commonly used format - fasta
## swissprot and embl

sub new_SWISSPROT_index {
    my ( $self, $index_dir, @files ) = @_;

    my %secondary_patterns;

    my $start_pattern   = "^ID   (\\S+)";
    my $primary_pattern = "^AC   (\\S+)\\;";

    $secondary_patterns{"ID"} = $start_pattern;

    my $index = Bio::DB::Flat::BinarySearch->new(
        -index_dir          => $index_dir,
        -format             => 'swissprot',
        -primary_pattern    => $primary_pattern,
        -primary_namespace  => "ACC",
        -start_pattern      => $start_pattern,
        -secondary_patterns => \%secondary_patterns
    );

    $index->build_index(@files);
}

sub new_EMBL_index {
    my ( $self, $index_dir, @files ) = @_;

    my %secondary_patterns;

    my $start_pattern     = "^ID   (\\S+)";
    my $primary_pattern   = "^AC   (\\S+)\\;";
    my $primary_namespace = "ACC";

    $secondary_patterns{"ID"} = $start_pattern;

    my $index = Bio::DB::Flat::BinarySearch->new(
        -index_dir          => $index_dir,
        -format             => 'embl',
        -primary_pattern    => $primary_pattern,
        -primary_namespace  => "ACC",
        -start_pattern      => $start_pattern,
        -secondary_patterns => \%secondary_patterns
    );

    $index->build_index(@files);

    return $index;
}

sub new_FASTA_index {
    my ( $self, $index_dir, @files ) = @_;

    my %secondary_patterns;

    my $start_pattern     = "^>";
    my $primary_pattern   = "^>(\\S+)";
    my $primary_namespace = "ACC";

    $secondary_patterns{"ID"} = "^>\\S+ +(\\S+)";

    my $index = Bio::DB::Flat::BinarySearch->new(
        -index_dir          => $index_dir,
        -format             => 'fasta',
        -primary_pattern    => $primary_pattern,
        -primary_namespace  => "ACC",
        -start_pattern      => $start_pattern,
        -secondary_patterns => \%secondary_patterns
    );

    $index->build_index(@files);

    return $index;
}

# EVERYTHING THAT FOLLOWS THIS
# is an awful hack - in reality Michele's code needs to be rewritten
# to use Bio::SeqIO, but I have too little time to do this -- LS
sub guess_alphabet {
    my $self = shift;
    my $line = shift;

    my $format = $self->format;
    return 'protein' if $format eq 'swissprot';

    if ( $format eq 'genbank' ) {
        return unless $line =~ /^LOCUS/;
        return 'dna' if $line =~ /\s+\d+\s+bp/i;
        return 'protein';
    }

    if ( $format eq 'embl' ) {
        return unless $line =~ /^ID/;
        return 'dna' if $line =~ / DNA;/i;
        return 'rna' if $line =~ / RNA;/i;
        return 'protein';
    }

    return;
}

# return (namespace,primary_pattern,start_pattern,secondary_pattern)
sub _guess_patterns {
    my $self   = shift;
    my $format = shift;
    if ( $format =~ /swiss(prot)?/i ) {
        return ( 'ID', "^ID   (\\S+)", "^ID   (\\S+)",
            { ACC => "^AC   (\\S+);" } );
    }

    if ($format =~ /embl/i) {
        return ('ID',
            "^ID   (\\S+[^; ])",
            "^ID   (\\S+[^; ])",
            {
             ACC     => q/^AC   (\S+);/,
             VERSION => q/^SV\s+(\S+)/
            });
     }

    if ( $format =~ /genbank/i ) {
        return (
            'ID',
            q/^LOCUS\s+(\S+)/,
            q/^LOCUS/,
            {
                ACC     => q/^ACCESSION\s+(\S+)/,
                VERSION => q/^VERSION\s+(\S+)/
            }
        );
    }

    if ( $format =~ /fasta/i ) {
        return ( 'ACC', '^>(\S+)', '^>(\S+)', );
    }

    $self->throw("I can't handle format $format");

}

1;
