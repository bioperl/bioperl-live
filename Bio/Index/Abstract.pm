
#
# BioPerl module for Bio::Index::Abstract
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::Abstract - Abstract interface for indexing a flat file

=head1 SYNOPSIS

You should not be using this module directly

=head1 DESCRIPTION

This object provides the basic mechanism to associate positions
in files with names. The position and filenames are stored in DBM
which can then be accessed later on. It is the equivalent of flat
file indexing (eg, SRS or efetch).

This object is the guts to the mechanism, which will be used by the
specific objects inherieting from it.


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

=head1 AUTHOR - Ewan Birney

Email - birney@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with an "_" (underscore).

=cut


# Let the code begin...


package Bio::Index::Abstract;
use vars qw($TYPE_AND_VERSION_KEY $AUTOLOAD @ISA @EXPORT_OK);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;



@ISA = qw(Bio::Root::Object Exporter);
@EXPORT_OK = qw();

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

# set stuff in self from @args
 return $make; # success - we hope!
}


=head2 record_id

  Title   : record_id
  Usage   : $obj->record_id( STRING );
  Function: Parses the ID for an entry from the string
            supplied, using the code in $obj->{'_id_parser'}
  Example : 
  Returns : scalar or exception
  Args    : STRING


=cut

sub record_id {
    my ($self, $line) = @_;

    if (my $id = $self->{'_id_parser'}->( $line )) {
        return $id;
    } else {
        $self->throw("Can't parse ID from line : $line");
    }
}


=head2 _type_and_version

  Title   : _type_and_version
  Usage   : Should be called by _initalize in the sub class
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
    my $type    = $self->_type_stamp();
    
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

sub _file_record {
    
}


=head2 id_parser

  Title   : id_parser
  Usage   : $obj->id_parser( CODE )
  Function: Stores or returns the code used by record_id
            to parse the ID for record from a string.  Useful
            for (for instance) specifying a different parser
            for different flavours of FASTA file.
  Example : $obj->id_parser( \&my_id_parser )
  Returns : ref to CODE if called without arguments
  Args    : CODE

=cut


sub id_parser {
    my( $self, $code ) = @_;
    
    if ($code) {
        $self->{'_id_parser'} = $code;
    } else {
        return $self->{'_id_parser'};
    }
}


=head2 _file_count

  Title   : _file_count
  Usage   : $obj->_file_count( INT )
  Function: Used by the index building sub in a sub class to
            track the number of files indexed.  Sets or gets
            the number of files indexed if called with or
            without an argument.
  Example : 
  Returns : 
  Args    : INT

=cut

# Sets or stores the number of files indexed
sub _file_count {
    my $self = shift;
    if (@_) {
        $self->db->{'__FILE_COUNT'} = shift;
    } else {
        return $self->db->{'__FILE_COUNT'};
    }    
}


=head2 

  Title   : 
  Usage   : $obj->_check_file_sizes()
  Function: Verifies that the files listed in the database
            are the same size as when the database was built,
            or throws an exception.  Should be called by the
            _initialize function in the sub class.
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
            $self->throw("file [ $file ] has changed size $stored_size -> $size");
        }
    }
    return 1;
}


=head2 add_record

  Title   : add_record
  Usage   : $obj->add_record( $id, @stuff );
  Function: Calls pack_record on @stuff, and adds the result
            of pack_record to the index database under key $id.
  Example : $obj->add_record( $id, $fileNumber, $begin, $end )
  Returns : TRUE on success or FALSE on failure
  Args    : ID LIST

=cut

sub add_record {
    my( $self, $id, @rec ) = @_;
    $self->db->{$id} = $self->pack_record( @rec );
}


=head2 pack_record

  Title   : pack_record
  Usage   : $packed_string = $obj->pack_record( LIST )
  Function: Packs an array of scalars into a single string
            joined by ASCII 034 (which is unlikely to be used
            in any of the strings), and returns it. 
  Example : $packed_string = $obj->pack_record( $fileNumber, $begin, $end )
  Returns : STRING or undef
  Args    : LIST

=cut

sub pack_record {
    my( $self, @args ) = @_;
    return join "\034", @args;
}

=head2 unpack_record

  Title   : unpack_record
  Usage   : $obj->unpack_record( STRING )
  Function: Splits the sting provided into a 3 element array,
            splitting on ASCII 034.
  Example : ( $fileNumber, $begin, $end ) = $obj->unpack_record( $self->db->{$id} )
  Returns : A 3 element ARRAY
  Args    : STRING containing ASCII 034

=cut

sub unpack_record {
    my( $self, @args ) = @_;
    return split /\034/, $args[0], 3;
}

=head2 empty

  Title   : empty
  Usage   : $obj->empty()
  Function: Empties an existing index database.  You might
            want to use this in a script which rebuilds an
            index file.  (On the other hand, you might like
            to build your new index under another name, then
            rename it, in a production environment.)
  Example : 
  Returns : TRUE on success (probably)
  Args    : 

=cut

# Empty the database
sub empty {
    my $self = shift;
    %{$self->db} = ();
}

=head2 _version

  Title   : _version
  Usage   : $type = $obj->_version()
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


=head2 _type_stamp

  Title   : _type_stamp
  Usage   : $type = $obj->_type_stamp()
  Function: Returns a string which identifes the type of index
            module.  Used to permanently identify an index as
            belonging to a particular indexing module.  Must be
            provided in sub class.
  Example : 
  Returns : 
  Args    : none

=cut

sub _type_stamp {
    my $self = shift;
    
    $self->throw("In Bio::Index::Abstract, no _type_stamp method in sub class");
}

1;
