
#
# BioPerl module for Bio::Index::Abstract
#
# Cared for by James Gilbert <jgrg@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::Fasta - Interface for indexing (multiple) fasta files

=head1 SYNOPSIS

    # Complete code for making an index for several
    # fasta files
    use Bio::Index::Fasta;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::Fasta->new($Index_File_Name, 'WRITE');
    $inx->make_index(@ARGV);

    # Print out several sequences present in the index
    # in gcg format
    use Bio::Index::Fasta;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::Fasta->new($Index_File_Name);

    foreach my $id (@ARGV) {
        my $seq = $inx->fetch($id); # Returns Bio::Seq object
        print $seq->layout('GCG');
    }

    # or, alternatively

    my $seq = $inx->get_Seq_by_id($id); #identical to fetch   

=head1 DESCRIPTION

Inherits functions for managing dbm files from Bio::Index::Abstract.pm,
and provides the basic funtionallity for indexing fasta files, and
retrieving the sequence from them. 

Bio::Index::Fasta supports the Bio::DB::BioSeqI interface, meaning
it can be used a a Sequence database for other parts of bioperl

=head1 FEED_BACK

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

=head1 AUTHOR - James Gilbert

Email - jgrg@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Index::Fasta;

use vars qw($VERSION @ISA @EXPORT_OK);
use strict;

use Bio::Index::Abstract;
use Bio::Seq;

@ISA = qw(Bio::Index::Abstract Bio::DB::BioSeqI Exporter);
@EXPORT_OK = qw();

sub _type_stamp {
    return '__FASTA__'; # What kind of index are we?
}

sub _version {
    return 0.2;
}
$VERSION = _version();



=head2 _initialize

  Title   : _initialize
  Usage   : $index->_initialize
  Function: Calls $index->SUPER::_initialize(), and then adds
            the default id parser for fasta files.
  Example : 
  Returns : 
  Args    : 

=cut

sub _initialize {
    my($self, $index_file, $write_flag) = @_;
    
    $self->SUPER::_initialize($index_file, $write_flag);
    $self->id_parser( \&default_id_parser );
}


=head2 _index_file

  Title   : _index_file
  Usage   : $index->_index_file( $file_name, $i )
  Function: Specialist function to index FASTA format files.
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
        );

    $begin = 0;
    $end   = 0;

    open FASTA, $file or $self->throw("Can't open file for read : $file");

    # Main indexing loop
    while (<FASTA>) {
        if (/^>/) {
            my $new_begin = tell(FASTA) - length( $_ );
            $end = $new_begin - 1;

            $self->add_record($id, $i, $begin, $end) if $id;

            $begin = $new_begin;
            ($id) = $self->record_id( $_ );
        }
    }
    # Don't forget to add the last record
    $end = tell(FASTA);
    $self->add_record($id, $i, $begin, $end) if $id;

    close FASTA;
    return 1;
}


# Should there be a prototype for this method in Index::Abstract.pm?
=head2 record_id

  Title   : record_id
  Usage   : $index->record_id( STRING );
  Function: Parses the ID for an entry from the string
            supplied, using the code in $index->{'_id_parser'}
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


=head2 id_parser

  Title   : id_parser
  Usage   : $index->id_parser( CODE )
  Function: Stores or returns the code used by record_id
            to parse the ID for record from a string.  Useful
            for (for instance) specifying a different parser
            for different flavours of FASTA file.
  Example : $index->id_parser( \&my_id_parser )
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



=head2 default_id_parser

  Title   : default_id_parser
  Usage   : $id = default_id_parser( $header )
  Function: The default Fasta ID parser for Fasta.pm
            Returns $1 from applying the regexp /^>\s*(\S+)/
            to $header.
  Example : 
  Returns : ID string
  Args    : a fasta header line string

=cut

sub default_id_parser {
    my $line = shift;
    $line =~ /^>\s*(\S+)/;
    return $1;
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
    
    my $db = $self->db();
    if (my $rec = $db->{ $id }) {
        my( @record );
        
        my ($file, $begin, $end) = $self->unpack_record( $rec );
        
        # Get the (possibly cached) filehandle
        my $fh = $self->_file_handle( $file );

        # Accumulate lines in @record until beyond end
        seek($fh, $begin, 0);
        while (defined(my $line = <$fh>)) {
            push(@record, $line);
            last if tell($fh) > $end;
        }
        
        $self->throw("Can't fetch sequence for record : $id")
            unless @record;
        
        # Parse record
        my $firstLine = shift @record;
        my ($name, $desc) = $firstLine =~ /^>\s*(\S+)\s*(.*?)\s*$/;
        chomp( @record );
        
        # Return a shiny Bio::Seq object
        return Bio::Seq->new( -ID   => $name,
                              -DESC => $desc,
                              -SEQ  => uc(join('', @record)) );
    } else {
	$self->throw("Unable to find a record for $id in Fasta index");
	return;
    }
}

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id()
 Function: retrieves a sequence object, identically to
           ->fetch, but here behaving as a Bio::DB::BioSeqI
 Returns : new Bio::Seq object
 Args    : string represents the id


=cut

sub get_Seq_by_id{
   my ($self,$id) = @_;

   return $self->fetch($id);
}


1;



