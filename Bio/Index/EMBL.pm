
#
# BioPerl module for Bio::Index::Abstract
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::EMBL - Interface for indexing (multiple) EMBL/Swissprot
.dat files (ie flat file embl/swissprot format).

=head1 SYNOPSIS

    # Complete code for making an index for several
    # EMBL files
    use Bio::Index::EMBL;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::EMBL->new($Index_File_Name, 'WRITE');
    $inx->make_index(@ARGV);

    # Print out several sequences present in the index
    # in gcg format
    use Bio::Index::EMBL;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::EMBL->new($Index_File_Name);

    foreach my $id (@ARGV) {
        my $seq = $inx->fetch($id); # Returns Bio::Seq object
        print $seq->layout('GCG');
    }

    # alternatively

    my $seq1 = $inx->get_Seq_by_id($id);
    my $seq2 = $inx->get_Seq_by_acc($acc);   

=head1 DESCRIPTION

Inherits functions for managing dbm files from Bio::Index::Abstract.pm,
and provides the basic funtionallity for indexing EMBL files, and
retrieving the sequence from them. Heavily snaffled from James Gilbert's
Fasta system.

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

=head1 AUTHOR - Ewan Birney

Email - birney@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Index::EMBL;

use vars qw($VERSION @ISA @EXPORT_OK);
use strict;

use Bio::Index::Abstract;
use Bio::Seq;

@ISA = qw(Bio::Index::Abstract Bio::DB::BioSeqI Exporter);
@EXPORT_OK = qw();

sub _type_stamp {
    return '__EMBL_FLAT__'; # What kind of index are we?
}

sub _version {
    return 0.1;
}
$VERSION = _version();



=head2 _initialize

  Title   : _initialize
  Usage   : $index->_initialize
  Function: Calls $index->SUPER::_initialize(), and then adds
            the default id parser for EMBL files.
  Example : 
  Returns : 
  Args    : 

=cut

sub _initialize {
    my($self, $index_file, $write_flag) = @_;
    
    $self->SUPER::_initialize($index_file, $write_flag);
}


=head2 _index_file

  Title   : _index_file
  Usage   : $index->_index_file( $file_name, $i )
  Function: Specialist function to index EMBL format files.
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
        );

    $begin = 0;
    $end   = 0;

    open EMBL, $file or $self->throw("Can't open file for read : $file");

    # Main indexing loop
    $id = $acc = undef;
    while (<EMBL>) {
	if( /^\/\// ) {
	    $end = tell(EMBL);
	    if( ! defined $id ) {
		$self->throw("Got to a end of entry line for an EMBL flat file with no parsed ID. Considering this a problem!");
		next;
	    }
	    if( ! defined $acc ) {
		$self->warn("For id [$id] in embl flat file, got no accession number. Storing id index anyway");
	    }

	    $self->add_record($id, $i, $begin, $end);

	    if( $acc ne $id ) {
		$self->add_record($acc, $i, $begin, $end);
	    }
	} elsif (/^ID\s+(\S+)/) {
	    $id = $1;
	    # not sure if I like this. Assummes tell is in bytes.
	    # we could tell before each line and save it.
            $begin = tell(EMBL) - length( $_ ); 
	    
	} elsif (/^AC\s+(\S+);/) { # ignore ? if there.
	    $acc =$1;
	} else {
	    # do nothing
	}
    }

    close EMBL;
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
    my ($desc,$acc,$out);
    my $db = $self->db();

    if (my $rec = $db->{ $id }) {
        my( @record );
        
        my ($file, $begin, $end) = $self->unpack_record( $rec );
        
        # Get the (possibly cached) filehandle
        my $fh = $self->_file_handle( $file );

        # move to start
        seek($fh, $begin, 0);


	#get id from file, and then loop to SQ line
        while (<$fh>) {
	    #print STDERR "Got $_";
	    /^SQ\s/ && last;
	    /^ID\s+(\S+)/ && do { $id = $1; };
	    /^DE\s+(.*?)\s+$/ && do { $desc .= $1; }; 
	    /^AC\s+(\S+?);/ && do { $acc .= $1; }; 
	    # accession numbers???
        }

        while (<$fh>) {
	    /^\/\// && last;
	    #print STDERR "Got $_";
	    s/[\W0-9]//g;
            push(@record, $_);
            last if tell($fh) > $end;
	}

        $self->throw("Can't fetch sequence for record : $id")
            unless @record;
        
        # Return a shiny Bio::Seq object
        $out = Bio::Seq->new( -ID   => $id,
                              -DESC => $desc,
			      -SEQ  => uc(join('', @record)) );
	$out->names()->{'acc'} = $acc if $acc;
	return $out;
    } else {
	$self->throw("Unable to find a record for $id in EMBL flat file index");
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

=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   : $seq = $db->get_Seq_by_acc()
 Function: retrieves a sequence object, identically to
           ->fetch, but here behaving as a Bio::DB::BioSeqI
 Returns : new Bio::Seq object
 Args    : string represents the accession number


=cut

sub get_Seq_by_acc {
   my ($self,$id) = @_;

   return $self->fetch($id);
}


1;



