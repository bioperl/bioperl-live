
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
    my $inx = Bio::Index::Fasta->new(
        -filename => $Index_File_Name,
        -write_flag => 1);
    $inx->make_index(@ARGV);

    # Print out several sequences present in the index
    # in Fasta format
    use Bio::Index::Fasta;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::Fasta->new($Index_File_Name);
    my $out = Bio::SeqIO->new('-format' => 'Fasta','-fh' => \*STDOUT);

    foreach my $id (@ARGV) {
        my $seq = $inx->fetch($id); # Returns Bio::Seq object
	$out->write_seq($seq);
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

use vars qw($VERSION @ISA);
use strict;

use Bio::Index::AbstractSeq;
use Bio::Seq;

@ISA = qw(Bio::Index::AbstractSeq);

#
# Suggested fix by Michael G Schwern <schwern@pobox.com> to
# get around a clash with CPAN shell...
#

BEGIN { 
    $VERSION = 0.2;
}

sub _version {
    return $VERSION;
}


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
    my $self = shift;
    
    $self->SUPER::_initialize(@_);
}


=head2 _file_format

 Title   : _file_format
 Function: The file format for this package, which is needed
           by the SeqIO system when reading the sequence.
 Returns : 'Fasta'

=cut

sub _file_format {
    return 'Fasta';
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
        $i,    # Index-number of file being indexed
        ) = @_;
    
    my( $begin,     # Offset from start of file of the start
                    # of the last found record.
        @id_list,   # List of ids from last record
        );

    $begin = 0;

    my $id_parser = $self->id_parser;

    open FASTA, $file or $self->throw("Can't open file for read : $file");

    # Main indexing loop
    while (<FASTA>) {
        if (/^>/) {
            # $begin is the position of the first character after the '>'
            my $begin = tell(FASTA) - length( $_ ) + 1;

            foreach my $id (&$id_parser($_)) {
                $self->add_record($id, $i, $begin);
            }
        }
    }

    close FASTA;
    return 1;
}

=head2 id_parser

  Title   : id_parser
  Usage   : $index->id_parser( CODE )
  Function: Stores or returns the code used by record_id to
            parse the ID for record from a string.  Useful
            for (for instance) specifying a different
            parser for different flavours of FASTA file. 
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
    my( $self, $code ) = @_;
    
    if ($code) {
        $self->{'_id_parser'} = $code;
    }
    return $self->{'_id_parser'} || \&default_id_parser;
}



=head2 default_id_parser

  Title   : default_id_parser
  Usage   : $id = default_id_parser( $header )
  Function: The default Fasta ID parser for Fasta.pm
            Returns $1 from applying the regexp /^>\s*(\S+)/
            to $header.
  Returns : ID string
  Args    : a fasta header line string

=cut

sub default_id_parser {    
    if ($_[0] =~ /^>\s*(\S+)/) {
        return $1;
    } else {
        return;
    }
}

1;



