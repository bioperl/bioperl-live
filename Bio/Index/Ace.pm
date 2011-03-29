=head1 NAME

Bio::Index::Ace - Interface for indexing the contigs and singletons in
(multiple) ACE-format assembly files.

=head1 SYNOPSIS


=head1 DESCRIPTION

Inherits functions for managing dbm files from Bio::Index::Abstract.pm,
and provides the basic funtionallity for indexing fasta files, and
retrieving the sequence from them. For best results 'use strict'.

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
reponsive experts will be able look at the problem and quickly address
it. Please include a thorough description of the problem with code and
data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Robert Buels

Email - rmb32@cornell.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Index::Ace;

use strict;
use warnings;

use base qw(Bio::Index::Abstract);

#
# Suggested fix by Michael G Schwern <schwern@pobox.com> to
# get around a clash with CPAN shell...
#

sub _version {
    return 0.2;
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
    my( $self, $file, $file_index_num ) = @_;

    my $id_parser = $self->id_parser;
    open my $ace_fh, '<', $file or $self->throw("Can't open file for read : $file");
    my $tell = 0;
    while ( my $line = <$ace_fh> ) {
        if( $line =~ /^CO / ) {
            for my $id ( $id_parser->($line) ) {
                $self->add_record($id, $file_index_num, $tell);
            }
        }
        $tell = tell( $ace_fh );
    }

    return 1;
}

=head2 id_parser

  Title   : id_parser
  Usage   : $index->id_parser( CODE )
  Function: Stores or returns the code used by record_id to
            parse the ID for record from a string.  Useful
            for (for instance) specifying a different
            parser for different flavours of ACE file. 
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
    if( $_[0] =~ /^CO\s+(\S+)/ ) {
        return $1;
    } else {
        return;
    }
}

=head2 fetch

  Title   : fetch
  Usage   : $index->fetch( $id )
  Function: Returns a Bio::Assembly::Contig or Bio::Assemly::Single from the index
  Example : $seq = $index->fetch( 'dJ67B12' )
  Returns : Bio::Seq object
  Args    : ID

=cut

sub fetch {
	my( $self, $id ) = @_;
	my $db = $self->db();

	my $rec = $db->{ $id }
          or return;

        my ($file, $begin) = $self->unpack_record( $rec );

        # Get the (possibly cached) SeqIO object
        my $assembly_io = $self->_assembly_io( $file );
        my $fh = $seqio->_fh();

        # move to start of record
        # $begin-- if( $^O =~ /mswin/i); # workaround for Win DB_File bug
        seek($fh, $begin, 0);

        my $contig_or_singleton = $assembly_io-># TODO figure out how to get this

	return $seq;
}
sub _assembly_io {
    my ( $self, $file ) = @_;

    return $self->{_cached_assemblyio}{$file} ||= do {
        # TODO: make an assemblyio
    };
}

1;
