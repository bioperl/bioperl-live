#
#
# BioPerl module for Bio::Index::Fasta
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by James Gilbert <jgrg@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::Fasta - Interface for indexing (multiple) fasta files

=head1 SYNOPSIS

    # Make an index for one or more fasta files
    use Bio::Index::Fasta;
    use strict;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::Fasta->new(-filename => $Index_File_Name,
                                     -write_flag => 1);
    $inx->make_index(@ARGV);


    # Once the index is made it can accessed, either in the
    # same script or a different one
    use Bio::Index::Fasta;
    use strict;

    my $Index_File_Name = shift;
    my $inx = Bio::Index::Fasta->new(-filename => $Index_File_Name);
    my $out = Bio::SeqIO->new(-format => 'Fasta',
                              -fh => \*STDOUT);

    foreach my $id (@ARGV) {
        my $seq = $inx->fetch($id); # Returns Bio::Seq object
         $out->write_seq($seq);
    }

    # or, alternatively
    my $id;
    my $seq = $inx->get_Seq_by_id($id); # identical to fetch()

=head1 DESCRIPTION

Inherits functions for managing dbm files from Bio::Index::Abstract.pm,
and provides the basic funtionallity for indexing fasta files, and
retrieving the sequence from them. For best results 'use strict'.

Bio::Index::Fasta supports the Bio::DB::BioSeqI interface, meaning
it can be used as a Sequence database for other parts of bioperl

Additional example code is available in scripts/index/*PLS and in 
the Bioperl Tutorial (L<http://www.bioperl.org/wiki/Bptutorial.pl>)

Note that by default the key for the sequence will be the first continuous
string after the 'E<gt>' in the fasta header. If you want to use a specific
substring of the fasta header you must use the id_parser() method.

You can also set or customize the unique key used to retrieve by 
writing your own function and calling the id_parser() method.
For example:

   $inx->id_parser(\&get_id);
   # make the index
   $inx->make_index($file_name);

   # here is where the retrieval key is specified
   sub get_id {
      my $line = shift;
      $line =~ /^>.+gi\|(\d+)/;
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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - James Gilbert

Email - jgrg@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Index::Fasta;

use strict;
use warnings;

use Bio::Seq;

use base qw(Bio::Index::AbstractSeq);

#
# Suggested fix by Michael G Schwern <schwern@pobox.com> to
# get around a clash with CPAN shell...
#

sub _version {
    return 0.2;
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
    );

    my $id_parser = $self->id_parser;

    open my $FASTA, '<', $file or $self->throw("Can't open file for read : $file");

    my $offset = ( $^O =~ /mswin/i ) ? 1 : 0;

    # Main indexing loop
    while (<$FASTA>) {
        if (/^>/) {
            
            # the following was fixed to allow validation - cjfields
            
            # $begin is the position of the first character after the '>'
            $begin = tell($FASTA) - length( $_ ) - $offset;
            
            foreach my $id (&$id_parser($_)) {
                $self->add_record($id, $i, $begin);
            }
        }
    }
    close $FASTA;
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
