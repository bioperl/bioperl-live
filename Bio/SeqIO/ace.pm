#
# BioPerl module for Bio::SeqIO::ace
#
# Cared for by James Gilbert <jgrg@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::ace - ace sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::SeqIO class.

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and
from ace file format.  It only parses a DNA or
Peptide objects contained in the ace file,
producing PrimarySeq objects from them.  All
other objects in the files will be ignored.  It
doesn't attempt to parse any annotation attatched
to the containing Sequence or Protein objects,
which would probably be impossible, since
everyone's ACeDB schema can be different.

It won't parse ace files containing Timestamps
correctly either.  This can easily be added if
considered necessary.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  bioperl-guts-l@bioperl.org             - Technically-oriented discussion
  http://www.bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS - James Gilbert

Email: jgrg@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::ace;
use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::Object

use Bio::SeqIO;

@ISA = qw(Bio::SeqIO);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  return unless my $make = $self->SUPER::_initialize(@args);
}

=head2 next_prmary_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::PrimarySeq object
 Args    : NONE

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    : NONE

=cut

{
    my %bio_mol_type = (
        'dna'       => 'dna',
        'peptide'   => 'protein',
    );
    
    sub next_primary_seq {
        my( $self, $as_next_seq ) = @_;
        local $/ = "";  # Split input on blank lines

        my $fh = $self->_filehandle;
        my( $type, $id );
        while (<$fh>) {
            if (($type, $id) = /^(DNA|Peptide)[\s:]+(.+?)\s*\n/si) {
                s/^.+$//m;  # Remove first line
                s/\s+//g;   # Remove whitespace
                last;
            }
        }
        # Return if there weren't any DNA or peptide objects
        return unless $type;
        
        # Choose the molecule type
        my $mol_type = $bio_mol_type{lc $type}
            or $self->throw("Can't get Bio::Seq molecule type for '$type'");

        # Remove quotes from $id
        $id =~ s/^"|"$//g;
        
        # Un-escape forward slashes, double quotes, percent signs,
        # semi-colons, tabs, and backslashes (if you're mad enough
        # to have any of these as part of object names in your acedb
        # database).
        $id =~ s/\\([\/"%;\t\\])/$1/g;

        if ($as_next_seq) {
            # Called as next_seq(), so give back a Bio::Seq
            return Bio::Seq->new(
                -seq        => $_,
	        -primary_id => $id,
	        -display_id => $id,
                -moltype    => $mol_type,
	        );
        } else {
            return Bio::PrimarySeq->new(
                -seq        => $_,
	        -primary_id => $id,
	        -display_id => $id,
                -moltype    => $mol_type,
	        );
        }
    }
    
}

sub next_seq {
    return next_primary_seq($_[0], 1);
}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq(@seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object(s)


=cut

sub write_seq {
    my ($self, @seq) = @_;
    
    foreach my $seq (@seq) {

        my $mol_type = $seq->moltype;
        my $id = $seq->display_id;
        
        # Escape special charachers in id
        $id =~ s/([\/"%;\t\\])/\\$1/g;
        
        # Print header for DNA or Protein object
        if ($mol_type eq 'dna') {
            $self->_print( 
                qq{\nSequence : "$id"\nDNA "$id"\n},
                qq{\nDNA : "$id"\n},
            );
        }
        elsif ($mol_type eq 'protein') {
            $self->_print(
                qq{\nProtein : "$id"\nPeptide "$id"\n},
                qq{\nPeptide : "$id"\n},
            );
        }
        else {
            $self->throw("Don't know how to produce ACeDB output for '$mol_type'");
        }

        # Print the sequence
        my $str = $seq->seq;
        my( $formatted_seq );
        while ($str =~ /(.{1,60})/g) {
            $formatted_seq .= "$1\n";
        }
        $self->_print($formatted_seq, "\n");
    }
    return 1;
}

1;
