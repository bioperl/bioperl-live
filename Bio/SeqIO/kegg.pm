# $Id$
#
# BioPerl module for Bio::SeqIO::kegg
#
# Cared for by Allen Day <allenday@ucla.edu>
#
# Copyright Allen Day
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::kegg - KEGG sequence input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the SeqIO handler system. Go:

    $stream = Bio::SeqIO->new(-file => $filename, -format => 'KEGG');

    while ( my $seq = $stream->next_seq() ) {
	# do something with $seq
    }


=head1 DESCRIPTION

This class transforms KEGG gene records into Bio::Seq objects.

=head2 Mapping of record properties to object properties

This section is supposed to document which sections and properties of
a GenBank databank record end up where in the Bioperl object model. It
is far from complete and presently focuses only on those mappings
which may be non-obvious. $seq in the text refers to the
Bio::Seq::RichSeqI implementing object returned by the parser for each
record.

=over 4

=item ENTRY

$seq-E<gt>primary_id

=item NAME

$seq-E<gt>display_id

=item DEFINITION

$seq-E<gt>annotation-E<gt>get_Annotations('description');

=item ORTHOLOG

grep {$_-E<gt>database eq 'KO'} $seq-E<gt>annotation-E<gt>get_Annotations('dblink')

=item CLASS

grep {$_-E<gt>database eq 'PATH'} $seq-E<gt>annotation-E<gt>get_Annotations('dblink')

=item POSITION

FIXME, NOT IMPLEMENTED

=item DBLINKS

$seq-E<gt>annotation-E<gt>get_Annotations('dblink')

=item CODON_USAGE

FIXME, NOT IMPLEMENTED

=item AASEQ

$seq-E<gt>translate-E<gt>seq

=item NTSEQ

$seq-E<gt>seq

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://www.bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Allen Day

Email allenday@ucla.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::kegg;
use vars qw(@ISA);
use strict;

use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::Species;
use Bio::Seq::SeqFactory;
use Bio::Annotation::Collection;
use Bio::Annotation::Comment;
use Bio::Annotation::DBLink;

@ISA = qw(Bio::SeqIO);
 
sub _initialize {
    my($self,@args) = @_;
    
    $self->SUPER::_initialize(@args); 
    # hash for functions for decoding keys.
    $self->{'_func_ftunit_hash'} = {}; 
    if( ! defined $self->sequence_factory ) {
	$self->sequence_factory(new Bio::Seq::SeqFactory
				(-verbose => $self->verbose(), 
				 -type => 'Bio::Seq::RichSeq'));
    }
}

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq::RichSeq object
 Args    :

=cut

sub next_seq {
  my ($self,@args) = @_;
  my $builder = $self->sequence_builder();
  my $seq;
  my %params;

  my $buffer;
  my (@acc, @features);
  my ($display_id, $annotation);
  my $species;

  # initialize; we may come here because of starting over
  @features = ();
  $annotation = undef;
  @acc = ();
  $species = undef;
  %params = (-verbose => $self->verbose); # reset hash
  local($/) = "///\n";

  $buffer = $self->_readline();

  return undef if( !defined $buffer ); # end of file
  $buffer =~ /^ENTRY/ ||
	$self->throw("KEGG stream with bad ENTRY line. Not KEGG in my book. Got '$buffer'");

  my %FIELDS;
  my @chunks = split /\n(?=\S)/, $buffer;

  foreach my $chunk (@chunks){
	my($key) = $chunk =~ /^(\S+)/;
	$FIELDS{$key} = $chunk;
  }

  my($entry_id,$entry_seqtype,$entry_species) = $FIELDS{ENTRY} =~ /^ENTRY\s+(\d+)\s+(\S+)\s+(\S+)\s*$/;

  my($name) = $FIELDS{NAME} =~ /^NAME\s+(.+)$/;

  my($definition) = $FIELDS{DEFINITION} =~ /^DEFINITION\s+(.+)$/s;
  $definition =~ s/\s+/ /gs;

  my($aa_length,$aa_seq) = $FIELDS{AASEQ} =~ /^AASEQ\s+(\d+)\n(.+)$/s;
  $aa_seq =~ s/\s+//g;
  my($nt_length,$nt_seq) = $FIELDS{NTSEQ} =~ /^NTSEQ\s+(\d+)\n(.+)$/s;
  $nt_seq =~ s/\s+//g;

  $annotation = Bio::Annotation::Collection->new();
  $annotation->add_Annotation('description',Bio::Annotation::Comment->new(-text => $definition));

  my($ortholog_db,$ortholog_id,$ortholog_desc) = $FIELDS{ORTHOLOG} =~ /^ORTHOLOG\s+(\S+):\s+(\S+)\s+(\S*)\s*$/;
  $annotation->add_Annotation('dblink',Bio::Annotation::DBLink->new(-database => $ortholog_db,
																	-primary_id => $ortholog_id,
																	-comment => $ortholog_desc
																   )
							 );

  $FIELDS{CLASS} =~ s/^CLASS\s+//;
  while($FIELDS{CLASS} =~ /.+?\[(\S+):(\S+)\]/gs){
	$annotation->add_Annotation('dblink',Bio::Annotation::DBLink->new(-database => $1, -primary_id => $2));
  }

  $FIELDS{DBLINKS} =~ s/^DBLINKS/       /;
  while($FIELDS{DBLINKS} =~ /\s+(\S+):\s+(\S+)\n/gs){
	$annotation->add_Annotation('dblink',Bio::Annotation::DBLink->new(-database => $1, -primary_id => $2));
  }

  $params{'-alphabet'}         = 'dna';
  $params{'-seq'}              = $nt_seq;
  $params{'-display_id'}       = $name;
  $params{'-accession_number'} = $entry_id;
  $params{'-species'}          = Bio::Species->new(-common_name => $entry_species);
  $params{'-annotation'}       = $annotation;

  $builder->add_slot_value(%params);
  $seq = $builder->make_object();

  return $seq;
}

1;
