# $Id$
#
# BioPerl module for Bio::DB::GenBank
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# 
# Added LWP support - Jason Stajich 2000-11-6
# completely reworked by Jason Stajich 2000-12-8
# to use WebDBSeqI

=head1 NAME

Bio::DB::GenBank - Database object interface to GenBank

=head1 SYNOPSIS

    use Bio::DB::GenBank;
    $gb = new Bio::DB::GenBank;

    $seq = $gb->get_Seq_by_id('MUSIGHBA1'); # Unique ID

    # or ...

    $seq = $gb->get_Seq_by_acc('J00522'); # Accession Number
    $seq = $gb->get_Seq_by_version('J00522.1'); # Accession.version
    $seq = $gb->get_Seq_by_gi('405830'); # GI Number

    # or ... best when downloading very large files, prevents
    # keeping all of the file in memory

    # also don't want features, just sequence so let's save bandwith
    # and request Fasta sequence
    $gb = new Bio::DB::GenBank(-retrievaltype => 'tempfile' , 
			       -format => 'Fasta');
    my $seqio = $gb->get_Stream_by_acc(['AC013798', 'AC021953'] );
    while( my $clone =  $seqio->next_seq ) {
      print "cloneid is ", $clone->display_id, " ", 
             $clone->accession_number, "\n";
    }
    # note that get_Stream_by_version is not implemented

=head1 DESCRIPTION

Allows the dynamic retrieval of Sequence objects (Bio::Seq) from the GenBank
database at NCBI, via an Entrez query.

WARNING: Please do NOT spam the Entrez web server with multiple requests.
NCBI offers Batch Entrez for this purpose.  Batch Entrez support will likely
be supported in a future version of DB::GenBank.

Note that when querying for GenBank accessions starting with 'NT_' you
will need to call $gb-E<gt>request_format('fasta') beforehand, because in
GenBank format (the default) the sequence part will be left out (the
reason is that NT contigs are rather annotation with references to
clones).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Aaron Mackey, Jason Stajich

Email amackey@virginia.edu
Email jason@chg.mc.duke.edu

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::GenBank;
use strict;
use vars qw(@ISA %PARAMSTRING);
use Bio::DB::NCBIHelper;

@ISA = qw(Bio::DB::NCBIHelper);
BEGIN {    
    %PARAMSTRING = ( 'batch'=>   {  'db'    => 'n',
				    'form'  => '1',			     
				    'title' => 'no', 
				},
#		     'batch'  => { 'cmd' => 'Retrieve',
#				   'db'  => 'n',
#				   
#			       },
		     'single'=>  { 'db'     => 'n',
				   'form'   => '1',			     
				   'title'  => 'no', 
			       },
		     'version'=> { 'pg'     => 'hist',
				   'type'   => 'acc',
			       },
		     'gi' =>     {  'db'    => 'n',
				    'form'  => '1',			     
				    'title' => 'no',
				}
		     );
}

# new is in NCBIHelper

# helper method to get db specific options

=head2 get_params

 Title   : get_params
 Usage   : my %params = $self->get_params($mode)
 Function: Returns key,value pairs to be passed to NCBI database
           for either 'batch' or 'single' sequence retrieval method
 Returns : a key,value pair hash
 Args    : 'single' or 'batch' mode for retrieval

=cut

sub get_params {
    my ($self, $mode) = @_;
    return %{$PARAMSTRING{$mode}};
}

# from Bio::DB::WebDBSeqI from Bio::DB::RandomAccessI

=head1 Routines Bio::DB::WebDBSeqI from Bio::DB::RandomAccessI

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id('ROA1_HUMAN')
 Function: Gets a Bio::Seq object by its name
 Returns : a Bio::Seq object
 Args    : the id (as a string) of a sequence
 Throws  : "id does not exist" exception

=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $seq = $db->get_Seq_by_acc($acc);
  Function: Gets a Seq object by accession numbers
  Returns : a Bio::Seq object
  Args    : the accession number as a string
  Note    : For GenBank, this just calls the same code for get_Seq_by_id()
  Throws  : "id does not exist" exception

=head2 get_Seq_by_gi

 Title   : get_Seq_by_gi
 Usage   : $seq = $db->get_Seq_by_gi('405830');
 Function: Gets a Bio::Seq object by gi number
 Returns : A Bio::Seq object
 Args    : gi number (as a string)
 Throws  : "gi does not exist" exception

=head2 get_Seq_by_version

 Title   : get_Seq_by_version
 Usage   : $seq = $db->get_Seq_by_version('X77802.1');
 Function: Gets a Bio::Seq object by sequence version
 Returns : A Bio::Seq object
 Args    : accession.version (as a string)
 Throws  : "acc.version does not exist" exception

=head1 Routines implemented by Bio::DB::NCBIHelper

=head2 get_request

 Title   : get_request
 Usage   : my $url = $self->get_request
 Function: HTTP::Request
 Returns : 
 Args    : %qualifiers = a hash of qualifiers (ids, format, etc)

=head2 get_Stream_by_batch

  Title   : get_Stream_by_batch
  Usage   : $seq = $db->get_Stream_by_batch($ref);
  Function: Retrieves Seq objects from Entrez 'en masse', rather than one
            at a time.  For large numbers of sequences, this is far superior
            than get_Stream_by_[id/acc]().
  Example :
  Returns : a Bio::SeqIO stream object
  Args    : $ref : either an array reference, a filename, or a filehandle
            from which to get the list of unique ids/accession numbers.

=head2 get_Stream_by_id

  Title   : get_Stream_by_id
  Usage   : $stream = $db->get_Stream_by_id( [$uid1, $uid2] );
  Function: Gets a series of Seq objects by unique identifiers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of unique identifiers for
                   the desired sequence entries

=head2 get_Stream_by_acc

  Title   : get_Stream_by_acc
  Usage   : $seq = $db->get_Stream_by_acc([$acc1, $acc2]);
  Function: Gets a series of Seq objects by accession numbers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of accession numbers for
                   the desired sequence entries
  Note    : For GenBank, this just calls the same code for get_Stream_by_id()

=cut

=head2 get_Stream_by_gi

  Title   : get_Stream_by_gi
  Usage   : $seq = $db->get_Seq_by_gi([$gi1, $gi2]);
  Function: Gets a series of Seq objects by gi numbers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of gi numbers for
                   the desired sequence entries
  Note    : For GenBank, this just calls the same code for get_Stream_by_id()

1;
__END__
