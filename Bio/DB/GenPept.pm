# $Id$
#
# BioPerl module for Bio::DB::GenPept
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

# completely reworked by Jason Stajich to use Bio::DB::WebDBSeqI 2000-12-12

=head1 NAME

Bio::DB::GenPept - Database object interface to GenPept

=head1 SYNOPSIS

    $gb = new Bio::DB::GenPept;

    $seq = $gb->get_Seq_by_id('195055'); # Unique ID

    # or ...

    $seq = $gb->get_Seq_by_acc('DEECTH'); # Accession Number
    
    my $seqio = $gb->get_Stream_by_id(['195055', 'DEECTH']);
    while( my $seq = $seqio->next_seq ) {
	print "seq is is ", $seq->display_id, "\n";
    }

=head1 DESCRIPTION

Allows the dynamic retrieval of Sequence objects (Bio::Seq) from the GenPept
database at NCBI, via an Entrez query.

WARNING: Please do NOT spam the Entrez web server with multiple requests.
NCBI offers Batch Entrez for this purpose.  Batch Entrez support will likely
be supported in a future version of DB::GenPept.

Currently the only return format supported by NCBI Entrez for GenPept
database is GenPept format, so any format specification passed to
GenPept will be ignored still be forced to GenPept format (which is
just GenBank format).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Aaron Mackey, Jason Stajich

Email amackey@virginia.edu
Email jason@chg.mc.duke.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::GenPept;
use strict;
use vars qw(@ISA $DEFAULTFORMAT %PARAMSTRING );
use Bio::DB::NCBIHelper;

@ISA = qw(Bio::DB::NCBIHelper);
BEGIN { 
    $DEFAULTFORMAT = 'genpept';	    
    %PARAMSTRING = ( batch  => { 'DB'          => 'p',
				 'REQUEST_TYPE'=> 'LIST_OF_GIS',
				 'HTML'        => 'FALSE',
				 'SAVETO'      => 'FALSE',
				 'NOHEADER'    => 'TRUE' },
		     single => { 'db'    => 'p',
				 'form'  => '6',			     
				 'title' => 'no',			     
			     }
		     );
}

# the new way to make modules a little more lightweight
sub new {
  my($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->request_format($self->default_format);
  return $self;
}

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

=head2 default_format

 Title   : default_format
 Usage   : my $format = $self->default_format
 Function: Returns default sequence format for this module
 Returns : string
 Args    : none

=cut

sub default_format {
    return $DEFAULTFORMAT;
}

# from Bio::DB::WebDBSeqI from Bio::DB::RandomAccessI

=head2 Routines Bio::DB::WebDBSeqI from Bio::DB::RandomAccessI

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id('ROA1_HUMAN')
 Function: Gets a Bio::Seq object by its name
 Returns : a Bio::Seq object
 Args    : the id (as a string) of a sequence
 Throws  : "id does not exist" exception

=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $seq = $db->get_Seq_by_acc('AAC73346');
  Function: Gets a Seq objects by accession number
  Returns : Bio::Seq object
  Args    : accession number to retrive by

=head2 Routines implemented by Bio::DB::NCBIHelper

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
  Usage   : $seq = $db->get_Stream_by_acc($acc);
  Function: Gets a series of Seq objects by accession numbers
  Returns : a Bio::SeqIO stream object
  Args    : $ref : a reference to an array of accession numbers for
                   the desired sequence entries
  Note    : For GenBank, this just calls the same code for get_Stream_by_id()

=head2 request_format

 Title   : request_format
 Usage   : my $format = $self->request_format;
           $self->request_format($format);
 Function: Get/Set sequence format retrieval
 Returns : string representing format
 Args    : $format = sequence format

=cut

# oberride to force format to be GenPept regardless
sub request_format {
    my ($self) = @_;
    return $self->SUPER::request_format($self->default_format());
}

1;
__END__
















