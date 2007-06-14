#
# $Id$
#
# BioPerl module for Bio::DB::XEMBL
#
# Cared for by Lincoln Stein
#
# Copyright Lincoln Stein
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::XEMBL - Database object interface for XEMBL entry retrieval

=head1 SYNOPSIS

  use Bio::DB::XEMBL;

  $embl = Bio::DB::XEMBL->new();

  # remember that XEMBL_ID does not equal GenBank_ID!
  $seq = $embl->get_Seq_by_id('BUM'); # EMBL ID
 	print "cloneid is ", $seq->id, "\n";

  # or changeing to accession number and Fasta format ...
  $seq = $embl->get_Seq_by_acc('J02231'); # XEMBL ACC
 	print "cloneid is ", $seq->id, "\n";

  # especially when using versions, you better be prepared
  # in not getting what what want
  eval {
      $seq = $embl->get_Seq_by_version('J02231.1'); # XEMBL VERSION
  };
  print "cloneid is ", $seq->id, "\n" unless $@;

  my $seqio = $embl->get_Stream_by_batch(['U83300','U83301','U83302']);
  while( my $clone =  $seqio->next_seq ) {
 	print "cloneid is ", $clone->id, "\n";
  }

=head1 DESCRIPTION

Allows the dynamic retrieval of Bio::Seq objects from the XEMBL
database. See L<Bio::Seq> for details.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Lincoln Stein

Email Lincoln Stein E<lt>lstein@cshl.orgE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::XEMBL;
use strict;
use SOAP::Lite;
# bsml parser appears broken...
use Bio::SeqIO::bsml;
use File::Temp 'tempfile';
use vars qw($MODVERSION);

use base qw(Bio::DB::RandomAccessI);
$MODVERSION = '0.2';

use constant DEFAULT_ENDPOINT => 'http://www.ebi.ac.uk:80/cgi-bin/xembl/XEMBL-SOAP.pl';

sub new {
    my ($class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
	$self->warn("The Bio::DB::XEMBL-related modules are deprecated as \n".
				"the XEMBL services are no longer available.  \nUse Bio::DB::DBFetch instead.\n");
    my $endpoint = $self->_rearrange([qw(ENDPOINT)]);
    $endpoint ||= DEFAULT_ENDPOINT;
    $self->endpoint($endpoint);
    return $self;
}

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id('ROA1_HUMAN')
 Function: Gets a Bio::Seq object by its name
 Returns : a Bio::Seq object
 Args    : the id (as a string) of a sequence
 Throws  : "id does not exist" exception


=cut

sub get_Seq_by_id {
   my ($self,@args) = @_;
   my $seqio = $self->get_Stream_by_batch([@args]);
   return $seqio->next_seq;
}

=head2 get_Stream_by_batch

  Title   : get_Stream_by_batch
  Usage   : $seq = $db->get_Stream_by_batch($ref);
  Function: Retrieves Seq objects from XEMBL 'en masse', rather than one
            at a time. Currently this is not particularly efficient, as
            it loads the entire result into memory and parses it.
  Example :
  Returns : a Bio::SeqIO stream object
  Args    : $ref : an array reference containing a list of unique 
            ids/accession numbers.

=cut

sub get_Stream_by_batch {
  my ($self, $ids) = @_;
  $self->throw("expected an array ref, but got $ids")
    unless ref($ids) eq 'ARRAY';
  my @args = @$ids;

  my $endpoint = $self->endpoint;
  my $som = SOAP::Lite
    ->uri('http://www.ebi.ac.uk/XEMBL')
    ->proxy($endpoint)
    ->getNucSeq(SOAP::Data->name(format=>'bsml'),
		SOAP::Data->name(ids=>"@args"));
  if ($som->fault) {
    $self->throw($som->faultstring);
  }
  my $result = $som->result;
  my($fh,$filename) = tempfile(File::Spec->tmpdir . '/bsmlXXXXXX',SUFFIX=>'.bsml');
  print $fh $result;
  close $fh;
  my $seqio = Bio::SeqIO->new(-file=>$filename,-format=>'bsml');
  unlink $filename;
  $seqio;
}

*get_Stream_by_id = \&get_Stream_by_batch;

=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   : $seq = $db->get_Seq_by_acc('X77802');
 Function: Gets a Bio::Seq object by accession number
 Returns : A Bio::Seq object
 Args    : accession number (as a string)
 Throws  : "acc does not exist" exception


=cut

sub get_Seq_by_acc{
   my ($self,@args) = @_;
   return $self->get_Seq_by_id(@args);
}

=head2 get_Seq_by_version

 Title   : get_Seq_by_version
 Usage   : $seq = $db->get_Seq_by_version('X77802.1');
 Function: Gets a Bio::Seq object by sequence version
 Returns : A Bio::Seq object
 Args    : accession.version (as a string)
 Throws  : "acc.version does not exist" exception

=cut

sub get_Seq_by_version{
   my ($self,@args) = @_;
   return $self->get_Seq_by_id(@args);
}

=head2 endpoint

 Title   : endpoint
 Usage   : $endpoint = $db->endpoint([$endpoint])
 Function: Gets/sets endpoint for SOAP connection
 Returns : old endpoint
 Args    : new endpoint(optional)

=cut

sub endpoint {
  my $self = shift;
  my $d = $self->{endpoint};
  $self->{endpoint} = shift if @_;
  $d;
}

=head2 new_from_registry

 Title   : new_from_registry
 Usage   : $db = Bio::DB::XEMBL->new_from_registry(%config)
 Function: creates a new Bio::DB::XEMBL object in a Bio::DB::Registry-
           compatible fashion
 Returns : new Bio::DB::XEMBL
 Args    : provided by the registry, see below
 Status  : Public

The following registry-configuration tags are recognized:

  location     Endpoint for the XEMBL service.  Currently the only
               known valid endpoint is 
               http://www.ebi.ac.uk:80/cgi-bin/xembl/XEMBL-SOAP.pl

NOTE: Since this info is supposed to be coming from WSDL, the location
is currently ignored.

=cut

sub new_from_registry {
    my ($self,%config) =  @_;
    my $location = $config{'location'} or $self->throw('Location must be specified.');
    my $index    = $self->new(-endpoint => $location);
}

1;
