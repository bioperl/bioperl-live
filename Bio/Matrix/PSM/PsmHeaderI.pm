#---------------------------------------------------------

=head1 NAME

Bio::Matrix::PSM::PsmHeaderI - handles the header data from a PSM file

=head1 SYNOPSIS

 use Bio::Matrix::PSM::IO;
 #Obtain an Bio::Matrix::PSM::IO object:
 my $psmIO= Bio::Matrix::PSM::IO->new(-file=>$file, -format=>'mast');

 #Get some general data about the file you are parsing:
 my $release=$psmIO->release;
 my $version=$psmIO->version;

 print "This analysis was performed using MAST version $version, release $release\n";

 #Now let's see what are the consensus sequences of the motifs fed as an input:
 my %seq=$psmIO->seq;

 #let's cycle through all consensus sequences now:

 foreach my $id ($psmIO->hid) {
   print "Motif $id is \t",$seq{$id},"\n";
 }

  #Finally look at the stuff we do not parse:
  my @inputfile=grep(/datafile/i,$psmIO->unstructured);

=head1 DESCRIPTION

Generally you should not use this object directly, you can access the
information through a PSM driver (See Bio::Matrix::PSM::IO). It is
handling the header data from a PSM file which may be very
different. This means that some of the methods will return undef
naturally, because this information is not present in the file which
is parsed. Some important data might be left over in the unstructured
part, and you might have to parse it yourself. I will try to
'structure' this header more in the near future.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

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

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Stefan Kirov

Email skirov@utk.edu

=head1 APPENDIX

=cut


# Let the code begin...
package Bio::Matrix::PSM::PsmHeaderI;
use Bio::Matrix::PSM::InstanceSite;
use Bio::Matrix::PSM::Psm;
use Bio::Matrix::PSM::IO;
use strict;
use base qw(Bio::Matrix::PSM::PsmI);

#Accessor methods, based on the driver
@Bio::Matrix::PSM::PsmHeader::MASTHEADER=qw(html version release 
					    seq hid length instances 
					    unstructured);
@Bio::Matrix::PSM::PsmHeader::MEMEHEADER=qw(html version release hid 
					    weight length unstructured);
@Bio::Matrix::PSM::PsmHeader::TRANSFACHEADER=qw(unstructured version release);
@Bio::Matrix::PSM::PsmHeader::ALLHEADER=qw(header release type version html 
					   release weight length hid 
					   seq instances unstructured);

=head2 new

 Title   : new
 Usage   : my $header= Bio::Matrix::PSM::PsmHeader->new
            ( -seq=>\%seq, -mid=>\%mid, -width=>\%width,
              -instances=>\%instances, -header=>\@header, -type=>'mast');
 Function: Creates a new Bio::Matrix::PSM::PsmHeader object
 Throws  :
 Example :
 Returns :  Bio::Matrix::PSM::PsmHeaderI object
 Args    :  hash


=cut

=head2 seq

 Title   : seq
 Usage   : my %seq= $header->seq();
 Function: Returns the sequence data as a hash, indexed by a 
           sequence ID (motif id or accession number)
           In case the input data is a motif it would return the 
           consenus seq for each of them (mast).
 Throws  :
 Example :
 Returns :  hash
 Args    :


=cut

sub seq {
     my $self = shift;
    $self->throw_not_implemented();
}


=head2 hid

 Title   : hid
 Usage   : my @ids= $header->hid();
 Function: Returns array with the motif/instance ids
 Throws  :
 Example :
 Returns :  array
 Args    :


=cut

sub hid {
     my $self = shift;
    $self->throw_not_implemented();
}

=head2 length

 Title   : length
 Usage   : my %length= $header->length();
 Function: Returns the length of the input sequence or motifs as a hash, indexed
           by a sequence ID (motif id or accession number)
 Throws  :
 Example :
 Returns :  hash
 Args    :


=cut

sub length {
     my $self = shift;
    $self->throw_not_implemented();
}

=head2 instances

 Title   : instances
 Usage   : my %instances= $header->length();
 Function: Returns the instance, used  as a hash, indexed
           by a sequence ID (motif id or accession number)
 Throws  :
 Example :
 Returns :  hash of Bio::Matrix::PSM::InstanceSite objects
 Args    :


=cut

sub instances {
     my $self = shift;
    $self->throw_not_implemented();
}

=head2 weights

 Title   : weights
 Usage   : my %weights= $header->weights();
 Function: Returns the weights of the input sequence as a hash, indexed
           by a sequence ID
 Throws  :
 Example :
 Returns :  hash
 Args    :


=cut

sub weights {
     my $self = shift;
    $self->throw_not_implemented();
}


=head2 unstuctured

 Title   : unstuctured
 Usage   : my @unstructured= $header->unstuctured();
 Function: Returns the unstructured data in the header as an array, one line per
           array element, all control symbols are removed with \W
 Throws  :
 Example :
 Returns :   array
 Args    :


=cut

sub unstructured {
     my $self = shift;
    $self->throw_not_implemented();
}

=head2 version

 Title   : version
 Usage   : my $version= $header->version;
 Function: Returns the version of the file being parsed if such exists
 Throws  :
 Example :
 Returns :  string
 Args    :


=cut

sub version {
     my $self = shift;
    $self->throw_not_implemented();
}

=head2 revision

 Title   : revision
 Usage   : my $revision= $header->revision;
 Function: Returns the revision of the file being parsed if such exists
 Throws  :
 Example :
 Returns :  string
 Args    :


=cut

sub revision {
     my $self = shift;
    $self->throw_not_implemented();
}

=head2 _check

 Title   : _check
 Usage   : if ($self->_check('weights') { #do something} else {return 0;}
 Function: Checks if the method called is aplicable to the file format
 Throws  :
 Example :
 Returns :  boolean
 Args    :  string


=cut

sub _check {
     my $self = shift;
    $self->throw_not_implemented();
}

1;
