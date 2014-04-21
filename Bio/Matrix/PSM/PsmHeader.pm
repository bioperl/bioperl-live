
=head1 NAME

Bio::Matrix::PSM::PsmHeader - PSM mast parser implementation

=head1 SYNOPSIS

  # See Bio::Matrix::PSM::IO for detailed documentation on how to use
  # PSM parsers

=head1 DESCRIPTION

Parser for mast. This driver unlike meme or transfac for example is
dedicated more to PSM sequence matches

=head1 FEEDBACK

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
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Stefan Kirov

Email skirov@utk.edu

=head1 APPENDIX

=cut


# Let the code begin...
package Bio::Matrix::PSM::PsmHeader;

use Bio::Matrix::PSM::InstanceSite;

use strict;
use base qw(Bio::Root::Root Bio::Matrix::PSM::PsmHeaderI);

#These define what structures within the
@Bio::Matrix::PSM::PsmHeader::MASTHEADER=qw(html version release seq hid 
					    length instances unstructured);
@Bio::Matrix::PSM::PsmHeader::MEMEHEADER=qw(html version release hid weight length unstructured);
@Bio::Matrix::PSM::PsmHeader::TRANSFACHEADER=qw(unstructured version release);
@Bio::Matrix::PSM::PsmHeader::PSIBLASTHEADER=qw(seq width ic);
@Bio::Matrix::PSM::PsmHeader::ALLHEADER=qw(header release type version html 
					   release weight length id 
					   seq instances unstructured);

=head2 new

 Title   : new
 Usage   : my $header= Bio::Matrix::PSM::PsmHeader->new(-seq=>\%seq, 
						       -mid=>\%mid, 
						       -width=>\%width,
                                                       -instances=>\%instances,
						       -header=>\@header,
						       -type=>'mast');
 Function: Creates a new Bio::Matrix::PSM::PsmHeader object
 Throws  :
 Example :
 Returns :  Bio::Matrix::PSM::PsmHeader object
 Args    :  hash


=cut

sub new {
    my ($class,@args)=@_;
    my $self = $class->SUPER::new(@args);
    return $self;
}

#parse version/release info here from the unstructured array
sub _initialize {
    my $self = shift;
    my $type=ref($self);
    $type=~s/\w+:://g;
    $self->{_type} = $type;
    my $dat=join(" ",grep(/version|release/i,@{$self->{unstructured}}));
    if ($dat && ($dat=~/version\b/i)) {
	$self->{version}=substr($dat,$+[0]+1);
	$self->{version}=~s/\s.+[^\d\.\:\/]//g;
	$self->{version}=~s/^\D//;
    }
    if ($dat && ($dat=~/release\b/i)) {
	my $rel=substr($dat,$+[0]+1);
	$rel=~s/[^\d\.\:\/\-]//g;
	$rel=~s/^\D//;
	if ($rel=~/\d\d:\d\d:\d\d/) { #Reformat if time is available too
	    my $time=substr($rel,$-[0]+1);
	    my $dat= substr($rel,0,$-[0]);
	    $self->{release}="$dat $time";
	}
	else {  $self->{release}=$rel; }
    }
    return $self;
}

=head2 seq

 Title   : seq
 Usage   : my %seq= $header->seq();
 Function: Returns the sequence data as a hash, indexed by a sequence ID (motif id or accession number)
           In case the input data is a motif it would return the consenus seq for each of them (mast).
 Throws  :
 Example :
 Returns :   hash
 Args    :


=cut

sub seq {
    my $self = shift;
    return () unless ($self->_check('seq'));
    return %{$self->{seq}};
}

=head2 hid

 Title   : hid
 Usage   : my @hid= $header->hid();
 Function: Returns array with the motif ids
 Throws  :
 Example :
 Returns :   array
 Args    :


=cut

sub hid {
    my $self = shift;
    return unless ($self->_check('hid'));
    my @header=@{$self->{hid}};
    return @header;
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
     return unless ($self->_check('length'));
    return $self->{length};
}

=head2 instances

 Title   : instances
 Usage   : my %instances= $header->instances();
 Function: Returns the info about the input data, contained in the header
 Throws  :
 Example :
 Returns : hash
 Args    :


=cut

sub instances {
      my $self = shift;
      return unless ($self->_check('instances'));
      return %{$self->{instances}};
}

=head2 weight

 Title   : weight
 Usage   : my %weights= $header->weight();
 Function: Returns the weights of the input sequence as a hash, indexed
           by a sequence ID
 Throws  :
 Example :
 Returns :  hash
 Args    :


=cut

sub weight {
    my $self = shift;
    return () unless ($self->_check('weight'));
    return %{$self->{weight}};
}


=head2 unstuctured

 Title   : unstuctured
 Usage   : my @unstructured= $header->unstuctured();
 Function: Returns the unstructured data in the header as an array, one line per
           array element, all control symbols are removed with \W
 Throws  :
 Example :
 Returns :  array
 Args    :


=cut

sub unstructured {
    my $self = shift;
    return @{$self->{unstructured}};
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
    return $self->{version};
}

=head2 release

 Title   : release
 Usage   : my $release= $header->release;
 Function: Returns the release of the file being parsed if such exists
 Throws  :
 Example :
 Returns :  string
 Args    :


=cut

sub release {
    my $self = shift;
    return $self->{release};
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
    my ($self,$method) = @_;
    my $type= $self->{'_type'};
    if ($type eq 'meme') { 
	return 0 unless (grep(/$method/,
				  @Bio::Matrix::PSM::PsmHeader::MEMEHEADER)); 
    } elsif ($type eq 'mast') { 
	return 0 unless (grep(/$method/,
				  @Bio::Matrix::PSM::PsmHeader::MASTHEADER));
    } elsif ($type eq 'transfac') { 
	return 0 unless (grep(/$method/,
				  @Bio::Matrix::PSM::PsmHeader::TRANSFACHEADER)); 
    } elsif ($type eq 'psiblast') { 
	return 0 unless (grep(/$method/,
				  @Bio::Matrix::PSM::PsmHeader::PSIBLASTHEADER)); 
    }
    return 1;
}

1;
