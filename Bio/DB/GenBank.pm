
#
# BioPerl module for Bio::DB::GenBank
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::GenBank - Database object interface to GenBank

=head1 SYNOPSIS

    $gb = new Bio::DB::GenBank;

    $seq = $gb->get_Seq_by_id('ROA1_HUMAN'); # Unique ID

    # or ...

    $seq = $gb->get_Seq_by_acc('X77802'); # Accession Number

=head1 DESCRIPTION

Allows the dynamic retrieval of Sequence objects (Bio::Seq) from the GenBank
database at NCBI, via an Entrez query.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::GenBank;
use vars qw($AUTOLOAD @ISA @EXPORT_OK);
use strict;

# Object preamble - inherits from Bio::DB::Abstract

use Bio::DB::Abstract;

use Bio::SeqIO;
use IO::Socket;

@ISA = qw(Bio::DB::Abstract Exporter);
@EXPORT_OK = qw();

# new() is inherited from Bio::DB::Abstract

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id($uid, $type);
           $seq = $db->get_Seq_by_id(-uid  => $uid,
                                     -type => $type );
 Function: Gets a Bio::Seq object by its unique identifier/name
 Returns : a Bio::Seq object
 Args    : $uid : the id (as a string) of the desired sequence entry
           $type : 'nucleotide' or 'protein' - specifies whether to
                   search GenBank or GenPept for the sequence
                   (defaults to 'nucleotide').  Optional.

=cut

sub get_Seq_by_id {

  my($self, @params) = @_;
  my($uid, $type) = $self->_rearrange( [qw(UID TYPE)], @params);

  $type ||= 'n'; # default of nucleotide
  my $entrez = ($type =~ m/^p/i ? 'db=p&' : 'db=n&') .
    "form=6&dopt=f&html=no&title=no&uid=$uid" ;

# most of this socket stuff is borrowed heavily from LWP::Simple, by
# Gisle Aas and Martijn Koster.  They copyleft'ed it, but we should give
# them full credit for this little diddy.

  my $sock = IO::Socket::INET->new(PeerAddr => 'www3.ncbi.nlm.nih.gov',
				   PeerPort => 80,
				   Proto    => 'tcp',
				   Timeout  => 60
				  );
  unless ($sock) {
    $@ =~ s/^.*?: //;
    $self->throw("Can't connect to GenBank ($@)");
  }

  $sock->autoflush(); # just for safety's sake

  print $sock join("\015\012" =>
		   "GET /htbin-post/Entrez/query?$entrez HTTP/1.0",
		   "Host: www3.ncbi.nlm.nih.gov",
		   "User-Agent: $0::Bio::DB::GenBank",
		   "", "");

  while(<$sock>) {
    if ( m,^HTTP/\d+\.\d+\s+(\d+)[^\012]\012, ) {
      my $code = $1;
      return undef unless $code =~ /^2/;
    }
    return undef if /^ERROR/;
    last if m/^------/; # Kludgy, but it's how L. Stein does Boulder too
  }

  return Bio::SeqIO->new(-fh => $sock, -format => 'Fasta')->next_seq()
}

=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $seq = $db->get_Seq_by_acc($uid, $type);
           $seq = $db->get_Seq_by_acc(-acc  => $acc,
                                      -type => $type );
  Function: Gets a Bio::Seq object by its accession number
  Returns : a Bio::Seq object
  Args    : $acc : the accession number of the desired sequence entry
            $type : 'nucleotide' or 'protein' - specifies whether to
                    search GenBank or GenPept for the sequence
                    (defaults to 'nucleotide').  Optional.
  Note    : For GenBank, this just calls the same code for get_Seq_by_id()

=cut

sub get_Seq_by_acc {

  my($self, @params) = @_;
  my($acc, $type) = $self->_rearrange( [qw(ACC TYPE)], @params);

  $type ||= 'n';

  return $self->get_Seq_by_id(-uid => $acc, -type => $type);
}

1;
__END__











