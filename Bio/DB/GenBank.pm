
package Bio::DB::GenBank;
use vars qw($AUTOLOAD @ISA @EXPORT_OK);
use strict;

use Bio::DB::Abstract;     # to inherit from

use Bio::SeqIO;                     # to parse and return Seq object
use IO::Socket;                     # to connect with GenBank

@ISA = qw(Bio::DB::Abstract Exporter);
@EXPORT_OK = qw();

sub get_Seq_by_id {

  my($self, @params) = @_;
  my($uid, $alpha) = $self->_rearrange( [qw(UID ALPHA)], @params);

  $alpha ||= 'n'; # default of nucleotide
  my $entrez = ($alpha =~ m/^p/i ? 'db=p&' : 'db=n&') .
    "form=6&dopt=g&html=no&title=no&uid=$uid" ;

  my $sock = IO::Socket::INET->new(PeerAddr => 'www3.ncbi.nlm.nih.gov',
				   PeerPort => 'http(80)',
				   Proto    => 'tcp' );
  select $sock;
  $sock->autoflush(1);
  print 'GET ' . "htbin-post/Entrez/query?$entrez";

  return Bio::SeqIO->new(-file => $sock, -format => 'EMBL')->next_seq()
}

1;

__END__
