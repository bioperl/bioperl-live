
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
    "form=6&dopt=f&html=no&title=no&uid=$uid" ;

  my $hostname = `hostname`;
  chop $hostname;
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
		   "GET /htbin-post/Entrez/query?$entrez HTTP/1.0\r\n",
		   "Host: www3.ncbi.nlm.nih.gov",
		   "User-Agent: $0::Bio::DB::GenBank",
		   "", "");

  while( <$sock> ) {
    print STDOUT;
  }

#  return Bio::SeqIO->new(-fh => $sock, -format => 'Fasta')->next_seq()
}

1;

__END__
