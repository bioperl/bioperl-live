
package Bio::DB::GenBank;
use vars qw($AUTOLOAD @ISA @EXPORT_OK);
use strict;

use Bio::DB::Abstract;     # to inherit from

use Bio::Seq;              # to return a Seq object
use Bio::Tools::WWW qw(:obj);       # to get url information
use LWP::Simple qw(get);           # to access the internet

@ISA = qw(Bio::DB::Abstract Exporter);
@EXPORT_OK = qw();

sub get_Seq_by_id {

  my($self, @params) = @_;
  my($uid, $alpha) = $self->_rearrange( [qw(UID ALPHA)], @params);

  $alpha ||= 'n'; # default of nucleotide
  my $entrez = ($alpha =~ m/^p$/i ?
		$BioWWW->search_url('gb_p') :
		$BioWWW->search_url('gb_n')
	       ) . $uid . '&html=no&title=no';

  my $seq = new Bio::Seq;
  my $data = get($entrez);
  print STDOUT $data; # debugging only AJM
  $seq->parse($data, 'genbank');

  return $seq;
}

1;

__END__
