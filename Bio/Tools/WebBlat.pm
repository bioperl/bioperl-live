# $Id$
#
# Perl module for Bio::Tools::WebBlat
#
# Cared for by Allen Day <allenday@ucla.edu>
#
# Copyright Allen Day
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::WebBlat - Run BLAT on UCSC CGI script

=head1 SYNOPSIS

  my $webblat = Bio::Tools::WebBlat->new();
  my $seq = Bio::Seq->new(display_id => 'foo' , seq => 'aataataat');
  my $searchio = $webblat->create_searchio($seq);

  while(my $result = $searchio->next_result){
    #process Bio::SearchIO::ResultI
  }

=head1 DESCRIPTION

Run BLAT at UCSC.  The useful method in this class after instantiation
is create_searchio.  This factory method takes a Bio::Seq object as
input, reformats it as Fasta, and sends it off to a CGI script running
at www.genome.ucsc.edu.  The resultant PSL embedded in HTML is
HTML-stripped and used to construct a Bio::SearchIO, which is passed
back to the user.

=head1 AUTHOR

  Allen Day E<lt>allenday@ucla.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::WebBlat;
use strict;
use base qw(Bio::Root::Root);

use Bio::SearchIO;
use Bio::SeqIO;
use Data::Dumper;
use IO::String;
use LWP::UserAgent;
use URI::Escape;

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::WebBlat->new();
 Function: Builds a new Bio::Tools::WebBlat object
 Returns : an instance of Bio::Tools::WebBlat
 Args    :


=cut

sub new {
  my($class,%arg) = @_;

  my $self = $class->SUPER::new(%arg);
  $self->init(%arg);

  return $self;
}

=head2 init

 Title   : init
 Usage   : $obj->init(%arg);
 Function: internal method.  initializes a new Bio::Tools::WebBlat object
 Returns : true on success
 Args    : args passed to new()


=cut

sub init {
  my($self,%arg) = @_;

  foreach my $arg (keys %arg){
    $self->$arg($arg{$arg}) if $self->can($arg);
  }

  return 1;
}

=head2 create_searchio

  Title   : create_searchio
  Usage   : my $searchio = $obj->create_searchio(%args);
  Function: connects to genome.ucsc.edu and runs a BLAT search against the
            organism/database you specify (see args).  search defaults to
            Human/Hg16
  Returns : returns a Bio::SearchIO object
  Args    : database  (optional) : database to search.  defaults to Hg16.
            organism  (optional) : organism to search.  defaults to Human.
            sequence  (required) : Bio::SeqI object to use in the search.

=cut

sub create_searchio {
  my ($self,%arg) = @_;

  $arg{database} ||= 'Hg16';
  $arg{organism} ||= 'Human';
  $self->throw("no Bio::Seq provided") unless ref($arg{sequence}) and $arg{sequence}->isa('Bio::Seq');

  my $fasta;
  my $out = Bio::SeqIO->new(-format => 'Fasta' , -fh => IO::String->new($fasta));
  $out->write_seq($arg{sequence});

  my $url = sprintf('http://www.genome.ucsc.edu/cgi-bin/hgBlat?org=%s&db=%s&type=%s&output=psl&userSeq=%s',
                    uri_escape($arg{organism}),
                    uri_escape($arg{database}),
                    uri_escape("BLAT's guess"),
                    uri_escape($fasta),
                   );

  my $ua = LWP::UserAgent->new(env_proxy => 1);
  my $response = $ua->get($url);
  if($response->is_success){

    my $html = $response->content;

    return unless $html =~ /psLayout/s;

    $html =~ s!^.+<pre>(.+?)</pre>.+$!$1!si;
    $html =~ s!</tt>$!!si; #yes, there is really bad tag nesting in the page

    my $psl = IO::String->new($html);

    my $searchio = Bio::SearchIO->new( -format => 'psl' , -fh => $psl );

    return $searchio;

  } else {
    $self->throw($response->status_line); }
}



#http://www.genome.ucsc.edu/cgi-bin/hgBlat?hgsid=30902160&org=Human&db=hg16&type=BLAT%27s+guess&output=psl&userSeq=%3EHS.22.q.1++++++%5B1-453%5D%0D%0AGATCTGATAAGTCCCAGGACTTCAGAAGAGCTGTGAGACCTTGGCCAAGTCACTTCCTCCTTCAGGAACAT

1;
