#
# $Id$
#
# BioPerl module for Bio::DB::SRS
#
# Written by Lorenz Pollak (lorenz.pollak@pharma.novartis.com)
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::SRS - Database object interface to SRS

=head1 SYNOPSIS

    $srs = new Bio::DB::SRS(-databases=>{'GENBANKNEW'=>'GenBank'},
                            -location=>'/bioinf/srs6/bin/getz');

    $seq = $srs->get_Seq_by_id('MUSIGHBA1'); # Unique ID

    # or ...

    $srs = new Bio::DB::SRS(-databases=>{'GENBANK'=>'GenBank',
                                         'EMBL'=>'EMBL'},
                            -location=>'http://srs.ebi.ac.uk'.
                                       '/srs6bin/cgi-bin/wgetz');

    $srs->ids(["AF010000","AE032*"]);
    $srs->min_seqlength(400);
    $ids_arrayref = $srs->get_ids();   # get array of found entry ids
    $seqs_arrayref = $srs->get_seqs(); # get array of found entries

=head1 DESCRIPTION

Allows the dynamic retrieval of Sequence objects (Bio::Seq) from a SRS
(Sequence Retrieval System) database that can be accessed locally or
via the http protocol (using the LWP modules). If you are behind a
firewall you just have to set the environment variable "http_proxy"
to fit your needs.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  bioperl-guts-l@bioperl.org         - Technically-oriented discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via email or the
web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Lorenz Pollak

Email lorenz.pollak@pharma.novartis.com

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::SRS;
use strict;
use vars qw(@ISA);

# Object preamble - inherits from Bio::DB::RandomAccessI

use Bio::DB::RandomAccessI;
use Bio::SeqIO;
use LWP::UserAgent;
use HTTP::Request;

@ISA = qw(Bio::Root::Object Bio::DB::RandomAccessI);

# this code is taken from the webblast module, thanks to the authors!
BEGIN {
  unless( eval "require HTTP::Request" and
	  eval "require LWP::UserAgent") {
    warn "\a\n".'='x50, "\n".
      "WARNING: COULDN'T LOAD THE LWP MODULE.\n\n".
      "   Download it from CPAN: http://www.perl.com/CPAN/.".
	"\n".'='x50, "\n\n";
  }
  unless( eval "require IO::Scalar") {
    warn "\a\n".'='x50, "\n".
      "WARNING: COULDN'T LOAD THE IO::Scalar MODULE.\n\n".
      "   This module is included in the IO-stringy collection\n".
      "   from CPAN: http://www.perl.com/CPAN/.".
	"\n".'='x50, "\n\n";
  }
}

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);
  
  my ($location,$is_local,$format,$databases) = 
    $self->_rearrange([qw(LOCATION
                          DATABASES
		     )],
		      @args);
  $location = "/bioinf/srs/srs/bin/irix64/getz" unless $location;
  $self->location($location);
  $self->databases($databases);
  # set stuff in self from @args  
  return $make; # success - we hope!
}

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $srs->get_Seq_by_id($id);
 Function: Gets a Bio::Seq object by its unique identifier/name
 Returns : a Bio::Seq object
 Args    : $id : the id (as a string) of the desired sequence entry

=cut

sub get_Seq_by_id {

  my $self = shift;
  my $id = shift or $self->throw("Must supply an identifier!\n");

  $self->ids([$id]);

  my $seqsref = $self->get_seqs();
  my $seq = $$seqsref[0];
  $self->throw("Unable to get sequence for id $id!\n") 
      if ( !defined $seq );
  return $seq;
}

=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $seq = $srs->get_Seq_by_acc($acc);
  Function: Gets a Bio::Seq object by its accession number
  Returns : a Bio::Seq object
  Args    : $acc : the accession number of the desired sequence entry

=cut

sub get_Seq_by_acc {

  my $self = shift;
  my $acc = shift or $self->throw("Must supply an accesion number!\n");
  
  $self->accs([$acc]);

  my $seqsref = $self->get_seqs();
  my $seq = $$seqsref[0];
  $self->throw("Unable to get sequence for id $acc!\n") 
      if ( !defined $seq );
  return $seq;
}

=head2 location

 Title   : location
 Usage   : $location = $srs->location($newloc);
 Function: Returns or sets the location of the srs binary (getz/wgetz).
           This location can either be a local path or an URL.
           If the string starts with "http://" then it is assumed
           that the following is an URL and is_local(0) will be set.
           Otherwise is_local(1) will be set automatically.
           Examples:
           "/bioinf/srs6/bin/getz"
           "http://srs.ebi.ac.uk/srs6bin/cgi-bin/wgetz"
 Returns : a string
 Args    : optional new location

=cut

sub location {
  my ($self, $location) = @_;

  if(defined $location) {
    # try to guess _is_local if not set
    if (! exists $self->{'_is_local'}) {
      $self->{'_is_local'} = 1;
      $self->{'_is_local'} = 0 if ($location =~ /^http:\/\//);
    }
    # do some testings
    if ($self->{'_is_local'}) {
      (-e $location) || 
	$self->throw($location." does not exist!");
      (-x $location) || 
	$self->throw($location." can not be executed!");
    }
    $self->{'_location'} = $location;
  }
  return ($self->{'_location'});
}

=head2 databases

 Title   : databases
 Usage   : $dbs_hashref = $srs->databases($newdbs_hashref);
 Function: Returns or sets the current database settings for
           retrieving sequence entries with SRS. This must be
           a reference to a hash with the following format:
           'DBNAME' => 'DBFORMAT'
           Example:
           $srs->databases({'EMBL' => 'embl',
                            'GENBANK' => 'GenBank',
                            'GENBANKNEW' => 'GenBank',
                            'SWISSPROT' => 'swiss',
                            'PIR' => 'pir',
                            'SPTREMBL' => 'swiss',
                            'REMTREMBL' => 'swiss'})
           Only specify those databases in which you actually want
           to perform a search / sequence retrieval.
 Returns : a hash reference
 Args    : optional new databases (as hash reference)

=cut

sub databases {
  my ($self, $dbsref) = @_;

  if(defined $dbsref) {
    if (ref($dbsref) ne "HASH") {
      $self->throw("databases must be a hash reference!\n".
		   "(format: 'dbname'=>'format')\n");
    }
    my $newdbsref = {};
    foreach my $dbname (keys %$dbsref) {
      $_ = $$dbsref{$dbname};
      $dbname = uc $dbname;
      if (/fasta/i) {
	$$newdbsref{$dbname} = "Fasta";
      }
      elsif (/embl/i) {
	$$newdbsref{$dbname} = "EMBL";
      }
      elsif (/genbank/i) {
	$$newdbsref{$dbname} = "GenBank";
      }
      elsif (/swiss/i) {
	$$newdbsref{$dbname} = "swiss";
      }
      elsif (/scf/i) {
	$$newdbsref{$dbname} = "SCF";
      }
      elsif (/pir/i) {
	$$newdbsref{$dbname} = "PIR";
      }
      elsif (/gcg/i) {
	$$newdbsref{$dbname} = "GCG";
      }
      elsif (/raw/i) {
	$$newdbsref{$dbname} = "raw";
      }
      elsif (/ace/i) {
	$$newdbsref{$dbname} = "ace";
      }
      else {
	$self->throw("Database format $_ unknown!\n");
      }
    }
    $self->{'_dbsref'} = $newdbsref;
  }
  return ($self->{'_dbsref'});
}

=head2 get_seqs

 Title   : get_seqs
 Usage   : $seqs_arrayref = $srs->get_seqs();
 Function: Returns an array of Bio::Seq objects for all entries
           that resulted from your actual retrieval settings
 Returns : a reference to an array of Bio::Seq objects
 Args    :

=cut

sub get_seqs {
  my ($self) = @_;
  my $ids_arrayref = $self->get_ids();
  $self->{'_full_records'} = 1;
  my $fh = $self->_run_srs($self->_build_query());
  my (@seqs_array,$streamfmt,$oldstreamfmt,$stream);
  $oldstreamfmt = 0;
  for(my $i=0; $i<@$ids_arrayref ; $i++) {
    next if ($ids_arrayref->[$i][2]); # double id
    $streamfmt = $self->{'_dbsref'}{$ids_arrayref->[$i][0]};
    # if we are using wgetz then skip first line because
    # wgetz prints "db:id" again
    if (!$self->{'_is_local'}) {
      my $dummy = <$fh>;
    }
    if ($streamfmt ne $oldstreamfmt) {
      $stream = Bio::SeqIO->new('-fh' => $fh, '-format' => $streamfmt);
      $oldstreamfmt = $streamfmt;
    }
    push @seqs_array, $stream->next_seq();
  }
  close $fh;
  return \@seqs_array;
}

=head2 get_ids

 Title   : get_ids
 Usage   : $ids_arrayref = $srs->get_ids();
 Function: Returns an array of database/id pairs for all entries
           that resulted from your actual retrieval settings
 Returns : a reference to a 2D array
           example: $dbname   = $ids_arrayref->[$i][0]
                    $entry_id = $ids_arrayref->[$i][1]
                    $double   = $ids_arrayref->[$i][2]
           this example would get values from the database/id pair
           at index number $i. The third value ($double) is normally
           set to zero. Only if the same id has already appeared in
           another database from your result $double will be set
           "true".
 Args    :
 Note    : Normally you should only need the function get_seqs

=cut

sub get_ids {
  my ($self) = @_;
  $self->{'_full_records'} = 0;
  my $fh = $self->_run_srs($self->_build_query());
  my (@ids_array,$ids_string,$double,$db,$id);
  $ids_string = "";
  while (<$fh>) {
    if (/([\w\d\-_]*):([\w\d\-_]*)/) {
      $db = $1; $id = $2;
      $double = 0;
      # check if id is identical with previous one
      if ($ids_string =~ /\s$id\s/) {$double = 1};
      push @ids_array, [$db,$id,$double];
      $ids_string .= " $id ";
    }
  }
  close $fh;
  return \@ids_array;
}

sub _strip_html {
  my ($self,$text) = @_;
  $text =~ s/<[^>]*>//gs;
  $text =~ s/&nbsp;/ /gs;
  $text =~ s/\s+\n/\n/gs;
  $text =~ s/\n+/\n/gs;
  $text =~ s/^\n//gs;
  return $text;
}

sub _run_srs {
  my ($self, $param) = @_;
  my $fh = Symbol::gensym();
  if ($self->{'_is_local'}) {
    open($fh, $self->{'_location'}." ".$param." |") || 
      $self->throw("Could not fork: $!\n");
  } else {
    my $ua = new LWP::UserAgent;
    my $req = new HTTP::Request 'GET',$self->{'_location'}."?".$param;
    $req->header('Accept' => 'text/html');
    $ua->env_proxy;
    my $res = $ua->request($req);
    if ($res->is_error) {
      $self->throw("An error occured while posting to ".$self->{'_location'}.
		  "\n".$res->status_line."\n");
    }
    my $content = $self->_strip_html($res->content);
    tie *$fh, 'IO::Scalar', \$content;
  }
  return $fh;
}

sub _build_query {
  my ($self) = @_;
  my $query = "";
  # 1) build query string
  if ((exists $self->{'_idsref'}) && (@{$self->{'_idsref'}})) {
    my $ids = join "|", @{$self->{'_idsref'}};
    $query .= "|[libs-id:".$ids."]";
  }
  if ((exists $self->{'_accsref'}) && (@{$self->{'_accsref'}})) {
    my $accs = join "|", @{$self->{'_accsref'}};
    $query .= "|[libs-acc:".$accs."]";
  }
  if ($self->{'_description'}) {
    $query .= "&[libs-description:".$self->{'_description'}."]";
  }
  if ($self->{'_min_seqlen'} || $self->{'_min_seqlen'}) {
    my $min_seqlen = "";
    my $max_seqlen = "";
    $min_seqlen = $self->{'_min_seqlen'} if ($self->{'_min_seqlen'});
    $max_seqlen = $self->{'_max_seqlen'} if ($self->{'_max_seqlen'});
    $query .= "&[libs-SeqLength#".$min_seqlen.":".$max_seqlen."]";
  }
  $query =~ s/^(\||&)//;
  my $libs = join " ",(keys %{$self->{'_dbsref'}});
  if ($self->{'_is_local'}) { # normal commandline ? (getz)
    $query =~ s/\[libs/[libs={$libs}/;
    $query = "\"$query\"";
  } else {                   # or url query string... (wgetz)
    $query =~ s/libs/{$libs}/g;
    $query =~ s/\s/_SP_/g;
    $query =~ s/#/%23/g;
  }
  # 2) build params
  $query = "-e ".$query if ($self->{'_full_records'});
  if (!$self->{'_is_local'}) { # url query string ? (wgetz)
    $query =~ s/\s/+/g;
  }
  return $query;
}

=head2 ids

 Title   : ids
 Usage   : $ids_arrayref = $srs->ids($newids_arrayref);
 Function: Returns or sets id strings for retrieving sequence
           entries with SRS (as reference to an array of strings).
           These ids will be joint with "|" (OR)
 Returns : an array reference
 Args    : optional new id search strings (as array reference)

=cut

sub ids {
  my ($self, $idsref) = @_;

  if(defined $idsref) {
    if (ref($idsref) ne "ARRAY") {
      $self->throw("ids must be an array reference!\n");
    }
    $self->{'_idsref'} = $idsref;
  }
  return ($self->{'_idsref'});
}

=head2 accs

 Title   : accs
 Usage   : $accs_arrayref = $srs->accs($newaccs_arrayref);
 Function: Returns or sets accession strings for retrieving sequence
           entries with SRS (as reference to an array of strings).
           These accessions will be joint with "|" (OR)
 Returns : an array reference
 Args    : optional new accession search strings (as array reference)

=cut

sub accs {
  my ($self, $accsref) = @_;

  if(defined $accsref) {
    if (ref($accsref) ne "ARRAY") {
      $self->throw("accs must be an array reference!\n");
    }
    $self->{'_accsref'} = $accsref;
  }
  return ($self->{'_accsref'});
}

=head2 description

 Title   : description
 Usage   : $desc = $srs->description($newdesc);
 Function: Returns or sets the string that should be used to
           search the Description field of sequence entries
 Returns : a string
 Args    : optional new description (as string)

=cut

sub description {
  my ($self, $description) = @_;

  if(defined $description) {
    $description =~ s/\s(\||&)\s/$1/;
    $description =~ s/^\s*(.*)\s*$/$1/;
    $self->{'_description'} = $description;
  }
  return ($self->{'_description'});
}

=head2 min_seqlen

 Title   : min_seqlen
 Usage   : $min = $srs->min_seqlen($newmin);
 Function: Returns or sets the minimum value for the length of
           retrieved sequences
 Returns : a scalar
 Args    : optional new value (scalar)

=cut

sub min_seqlen {
  my ($self, $min_seqlen) = @_;

  $self->{'_min_seqlen'} = $min_seqlen if(defined $min_seqlen);
  return ($self->{'_min_seqlen'});
}

=head2 max_seqlen

 Title   : max_seqlen
 Usage   : $max = $srs->max_seqlen($newmax);
 Function: Returns or sets the maximum value for the length of
           retrieved sequences
 Returns : a scalar
 Args    : optional new value (scalar)

=cut

sub max_seqlen {
  my ($self, $max_seqlen) = @_;

  $self->{'_max_seqlen'} = $max_seqlen if(defined $max_seqlen);
  return ($self->{'_max_seqlen'});
}

=head2 is_local

 Title   : is_local
 Usage   : $is_local = $srs->is_local($new_is_local);
 Function: Returns or sets whether the binary location that has been
           set with the function "location" should be regarded as
           local (= run local binary) or remote (= access location
           via http)
 Returns : a boolean value
 Args    : optional new value (boolean)
 Note    : Normally the function "location" will call "is_local" 
           with proper values automatically!

=cut

sub is_local {
  my ($self, $is_local) = @_;

  $self->{'_is_local'} = $is_local if(defined $is_local);
  return ($self->{'_is_local'});
}


1;
__END__
