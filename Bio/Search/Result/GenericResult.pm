# $Id$
#
# BioPerl module for Bio::Search::Result::GenericResult
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Result::GenericResult - Generic Implementation of Bio::Search::Result::ResultI interface applicable to most search results.

=head1 SYNOPSIS


    # typically one gets Results from a SearchIO stream
    use Bio::SearchIO;
    my $io = new Bio::SearchIO(-format => 'blast',
 			       -file   => 't/data/HUMBETGLOA.tblastx');
    my $result = $io->next_result;
    while( $hit = $result->next_hits()) {  
    # insert code here for hit processing
    }

    use Bio::Search::Result::GenericResult;
    my @hits = (); # would be a list of Bio::Search::Hit::HitI objects
    # typically these are created from a Bio::SearchIO stream
    my $result = new Bio::Search::Result::GenericResult
	( -query_name        => 'HUMBETGLOA',
	  -query_accession   => ''
	  -query_description => 'Human haplotype C4 beta-globin gene, complete cds.'
	  -query_length      => 3002
	  -database_name     => 'ecoli.aa'
	  -database_letters  => 4662239,
	  -database_entries  => 400,
	  -parameters        => { 'e' => '0.001' },
	  -statistics        => { 'kappa' => 0.731 },
	  -algorithm         => 'blastp',
           -algorithm_version => '2.1.2',
	  );

    my $id = $result->query_name();

    my $desc = $result->query_description();

    my $name = $result->database_name();

    my $size = $result->database_letters();

    my $num_entries = $result->database_entries();

    my $gap_ext = $result->get_parameter('e');

    my @params = $result->available_parameters;

    my $kappa = $result->get_statistic('kappa');

    my @statnames = $result->available_statistics;



=head1 DESCRIPTION

This object is an implementation of the Bio::Search::Result::ResultI
interface and provides a generic place to store results from a
sequence database search.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich and Steve Chervitz

Email jason@bioperl.org
Email sac@bioperl.org

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::Result::GenericResult;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::Search::Result::ResultI;

@ISA = qw(Bio::Root::Root Bio::Search::Result::ResultI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Search::Result::GenericResult();
 Function: Builds a new Bio::Search::Result::GenericResult object 
 Returns : Bio::Search::Result::GenericResult
 Args    : -query_name        => Name of query Sequence
           -query_accession   => Query accession number (if available)
           -query_description => Description of query sequence
           -query_length      => Length of query sequence
           -database_name     => Name of database
           -database_letters  => Number of residues in database
           -database_entries  => Number of entries in database
           -parameters        => hash ref of search parameters (key => value)
           -statistics        => hash ref of search statistics (key => value)
           -algorithm         => program name (blastx)
           -algorithm_version => version of the algorithm (2.1.2)
=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  $self->{'_hits'} = [];
  $self->{'_hitindex'} = 0;
  $self->{'_statistics'} = {};
  $self->{'_parameters'} = {};

  my ($qname,$qacc,$qdesc,$qlen,
      $dbname,$dblet,$dbent,$params,   
      $stats, $hits, $algo, $algo_v) = $self->_rearrange([qw(QUERY_NAME
							     QUERY_ACCESSION
							     QUERY_DESCRIPTION
							     QUERY_LENGTH
							     DATABASE_NAME
							     DATABASE_LETTERS
							     DATABASE_ENTRIES
							     PARAMETERS
							     STATISTICS
							     HITS
							     ALGORITHM
							     ALGORITHM_VERSION
							     )],@args);

  defined $algo   && $self->algorithm($algo);
  defined $algo_v && $self->algorithm_version($algo_v);

  defined $qname && $self->query_name($qname);
  defined $qacc  && $self->query_accession($qacc);
  defined $qdesc && $self->query_description($qdesc);
  defined $qlen  && $self->query_length($qlen);
  defined $dbname && $self->database_name($dbname);
  defined $dblet  && $self->database_letters($dblet);
  defined $dbent  && $self->database_entries($dbent);

  if( defined $params ) {
      if( ref($params) !~ /hash/i ) {
	  $self->throw("Must specify a hash reference with the the parameter '-parameters");
      }
      while( my ($key,$value) = each %{$params} ) {
	  $self->add_parameter($key,$value);
      }
  }
  if( defined $stats ) {
      if( ref($stats) !~ /hash/i ) {
	  $self->throw("Must specify a hash reference with the the parameter '-statistics");
      }
      while( my ($key,$value) = each %{$stats} ) {
	  $self->add_statistic($key,$value);
      }
  }

  if( defined $hits  ) { 
      $self->throw("Must define arrayref of Hits when initializing a $class\n") unless ref($hits) =~ /array/i;

      foreach my $s ( @$hits ) {
	  $self->add_hit($s);
      }
  }
  return $self;
}

=head2 algorithm

 Title   : algorithm
 Usage   : my $r_type = $hsp->algorithm
 Function: Obtain the name of the algorithm used to obtain the Result
 Returns : string (e.g., BLASTP)
 Args    : [optional] scalar string to set value

=cut

sub algorithm{
    my ($self,$value) = @_;
    my $previous = $self->{'_algorithm'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_algorithm'} = $value;
    } 
    return $previous;   
}

=head2 algorithm_version

 Title   : algorithm_version
 Usage   : my $r_version = $hsp->algorithm_version
 Function: Obtain the version of the algorithm used to obtain the Result
 Returns : string (e.g., 2.1.2)
 Args    : [optional] scalar string to set algorithm version value

=cut

sub algorithm_version{
    my ($self,$value) = @_;
    my $previous = $self->{'_algorithm_version'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_algorithm_version'} = $value;
    } 

    return $previous;   
}

=head2 Bio::Search::Result::ResultI interface methods

Bio::Search::Result::ResultI implementation

=head2 next_hit

 Title   : next_hit
 Usage   : while( $hit = $result->next_hit()) { ... }
 Function: Returns the next available Hit object, representing potential
           matches between the query and various entities from the database.
 Returns : a Bio::Search::Hit::HitI object or undef if there are no more.
 Args    : none


=cut

sub next_hit {
    my ($self,@args) = @_;
    my $index = $self->_nexthitindex;
    return undef if $index > scalar @{$self->{'_hits'}};
    return $self->{'_hits'}->[$index];    
}

=head2 query_name

 Title   : query_name
 Usage   : $id = $result->query_name();
 Function: Get the string identifier of the query used by the
           algorithm that performed the search.
 Returns : a string.
 Args    : [optional] new string value for query name

=cut

sub query_name {
    my ($self,$value) = @_;
    my $previous = $self->{'_queryname'};
    if( defined $value || ! defined $previous ) {
	$value = $previous = '' unless defined $value;
	$self->{'_queryname'} = $value;
    } 
    return $previous;
}

=head2 query_accession

 Title   : query_accession
 Usage   : $id = $result->query_accession();
 Function: Get the accession (if available) for the query sequence
 Returns : a string
 Args    : [optional] new string value for accession

=cut

sub query_accession {
    my ($self,$value) = @_;
    my $previous = $self->{'_queryacc'};
    if( defined $value || ! defined $previous ) {
	$value = $previous = '' unless defined $value;
	$self->{'_queryacc'} = $value;
    } 
    return $previous;
}

=head2 query_length

 Title   : query_length
 Usage   : $id = $result->query_length();
 Function: Get the length of the query sequence
           used in the search.
 Returns : a number
 Args    :  [optional] new integer value for query length

=cut

sub query_length {
    my ($self,$value) = @_;
    my $previous = $self->{'_querylength'};
    if( defined $value || ! defined $previous ) {
	$value = $previous = 0 unless defined $value;
	$self->{'_querylength'} = $value;
    } 
    return $previous;
}

=head2 query_description

 Title   : query_description
 Usage   : $id = $result->query_description();
 Function: Get the description of the query sequence
           used in the search.
 Returns : a string
 Args    : [optional] new string for the query description

=cut

sub query_description {
    my ($self,$value) = @_;
    my $previous = $self->{'_querydesc'};
    if( defined $value || ! defined $previous ) {
	$value = $previous = '' unless defined $value;
	$self->{'_querydesc'} = $value;
    } 
    return $previous;
}


=head2 database_name

 Title   : database_name
 Usage   : $name = $result->database_name()
 Function: Used to obtain the name of the database that the query was searched
           against by the algorithm.
 Returns : a scalar string
 Args    : [optional] new string for the db name

=cut

sub database_name {
    my ($self,$value) = @_;
    my $previous = $self->{'_dbname'};
    if( defined $value || ! defined $previous ) {
	$value = $previous = '' unless defined $value;
	$self->{'_dbname'} = $value;
    } 
    return $previous;
}

=head2 database_letters

 Title   : database_letters
 Usage   : $size = $result->database_letters()
 Function: Used to obtain the size of database that was searched against.
 Returns : a scalar integer (units specific to algorithm, but probably the
           total number of residues in the database, if available) or undef if
           the information was not available to the Processor object.
 Args    : [optional] new scalar integer for number of letters in db 


=cut

sub database_letters {
    my ($self,$value) = @_;
    my $previous = $self->{'_dbletters'};
    if( defined $value || ! defined $previous ) {
	$value = $previous = '' unless defined $value;
	$self->{'_dbletters'} = $value;
    } 
    return $previous;
}

=head2 database_entries

 Title   : database_entries
 Usage   : $num_entries = $result->database_entries()
 Function: Used to obtain the number of entries contained in the database.
 Returns : a scalar integer representing the number of entities in the database
           or undef if the information was not available.
 Args    : [optional] new integer for the number of sequence entries in the db


=cut

sub database_entries {
    my ($self,$value) = @_;
    my $previous = $self->{'_dbentries'};
    if( defined $value || ! defined $previous ) {
	$value = $previous = '' unless defined $value;
	$self->{'_dbentries'} = $value;
    } 
    return $previous;
}

=head2 get_parameter

 Title   : get_parameter
 Usage   : my $gap_ext = $report->get_parameter('gapext')
 Function: Returns the value for a specific parameter used
           when running this report
 Returns : string
 Args    : name of parameter (string)

=cut

sub get_parameter{
   my ($self,$name) = @_;
   return $self->{'_parameters'}->{$name};
}

=head2 available_parameters

 Title   : available_parameters
 Usage   : my @params = $report->available_paramters
 Function: Returns the names of the available parameters
 Returns : Return list of available parameters used for this report
 Args    : none

=cut

sub available_parameters{
   my ($self) = @_;
   return keys %{$self->{'_parameters'}};
}


=head2 get_statistic

 Title   : get_statistic
 Usage   : my $gap_ext = $report->get_statistic('kappa')
 Function: Returns the value for a specific statistic available 
           from this report
 Returns : string
 Args    : name of statistic (string)

=cut

sub get_statistic{
   my ($self,$key) = @_;
   return $self->{'_statistics'}->{$key};
}

=head2 available_statistics

 Title   : available_statistics
 Usage   : my @statnames = $report->available_statistics
 Function: Returns the names of the available statistics
 Returns : Return list of available statistics used for this report
 Args    : none

=cut

sub available_statistics{
   my ($self) = @_;
   return keys %{$self->{'_statistics'}};
}

=head2 Bio::Search::Report 

Bio::Search::Result::GenericResult specific methods

=head2 add_hit

 Title   : add_hit
 Usage   : $report->add_hit($hit)
 Function: Adds a HitI to the stored list of hits
 Returns : Number of HitI currently stored
 Args    : Bio::Search::Hit::HitI

=cut

sub add_hit {
    my ($self,$s) = @_;
    if( $s->isa('Bio::Search::Hit::HitI') ) { 
	push @{$self->{'_hits'}}, $s;
    } else { 
	$self->warn("Passed in " .ref($s). 
		    " as a Hit which is not a Bio::Search::HitI... skipping");
    }
    return scalar @{$self->{'_hits'}};
}


=head2 rewind

 Title   : rewind
 Usage   : $result->rewind;
 Function: Allow one to reset the Hit iteration to the beginning
           Since this is an in-memory implementation
 Returns : none
 Args    : none

=cut

sub rewind{
   my ($self) = @_;
   $self->{'_hitindex'} = 0;
}


=head2 _nexthitindex

 Title   : _nexthitindex
 Usage   : private

=cut

sub _nexthitindex{
   my ($self,@args) = @_;
   return $self->{'_hitindex'}++;
}



=head2 add_parameter

 Title   : add_parameter
 Usage   : $report->add_parameter('gapext', 11);
 Function: Adds a parameter
 Returns : none
 Args    : key  - key value name for this parama
           value - value for this parameter

=cut

sub add_parameter{
   my ($self,$key,$value) = @_;
   $self->{'_parameters'}->{$key} = $value;
}


=head2 add_statistic

 Title   : add_statistic
 Usage   : $report->add_statistic('lambda', 2.3);
 Function: Adds a parameter
 Returns : none
 Args    : key  - key value name for this parama
           value - value for this parameter

=cut

sub add_statistic {
   my ($self,$key,$value) = @_;
   $self->{'_statistics'}->{$key} = $value;
   return;
}

1;
