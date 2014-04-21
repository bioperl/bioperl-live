#
# BioPerl module for Bio::Search::Result::HMMERResult
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Result::HMMERResult - A Result object for HMMER results

=head1 SYNOPSIS

    use Bio::Search::Result::HMMERResult;
    my $result = Bio::Search::Result::HMMERResult->new
    ( -hmm_name => 'pfam',
      -sequence_file => 'roa1.pep',
      -hits => \@hits);

    # generally we use Bio::SearchIO to build these objects
    use Bio::SearchIO;
    my $in = Bio::SearchIO->new(-format => 'hmmer',
			       -file   => 'result.hmmer');
    while( my $result = $in->next_result ) {
	print $result->query_name, " ", $result->algorithm, " ", $result->num_hits(), " hits\n";
    }

=head1 DESCRIPTION

This is a specialization of L<Bio::Search::Result::GenericResult>.
There are a few extra methods, specifically L<sequence_file>,
L<hmm_name>, L<next_models>, and L<models>.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::Result::HMMERResult;
use strict;



use base qw(Bio::Search::Result::GenericResult);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::Result::HMMERResult->new();
 Function: Builds a new Bio::Search::Result::HMMERResult object 
 Returns : Bio::Search::Result::HMMERResult
 Args    : -hmm_name => string, name of hmm file
           -sequence_file => name of the sequence file

plus Bio::Search::Result::GenericResult parameters

           -query_name        => Name of query Sequence
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
           -program_reference => literature reference string for this algorithm

=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($hmm,$seqfile) = $self->_rearrange([qw(HMM_NAME SEQUENCE_FILE)],
					 @args);
  
  defined( $seqfile) && $self->sequence_file($seqfile);
  defined( $hmm) && $self->hmm_name($hmm);

  return $self;
}


=head2 hmm_name

 Title   : hmm_name
 Usage   : $obj->hmm_name($newval)
 Function: Get/Set the value of hmm_name
 Returns : value of hmm_name
 Args    : newvalue (optional)


=cut

sub hmm_name{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_hmm_name'} = $value;
    }
    return $self->{'_hmm_name'};
}


=head2 sequence_file

 Title   : sequence_file
 Usage   : $obj->sequence_file($newval)
 Function: Get/Set the value of sequence_file
 Returns : value of sequence_file
 Args    : newvalue (optional)


=cut

sub sequence_file{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_sequence_file'} = $value;
    }
    return $self->{'_sequence_file'};

}


=head2 next_model

 Title   : next_model
 Usage   : my $domain = $result->next_model
 Function: Returns the next domain - this
           is an alias for next_hit
 Returns : L<Bio::Search::Hit::HitI> object
 Args    : none


=cut

sub next_model{ shift->next_hit }

=head2 models

 Title   : models
 Usage   : my @domains = $result->models;
 Function: Returns the list of HMM models seen - this
           is an alias for hits()
 Returns : Array of L<Bio::Search::Hit::HitI> objects
 Args    : none


=cut

sub models{ shift->hits }

=head2 Bio::Search::Result::GenericResult inherited methods

=cut

=head2 algorithm

 Title   : algorithm
 Usage   : my $r_type = $hsp->algorithm
 Function: Obtain the name of the algorithm used to obtain the Result
 Returns : string (e.g., BLASTP)
 Args    : [optional] scalar string to set value

=cut

=head2 algorithm_version

 Title   : algorithm_version
 Usage   : my $r_version = $hsp->algorithm_version
 Function: Obtain the version of the algorithm used to obtain the Result
 Returns : string (e.g., 2.1.2)
 Args    : [optional] scalar string to set algorithm version value

=cut

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

=head2 query_name

 Title   : query_name
 Usage   : $id = $result->query_name();
 Function: Get the string identifier of the query used by the
           algorithm that performed the search.
 Returns : a string.
 Args    : [optional] new string value for query name

=cut

=head2 query_accession

 Title   : query_accession
 Usage   : $id = $result->query_accession();
 Function: Get the accession (if available) for the query sequence
 Returns : a string
 Args    : [optional] new string value for accession

=cut

=head2 query_length

 Title   : query_length
 Usage   : $id = $result->query_length();
 Function: Get the length of the query sequence
           used in the search.
 Returns : a number
 Args    :  [optional] new integer value for query length

=cut

=head2 query_description

 Title   : query_description
 Usage   : $id = $result->query_description();
 Function: Get the description of the query sequence
           used in the search.
 Returns : a string
 Args    : [optional] new string for the query description

=cut

=head2 database_name

 Title   : database_name
 Usage   : $name = $result->database_name()
 Function: Used to obtain the name of the database that the query was searched
           against by the algorithm.
 Returns : a scalar string
 Args    : [optional] new string for the db name

=cut

=head2 database_letters

 Title   : database_letters
 Usage   : $size = $result->database_letters()
 Function: Used to obtain the size of database that was searched against.
 Returns : a scalar integer (units specific to algorithm, but probably the
           total number of residues in the database, if available) or undef if
           the information was not available to the Processor object.
 Args    : [optional] new scalar integer for number of letters in db 


=cut

=head2 database_entries

 Title   : database_entries
 Usage   : $num_entries = $result->database_entries()
 Function: Used to obtain the number of entries contained in the database.
 Returns : a scalar integer representing the number of entities in the database
           or undef if the information was not available.
 Args    : [optional] new integer for the number of sequence entries in the db


=cut

=head2 get_parameter

 Title   : get_parameter
 Usage   : my $gap_ext = $report->get_parameter('gapext')
 Function: Returns the value for a specific parameter used
           when running this report
 Returns : string
 Args    : name of parameter (string)

=cut

=head2 available_parameters

 Title   : available_parameters
 Usage   : my @params = $report->available_paramters
 Function: Returns the names of the available parameters
 Returns : Return list of available parameters used for this report
 Args    : none

=cut

=head2 get_statistic

 Title   : get_statistic
 Usage   : my $gap_ext = $report->get_statistic('kappa')
 Function: Returns the value for a specific statistic available 
           from this report
 Returns : string
 Args    : name of statistic (string)

=cut

=head2 available_statistics

 Title   : available_statistics
 Usage   : my @statnames = $report->available_statistics
 Function: Returns the names of the available statistics
 Returns : Return list of available statistics used for this report
 Args    : none

=cut

=head2 Bio::Search::Result::GenericResult specific methods

=cut

=head2 add_hit

 Title   : add_hit
 Usage   : $report->add_hit($hit)
 Function: Adds a HitI to the stored list of hits
 Returns : Number of HitI currently stored
 Args    : Bio::Search::Hit::HitI

=cut

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


=head2 add_parameter

 Title   : add_parameter
 Usage   : $report->add_parameter('gapext', 11);
 Function: Adds a parameter
 Returns : none
 Args    : key  - key value name for this parama
           value - value for this parameter

=cut

=head2 add_statistic

 Title   : add_statistic
 Usage   : $report->add_statistic('lambda', 2.3);
 Function: Adds a parameter
 Returns : none
 Args    : key  - key value name for this parama
           value - value for this parameter

=cut

=head2 num_hits

 Title   : num_hits
 Usage   : my $hitcount= $result->num_hits
 Function: returns the number of hits for this query result
 Returns : integer
 Args    : none


=cut

=head2 hits

 Title   : hits
 Usage   : my @hits = $result->hits
 Function: Returns the available hits for this Result
 Returns : Array of L<Bio::Search::Hit::HitI> objects
 Args    : none


=cut

=head2 program_reference

 Title   : program_reference
 Usage   : $obj->program_reference($newval)
 Function: 
 Returns : value of the literature reference for the algorithm
 Args    : newvalue (optional)


=cut

1;
