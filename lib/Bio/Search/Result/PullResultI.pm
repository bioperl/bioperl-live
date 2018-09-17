#
# BioPerl module Bio::Search::Result::PullResultI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Result::PullResultI - Bio::Search::Result::ResultI interface for
                                  'pull' parsers

=head1 SYNOPSIS

    # This is an interface and cannot be instantiated

    # typically one gets Results from a SearchIO stream
    use Bio::SearchIO;
    my $io = Bio::SearchIO->new(-format => 'hmmer_pull',
                                -file   => 't/data/hmmpfam.out');

    my $result = $io->next_result;

    while( $hit = $result->next_hit()) {
        # enter code here for hit processing
    }

    my $id = $result->query_name();
    my $desc = $result->query_description();
    my $dbname = $result->database_name();
    my $size = $result->database_letters();
    my $num_entries = $result->database_entries();
    my $gap_ext = $result->get_parameter('gapext');
    my @params = $result->available_parameters;
    my $kappa = $result->get_statistic('kappa');
    my @statnames = $result->available_statistics;

=head1 DESCRIPTION

Bio::Search::Result::ResultI objects are data structures containing
the results from the execution of a search algorithm.  As such, it may
contain various algorithm specific information as well as details of
the execution, but will contain a few fundamental elements, including
the ability to return Bio::Search::Hit::HitI objects.

PullResultI is for fast implementations that only do parsing work on the result
data when you actually request information by calling one of the ResultI
methods.

Many methods of ResultI are implemented in a way suitable for inheriting classes
that use Bio::PullParserI. It only really makes sense for PullResult modules to
be created by (and have as a -parent) SearchIO modules written using
PullParserI.

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

=head1 AUTHOR Sendu Bala

Email bix@sendu.me.uk

=head1 COPYRIGHT

Copyright (c) 2006 Sendu Bala.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::Result::PullResultI;

use strict;

use Bio::Search::GenericStatistics;
use Bio::Tools::Run::GenericParameters;

use base qw(Bio::PullParserI Bio::Search::Result::ResultI);

=head2 _setup

 Title   : _setup
 Usage   : $self->_setup(@args)
 Function: Implementers should call this to setup common fields and deal with
           common arguments to new().
 Returns : n/a
 Args    : @args received in new().

=cut

sub _setup {
    my ($self, @args) = @_;
    
    # fields most subclasses probably will want
    $self->_fields( { ( next_hit => undef,
                        num_hits => undef,
                        hits => undef,
                        no_hits_found => undef,
                        query_name => undef,
                        query_accession => undef,
                        query_length => undef,
                        query_description => undef  ) } );
    
    my ($parent, $chunk, $params, $stats) = $self->_rearrange([qw(PARENT
                                                                  CHUNK
                                                                  PARAMETERS
                                                                  STATISTICS)],
                                                              @args);
    $self->throw("Need -parent or -chunk to be defined") unless $parent || $chunk;
    
	$self->parent($parent) if $parent;
    
    if ($chunk) {
        my ($io, $start, $end) = (undef, 0, undef);
        if (ref($chunk) eq 'ARRAY') {
            ($io, $start, $end) = @{$chunk};
        }
        else {
            $io = $chunk;
        }
        $self->chunk($io, -start => $start, -end => $end);
    }
    
    if (defined $params) {
        if (ref($params) !~ /hash/i) {
            $self->throw("Must specify a hash reference with the the parameter '-parameters");
        }
        while (my ($key,$value) = each %{$params}) {
            $self->add_parameter($key, $value);
        }
    }
    if (defined $stats) {
        if (ref($stats) !~ /hash/i) {
            $self->throw("Must specify a hash reference with the the parameter '-statistics");
        }
        while (my ($key,$value) = each %{$stats}) {
            $self->add_statistic($key, $value);
        }
    }
}

#
# Some of these methods are written explitely to avoid ResultI throwing not
# implemented; if it didn't do that then PullParserI AUTOLOAD would have
# cought all them.
#

=head2 next_hit

 Title   : next_hit
 Usage   : while( $hit = $result->next_hit()) { ... }
 Function: Returns the next available Hit object, representing potential
           matches between the query and various entities from the database.
 Returns : a Bio::Search::Hit::HitI object or undef if there are no more.
 Args    : none

=cut

sub next_hit {
    return shift->get_field('next_hit');
}

=head2 sort_hits

 Title		: sort_hits
 Usage		: $result->sort_hits(\&sort_function)
 Function	: Sorts the available hit objects by a user-supplied function.
              Defaults to sort by descending score.
 Returns	: n/a
 Args		: A coderef for the sort function. See the documentation on the Perl
              sort() function for guidelines on writing sort functions.  
 Note		: To access the special variables $a and $b used by the Perl sort()
              function the user function must access Bio::Search::Result::ResultI namespace. 
              For example, use : 
              $result->sort_hits(sub{$Bio::Search::Result::ResultI::a->length <=> 
					                 $Bio::Search::Result::ResultI::b->length});
              NOT $result->sort_hits($a->length <=>$b->length);

=cut

# In ResultI. subclasses will probably want to override since sort_hits normally
# calls hits().

=head2 query_name

 Title   : query_name
 Usage   : $id = $result->query_name();
 Function: Get the string identifier of the query used by the
           algorithm that performed the search.
 Returns : a string.
 Args    : none

=cut

sub query_name {
    return shift->get_field('query_name');
}

=head2 query_accession

 Title   : query_accession
 Usage   : $id = $result->query_accession();
 Function: Get the accession (if available) for the query sequence
 Returns : a string
 Args    : none

=cut

sub query_accession {
    return shift->get_field('query_accession');
}

=head2 query_length

 Title   : query_length
 Usage   : $id = $result->query_length();
 Function: Get the length of the query sequence used in the search.
 Returns : a number
 Args    : none

=cut

sub query_length {
    return shift->get_field('query_length');
}

=head2 query_description

 Title   : query_description
 Usage   : $id = $result->query_description();
 Function: Get the description of the query sequence
           used in the search.
 Returns : a string
 Args    : none

=cut

sub query_description {
    return shift->get_field('query_description');
}

=head2 database_name

 Title   : database_name
 Usage   : $name = $result->database_name()
 Function: Used to obtain the name of the database that the query was searched
           against by the algorithm.
 Returns : a scalar string
 Args    : none

=cut

sub database_name {
    return shift->get_field('database_name');
}

=head2 database_letters

 Title   : database_letters
 Usage   : $size = $result->database_letters()
 Function: Used to obtain the size of database that was searched against.
 Returns : a scalar integer (units specific to algorithm, but probably the
           total number of residues in the database, if available) or undef if
           the information was not available to the Processor object.
 Args    : none

=cut

sub database_letters {
    return shift->get_field('database_letters');
}

=head2 database_entries

 Title   : database_entries
 Usage   : $num_entries = $result->database_entries()
 Function: Used to obtain the number of entries contained in the database.
 Returns : a scalar integer representing the number of entities in the database
           or undef if the information was not available.
 Args    : none

=cut

sub database_entries {
    return shift->get_field('database_entries');
}

=head2 algorithm

 Title   : algorithm
 Usage   : my $r_type = $result->algorithm
 Function: Obtain the name of the algorithm used to obtain the Result
 Returns : string (e.g., BLASTP)
 Args    : [optional] scalar string to set value

=cut

sub algorithm {
   return shift->get_field('algorithm');
}

=head2 algorithm_version

 Title   : algorithm_version
 Usage   : my $r_version = $result->algorithm_version
 Function: Obtain the version of the algorithm used to obtain the Result
 Returns : string (e.g., 2.1.2)
 Args    : [optional] scalar string to set algorithm version value

=cut

sub algorithm_version {
   return shift->get_field('algorithm_version');
}

=head2 algorithm_reference

 Title   : algorithm_reference
 Usage   : $obj->algorithm_reference($newval)
 Function: 
 Returns : value of the literature reference for the algorithm
 Args    : newvalue (optional)
 Comments: The default implementation in ResultI returns an empty string
           rather than throwing a NotImplemented exception, since
           the ref may not always be available and is not critical.

=cut

sub algorithm_reference {
   my ($self) = @_;
   return '';
}

=head2 num_hits

 Title   : num_hits
 Usage   : my $hitcount= $result->num_hits
 Function: returns the number of hits for this query result
 Returns : integer
 Args    : none

=cut

sub num_hits {
   return shift->get_field('num_hits');
}

=head2 hits

 Title   : hits
 Usage   : my @hits = $result->hits
 Function: Returns the HitI objects contained within this Result
 Returns : Array of Bio::Search::Hit::HitI objects
 Args    : none

See Also: L<Bio::Search::Hit::HitI>

=cut

sub hits {
   return shift->get_field('hits');
}

=head2 no_hits_found

 Usage     : $nohits = $blast->no_hits_found();
 Function  : Get boolean indicator indicating whether or not any hits
             were present in the report.

             This is NOT the same as determining the number of hits via
             the hits() method, which will return zero hits if there were no
             hits in the report or if all hits were filtered out during the
             parse.

             Thus, this method can be used to distinguish these possibilities
             for hitless reports generated when filtering.

 Returns   : Boolean
 Args      : none

=cut

sub no_hits_found {
    return shift->get_field('no_hits_found');
}

=head2 rewind

 Title   : rewind
 Usage   : $result->rewind;
 Function: Allow one to reset the Hit iterator to the beginning
           Since this is an in-memory implementation
 Returns : none
 Args    : none

=cut

sub rewind {
   shift->throw_not_implemented();
}

=head2 get_parameter

 Title   : get_parameter
 Usage   : my $gap_ext = $result->get_parameter('gapext')
 Function: Returns the value for a specific parameter used
           when running this result
 Returns : string
 Args    : name of parameter (string)

=cut

sub get_parameter {
    my ($self, $param) = @_;
    $param || return;
    return unless defined $self->{_parameters};
    return $self->{_parameters}->get_parameter($param);
}

=head2 available_parameters

 Title   : available_parameters
 Usage   : my @params = $result->available_parameters
 Function: Returns the names of the available parameters
 Returns : Return list of available parameters used for this result
 Args    : none

=cut

sub available_parameters {
    my $self = shift;
    return () unless defined $self->{_parameters};
    return $self->{_parameters}->available_parameters;
}

=head2 add_parameter

 Title   : add_parameter
 Usage   : $result->add_parameter('gapext', 11);
 Function: Adds a parameter
 Returns : none
 Args    : key  - key value name for this parama
           value - value for this parameter

=cut

sub add_parameter {
    my ($self, $key, $value) = @_;
    unless (exists $self->{_parameters}) {
        $self->{_parameters} = Bio::Tools::Run::GenericParameters->new();
    }
    $self->{_parameters}->set_parameter($key => $value);
}

=head2 get_statistic

 Title   : get_statistic
 Usage   : my $gap_ext = $result->get_statistic('kappa')
 Function: Returns the value for a specific statistic available 
           from this result
 Returns : string
 Args    : name of statistic (string)

=cut

sub get_statistic {
    my ($self, $stat) = @_;
    $stat || return;
    return unless defined $self->{_statistics};
    return $self->{_statistics}->get_statistic($stat);
}

=head2 available_statistics

 Title   : available_statistics
 Usage   : my @statnames = $result->available_statistics
 Function: Returns the names of the available statistics
 Returns : Return list of available statistics used for this result
 Args    : none

=cut

sub available_statistics {
    my $self = shift;
    return () unless defined $self->{_statistics};
    return $self->{_statistics}->available_statistics;
}

=head2 add_statistic

 Title   : add_statistic
 Usage   : $result->add_statistic('lambda', 2.3);
 Function: Adds a statistic
 Returns : none
 Args    : key  - key value name for this statistic
           value - value for this statistic

=cut

sub add_statistic {
    my ($self, $key, $value) = @_;
    unless (exists $self->{_statistics}) {
        $self->{_statistics} = Bio::Search::GenericStatistics->new();
    }
    $self->{_statistics}->set_statistic($key => $value);
}

1;
