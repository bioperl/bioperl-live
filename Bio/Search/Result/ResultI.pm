#-----------------------------------------------------------------
#
# BioPerl module Bio::Search::Result::ResultI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# Originally created by Aaron Mackey <amackey@virginia.edu>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Result::ResultI - Abstract interface to Search Result objects

=head1 SYNOPSIS

# Bio::Search::Result::ResultI objects cannot be instantiated since this
# module defines a pure interface.

# Given an object that implements the Bio::Search::Result::ResultI  interface,
# you can do the following things with it:

    use Bio::SearchIO;
    my $io = Bio::SearchIO->new(-format => 'blast',
                                -file   => 't/data/HUMBETGLOA.tblastx');
    my $result = $io->next_result;
    while( $hit = $result->next_hit()) { # enter code here for hit processing
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

=head1 AUTHOR 

Aaron Mackey E<lt>amackey@virginia.eduE<gt>  (original author)

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports and comments.

=head1 COPYRIGHT

Copyright (c) 1999-2001 Aaron Mackey, Steve Chervitz. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

#'
# Let the code begin...


package Bio::Search::Result::ResultI;

use strict;


use base qw(Bio::AnalysisResultI);


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
    $self->throw_not_implemented;
}

=head2 sort_hits

 Title		: sort_hits
 Usage		: $result->sort_hits(\&sort_function)
 Function	: Sorts the available hit objects by a user-supplied function. Defaults to sort
                  by descending score.
 Returns	: n/a
 Args		: A coderef for the sort function.  See the documentation on the Perl sort() 
                  function for guidelines on writing sort functions.  
 Note		: To access the special variables $a and $b used by the Perl sort() function 
                  the user function must access Bio::Search::Result::ResultI namespace. 
                  For example, use : 
                  $result->sort_hits( sub{$Bio::Search::Result::ResultI::a->length <=> 
					      $Bio::Search::Result::ResultI::b->length});
                   NOT $result->sort_hits($a->length <=>$b->length);

=cut

sub sort_hits {
    my ($self, $coderef) = @_;
    my @sorted_hits;

    if ($coderef)  {
	$self->throw('sort_hits requires a sort function passed as a subroutine reference')
	    unless (ref($coderef) eq 'CODE');
    }
    else {
	$coderef = \&_default_sort_hits;
	# throw a warning?
    }

    my @hits = $self->hits();
    
    eval {@sorted_hits = sort $coderef @hits };

   if ($@) {
       $self->throw("Unable to sort hits: $@");
   }
   else {
       $self->{'_hits'} = \@sorted_hits;
       $self->{'_no_iterations'} = 1; # to bypass iteration checking in hits() method
       1;
   }
}

=head2 _default sort_hits

  Title	: _default_sort_hits
  Usage	: Do not call directly.
  Function: Sort hits in descending order by score
  Args	: None
  Returns: 1 on success
  Note	: Used by $result->sort_hits()

=cut

sub _default_sort_hits {
    $Bio::Search::Result::ResultI::b->score <=> 
	    $Bio::Search::Result::ResultI::a->score;

}

=head2 query_name

 Title   : query_name
 Usage   : $id = $result->query_name();
 Function: Get the string identifier of the query used by the
           algorithm that performed the search.
 Returns : a string.
 Args    : none

=cut

sub query_name {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}

=head2 query_accession

 Title   : query_accession
 Usage   : $id = $result->query_accession();
 Function: Get the accession (if available) for the query sequence
 Returns : a string
 Args    : none

=cut

sub query_accession {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
}


=head2 query_length

 Title   : query_length
 Usage   : $id = $result->query_length();
 Function: Get the length of the query sequence
           used in the search.
 Returns : a number
 Args    : none

=cut

sub query_length {
    my ($self,@args) = @_;
    $self->throw_not_implemented;
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
    my ($self,@args) = @_;
    $self->throw_not_implemented;
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
    my ($self,@args) = @_;

    $self->throw_not_implemented;
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
    my ($self,@args) = @_;
    $self->throw_not_implemented();
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
    my ($self,@args) = @_;

    $self->throw_not_implemented();
}

=head2 get_parameter

 Title   : get_parameter
 Usage   : my $gap_ext = $result->get_parameter('gapext')
 Function: Returns the value for a specific parameter used
           when running this result
 Returns : string
 Args    : name of parameter (string)

=cut

sub get_parameter{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 available_parameters

 Title   : available_parameters
 Usage   : my @params = $result->available_parameters
 Function: Returns the names of the available parameters
 Returns : Return list of available parameters used for this result
 Args    : none

=cut

sub available_parameters{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 get_statistic

 Title   : get_statistic
 Usage   : my $gap_ext = $result->get_statistic('kappa')
 Function: Returns the value for a specific statistic available 
           from this result
 Returns : string
 Args    : name of statistic (string)

=cut

sub get_statistic{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 available_statistics

 Title   : available_statistics
 Usage   : my @statnames = $result->available_statistics
 Function: Returns the names of the available statistics
 Returns : Return list of available statistics used for this result
 Args    : none

=cut

sub available_statistics{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 algorithm

 Title   : algorithm
 Usage   : my $r_type = $result->algorithm
 Function: Obtain the name of the algorithm used to obtain the Result
 Returns : string (e.g., BLASTP)
 Args    : [optional] scalar string to set value

=cut

sub algorithm{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 algorithm_version

 Title   : algorithm_version
 Usage   : my $r_version = $result->algorithm_version
 Function: Obtain the version of the algorithm used to obtain the Result
 Returns : string (e.g., 2.1.2)
 Args    : [optional] scalar string to set algorithm version value

=cut

sub algorithm_version{
   my ($self) = @_;
   $self->throw_not_implemented();
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

sub algorithm_reference{
   my ($self) = @_;
   return '';
}

=head2 rid

 Title   : rid
 Usage   : $obj->rid($newval)
 Function:
 Returns : value of the BLAST Request ID (eg. RID: ZABJ4EA7014)
 Args    : newvalue (optional)
 Comments: The default implementation in ResultI returns an empty string
           rather than throwing a NotImplemented exception, since
           the RID may not always be available and is not critical.
           See: (1) http://www.ncbi.nlm.nih.gov/Class/MLACourse/Modules/BLAST/rid.html
                (2) http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node63.html
=cut

sub rid{
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

sub num_hits{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 hits

 Title   : hits
 Usage   : my @hits = $result->hits
 Function: Returns the HitI objects contained within this Result
 Returns : Array of Bio::Search::Hit::HitI objects
 Args    : none

See Also: L<Bio::Search::Hit::HitI>

=cut

sub hits{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}


=head2 no_hits_found

 Usage     : $nohits = $blast->no_hits_found();
 Purpose   : Get boolean indicator indicating whether or not any hits
             were present in the report.

             This is NOT the same as determining the number of hits via
             the hits() method, which will return zero hits if there were no
             hits in the report or if all hits were filtered out during the parse.

             Thus, this method can be used to distinguish these possibilities
             for hitless reports generated when filtering.

 Returns   : Boolean
 Argument  : none

=cut

#-----------
sub no_hits_found { shift->throw_not_implemented }



=head2 set_no_hits_found

 Usage     : $blast->set_no_hits_found(); 
 Purpose   : Set boolean indicator indicating whether or not any hits
             were present in the report.
 Returns   : n/a
 Argument  : none

=cut

sub set_no_hits_found { shift->throw_not_implemented }

1;


