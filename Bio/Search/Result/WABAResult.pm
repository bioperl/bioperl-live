#
# BioPerl module for Bio::Search::Result::WABAResult
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

Bio::Search::Result::WABAResult - Result object for WABA alignment output

=head1 SYNOPSIS

# use this object exactly as you would a GenericResult
# the only extra method is query_database which is the 
# name of the file where the query sequence came from

=head1 DESCRIPTION

This object is for WABA result output, there is little difference
between this object and a GenericResult save the addition of one
method query_database.  Expect many of the fields for GenericResult to
be empty however as WABA was not intended to provide a lot of extra
information other than the alignment.

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


package Bio::Search::Result::WABAResult;
use strict;


use base qw(Bio::Search::Result::GenericResult);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::Result::WABAResult->new();
 Function: Builds a new Bio::Search::Result::WABAResult object 
 Returns : Bio::Search::Result::WABAResult
 Args    : -query_database => "name of the database where the query came from"


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($db) = $self->_rearrange([qw(QUERY_DATABASE)], @args);
  defined $db && $self->query_database($db);
  return $self;
}

=head2 query_database

 Title   : query_database
 Usage   : $obj->query_database($newval)
 Function: Data field for the database filename where the 
           query sequence came from
 Returns : value of query_database
 Args    : newvalue (optional)


=cut

sub query_database{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'query_database'} = $value;
    }
    return $self->{'query_database'};
}


=head2 All other methods are inherited from Bio::Search::Result::GenericResult

See the L<Bio::Search::Result::GenericResult> for complete
documentation of the rest of the methods that are available for this
module.

=cut

1;
