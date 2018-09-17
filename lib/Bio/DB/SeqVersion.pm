#
# BioPerl module for Bio::DB::SeqVersion
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Brian Osborne
#
# Copyright Brian Osborne 2006
#
# You may distribute this module under the same terms as Perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::DB::SeqVersion - front end to querying databases for identifier 
versions

=head1 SYNOPSIS

  use Bio::DB::SeqVersion;

  my $query = Bio::DB::SeqVersion->new(-type => 'gi');

  my @all_gis = $query->get_all(2);

  my $live_gi = $query->get_recent(2);

=head1 DESCRIPTION

The default type is 'gi'.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Brian Osborne

Email bosborne at alum.mit.edu

=head1 CONTRIBUTORS

Torsten Seemann - torsten.seemann AT infotech.monash.edu.au

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::SeqVersion;
use strict;

use base qw(Bio::WebAgent Bio::Root::Root);

# Private class variable

my $DEFAULTIDTYPE = 'gi'; # sub default_id_type()

=head2 new()

 Usage   : my $obj = Bio::DB::SeqVersion->new();
 Function: Create a Bio::DB::SeqVersion object 
 Returns : An instance of Bio::DB::SeqVersion
 Args    : -type      Identifier namespace, default is 'gi' 

=cut

sub new {
  my($class,@args) = @_;

  if( $class =~ /Bio::DB::SeqVersion::\S+/ ) {
    my ($self) = $class->SUPER::new(@args);
    $self->_initialize(@args);
    return $self;
  } 
  else {
    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param; # lowercase keys

    # we delete '-type' so it doesn't get passed to the sub-class constructor
    # note: delete() returns the value of the item deleted (undef if non-existent)
    my $type = lc( delete($param{'-type'}) || $DEFAULTIDTYPE );

    return unless( $class->_load_seqversion_module($type) );
    
    # we pass %param here, not @args, as we have filtered out -type
    return "Bio::DB::SeqVersion::$type"->new(%param);
  }
}

=head2 get_recent()

 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_recent {
  my ($self,@args) = @_;
  $self->throw_not_implemented();
}

=head2 get_all()

 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_all {
	my ($self,@args) = @_;
	$self->throw_not_implemented();
}

=head2 _load_seqversion_module

 Title   : _load_seqversion_module
 Usage   : Used internally
 Function: Loads up a module at run time on demand
 Example :
 Returns :
 Args    : Name of identifier type

=cut

sub _load_seqversion_module {
	my ($self,$db) = @_;
	my $module = "Bio::DB::SeqVersion::" . $db;
	my $ok;

	eval { $ok = $self->_load_module($module) };
	if ( $@ ) {
		print STDERR $@;
		print STDERR <<END;
$self: $module cannot be found
Exception $@
For more information about the Bio::DB::SeqVersion system please see
the Bio::DB::SeqVersion docs.
END
		;
	}
	return $ok;
}

=head2 default_id_type

 Title   : default_id_type
 Usage   : my $type = $self->default_id_type
 Function: Returns default identifier type for this module
 Returns : string
 Args    : none

=cut

sub default_id_type {
    return $DEFAULTIDTYPE;
}

1;
