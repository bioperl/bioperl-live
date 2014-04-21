#
# BioPerl module for Bio::DB::Expression
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Allen Day <allenday@ucla.edu>
#
# Copyright Allen Day
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Expression - DESCRIPTION of Object

=head1 SYNOPSIS

  use Bio::DB::Expression;
  my $db = Bio::DB::Expression->new( -source => 'geo' );
  my @platforms = $db->get_platforms();
  foreach my $platform ( @platforms ) {
    my @datasets = $platform->get_datasets();
    foreach my $dataset ( @datasets ) {
      my @samples = $dataset->get_samples();
      foreach my $sample ( @samples ) {
        #...
      }
    }
  }

=head1 DESCRIPTION

Describe the object here

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

=head1 AUTHOR - Allen Day

Email allenday@ucla.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::DB::Expression;
use strict;
use base qw(Bio::Root::HTTPget Bio::Root::Root);
use Bio::Root::HTTPget;
our $DefaultSource = 'geo';

=head2 new()

 Usage   : my $obj = Bio::DB::Expression->new();
 Function: Builds a new Bio::DB::Expression object 
 Returns : an instance of Bio::DB::Expression
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  if( $class =~ /Bio::DB::Expression::(\S+)/ ) {
    my ($self) = $class->SUPER::new(@args);
    $self->_initialize(@args);
    return $self;
  } else {
    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
    my $source = $param{'-source'} || $DefaultSource;

    $source = "\L$source";	# normalize capitalization to lower case

    # normalize capitalization
    return unless( $class->_load_expression_module($source) );
    return "Bio::DB::Expression::$source"->new(@args);
  }
}

=head2 get_platforms()

 Usage   :
 Function:
 Example :
 Returns : a list of Bio::Expression::Platform objects
 Args    :

=cut

sub get_platforms {
  my ($self,@args) = @_;
  $self->throw_not_implemented();
}

=head2 get_samples()

 Usage   :
 Function:
 Example :
 Returns : a list of Bio::Expression::Sample objects
 Args    :

=cut

sub get_samples {
  my ($self,@args) = @_;
  $self->throw_not_implemented();
}

=head2 get_contacts()

 Usage   :
 Function:
 Example :
 Returns : a list of Bio::Expression::Contact objects
 Args    :

=cut

sub get_contacts {
  my ($self,@args) = @_;
  $self->throw_not_implemented();
}

=head2 get_datasets()

 Usage   :
 Function:
 Example :
 Returns : a list of Bio::Expression::DataSet objects
 Args    :

=cut

sub get_datasets {
  my ($self,@args) = @_;
  $self->throw_not_implemented();
}




=head2 _load_expression_module

 Title   : _load_expression_module
 Usage   : *INTERNAL Bio::DB::Expression stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_expression_module {
    my ($self, $source) = @_;
    my $module = "Bio::DB::Expression::" . $source;
    my $ok;

    eval { $ok = $self->_load_module($module) };
    if ( $@ ) {
	print STDERR $@;
	print STDERR <<END;
$self: $source cannot be found
Exception $@
For more information about the Bio::DB::Expression system please see
the Bio::DB::Expression docs.  This includes ways of checking for 
formats at compile time, not run time.
END
  ;
    }
    return $ok;
}


1;
