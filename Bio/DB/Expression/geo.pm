=head1 NAME

Bio::DB::Expression::geo - *** DESCRIPTION of Class

=head1 SYNOPSIS

*** Give standard usage here

=head1 DESCRIPTION

*** Describe the object here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

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

=head1 AUTHOR

Allen Day E<lt>allenday@ucla.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a '_'.  Methods are
in alphabetical order for the most part.

=cut


# Let the code begin...

package Bio::DB::Expression::geo;
use strict;
use base qw(Bio::DB::Expression);

use Bio::Expression::Contact;
use Bio::Expression::DataSet;
use Bio::Expression::Platform;
use Bio::Expression::Sample;

use constant URL_PLATFORMS => 'http://www.ncbi.nlm.nih.gov/geo/query/browse.cgi?pgsize=100000&mode=platforms&submitter=-1&filteron=0&filtervalue=-1&private=1&sorton=pub_date&sortdir=1&start=1';
use constant URL_PLATFORM => 'http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?form=text&view=full&acc=';
use constant URL_DATASET => 'http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?form=text&view=full&acc=';
use constant URL_SAMPLE => 'http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?form=text&view=full&acc=';

=head2 _initialize()

 Usage   : $obj->_initialize(%arg);
 Function: Internal method to initialize a new Bio::DB::Expression::geo object
 Returns : true on success
 Args    : Arguments passed to new()

=cut

sub _initialize {
  my($self,%arg) = @_;

  foreach my $arg (keys %arg){
    my $marg = $arg;
    $marg =~ s/^-//;
    $self->$marg($arg{$arg}) if $self->can($marg);
  }

  return 1;
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

  my $doc = $self->_get_url( URL_PLATFORMS );
  $doc =~ s!^.+?>Release date<.+?</tr>(.+)</table>!$1!gs;

  my @platforms = ();
  my @records = split m!</tr>\s+<tr>!, $doc;

  foreach my $record ( @records ) {
    my ($platform_acc,$name,$tax_acc,$contact_acc,$contact_name) =
      $record =~ m!acc\.cgi\?acc=(.+?)".+?<td.+?>(.+?)<.+?<td.+?>.+?<.+?<td.+?>.+?href=".+?id=(.+?)".+?<td.+?OpenSubmitter\((\d+?)\).+?>(.+?)<!s;
    next unless $platform_acc;

    my $platform = Bio::Expression::Platform->new(
                                                  -accession => $platform_acc,
                                                  -name => $name,
                                                  -_taxon_id => $tax_acc,
                                                  -contact => Bio::Expression::Contact->new(
                                                                                            -source => 'geo',
                                                                                            -accession => $contact_acc,
                                                                                            -name => $contact_name,
                                                                                            -db => $self
                                                                                           ),
                                                  -db => $self,
                                                 );
    push @platforms, $platform;
  }

  return @platforms;
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

 Usage   : $db->get_datasets('accession');
 Function:
 Example :
 Returns : a list of Bio::Expression::DataSet objects
 Args    :

=cut

sub get_datasets {
  my ($self,$platform) = @_;

  my @lines = split /\n/, $self->_get_url( URL_PLATFORM . $platform->accession );

  my @datasets = ();

  foreach my $line ( @lines ) {
    my ($dataset_acc) = $line =~ /^\!Platform_series_id = (\S+?)\s*$/;
    next unless $dataset_acc;

    my $dataset = Bio::Expression::DataSet->new(
                                                -accession => $dataset_acc,
                                                -platform => $platform,
                                                -db => $self,
                                               );

    push @datasets, $dataset;
  }

  return @datasets;
}

sub fill_sample {
  my ( $self, $sample ) = @_;

  my @lines = split /\n/, $self->_get_url( URL_SAMPLE. $sample->accession );

  foreach my $line ( @lines ) {
    if ( my ($name) = $line =~ /^\!Sample_title = (.+?)\s*$/ ) {
      $sample->name( $name );
    }
    elsif ( my ($desc) = $line =~ /^\!Sample_characteristics.*? = (.+?)\s*$/ ) {
      $sample->description( $desc );
    }
    elsif ( my ($source_name) = $line =~ /^\!Sample_source_name.*? = (.+?)\s*$/ ) {
      $sample->source_name( $source_name );
    }
    elsif ( my ($treatment_desc) = $line =~ /^\!Sample_treatment_protocol.*? = (.+?)\s*$/ ) {
      $sample->treatment_description( $treatment_desc );
    }
  }
  return 1;
}

sub fill_dataset {
  my ( $self, $dataset ) = @_;

  my @lines = split /\n/, $self->_get_url( URL_DATASET . $dataset->accession );

  my @samples = ();

  foreach my $line ( @lines ) {
    if ( my ($sample_acc) = $line =~ /^\!Series_sample_id = (\S+?)\s*$/ ) {
      my $sample = Bio::Expression::Sample->new(
                                                -accession => $sample_acc,
                                                -dataset => $dataset,
                                                -db => $self,
                                               );
      push @samples, $sample;
    }
    elsif ( my ($pubmed_acc) = $line =~ /^\!Series_pubmed_id = (.+?)\s*$/ ) {
      $dataset->pubmed_id( $pubmed_acc );
    }
    elsif ( my ($web_link) = $line =~ /^\!Series_web_link = (.+?)\s*$/ ) {
      $dataset->web_link( $web_link );
    }
    elsif ( my ($contact) = $line =~ /^\!Series_contact_name = (.+?)\s*$/ ) {
      $dataset->contact( $contact );
    }
    elsif ( my ($name) = $line =~ /^\!Series_title = (.+?)\s*$/ ) {
      $dataset->name( $name );
    }
    elsif ( my ($desc) = $line =~ /^\!Series_summary = (.+?)\s*$/ ) {
      $dataset->description( $desc );
    }
    elsif ( my ($design) = $line =~ /^\!Series_type = (.+?)\s*$/ ) {
      $dataset->design( $design );
    }
    elsif ( my ($design_desc) = $line =~ /^\!Series_overall_design = (.+?)\s*$/ ) {
      $dataset->design_description( $design_desc );
    }
  }

  $dataset->samples(\@samples);
}

#################################################

=head2 _platforms_doc()

 Usage   :
 Function:
 Example :
 Returns : an HTML document containing a table of all platforms
 Args    :


=cut

sub _get_url {
  my ($self,$url) = @_;

  my $response;
  eval {
    $response = $self->get( $url );
  };
  if( $@ ) {
    $self->warn("Can't query website: $@");
    return;
  }
  $self->debug( "resp is $response\n"); 

  return $response;
}


1;
