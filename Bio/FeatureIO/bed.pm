=pod

=head1 NAME

Bio::FeatureIO::bed - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

http://www.genome.ucsc.edu/goldenPath/help/customTrack.html#BED

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
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Allen Day

Email allenday@ucla.edu

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::FeatureIO::bed;

use strict;
use base qw(Bio::FeatureIO);
use Bio::SeqFeature::Annotated;
use Bio::OntologyIO;

sub _initialize {
  my($self,%arg) = @_;

  $self->SUPER::_initialize(%arg);

  #read headers
  my $directive;
  while(($directive = $self->_readline()) && $directive =~ /^##/){
    $self->_handle_directive($directive)
  }
  $self->_pushback($directive);

  $self->so(
            Bio::Ontology::OntologyStore->get_ontology('Sequence Ontology')
           );
}

sub write_feature {
  my($self,$feature) = @_;
  $self->throw("only Bio::SeqFeature::Annotated objects are writeable") unless $feature->isa('Bio::SeqFeature::Annotated');

  my $chrom       = $feature->seq_id || '';
  my $chrom_start = $feature->start  || 0;
  my $chrom_end   = $feature->stop   || 0;

  #try to make a reasonable name
  my $name        = undef;
  if(my @v = ($feature->annotation->get_Annotations('Name'))){
    $name = $v[0];
    $self->warn("only using first of feature's multiple names: ".join ',', map {$_->value} @v) if scalar(@v) > 1;
  } elsif(my @v = ($feature->annotation->get_Annotations('ID'))){
    $name = $v[0];
    $self->warn("only using first of feature's multiple IDs: ".join ',', map {$_->value} @v) if scalar(@v) > 1;
  } else {
    $name = 'anonymous';
  }

  my $score = $feature->score || 0;
  my $strand = $feature->strand == 0 ? '-' : '+'; #default to +
  my $thick_start = '';  #not implemented, used for CDS
  my $thick_end = '';    #not implemented, used for CDS
  my $reserved = 0;
  my $block_count = '';  #not implemented, used for sub features
  my $block_sizes = '';  #not implemented, used for sub features
  my $block_starts = ''; #not implemented, used for sub features

  print join("\t",($chrom,$chrom_start,$chrom_end,$name,$score,$strand,$thick_start,$thick_end,$reserved,$block_count,$block_sizes, $block_starts)),"\n";
  $self->write_feature($_) foreach $feature->get_SeqFeatures();
}

sub next_feature {
  shift->throw_not_implemented();
}

1;
