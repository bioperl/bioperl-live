
#
# BioPerl module for Bio::Tools::EPCR
#
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::EPCR - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org            - General discussion
http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@chg.mc.duke.edu

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::EPCR;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::SeqAnalysisParserI;
use Bio::SeqFeature::FeaturePair;
use Bio::SeqFeature::Generic;

@ISA = qw(Bio::SeqAnalysisParserI Bio::Root::IO Bio::Root::RootI);

=head2 new

 Title   : new
 Usage   : my $epcr = new Bio::Tools::EPCR(-file => $file);
 Function: Initializes a new EPCR parser
 Returns : Bio::Tools::EPCR
 Args    : -fh   => filehandle
           OR
           -file => filename

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->_initialize_io(@args);
  
  return $self;
}

=head2 next_feature

 Title   : next_feature
 Usage   : $seqfeature = $obj->next_feature();
 Function: Returns the next feature available in the analysis result, or
           undef if there are no more features.
 Example :
 Returns : A Bio::SeqFeatureI implementing object, or undef if there are no
           more features.
 Args    : none    

=cut

sub next_feature {
    my ($self) = @_;
    my $line = $self->_readline;
    return undef unless defined($line);
    chomp($line);
    my($seqname,$location,$mkrname, $rest) = split(/\s+/,$line,4);
    
    my ($start,$end) = ($location =~ /(\S+)\.\.(\S+)/);

    # If we require that e-PCR is run with D=1 we can detect a strand
    # for now hardcoded to 0

    my $strand = 0;
    my $markerfeature = new Bio::SeqFeature::Generic ( '-start'   => $start,
						       '-end'     => $end,
						       '-strand'  => $strand,
						       '-source'  => 'e-PCR',
						       '-primary' => 'sts',
						       '-seqname' => $seqname,
						       '-tag'     => {
							   'name'=> $mkrname,
							   'note' => $rest,
						       });
    return $markerfeature;
}

1;
