# $Id$
#
# BioPerl module for Bio::Tools::Genewise
#
# Copyright Fugu Team 
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Genewise - Results of one Genewise run

=head1 SYNOPSIS

  use Bio::Tools::Genewise;
  my $gw = Bio::Tools::Genewise(-file=>"genewise.out");

  while (my $gene = $gw->next_prediction){
      my @transcripts = $gw->transcripts;
      foreach my $t(@transcripts){
        my @exons =  $t->exons;
        foreach my $e(@exons){
            print $e->start." ".$e->end."\n";
        }
      }
  }

=head1 DESCRIPTION

This is the parser for the output of Genewise. It takes either a file handle or
a file name and returns a Bio::SeqFeature::Gene::GeneStructure object.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Fugu Team 

 Email: fugui@worf.fugu-sg.org

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Genewise;
use vars qw(@ISA $Srctag);
use strict;
use Symbol;

use Bio::Root::Root;
use Bio::Root::IO;
use Bio::Tools::AnalysisResult;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Gene::Exon;
use Bio::SeqFeature::FeaturePair;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::GeneStructure;

@ISA = qw(Bio::Root::Root Bio::Root::IO);
$Srctag = 'genewise';

=head2 new

 Title   : new
 Usage   : $obj->new(-file=>"genewise.out");
           $obj->new(-fh=>\*GW);
 Function: Constructor for genewise wrapper. Takes either a file or filehandle
 Example :
 Returns : L<Bio::Tools::Genewise>

=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->_initialize_io(@args);
  return $self;
}

=head2 _get_strand

 Title   : _get_strand
 Usage   : $obj->_get_strand
 Function: takes start and end values, swap them if start>end and returns end
 Example :
 Returns :$start,$end,$strand

=cut

sub _get_strand {
  my ($self,$start,$end) = @_;
  $start || $self->throw("Need a start");
  $end   || $self->throw("Need an end");
  my $strand;
  if ($start > $end) {
    my $tmp = $start;
    $start = $end;
    $end = $tmp;
    $strand = -1;
  }
  else {
    $strand = 1;
  }
  return ($start,$end,$strand);
}

=head2 score

 Title   : score
 Usage   : $obj->score
 Function: get/set for score info
 Example :
 Returns : a score value

=cut

sub _score {
  my ($self,$val) = @_;
  if($val){
    $self->{'_score'} = $val;
  }
  return $self->{'_score'};
}

=head2 _prot_id

 Title   : _prot_id
 Usage   : $obj->_prot_id
 Function: get/set for protein id 
 Example :
 Returns :a protein id

=cut

sub _prot_id {
  my ($self,$val) = @_;
  if($val){
    $self->{'_prot_id'} = $val;
  }
  return $self->{'_prot_id'};
}

=head2 _target_id

 Title   : _target_id
 Usage   : $obj->_target_id
 Function: get/set for genomic sequence id
 Example :
 Returns :a target id

=cut

sub _target_id {
  my ($self,$val) = @_;
  if($val){
    $self->{'_target_id'} = $val;
  }
  return $self->{'_target_id'};
}


=head2 next_prediction

 Title   : next_prediction
 Usage   : while($gene = $genewise->next_prediction()) {
                  # do something
           }
 Function: Returns the gene structure prediction of the Genewise result
           file. Call this method repeatedly until FALSE is returned.

 Example :
 Returns : a Bio::SeqFeature::Gene::GeneStructure object
 Args    :

=cut


sub next_prediction {
    my ($self) = @_;
    my $genes = new Bio::SeqFeature::Gene::GeneStructure(-source => $Srctag);
    my $transcript = new Bio::SeqFeature::Gene::Transcript(-source => $Srctag);
    $/ = "//";
    my $score;
    my $prot_id;
    my $target_id;
    while ($_ = $self->_readline) {
	$self->debug( $_ ) if( $self->verbose > 0);
        ($score) = $_=~m/Score\s+(\d+[\.][\d]+)/;
        $self->_score($score) unless defined $self->_score;
        ($prot_id) = $_=~m/Query protein:\s+([\w\.]+)/;
        $self->_prot_id($prot_id) unless defined $self->_prot_id;
        ($target_id) = $_=~m/Target Sequence\s+([\w\.]+)/;	
        $self->_target_id($target_id) unless defined $self->_target_id;
        next unless /Gene\s+\d+\n/;

        #grab exon + supporting feature info
        my @exons;
	
	unless ( @exons = $_ =~ m/(Exon .+\s+Supporting .+)/g ) {
	    @exons = $_ =~ m/(Exon .+\s+)/g;

	}
        my $nbr = 1;

        #loop through each exon-supporting feature pair
        foreach my $e (@exons){
	    my ($e_start,$e_end,$phase) = $e =~ m/Exon\s+(\d+)\s+(\d+)\s+phase\s+(\d+)/;
	    my $e_strand;
	    ($e_start,$e_end,$e_strand) = $self->_get_strand($e_start,$e_end);
	    $transcript->strand($e_strand) unless $transcript->strand != 0;
	    
	    my $exon = new Bio::SeqFeature::Gene::Exon
		(-seq_id=>$self->_target_id,
		 -source => $Srctag,
		 -start=>$e_start, 
		 -end=>$e_end, 
		 #-frame => $phase,
		 -strand=>$e_strand);
	    $exon->add_tag_value('phase',$phase);
	    if( $self->_prot_id ) {
		$exon->add_tag_value('Sequence',"Protein:".$self->_prot_id);
	    }
	    $exon->add_tag_value("Exon",$nbr++);

	    if( $e =~ m/Supporting\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
		my ($geno_start,$geno_end,
		    $prot_start,
		    $prot_end) = ($1,$2,$3,$4);
          
		my $prot_strand;
		($prot_start,$prot_end,
		 $prot_strand) = $self->_get_strand($prot_start,$prot_end);
		
		my $pf = new Bio::SeqFeature::Generic
		    ( -start   => $prot_start,
		      -end     => $prot_end,
		      -seq_id  => $self->_prot_id,
		      -score   => $self->_score,
		      -strand  => $prot_strand,
		      -source  => $Srctag,
		      -primary=> 'supporting_protein_feature',
		      );
		my $geno_strand;
		($geno_start,$geno_end,
		 $geno_strand) = $self->_get_strand($geno_start,$geno_end);
		my $gf = new Bio::SeqFeature::Generic
		    ( -start   => $geno_start,
		      -end     => $geno_end,
		      -seq_id  => $self->_target_id,
		      -score   => $self->_score,
		      -strand  => $geno_strand,
		      -source  => $Srctag,
		      -primary => 'supporting_genomic_feature',
		      );
		my $fp = new Bio::SeqFeature::FeaturePair(-feature1=>$gf,
							  -feature2=>$pf);
	  
		$exon->add_tag_value( 'supporting_feature' => $fp );
	    }
          $transcript->add_exon($exon);
        }
	$transcript->seq_id($self->_target_id);
	$genes->add_transcript($transcript);
	$genes->seq_id($self->_target_id);
	return $genes;
    }
}
1;
