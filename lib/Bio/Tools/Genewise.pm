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
    my @transcripts = $gene->transcripts;
      foreach my $t(@transcripts){
        my @exons =  $t->exons;
        foreach my $e(@exons){
            print $e->start." ".$e->end."\n";
        }
      }
  }

=head1 DESCRIPTION

This is the parser for the output of Genewise. It takes either a file
handle or a file name and returns a 
Bio::SeqFeature::Gene::GeneStructure object.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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

=head1 AUTHOR - Fugu Team, Jason Stajich 

 Email: fugui@worf.fugu-sg.org
 Email: jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Genewise;
use vars qw($Srctag);
use strict;
use Symbol;

use Bio::Tools::AnalysisResult;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Gene::Exon;
use Bio::SeqFeature::FeaturePair;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::GeneStructure;

use base qw(Bio::Root::Root Bio::Root::IO);
$Srctag = 'genewise';

=head2 new

 Title   : new
 Usage   : $obj->new(-file=>"genewise.out");
           $obj->new(-fh=>\*GW);
 Function: Constructor for genewise wrapper. Takes either a file or filehandle
 Example :
 Returns : Bio::Tools::Genewise object

See L<Bio::Tools::Genewise>

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
 Function: takes start and end values, swap them if start>end and 
           returns end
 Example :
 Returns :$start,$end,$strand

=cut

sub _get_strand {
  my ($self,$start,$end) = @_;
  defined($start) || $self->throw("Need a start");
  defined($end)   || $self->throw("Need an end");
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

=head2 _score

 Title   : _score
 Usage   : $obj->_score
 Function: get/set for score info
 Returns : a score value

=cut

sub _score {
    my $self = shift;
    return $self->{'_score'} = shift if @_;
    return $self->{'_score'};
}

=head2 _prot_id

 Title   : _prot_id
 Usage   : $obj->_prot_id
 Function: get/set for protein id 
 Returns :a protein id

=cut

sub _prot_id {
    my $self = shift;
    return $self->{'_prot_id'} = shift if @_;
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
    my $self = shift;
    return $self->{'_target_id'} = shift if @_;
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

See L<Bio::SeqFeature::Gene::GeneStructure>

=cut


sub next_prediction {
    my ($self) = @_;

    unless ( $self->parsed ){
	$self->_parse_genes;
	$self->parsed(1);
    }
    return shift @{$self->{'_genes'}};
}

sub parsed {
    my $self = shift;
    return $self->{'_parsed'} = 1 if @_ && $_[0]; 
    return $self->{'_parsed'};
}
  
sub _parse_genes {
    my ($self) = @_;
    my (@alignments,@genes);
    local ($/) = "//";

    while ( defined($_ = $self->_readline) ) {
	$self->debug( $_ );
	while( /Alignment\s+(\d+)\s+Score\s+(\S+)\s+\(Bits\)/g ) {
	    $alignments[$1] = $2;
	}
	if( /Score\s+(\-?\d+(\.\d+)?)/ ) {
	    $self->_score($1);# unless defined $self->_score;    
	}

	if( /Query\s+(?:protein|model)\:\s+(\S+)/ ) {
	    $self->_prot_id($1); #unless defined $self->_prot_id;
	} 

	if( /Target Sequence\s+(\S+)/ ) {	
	    $self->_target_id($1);# unless defined $self->_target_id;
	}	
	
	next unless /Gene\s+\d+\n/;
	my @genes_txt = split(/Gene\s+\d+\n/,$_);
	shift @genes_txt;	#remove first empty entry
	my $gene_num = 1;
	foreach my $gene_txt (@genes_txt) {
	    my $score = $alignments[$gene_num];
	    # If genewise has assigned a strand to the gene as a whole
	    # overall gene start and end
	    my ($g_start, $g_end, $type) = 
		$gene_txt =~ m/Gene\s+
		(\d+)[\s-]+	# start (1-based)
		(\d+)\s+	# end
		(?:\[(\w+)\])?	# 
		/x;
	    my $g_strand;
	    my $source_tag = $type ? "$Srctag". "_$type" : $Srctag;
	    my $genes = Bio::SeqFeature::Gene::GeneStructure->new
		(-source => $source_tag,
		 -score  => $self->_score,
		 );
	    $genes->add_tag_value('ID', $self->_prot_id.".gene");
	    my $transcript = Bio::SeqFeature::Gene::Transcript->new
		(-source => $source_tag,
		 -score  => $score);
	    ($g_start, $g_end, $g_strand) = $self->_get_strand($g_start,$g_end);
	    $genes->strand($g_strand);

	    # grab exon + supporting feature info
	    my @exons;
	    unless ( @exons = $gene_txt =~ m/(Exon .+\s+Supporting .+)/g ) {
		@exons = $gene_txt =~ m/(Exon .+\s+)/g;
	    }
	    my $nbr = 1;

	    # loop through each exon-supporting feature pair
	    foreach my $e (@exons){
		my ($e_start,$e_end,$phase) = 
		    $e =~ m/Exon\s+
		    (\d+)[\s-]+	# start (1 based)
		    (\d+)\s+	# end
		    phase\s+(\d+) # phase
		    /x;
		my $e_strand;
		($e_start,$e_end,$e_strand) = $self->_get_strand($e_start,
								 $e_end);
		$transcript->strand($e_strand) unless $transcript->strand != 0;
		my $exon = Bio::SeqFeature::Gene::Exon->new
		    (-seq_id =>$self->_target_id,
		     -source => $source_tag,
		     -start  =>$e_start, 
		     -end    =>$e_end, 
		     -score  => $score,
		     #-frame => $phase,
		     -strand =>$e_strand);

		$exon->add_tag_value('phase',$phase);
		$exon->is_coding(1);
		if( $self->_prot_id ) {
		    $exon->add_tag_value('Parent',$self->_prot_id);
		}
		$exon->add_tag_value("Exon",$nbr++);
		if( $e =~ m/Supporting\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
		    my ($geno_start,$geno_end,
			$prot_start, $prot_end) = ($1,$2,$3,$4);
		    my $prot_strand;
		    ($prot_start,$prot_end,
		     $prot_strand) = $self->_get_strand($prot_start,$prot_end);
		    my $pf = Bio::SeqFeature::Generic->new
			( -start   => $prot_start, 
			  -end     => $prot_end,
			  -seq_id  => $self->_prot_id,
			  -score   => $score,
			  -strand  => $prot_strand,
			  -source  => $source_tag,
			  -primary_tag => 'supporting_protein_feature',);
		    my $geno_strand;
		    ($geno_start,$geno_end,
		     $geno_strand) = $self->_get_strand($geno_start,$geno_end);
		    my $gf = Bio::SeqFeature::Generic->new 
			( -start   => $geno_start,
			  -end     => $geno_end,
			  -seq_id  => $self->_target_id,
			  -score   => $score,
			  -strand  => $geno_strand,
			  -source  => $source_tag,
			  -primary_tag => 'supporting_genomic_feature',);
		    my $fp = Bio::SeqFeature::FeaturePair->new
			(-feature1 =>$gf,
			 -feature2 =>$pf);
		    $exon->add_tag_value( 'supporting_feature',$fp );
		    if( $self->_prot_id ) {
			$exon->add_tag_value('Target','Protein:'.$self->_prot_id);
			$exon->add_tag_value('Target',$prot_start);
			$exon->add_tag_value('Target',$prot_end);
		    }
		}
		$transcript->add_exon($exon);
	    }
	    $transcript->seq_id($self->_target_id);
	    $transcript->add_tag_value('ID', $self->_prot_id);
	    $transcript->add_tag_value('Parent', $self->_prot_id.".gene");
	    $genes->add_transcript($transcript);
	    $genes->seq_id($self->_target_id);
	    push @genes, $genes;
	    $gene_num++;
	}
    }
    $self->{'_genes'} = \@genes;
}

1;
