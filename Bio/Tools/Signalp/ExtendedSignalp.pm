#
# BioPerl module for Bio::Tools::Signalp::ExtendedSignalp
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>
#
# Copyright Emmanuel Quevillon
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Signalp::ExtendedSignalp - enhanced parser for Signalp output

=head1 SYNOPSIS

 use Bio::Tools::Signalp::ExtendedSignalp;
 my $params = [qw(maxC maxY maxS meanS D)];
 my $parser = new Bio::Tools::Signalp::ExtendedSignalp(
                                                       -fh      => $filehandle
                                                       -factors => $params
                                                      );

 $parser->factors($params);
 while( my $sp_feat = $parser->next_feature ) {
       #do something
       #eg
       push @sp_feat, $sp_feat;
 }

=head1 DESCRIPTION

# Please direct questions and support issues to I<bioperl-l@bioperl.org> 

Parser module for Signalp.

Based on the EnsEMBL module Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp
originally written by Marc Sohrmann (ms2 a sanger.ac.uk) Written in BioPipe by
Balamurugan Kumarasamy (savikalpa a fugu-sg.org) Cared for by the Fugu
Informatics team (fuguteam@fugu-sg.org)

You may distribute this module under the same terms as perl itself

Compared to the original SignalP, this method allow the user to filter results
out based on maxC maxY maxS meanS and D factor cutoff for the Neural Network (NN)
method only. The HMM method does not give any filters with 'YES' or 'NO' as result.

The user must be aware that the filters can only by applied on NN method.
Also, to ensure the compatibility with original Signalp parsing module, the user
must know that by default, if filters are empty, max Y and mean S filters are
automatically used to filter results.

If the used gives a list, then the parser will only report protein having 'YES'
for each factor.

This module supports parsing for full, summary and short output form signalp.
Actually, full and summary are equivalent in terms of filtering results.

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

=head1 AUTHOR

 Based on the Bio::Tools::Signalp module
 Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are usually preceded with a _

=cut

package Bio::Tools::Signalp::ExtendedSignalp;

use strict;
use Data::Dumper;
use Bio::SeqFeature::Generic;
# don't need Bio::Root::Root/IO (already in inheritance tree)
use base qw(Bio::Tools::Signalp Bio::Tools::AnalysisResult);

#Supported arguments
my $FACTS = {
	     'maxC'  => 1,
	     'maxS'  => 1,
	     'maxY'  => 1,
	     'meanS' => 1,
	     'D'     => 1,
	    };

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Signalp::ExtendedSignalp();
 Function: Builds a new Bio::Tools::Signalp::ExtendedSignalp object
 Returns : Bio::Tools::Signalp::ExtendedSignalp
 Args    : -fh/-file => $val, # for initing input, see Bio::Root::IO


=cut

sub new {
      my($class,@args) = @_;

      my $self = $class->SUPER::new(@args);
      $self->_initialize_io(@args);

      my $factors = $self->_rearrange([qw(FACTORS)], @args);
      #To behave like the parent module (Bio::Tools::Signalp) we default factors to these two factors
      if($factors && scalar(@$factors)){
	  $factors = $factors;
      }
      else{
	  $factors = [qw(maxY meanS)];
      }
      $factors && $self->factors($factors);
						
      return $self;
}

=head2 next_feature

 Title   : next_feature
 Usage   : my $feat = $signalp->next_feature
 Function: Get the next result feature from parser data
 Returns : Bio::SeqFeature::Generic
 Args    : none


=cut

sub next_feature {

    my ($self) = @_;

    if(!$self->_parsed()){
	$self->_parse();
    }

    return shift @{$self->{_features}} || undef;

}

=head2 _filterok

 Title   : _filterok
 Usage   : my $feat = $signalp->_filterok
 Function: Check if the factors required by the user are all ok.
 Returns : 1/0
 Args    : hash reference


=cut

sub _filterok {

    my($self, $hash) = @_;

    #We hope everything will be fine ;)
    my $bool = 1;

    #If the user did not give any filter, we keep eveything
    return $bool unless keys %{$self->{_factors}};

    #If only one of the factors parsed is equal to NO based on the user factors cutoff
    #Then the filter is not ok.
    foreach my $fact (keys %{$self->factors()}){
	if(exists($hash->{$fact}) && $hash->{$fact} =~ /^N/){
	    $bool = 0;
	}
    }

    return $bool;

}

=head2 factors

 Title   : factors
 Usage   : my $feat = $signalp->factors
 Function: Get/Set the filters required from the user
 Returns : hash
 Args    : array reference


=cut

sub factors {

    my($self, $array) = @_;

    if($array){
	$self->{_factors} = { };
	foreach my $f (@$array){
	    if(exists($FACTS->{$f})){
		$self->{_factors}->{$f} = 1;
	    }
	    else{
		$self->throw("[$f] incorrect factor. Supported:\n- ".join("\n- ", keys %$FACTS)."\n");
	    }
	}
    }

    return $self->{_factors};

}

=head2 _parsed

 Title   : _parsed
 Usage   : obj->_parsed()
 Function: Get/Set if the result is parsed or not
 Returns : 1/0 scalar
 Args    : On set 1


=cut

sub _parsed {

    my($self, $parsed) = @_;

    if(defined($parsed)){
	$self->{_parsed} = $parsed;
    }

    return $self->{_parsed};

}

=head2 _parse

 Title   : _parse
 Usage   : obj->_parse
 Function: Parse the SignalP result
 Returns :
 Args    :


=cut

sub _parse {

    my($self) = @_;

    #Let's read the file...
    while (my $line = $self->_readline()) {

	chomp $line;
	#We want to be sure to catch the first non empty line to be ablte to determine
	#which format we are working with...
	next unless ($line =~ /^>(\S+)|^# SignalP-[NHM]+ \S+ predictions/);

	if($line =~ /^>(\S+)/){
	    $self->_pushback($line);
	    $self->_parse_summary_format();
	    last;
	}
	elsif($line =~ /^# SignalP-[NHM]+ \S+ predictions/){
	    $self->_pushback($line);
	    $self->_parse_short_format();
	    last;
	}
	else{
	    $self->throw("Unable to determine the format type.");
	}
    }

    return;
}

=head2 _parse_summary_format

 Title   : _parse_summary_format
 Usage   : $self->_parse_summary_format
 Function: Method to parse summary/full format from signalp output
           It automatically fills filtered features.
 Returns :
 Args    :

=cut

sub _parse_summary_format {

    my($self) = @_;

    my $feature = undef;
    my $ok = 0;

    while(my $line = $self->_readline()){

	if($line =~ /^SignalP-NN result:/){
	    $self->_pushback($line);
	    $feature = $self->_parse_nn_result($feature);
	}
	if($line =~ /^SignalP-HMM result:/){
	    $self->_pushback($line);
	    $feature = $self->_parse_hmm_result($feature);
	}

	if($line =~ /^---------/ && $feature){
	    my $new_feature = $self->create_feature($feature);
	    push @{$self->{_features}}, $new_feature if $new_feature;
	    $feature = undef;
	}
    }

    return;
}


=head2 _parse_nn_result

 Title   : _parse_nn_result
 Usage   : obj->_parse_nn_result
 Function: Parses the Neuronal Network (NN) part of the result
 Returns : Hash reference
 Args    :


=cut

sub _parse_nn_result {

    my($self, $feature) = @_;

    my $ok   = 0;
    my %facts;

    #SignalP-NN result:
    #>MGG_11635.5           length = 100
    ## Measure  Position  Value  Cutoff  signal peptide?
    #  max. C    37       0.087   0.32   NO
    #  max. Y    37       0.042   0.33   NO
    #  max. S     3       0.062   0.87   NO
    #  mean S     1-36    0.024   0.48   NO
    #       D     1-36    0.033   0.43   NO

    while(my $line = $self->_readline()){

	chomp $line;

	if($line =~ /^SignalP-NN result:/){
	    $ok = 1;
	    next;
	}

	$self->throw("Wrong line for parsing NN results.") unless $ok;

	if ($line=~/^\>(\S+)\s+length/) {
	    $self->seqname($1);
	    %facts = ();
	    next;
	}
	elsif($line =~ /max\.\s+C\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
	    $feature->{maxCprob} = $1;
	    $facts{maxC} = $2;
	    next;
	}
	elsif ($line =~ /max\.\s+Y\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
	    $feature->{maxYprob} = $1;
	    $facts{maxY} = $2;
	    next;
	}
	elsif($line =~ /max\.\s+S\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
	    $feature->{maxSprob} = $1;
	    $facts{maxS} = $2;
	    next;
	}
	elsif ($line=~/mean\s+S\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
	    $feature->{meanSprob} = $1;
	    $facts{meanS} = $2;
	    next;
	}
	elsif ($line=~/\s+D\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
	    $feature->{Dprob} = $1;
	    $facts{D} = $2;
	    next;
	}
	#If we don't have this line it means that all the factors cutoff are equal to 'NO'
	elsif ($line =~ /Most likely cleavage site between pos\.\s+(\d+)/) {
	    #if($self->_filterok(\%facts)){
		#$feature->{name}       = $self->seqname();
		#$feature->{start}      = 1;
		$feature->{end}        = $1 + 1; #To be consistent with end given in short format
	    #}
	    #return $feature;
	}
	elsif($line =~ /^\s*$/){
	    last;
	}
    }

    if($self->_filterok(\%facts)){
	$feature->{name}  = $self->seqname();
	$feature->{start} = 1;
	$feature->{nnPrediction} = 'signal-peptide';
    }

    return $feature;
}


=head2 _parse_hmm_result

 Title   : _parse_hmm_result
 Usage   : obj->_parse_hmm_result
 Function: Parses the Hiden Markov Model (HMM) part of the result
 Returns : Hash reference
 Args    :

=cut

sub _parse_hmm_result {

    my ($self, $feature_hash) = @_;

    my $ok = 0;

    #SignalP-HMM result:
    #>MGG_11635.5
    #Prediction: Non-secretory protein
    #Signal peptide probability: 0.000
    #Signal anchor probability: 0.000
    #Max cleavage site probability: 0.000 between pos. -1 and  0

    while(my $line = $self->_readline()){

	chomp $line;
	next if $line =~ /^\s*$/o;

	if($line =~ /^SignalP-HMM result:/){
	    $ok = 1;
	    next;
	}

	$self->throw("Wrong line for parsing HMM result.") unless $ok;

	if($line =~ /^>(\S+)/){
	    #In case we already seen a name with NN results
	    $feature_hash->{name} = $1 unless $self->seqname();
	}
        elsif($line =~ /Prediction: (.+)$/){
            $feature_hash->{hmmPrediction} = $1;
        }
	elsif($line =~ /Signal peptide probability: ([0-9\.]+)/){
            $feature_hash->{peptideProb} = $1;
        }
	elsif($line =~ /Signal anchor probability: ([0-9\.]+)/){
            $feature_hash->{anchorProb} = $1;
        }
	elsif($line =~ /Max cleavage site probability: (\S+) between pos. \S+ and (\S+)/){
	    $feature_hash->{cleavageSiteProb} = $1;
	    #Strange case, if we don't have an end value in NN result (no nn method launched)
	    #We try anyway to get an end value, unless this value is lower than 1 which is
	    #the start
	    $feature_hash->{end}   = $2 if($2 > 1 && !$feature_hash->{end});
	    $feature_hash->{start} = 1 unless $feature_hash->{start};
	    last;
	}
    }

    return $feature_hash;
}

=head2 _parse_short_format

 Title   : _parse_short_format
 Usage   : $self->_parse_short_format
 Function: Method to parse short format from signalp output
           It automatically fills filtered features.
 Returns :
 Args    :

=cut

sub _parse_short_format {

                                my($self) = @_;

    my $ok = 0;
    my $method = undef;
    $self->{_oformat} = 'short';

    #Output example
    # SignalP-NN euk predictions                                   	                # SignalP-HMM euk predictions
    # name                Cmax  pos ?  Ymax  pos ?  Smax  pos ?  Smean ?  D     ? 	# name      !  Cmax  pos ?  Sprob ?
    #Q5A8M1_CANAL          0.085  27 N  0.190  35 N  0.936  27 Y  0.418 N  0.304 N	Q5A8M1_CANAL  Q  0.001  35 N  0.002 N
    #O74127_YARLI          0.121  21 N  0.284  21 N  0.953  11 Y  0.826 Y  0.555 Y	O74127_YARLI  S  0.485  23 N  0.668 Y
    #Q5VJ86_9PEZI          0.355  24 Y  0.375  24 Y  0.798  12 N  0.447 N  0.411 N	Q5VJ86_9PEZI  Q  0.180  23 N  0.339 N
    #Q5A8U5_CANAL          0.085  27 N  0.190  35 N  0.936  27 Y  0.418 N  0.304 N	Q5A8U5_CANAL  Q  0.001  35 N  0.002 N

    while(my $line = $self->_readline()){
	
	chomp $line;
	next if $line =~ /^\s*$|^# name/;

	if($line =~ /^#/){
	    $method = $line =~ /SignalP-NN .+ SignalP-HMM/ ?
	                                            'both' : $line =~ /SignalP-NN/ ?
							                      'nn' : 'hmm';
	    next;
	}

	#$self->throw("It looks like the format is not 'short' format.") unless($ok);

	my @data = split(/\s+/, $line);
	$self->seqname($data[0]);

	my $factors = { };
	my $feature = { };

	#NN results gives more fields than HMM
	if($method eq 'both' || $method eq 'nn'){

	    $feature->{maxCprob} = $data[1];
	    $factors->{maxC}     = $data[3];
	    $feature->{maxYprob} = $data[4];
	    $factors->{maxY}     = $data[6];
	    $feature->{maxSprob} = $data[7];
	    $factors->{maxS}     = $data[9];
	    $feature->{meanSprob}= $data[10];
	    $factors->{meanS}    = $data[11];
	    $feature->{Dprob}    = $data[12];
	    $factors->{D}        = $data[13];
	    #It looks like the max Y position is reported as the most likely cleavage position
	    $feature->{end}      = $data[5];
	    $feature->{nnPrediction} = 'signal-peptide';

	    if($method eq 'both'){
		$feature->{hmmPrediction}    = $data[15] eq 'Q' ? 'Non-secretory protein' : 'Signal peptide';
		$feature->{cleavageSiteProb} = $data[16];
		$feature->{peptideProb}      = $data[19];
	    }
	}
	elsif($method eq 'hmm'){
	    #In short output anchor probability is not given
	    $feature->{hmmPrediction}    = $data[1] eq 'Q' ? 'Non-secretory protein' : 'Signal peptide';
	    $feature->{cleavageSiteProb} = $data[2];
	    $feature->{peptideProb}      = $data[5];
	    #It looks like the max cleavage probability position is given by the Cmax proability
	    $feature->{end} = $data[3];
	}

	#Unfortunately, we cannot parse the filters for hmm method.
	if($self->_filterok($factors)){
	    $feature->{name}        = $self->seqname();
	    $feature->{start}       = 1;
	    $feature->{source}      = 'Signalp';
	    $feature->{primary}     = 'signal_peptide';
	    $feature->{program}     = 'Signalp';
	    $feature->{logic_name}  = 'signal_peptide';

	    my $new_feat = $self->create_feature($feature);
	    push @{$self->{_features}}, $new_feat if $new_feat;
	}
    }

    return;
}

=head2 create_feature

 Title   : create_feature
 Usage   : obj->create_feature(\%feature)
 Function: Internal(not to be used directly)
 Returns :
 Args    :


=cut

sub create_feature {

    my ($self, $feat) = @_;

    #If we don't have neither start nor end, we return.
    unless($feat->{name} && $feat->{start} && $feat->{end}){
	return;
    }

    # create feature object
    my $feature = Bio::SeqFeature::Generic->new(
						-seq_id     => $feat->{name},
						-start      => $feat->{start},
						-end        => $feat->{end},
						-score      => defined($feat->{peptideProb}) ? $feat->{peptideProb} : '',
						-source     => 'Signalp',
						-primary    => 'signal_peptide',
						-logic_name => 'signal_peptide',
					       );

    $feature->add_tag_value('peptideProb', $feat->{peptideProb});
    $feature->add_tag_value('anchorProb', $feat->{anchorProb});
    $feature->add_tag_value('evalue',$feat->{anchorProb});
    $feature->add_tag_value('percent_id','NULL');
    $feature->add_tag_value("hid",$feat->{primary});
    $feature->add_tag_value('signalpPrediction', $feat->{hmmPrediction});
    $feature->add_tag_value('cleavageSiteProb', $feat->{cleavageSiteProb}) if($feat->{cleavageSiteProb});
    $feature->add_tag_value('nnPrediction', $feat->{nnPrediction})         if($feat->{nnPrediction});
    $feature->add_tag_value('maxCprob', $feat->{maxCprob})   if(defined($feat->{maxCprob}));
    $feature->add_tag_value('maxSprob', $feat->{maxSprob})   if(defined($feat->{maxSprob}));
    $feature->add_tag_value('maxYprob', $feat->{maxYprob})   if(defined($feat->{maxYprob}));
    $feature->add_tag_value('meanSprob', $feat->{meanSprob}) if(defined($feat->{meanSprob}));
    $feature->add_tag_value('Dprob', $feat->{Dprob})         if(defined($feat->{Dprob}));

    return $feature;

}

=head2 seqname

 Title   : seqname
 Usage   : obj->seqname($name)
 Function: Internal(not to be used directly)
 Returns :
 Args    :


=cut

sub seqname{
    my ($self,$seqname)=@_;

    if (defined($seqname)){
        $self->{'seqname'} = $seqname;
    }

    return $self->{'seqname'};

}


1;


