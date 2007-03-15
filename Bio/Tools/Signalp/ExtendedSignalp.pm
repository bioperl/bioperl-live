# Parser module for Signalp Bio::Tools::Signalp::ExtendedSignalP
#
# 
# Based on the EnsEMBL module
# Bio::EnsEMBL::Pipeline::Runnable::Protein::Signalp originally
# written by Marc Sohrmann (ms2@sanger.ac.uk) Written in BioPipe by
# Balamurugan Kumarasamy <savikalpa@fugu-sg.org> Cared for by the Fugu
# Informatics team (fuguteam@fugu-sg.org)

# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Signalp::ExtendedSignalP

=head1 SYNOPSIS

 use Bio::Tools::Signalp::ExtendedSignalp;
 my $params = [qw(maxC maxY maxS meanS D)];
 my $parser = new Bio::Tools::Signalp::ExtendedSignalp(-fh =>$filehandle
                                                       -factors => $params);

 #Compare to the original SignalP, this method allow the user to filter results
 #out based on maxC maxY maxS meanS and D cutoff factors.
 #If the user gives a list, then the parser will only report protein having 'YES'
 #for each factor.

 $parser->factors($params);
 while( my $sp_feat = $parser->next_result ) {
       #do something
       #eg
       push @sp_feat, $sp_feat;
 }

=head1 DESCRIPTION

 Parser for Signalp output

=head1 FEEDBACK

=head2 Mailing Lists

 User feedback is an integral part of the evolution of this and other
 Bioperl modules. Send your comments and suggestions preferably to
 the Bioperl mailing list.  Your participation is much appreciated.

 bioperl-l@bioperl.org              - General discussion
 http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted va the
web:

 http://bugzilla.bioperl.org/

=head1 AUTHOR

 Based on the Bio::Tools::Signalp module
 Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are usually preceded with a _

=cut

package Bio::Tools::Signalp::ExtendedSignalp;
use vars qw(@ISA);
use strict;

use constant DEBUG => 0;
use Bio::Root::Root;
use Bio::Root::IO;
use Bio::SeqFeature::Generic;
use Bio::Tools::Signalp;
@ISA = qw(Bio::Root::Root Bio::Root::IO Bio::Tools::Signalp);

#Supported factors
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
      $factors && $self->factors($factors);
						
      return $self;
}

=head2 next_result

 Title   : next_result
 Usage   : my $feat = $signalp->next_result
 Function: Get the next result from parser data
 Returns : Bio::SeqFeature::Generic
 Args    : none


=cut

sub next_result {

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
	if($hash->{$fact} eq 'NO'){
	    $bool = 0;
	}
    }

    return $bool;

}

=head2 factors

 Title   : factors
 Usage   : my $factors = $signalp->factors
 Function: Get/Set the filters required from the user
 Returns : hash
 Args    : array reference

=cut

sub factors {

    my($self, $array) = @_;

    foreach my $f (@$array){
	if(exists($FACTS->{$f})){
	    $self->{_factors}->{$f} = 1;
	}
	else{
	    $self->throw("Factor [$f] incorrect. Supported:\n- ".join("\n- ", keys %$FACTS)."\n");
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

    my $idseen = 0;
    my $feature = undef;

    while (my $line = $self->_readline()) {

	chomp $line;
	if($line =~ /^>(\S+)/){
	    $idseen = 1;
	    next;
	}

	if($line =~ /^SignalP-NN result:/){
	    $self->_pushback($line);
	    $feature = $self->_parse_nn_result($feature);
	}
	elsif($line =~ /^SignalP-HMM result:/){
	    $self->_pushback($line);
	    $feature = $self->_parse_hmm_result($feature);
	}

	if($idseen && $feature){
	    my $new_feature = $self->create_feature($feature);
	    push @{$self->{_features}}, $new_feature if $new_feature;
	    $idseen = 0;
	    $feature = { };
	}
    }

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
    my $id;
    my %facts;
    my $end  = 0;

    #SignalP-NN result:
    #>MGG_11635.5           length = 100
    ## Measure  Position  Value  Cutoff  signal peptide?
    #  max. C    37       0.087   0.32   NO
    #  max. Y    37       0.042   0.33   NO
    #  max. S     3       0.062   0.87   NO
    #  mean S     1-36    0.024   0.48   NO
    #       D     1-36    0.033   0.43   NO

    while(my $line = $self->_readline()){

	if($line =~ /^SignalP-NN result:/){
	    $ok = 1;
	    next;
	}

	$self->throw("Wrong line for parsing NN results.") unless $ok;

	if ($line=~/^\>(\S+)\s+length\s+=\s+(\d+)/) {
	    $id = $1;
	    $self->seqname($id);
	    %facts = ();
	    $feature->{length} = $2;
	    next;
	}
	elsif($line =~ /max\.\s+C\s+\S+\s+(\S+)\s+\S+\s+(\S+)/) {
	    $facts{maxC} = $2;
	    $feature->{maxC} = $1;
	    next;
	}
	elsif ($line =~ /max\.\s+Y\s+\S+\s+(\S+)\s+\S+\s+(\S+)/) {
	    $facts{maxY} = $2;
	    $feature->{maxY} = $1;
	    next;
	}
	elsif($line =~ /max\.\s+S\s+\S+\s+(\S+)\s+\S+\s+(\S+)/) {
	    $facts{maxS} = $2;
	    $feature->{maxS} = $1;
	    next;
	}
	elsif ($line=~/mean\s+S\s+\S+\s+(\S+)\s+\S+\s+(\S+)/) {
	    $facts{meanS} = $2;
	    $feature->{meanS} = $1;
	    next;
	}
	elsif ($line=~/\s+D\s+\S+\s+(\S+)\s+\S+\s+(\S+)/) {
	    $facts{D} = $2;
	    $feature->{D} = $1;
	    #flag the end of the scores in case we don't have a signal peptide
	    #prediction following this current line...
	    $end = 1;
	    next;
	}
	#If we don't have this line it means that all the factors cutoff are equal to 'NO'
	elsif ($line =~ /Most likely cleavage site between pos\.\s+(\d+)/ || $end) {

	    if($self->_filterok(\%facts)){

		$feature->{name}       = $self->seqname();
		$feature->{start}      = 1;
		$feature->{end}        = $1;
		$feature->{source}     = 'Signalp';
		$feature->{primary}    = 'signal_peptide';
		$feature->{program}    = 'Signalp';
		$feature->{logic_name} = 'signal_peptide';

	    }else{
		print STDERR "FILTERS NOT OK FOR ",$self->seqname(), "\n" if DEBUG;
	    }
	    return $feature;
	}
    }

    return;
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

    while(my $line = $self->_readline){
        chomp $line;

	if($line =~ /^SignalP-HMM result:/){
	    $ok = 1;
	    next;
	}

	$self->throw("Wrong line for parsing HMM result.") unless $ok;

        if($line =~ /Prediction: (.+)$/){
            $feature_hash->{hmmPrediction} = $1;
        }elsif($line =~ /Signal peptide probability: ([0-9\.]+)/){
            $feature_hash->{peptideProb} = $1;
        }elsif($line =~ /Signal anchor probability: ([0-9\.]+)/){
            $feature_hash->{anchorProb} = $1;
            last;
        }
    }

    return $feature_hash;
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
    unless($feat->{name}){# && $feat->{start} && $feat->{end}){
	return;
    }

    # create feature object
    my $feature = Bio::SeqFeature::Generic->new(
						-seq_id      => $feat->{name},
						-start       => $feat->{start} || '',
						-end         => $feat->{end}   || '',
						-score       => $feat->{score},
						-source      => $feat->{source},
						-primary     => $feat->{primary},
						-logic_name  => $feat->{logic_name},
					       );

    $feature->score($feat->{peptideProb});
    $feature->add_tag_value('peptideProb', $feat->{peptideProb});
    $feature->add_tag_value('anchorProb', $feat->{anchorProb});
    $feature->add_tag_value('evalue',$feat->{anchorProb});
    $feature->add_tag_value('percent_id','NULL');
    $feature->add_tag_value("hid",$feat->{primary});
    $feature->add_tag_value('signalpPrediction', $feat->{hmmPrediction});
    $feature->add_tag_value('maxC_score', $feat->{maxC});
    $feature->add_tag_value('maxS_score', $feat->{maxS});
    $feature->add_tag_value('meanS_score', $feat->{meanS});
    $feature->add_tag_value('maxY_score', $feat->{maxY});
    $feature->add_tag_value('D_score', $feat->{D});

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

    my ($self, $seqname)=@_;

    if (defined($seqname)){
        $self->{'seqname'} = $seqname;
    }

    return $self->{'seqname'};

}

1;


