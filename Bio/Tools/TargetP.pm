#
# Bioperl module for TargetP
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::TargetP - Results of one TargetP run

=head1 SYNOPSIS

   use Bio::Tools::TargetP;

   #filename for  TargetP result :
   $targetp = Bio::Tools::TargetP->new(-file => 'targetp.out');

   # filehandle for TargetP :
   $targetp = Bio::Tools::TargetP->new( -fh  => \*INPUT );

   ### targetp v1.1 prediction results ##################################
   #Number of query sequences:  11
   #Cleavage site predictions included.
   #Using NON-PLANT networks.
   #
   #Name                  Len            mTP     SP  other  Loc  RC  TPlen
   #----------------------------------------------------------------------
   #swall|Q9LIP3|C72Y_AR  500          0.245  0.935  0.009   S    2     22
   #swall|Q52813|AAPQ_RH  400          0.170  0.462  0.577   _    5      -
   #swall|O86459|AAT_RHI  400          0.346  0.046  0.660   _    4      -



   # parse the results
   while($feature = $targetp->next_prediction()) {

           #$feature is a Bio::SeqFeature::Generic object
           my $method     = $targetp->analysis_method();
           my $vesion     = $targetp->analysis_method_version() || $feature->source();
           my $seqid      = $feature->seq_id();
           # ...
     }

   # essential if you gave a filename at initialization (otherwise the file
   # will stay open)
   $targetp->close();

=head1 DESCRIPTION

TargetP modules will provides parsed information about protein
localization.  It reads in a targetp output file.  It parses the
results, and returns a Bio::SeqFeature::Generic object for each
seqeunces found to have a subcellular localization

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

=head1 AUTHORS - Emmanuel Quevillon

Email emmanuel.quevillon@versailles.inra.fr

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::TargetP;
use strict;
use Bio::Tools::AnalysisResult;
use Bio::SeqFeature::Generic;
use Data::Dumper;

use base qw(Bio::Tools::AnalysisResult);


#Definition of 'Loc' field according to http://www.cbs.dtu.dk/services/TargetP/output.php
my $MAPLOC = {
	      'S' => 'Secretory pathway',
	      'M' => 'Mitochondrion',
	      'C' => 'Chloroplast',
	      '_' => 'Any other',
	      '*' => 'Unknown',
	      '?' => 'Unknown',
	     };


=head1 analysis_method

 Usage     : $self->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
 Returns   : String
 Argument  : n/a

=cut

sub analysis_method {

    my ($self, $method) = @_;

    if($method && ($method !~ /TargetP/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }

    return $self->SUPER::analysis_method($method);
}

=head1 network

  Title   : network
  Usage   : $self->network($network)
  Function: This method Get/Set the network used for the analysis (PLANT or NON-PLANT)
  Example :
  Returns : string
  Arguments: On set, the network used

=cut

sub network {

    my($self, $net) = @_;

    if(defined($net)){
	$self->{'_network'} = $net;
    }

    return $self->{'_network'};

}


=head1 cleavage

  Title    :  cleavage
  Usage    : $self->cleavage($cleavage)
  Function : This method Get/Set if SignalP program was used to run TargetP
  Example  :
  Returns  : 1 or 0
  Arguments: On set, the cleavage used or not

=cut

sub cleavage {

    my($self, $cleavage) = @_;

    if(defined($cleavage)){
	$self->{'_cleavage'} = $cleavage =~ /not included/ ? '0' : '1';
    }

    return $self->{'_cleavage'};

}


=head1 next_prediction

  Usage    : $targetp->next_prediction()
  Purpose  : Returns the next TargetP prediction
  Returns  : A Bio::SeqFeature::Generic object
  Arguments: n/a

=cut

sub next_prediction {

    my($self) = @_;

    unless($self->_parsed()){
	$self->_parse_results();
	$self->_parsed(1);
    }

    return shift @{$self->{'_features'}} || undef;
}

=head1 create_feature

  Title     : create_feature
  Usage     : $self->create_feature(\%hash);
  Function  : This method creates a new Bio::SeqFeature::Generic object
  Example   : 
  Returns   : Bio::SeqFeature::Generic
  Arguments : hash reference

=cut

sub create_feature {

    my($self, $feat) = @_;

    $self->throw("Need a reference to hash table") unless($feat && ref($feat) eq 'HASH');

    my $feature = Bio::SeqFeature::Generic->new(
						-seq_id      => $feat->{seqid},
						-source_tag  => $self->analysis_method(),
						-primary_tag => 'signal_peptide',        #Sequence Ontology compliant
						-strand      => '+',
					       );

    if(defined($feat->{seqlen})){
	$feature->start(1);
	$feature->end($feat->{seqlen});
    }
    $feature->add_tag_value('location',            $MAPLOC->{$feat->{loc}})   if(exists($MAPLOC->{$feat->{loc}}));
    $feature->add_tag_value('chloroplastCutOff',   $feat->{cTP})              if(defined($feat->{cTP}));
    $feature->add_tag_value('mitochondrionCutOff', $feat->{mTP})              if(defined($feat->{mTP}));
    $feature->add_tag_value('signalPeptideCutOff', $feat->{SP})               if(defined($feat->{SP}));
    $feature->add_tag_value('otherCutOff',         $feat->{other})            if(defined($feat->{other}));
    $feature->add_tag_value('reliabilityClass',    $feat->{RC})               if(defined($feat->{RC}));
    $feature->add_tag_value('signalPeptideLength', $feat->{TPLen})            if(defined($feat->{TPLen}));

    $feature->add_tag_value('network',             $self->network());

    return $feature;

}


=head2 PRIVATE METHODS

=cut

=head2 _initialize_state

 Title   : _initialize_state
 Usage   : n/a; usually called by _initialize() itself called by new()
 Function: This method is supposed to reset the state such that any 'history'
           is lost. State information that does not change during object
           lifetime is not considered as history, e.g. parent, name, etc shall
           not be reset. An inheriting object should only be concerned with
           state information it introduces itself, and for everything else
           call SUPER::_initialize_state(@args).

           The argument syntax is the same as for new() and _initialize(),
           i.e., named parameters following the -name=>$value convention.
           The following parameters are dealt with by the implementation
           provided here:
              -INPUT, -FH, -FILE
           (tags are case-insensitive).
 Example :
 Returns :
 Args    :

=cut

sub _initialize_state  {

  	my ($self,@args,) = @_;
  	# first call the inherited method!
  	$self->SUPER::_initialize_state(@args);

  	# our private state variables
  	$self->{'_features'}   = [ ];
  	$self->{'_parameters'} = undef;
	$self->{'_format'}     = undef;
	$self->{'_network'}    = undef;
	$self->{'_cleavage'}   = undef;
	$self->{'_parsed'}     = 0;

  	$self->analysis_method('TargetP');

	return 1;
}

=head2 _predictions

  Usage    : $targetp->_prediction()
  Purpose  : Returns the number of TargetP predictions
  Returns  : A scalar (number)
  Arguments: n/a

=cut

sub _predictions {

    my($self) = @_;

    return scalar(@{$self->{'_features'}}) || 0;
}


=head2 _parsed

 Title     : _parsed
 Usage     : $targetp->_parsed(1)
 Function  : This method is used to know if the output result is parsed or not
             For internal use only
 Example   :
 Returns   : 1/0
 Arguments : 1/0 for setting

=cut

sub _parsed {

    my($self, $value) = @_;

    if(defined($value)){
	$self->{'_parsed'} = $value;
    }

    return $self->{'_parsed'};
}



=head2 _parse_results

  Title    : _parse_results
  Usage    : $self->_parse_results()
  Function : This method parses a TargetP output
             For internal use only
  Example  :
  Returns  : n/a
  Arguments: none

=cut

sub _parse_results {

    my($self) = @_;


    ### targetp v1.1 prediction results ##################################
    #Number of query sequences:  11
    #Cleavage site predictions included.
    #Using NON-PLANT networks.
    #
    #Name                  Len            mTP     SP  other  Loc  RC  TPlen
    #----------------------------------------------------------------------
    #swall|Q9LIP3|C72Y_AR  500          0.245  0.935  0.009   S    2     22
    #swall|Q52813|AAPQ_RH  400          0.170  0.462  0.577   _    5      -
    #swall|O86459|AAT_RHI  400          0.346  0.046  0.660   _    4      -


    while(defined(my $line = $self->_readline())){

	if($line =~ /targetp (v[\d\.]+)/){

	    $self->analysis_method_version($1);

	}elsif($line =~ /Cleavage site predictions (.*)/){

	    $self->cleavage($1);

	}elsif($line =~ /Using (\S+) networks/){

	    $self->network($1);

	}elsif($line =~ /^Name/){

	    #We skip the next line which is '------------------'
	    $self->_readline();

	    my $hash = { };

	    while(defined(my $line = $self->_readline())){

		last if($line =~ /^----/);

		my $hash = $self->_parse_line($line);

		my $new_feature = $self->create_feature($hash);
		
		$self->_add_feature($new_feature);
	    }
	}
    }

    return;
}

=head2 _parse_line

 Title    : _parse_line
 Usage    : $self->_parse_line($line)
 Function : This method parses the line result
            For internal use only
 Example  :
 Returns  : Hash reference
 Arguemnts: line to parse

=cut

sub _parse_line {

    my($self, $line) = @_;

    $self->throw("No line to parse given") unless($line);

    my $hash = { };
    my ($seqid, $seqlen, $cTP, $mTP, $SP, $other, $loc, $RC, $TPlen);

    if($self->network() eq 'NON-PLANT'){

	($seqid, $seqlen, $mTP, $SP, $other, $loc, $RC, $TPlen) = split(/\s+/, $line);

    }else{

	($seqid, $seqlen, $cTP, $mTP, $SP, $other, $loc, $RC, $TPlen) = split(/\s+/, $line);

    }

    $hash->{seqid}  = $seqid;
    $hash->{seqlen} = $seqlen;
    $hash->{cTP}    = $cTP || undef;
    $hash->{mTP}    = $mTP;
    $hash->{SP}     = $SP;
    $hash->{other}  = $other;
    $hash->{loc}    = $loc;
    $hash->{RC}     = $RC;
    $hash->{TPLen}  = ($TPlen && $TPlen =~ /\d+/) ? $TPlen : undef;

    return $hash;

}

=head2 _add_feature

 Title    : _add_feature
 Usage    : $self->_add_feature($feature)
 Function : This method stores a feature object
            For internal use only
 Example  :
 Returns  : n/a
 Arguments: Bio::SeqFeature::Generic

=cut

sub _add_feature {

    my($self, $feature) = @_;

    $self->throw("Need a Bio::SeqFeature::Generic object") unless $feature->isa("Bio::SeqFeature::Generic");

    push(@{$self->{'_features'}}, $feature);

    return;

}

=head2 _toString_location

 Title    : _toString_location
 Usage    : $self->_toString_location($key)
 Function : This method convert the 'one letter code' location to 
            the corresponding definition
            For internal use only
 Example  :
 Returns  : Location or undef
 Arguments: String

=cut

sub _toString_location {

    my($self, $key) = @_;

    if($key && exists($MAPLOC->{$key})){
	return $MAPLOC->{$key};
    }

    return;
}



1;
