
#
# BioPerl module for Bio::Tools::MZEF
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::MZEF - Results of one MZEF run

=head1 SYNOPSIS

   $mzef = Bio::Tools::MZEF->new(-file => 'result.mzef');
   # filehandle:
   $mzef = Bio::Tools::MZEF->new( -fh  => \*INPUT );
   # to indicate that the sequence was reversed prior to feeding it to MZEF
   # and that you want to have this reflected in the strand() attribute of 
   # the exons, as well have the coordinates translated to the non-reversed
   # sequence
   $mzef = Bio::Tools::MZEF->new( -file   => 'result.mzef',
                                  -strand => -1 );

   # parse the results
   while($gene = $mzef->next_prediction()) {
       # $gene is an instance of Bio::Tools::Prediction::Gene
       
       # $gene->exons() returns an array of 
       # Bio::Tools::Prediction::Exon objects
       # all exons:
       @exon_arr = $gene->exons();

       # internal exons only
       @intrl_exons = $gene->exons('Internal');
       # note that presently MZEF predicts only internal exons!
   }

   # essential if you gave a filename at initialization (otherwise the file
   # will stay open)
   $mzef->close();

=head1 DESCRIPTION

The MZEF module provides a parser for MZEF gene structure prediction
output.


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
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Hilmar Lapp

Email hlapp@gmx.net (or hilmar.lapp@pharma.novartis.com)

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::MZEF;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;
use Bio::Tools::Prediction::Gene;
use Bio::Tools::Prediction::Exon;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

  my($fh,$file,$strand) =
      $self->_rearrange([qw(FH
			    FILE
			    STRAND
			    )],
			@args);

  $self->{'readbuffer'} = "";
  $strand = 1 unless defined($strand);
  $self->{'_strand'} = $strand;

  if( defined $fh && defined $file ) {
      $self->throw("You have defined both a filehandle and file to read from. Not good news!");
  }
  if((defined $file) && ($file ne '')) {
      $fh = Symbol::gensym();
      open ($fh,$file)
	  || $self->throw("Could not open $file for Fasta stream reading $!");
  }
  if((! defined($fh)) && ($file eq "")) {
      $fh = \*STDIN;
  }
  $self->_filehandle($fh) if defined $fh;

  # private state variables
  $self->{'_preds_parsed'} = 0;
  $self->{'_has_cds'} = 0;

  # set stuff in self from @args
  return $make; # success - we hope!
}


=head2 close

 Title   : close
 Usage   : $sim4_result->close()
 Function: Closes the file handle associated with this result file
 Example :
 Returns :
 Args    :

=cut

sub close {
   my ($self, @args) = @_;

   $self->{'_filehandle'} = undef;
}

=head2 _pushback

 Title   : _pushback
 Usage   : $obj->_pushback($newvalue)
 Function: puts a line previously read with _readline back into a buffer
 Example :
 Returns :
 Args    : newvalue

=cut

sub _pushback {
  my ($obj, $value) = @_;
  $obj->{'readbuffer'} .= $value;
}


=head2 _filehandle

 Title   : _filehandle
 Usage   : $obj->_filehandle($newval)
 Function:
 Example :
 Returns : value of _filehandle
 Args    : newvalue (optional)


=cut

sub _filehandle {
    my ($obj, $value) = @_;
    if(defined $value) {
	$obj->{'_filehandle'} = $value;
    }
    return $obj->{'_filehandle'};
}


=head2 _readline

 Title   : _readline
 Usage   : $obj->_readline
 Function:
 Example :
 Returns : reads a line of input

=cut

sub _readline {
  my $self = shift;
  my $fh = $self->_filehandle();
  my $line;

  # if the buffer been filled by _pushback then return the buffer
  # contents, rather than read from the filehandle
  if ( defined $self->{'readbuffer'} ) {
      $line = $self->{'readbuffer'};
      undef $self->{'readbuffer'};
  } else {
      $line = defined($fh) ? <$fh> : <>;
  }
  return $line;
}

sub DESTROY {
    my $self = shift;

    $self->close();
}

=head2 next_prediction

 Title   : next_prediction
 Usage   : while($gene = $mzef->next_prediction()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the MZEF result
           file. Call this method repeatedly until FALSE is returned.

           Note that with the present version of MZEF there will only be one
           object returned, because MZEF does not predict individual genes
           but just potential internal exons.
 Example :
 Returns : A Bio::Tools::Prediction::Gene object.
 Args    :

=cut

sub next_prediction {
    my ($self) = @_;
    my $gene;

    # if the prediction section hasn't been parsed yet, we do this now
    $self->_parse_predictions() unless $self->_predictions_parsed();

    # get next gene structure
    $gene = $self->_prediction();

    if($gene) {
	# there is no predicted peptide or CDS sequence which could be filled
	# in here
    }

    return $gene;
}

=head2 _parse_predictions

 Title   : _parse_predictions()
 Usage   : $obj->_parse_predictions()
 Function: Parses the prediction section. Automatically called by
           next_prediction() if not yet done.
 Example :
 Returns : 

=cut

sub _parse_predictions {
    my ($self) = @_;
    my ($method); # set but not used presently
    my $exon_tag = "InternalExon";
    my $gene;
    # my $seqname; # name given in output is poorly formatted
    my $seqlen;
    my $prednr = 1;

    while(defined($_ = $self->_readline())) {
	if(/^\s*(\d+)\s*-\s*(\d+)\s+/) {
	    # exon or signal
	    if(! defined($gene)) {
		$gene = Bio::Tools::Prediction::Gene->new(
                                       '-primary' => "GenePrediction$prednr",
				       '-source' => 'MZEF');
	    }
	    # we handle start-end first because may not be space delimited
	    # for large numbers
	    my ($start,$end) = ($1,$2);
	    s/^\s*(\d+)\s*-\s*(\d+)\s+//;
	    # split the rest into fields
	    chomp();
	    # format: Coordinates P Fr1 Fr2 Fr3 Orf 3ss Cds 5ss
	    # index:              0   1   2   3   4   5   6   7
	    my @flds = split(' ', $_);
	    # create the feature object depending on the type of signal --
	    # which is always an (internal) exon for MZEF
	    my $predobj = Bio::Tools::Prediction::Exon->new();
	    # set common fields
	    $predobj->source_tag('MZEF');
	    $predobj->significance($flds[0]);
	    $predobj->score($flds[0]); # what shall we set as overall score?
	    $predobj->strand($self->{'_strand'}); # MZEF searches only one
	    if($predobj->strand() == 1) {
		$predobj->start($start);
		$predobj->end($end);
	    } else {
		$predobj->start($seqlen-$end+1);
		$predobj->end($seqlen-$start+1);
	    }
	    # set scores
	    $predobj->start_signal_score($flds[5]);
	    $predobj->end_signal_score($flds[7]);
	    $predobj->coding_signal_score($flds[6]);
	    # frame -- we simply extract the one with highest score from the
	    # orf field, and store the individual scores for now
	    my $frm = index($flds[4], "1");
	    $predobj->frame(($frm < 0) ? undef : $frm);
	    $predobj->primary_tag($exon_tag);
	    # add to gene structure (should be done only when start and end
	    # are set, in order to allow for proper expansion of the range)
	    $gene->add_exon($predobj);		
	    next;
	}
	if(/^\s*Internal .*(MZEF)/) {
	    $method = $1;
	    next;
	}
	if(/^\s*File_Name:\s+(\S+)\s+Sequence_length:\s+(\d+)/) {
	    # $seqname = $1; # this is too poor currently (file name truncated
                             # to 10 chars) in order to be sensible enough
	    $seqlen = $2;
	    next;
	}
    }
    # $gene->seqname($seqname);
    $self->_add_prediction($gene) if defined($gene);
    $self->_predictions_parsed(1);
}

=head2 _prediction

 Title   : _prediction()
 Usage   : $gene = $obj->_prediction()
 Function: internal
 Example :
 Returns : 

=cut

sub _prediction {
    my ($self) = @_;

    return undef unless(exists($self->{'_preds'}) && @{$self->{'_preds'}});
    return shift(@{$self->{'_preds'}});
}

=head2 _add_prediction

 Title   : _add_prediction()
 Usage   : $obj->_add_prediction($gene)
 Function: internal
 Example :
 Returns : 

=cut

sub _add_prediction {
    my ($self, $gene) = @_;

    if(! exists($self->{'_preds'})) {
	$self->{'_preds'} = [];
    }
    push(@{$self->{'_preds'}}, $gene);
}

=head2 _predictions_parsed

 Title   : _predictions_parsed
 Usage   : $obj->_predictions_parsed
 Function: internal
 Example :
 Returns : TRUE or FALSE

=cut

sub _predictions_parsed {
    my ($self, $val) = @_;

    $self->{'_preds_parsed'} = $val if $val;
    if(! exists($self->{'_preds_parsed'})) {
	$self->{'_preds_parsed'} = 0;
    }
    return $self->{'_preds_parsed'};
}


1;
