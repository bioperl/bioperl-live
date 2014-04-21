#
# BioPerl module for Bio::Tools::Glimmer
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Glimmer - parser for Glimmer 2.X/3.X prokaryotic and 
GlimmerM/GlimmerHMM eukaryotic gene predictions

=head1 SYNOPSIS

   use Bio::Tools::Glimmer;

   # file
   my $parser = Bio::Tools::Glimmer->new(-file => $file);
   # filehandle:
   $parser = Bio::Tools::Glimmer->new( -fh  => \*INPUT );
   # provide a sequence identifier (Glimmer 2.X)
   my $parser = Bio::Tools::Glimmer->new(-file => $file, -seqname => seqname);
   # force format (override automatic detection)
   my $parser = Bio::Tools::Glimmer->new(-file => $file, -format => 'GlimmerM');

   # parse the results
   # note: this class is-a Bio::Tools::AnalysisResult which implements
   # Bio::SeqAnalysisParserI, i.e., $glimmer->next_feature() is the same

   while(my $gene = $parser->next_prediction()) {
       # For eukaryotic input (GlimmerM/GlimmerHMM), $gene will be an instance
       # of Bio::Tools::Prediction::Gene, which inherits off
       # Bio::SeqFeature::Gene::Transcript, and $gene->exons() will return an
       # array of Bio::Tools::Prediction::Exon objects.
       # For prokaryotic input (Glimmer2.X/Glimmer3.X), $gene will be an
       # instance of Bio::SeqFeature::Generic

       # all exons (eukaryotic only):
       @exon_arr = $gene->exons();
       # initial exons only
       @init_exons = $gene->exons('Initial');
       # internal exons only
       @intrl_exons = $gene->exons('Internal');
       # terminal exons only
       @term_exons = $gene->exons('Terminal');
   }

=head1 DESCRIPTION

This is a module for parsing Glimmer, GlimmerM and GlimmerHMM predictions.  
It will create gene objects from the prediction report which can 
be attached to a sequence using Bioperl objects, or output as GFF 
suitable for loading into Bio::DB::GFF for use with Gbrowse.

Glimmer is open source and available at
L<http://www.cbcb.umd.edu/software/glimmer/>.

GlimmerM is open source and available at 
L<http://www.tigr.org/software/glimmerm/>.

GlimmerHMM is open source and available at
L<http://www.cbcb.umd.edu/software/GlimmerHMM/>.

Note that Glimmer 2.X will only process the first
sequence in a fasta file, and the prediction report does not contain any
sort of sequence identifier

Note that Glimmer 3.X produces two output files.  This module only parses
the .predict file.


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
email or the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 CONTRIBUTORS

Torsten Seemann

Mark Johnson

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Glimmer;
use strict;

use Bio::Factory::FTLocationFactory;
use Bio::Tools::Prediction::Gene;
use Bio::Tools::Prediction::Exon;

use base qw(Bio::Tools::AnalysisResult);

sub _initialize_state {
    my($self,@args) = @_;

    # first call the inherited method!
    my $make = $self->SUPER::_initialize_state(@args);

    $self->{'_preds_parsed'} = 0;
    # array of pre-parsed predictions
    $self->{'_preds'} = [];
}

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Glimmer->new();
 Function: Builds a new Bio::Tools::Glimmer object 
 Returns : an instance of Bio::Tools::Glimmer
 Args    : format ('Glimmer', 'GlimmerM', 'GlimmerHMM'), seqname


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($format, $seqname, $seqlength, $detail) =
      $self->_rearrange([qw(FORMAT SEQNAME SEQLENGTH DETAIL)], @args);

  # override automagic format detection
  if (defined($format) &&
      (($format eq 'Glimmer')  ||
       ($format eq 'GlimmerM') ||
       ($format eq 'GlimmerHMM'))
  ) {
      $self->_format($format);
  }
  
  if (defined($detail)) {
      $self->_format('Glimmer');
      $self->_detail_file($detail);
  }
  
  # hardwire seq_id when creating gene and exon objects (Glimmer 2.X)
  $self->_seqname($seqname) if defined($seqname);
  
  # store the length of the input sequence (Glimmer 2.X)
  $self->_seqlength($seqlength) if defined($seqlength);
  
  return $self;
}

=head2 analysis_method

 Usage     : $glimmer->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
             /glimmer/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method { 
#-------------
    my ($self, $method) = @_;  
    if($method && ($method !~ /glimmer/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}

=head2 next_feature

 Title   : next_feature
 Usage   : while($gene = $glimmer->next_feature()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the Glimmer result
           file. Call this method repeatedly until FALSE is returned.

           The returned object is actually a SeqFeatureI implementing object.
           This method is required for classes implementing the
           SeqAnalysisParserI interface, and is merely an alias for 
           next_prediction() at present.

 Example :
 Returns : A Bio::Tools::Prediction::Gene object.
 Args    :

=cut

sub next_feature {
    my ($self,@args) = @_;
    # even though next_prediction doesn't expect any args (and this method
    # does neither), we pass on args in order to be prepared if this changes
    # ever
    return $self->next_prediction(@args);
}

=head2 next_prediction

 Title   : next_prediction
 Usage   : while($gene = $glimmer->next_prediction()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the Glimmer result
           file. Call this method repeatedly until FALSE is returned.

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

    
    my %method = (
                  'Glimmer'    => '_parse_prokaryotic',
                  'GlimmerM'   => '_parse_eukaryotic',
                  'GlimmerHMM' => '_parse_eukaryotic',
                  '_DEFAULT_'  => '_parse_eukaryotic',
              );
    
    my $format = $self->_format();
    
    if (!$format) {
        
        while (my $line = $self->_readline()) {

            if ( $line =~ /^Glimmer\S*\s+\(Version\s*\S+\)/ ) {
                $format = 'GlimmerM';
                $self->_pushback($line);
                last;
            }
            elsif ( $line =~ /^Glimmer\S*$/ ) {
                $format = 'GlimmerHMM';
                $self->_pushback($line);
                last;
            }
            elsif ($line =~ /^Putative Genes:$/) {
                $format = 'Glimmer';
                $self->_pushback($line);
                last;
            }
            elsif ($line =~ /^>(\S+)/) {
                $format = 'Glimmer';
                $self->_pushback($line);
                last;
            }
            
        }
        
    }

    my $method =
        (exists($method{$format})) ? $method{$format} : $method{'_DEFAULT_'};

    return $self->$method();
    
}


=head2 _parse_eukaryotic

 Title   : _parse_eukaryotic()
 Usage   : $obj->_parse_eukaryotic()
 Function: Parses the prediction section. Automatically called by
           next_prediction() if not yet done.
 Example :
 Returns : 

=cut

sub _parse_eukaryotic {
    my ($self) = @_;

    my ($gene,$seqname,$seqlen,$source,$lastgenenum);
    
    while(defined($_ = $self->_readline())) {
	if( /^(Glimmer\S*)\s+\(Version\s*(\S+)\)/ ) {
	    $source = "$1_$2";
	    next;
	} elsif( /^(GlimmerHMM\S*)$/ ) { # GlimmerHMM has no version
	    $source = $1;
	    next;
	} elsif(/^Sequence name:\s+(.+)$/ ) {
	    $seqname = $1;
	    next;
	} elsif( /^Sequence length:\s+(\S+)/ ) {
	    $seqlen = $1;
	    next;
	} elsif( m/^(Predicted genes)|(Gene)|\s+\#/ || /^\s+$/ ) { 
	    next;
	    
	} elsif( # GlimmerM/HMM gene-exon prediction line
		 /^\s+(\d+)\s+ # gene num
		 (\d+)\s+      # exon num
		 ([\+\-])\s+   # strand
		 (\S+)\s+      # exon type
		 (\d+)\s+(\d+) # exon start, end
		 \s+(\d+)      # exon length		 
		 /ox ) {
	    my ($genenum,$exonnum,$strand,$type,$start,$end,$len) = 
		( $1,$2,$3,$4,$5,$6,$7);
	    if( ! $lastgenenum || $lastgenenum != $genenum) {		
		$self->_add_prediction($gene) if ( $gene );
		$gene = Bio::Tools::Prediction::Gene->new
		    (
		     '-seq_id'      => $seqname,
		     '-primary_tag' => "gene",
		     '-source_tag'  => $source,
		     '-tag'         => { 'Group' => "GenePrediction$genenum"},
		     );
	    }
	    my $exon = Bio::Tools::Prediction::Exon->new
		('-seq_id'     => $seqname,
		 '-start'      => $start,
		 '-end'        => $end,
		 '-strand'     => $strand eq '-' ? '-1' : '1',
		 '-source_tag' => $source,
		 '-primary_tag'=> 'exon',
		 '-tag'         => { 'Group' => "GenePrediction$genenum"},
		 );
	    $gene->add_exon($exon,lc($type));
	    $lastgenenum = $genenum;
	}
    }
    $self->_add_prediction($gene) if( $gene );
    $self->_predictions_parsed(1);
}

=head2 _parse_prokaryotic

 Title   : _parse_prokaryotic()
 Usage   : $obj->_parse_prokaryotic()
 Function: Parses the prediction section. Automatically called by
           next_prediction() if not yet done.
 Example :
 Returns : 

=cut

sub _parse_prokaryotic {
    my ($self) = @_;

    # default value, possibly overriden later
    my $source = 'Glimmer';

    # Store the sequence length(s) here, either from the
    # seqlength arg to the constructor, or from the
    # Glimmer 3.X detail file
    my %seqlength = ( );
    
    # Glimmer 2.X does not provide a sequence identifer
    # in the prediction report (will default to unknown
    # if not specified in the seqname arg to the
    # constructor
    #
    # Glimmer 2.X does not report the length of the
    # input sequence, either (will default to undef
    # if not specified in the seqlength arg to the
    # constructor
    my $seqname = $self->_seqname();
    my $seqlength = $self->_seqlength();

    if (defined($seqlength)) {
        $seqlength{$seqname} = $seqlength
    }

    # Parse the detail file, if we have one (Glimmer 3.X)
    my $detail_file = $self->_detail_file();
    
    if (defined($detail_file)) {

        my $io = Bio::Root::IO->new(-file => $detail_file);
        my $seqname;
        
        while (defined($_ = $io->_readline())) {
            if ($_ =~ /^>(\S+)/) {
                $seqname = $1;
                next;
            }

            if (defined($seqname) && ($_ =~ /^Sequence length = (\d+)$/)) {
                $seqlength{$seqname} = $1;
                next;
            }
        }
    }
    
    my $location_factory = Bio::Factory::FTLocationFactory->new();
    
    while(defined($_ = $self->_readline())) {
        # Glimmer 3.X does provide a sequence identifier -
        # beware whitespace at the end (comes through from
        # the fasta file)
        if ($_ =~ /^Putative Genes:$/) {
            $source = 'Glimmer_2.X';
            next;
        }
        # Glimmer 3.X sequence identifier
        elsif ($_ =~ /^>(\S+)/) {
            $seqname = $1;
            $seqlength = $seqlength{$seqname};
            $source = 'Glimmer_3.X';
            next;
        }        
        elsif (
               # Glimmer 2.X prediction 
               (/^\s+(\d+)\s+      # gene num
                (\d+)\s+(\d+)\s+   # start, end
                \[([\+\-])(\d{1})\s+ # strand, frame
                /ox ) ||
               # Glimmer 3.X prediction
               (/^[^\d]+(\d+)\s+    # orf (numeric portion)
                (\d+)\s+(\d+)\s+   # start, end
                ([\+\-])(\d{1})\s+   # strand, frame
                ([\d\.]+)          # score
                /ox)) {
	    my ($genenum,$start,$end,$strand,$frame,$score) = 
		( $1,$2,$3,$4,$5,$6 );

            my $circular_prediction = 0;

            # Check for a circular prediction before we
            # start fiddling with the coordinates
            if ($strand eq '+') {
                if ($start > $end) {
                    $circular_prediction = 1;
                }
            }
            else {
                if ($start < $end) {
                    $circular_prediction = 1;
                }
            }

            if ($circular_prediction) {
                unless (defined($seqlength)) {
                    $self->throw("need to know the sequence length to handle wraparound genes");
                }
            }
            
            # Glimmer 2.X predictions do not include
            # the stop codon - this might extend the
            # prediction off either end of the sequence.
            # This works fine even on circular/wraparound
            # predictions.  
            if ($source eq 'Glimmer_2.X') {
                if ($strand eq '+') { 
                    $end += 3;
                }
                else {
                    $end -= 3;
                }
            }

            # We might have extended a Glimmer 2.X prediction
            # beyond the boundaries of the input sequence.
            # Also, Glimmer 3.X (with -X) will output predictions
            # with coordinates less than 1 or greater than the
            # length of the sequence.
            my ($fst, $fend);
            foreach my $coord ($start, $end) {
                if ($coord < 1) {
                    $coord = '<1';
                    $fst++;
                } elsif (defined($seqlength) && ($coord > $seqlength)) {
                    $coord = ">$seqlength";
                    $fend++;
                }
            }
            
            my $location_string;

            if ($circular_prediction) {
                if ($strand eq '+') {
                    $location_string = "join($start..$seqlength,1..$end)";
                }
                else {
                    $location_string = "join($start..1,$seqlength..$end)";
                }
                }
            else {
                # start must always be less than end for gene locations
                if ($strand eq '-' && !$fst && !$fend && $start > $end) {
                    ($start, $end) = ($end, $start);
                }
                $location_string = "$start..$end";
            }
            
            my $location_object =
                $location_factory->from_string($location_string);
            
            # convert glimmer's frame range from 1-3 to SeqFeature's 0-2.
            $frame--;
            
            my $gene = Bio::SeqFeature::Generic->new
                (
                 '-seq_id'      => $seqname,
                 '-location'   => $location_object,
                 '-strand'     => $strand eq '-' ? '-1' : '1',
                 '-frame'		=> $frame,
                 '-source_tag'  => $source,
                 '-display_name' => "orf$genenum",
                 '-primary_tag'=> 'gene',
                 '-tag'         => { 'Group' => "GenePrediction_$genenum"},
                 '-score'       => $score || undef
             );
            
            $self->_add_prediction($gene) 
	}
    }
    
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

    return unless(exists($self->{'_preds'}) && @{$self->{'_preds'}});
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

=head2 _seqname

 Title   : _seqname
 Usage   : $obj->_seqname($seqname)
 Function: internal (for Glimmer 2.X)
 Example :
 Returns : String

=cut

sub _seqname {
    my ($self, $val) = @_;

    $self->{'_seqname'} = $val if $val;
    if(! exists($self->{'_seqname'})) {
	$self->{'_seqname'} = 'unknown';
    }
    return $self->{'_seqname'};
}

=head2 _seqlength

 Title   : _seqlength
 Usage   : $obj->_seqlength($seqlength)
 Function: internal (for Glimmer 2.X)
 Example :
 Returns : String

=cut

sub _seqlength {
    my ($self, $val) = @_;

    $self->{'_seqlength'} = $val if $val;
    return $self->{'_seqlength'};
}

=head2 _format

 Title   : _format
 Usage   : $obj->_format($format)
 Function: internal
 Example :
 Returns : String

=cut

sub _format {
    my ($self, $val) = @_;

    $self->{'_format'} = $val if $val;

    return $self->{'_format'};
}

=head2 _detail_file

 Title   : _detail_file
 Usage   : $obj->_detail_file($filename)
 Function: internal (for Glimmer 3.X)
 Example :
 Returns : String

=cut

sub _detail_file {
    my ($self, $val) = @_;

    $self->{'_detail_file'} = $val if $val;
    return $self->{'_detail_file'};
}

1;
