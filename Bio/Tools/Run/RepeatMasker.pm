# $Id$
# BioPerl module for RepeatMasker
#
# Cared for by Shawn Hoon <shawnh@fugu-sg.org>
#
# Copyright Shawn Hoon
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::RepeatMasker -
Wrapper for RepeatMasker Program

=head1 SYNOPSIS

  use Bio::Tools::Run::RepeatMasker;

  my @params=("mam" => 1,"noint"=>1);
  my $factory = Bio::Tools::Run::RepeatMasker->new(@params);
  $in  = Bio::SeqIO->new(-file => "contig1.fa",
                         -format => 'fasta');
  my $seq = $in->next_seq();

  #return an array of Bio::SeqFeature::FeaturePair objects
  my @feats = $factory->mask($seq); 

or

  $factory->mask($seq);
  my @feats = $factory->repeat_features;

  #return the masked sequence, a Bio::SeqI object
  my $masked_seq = $factory->masked_seq;

=head1 DESCRIPTION

RepeatMasker is a program that screens DNA sequences for interspersed
repeats known to exist in mammalian genomes as well as for low
complexity DNA sequences. For more information, on the program and its
usage, please refer to http://repeatmasker.genome.washington.edu/.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bioperl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Shawn Hoon


Email shawnh@fugu-sg.org


=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut


package Bio::Tools::Run::RepeatMasker;

use vars qw($AUTOLOAD @ISA $PROGRAM $PROGRAMDIR
            $TMPDIR $TMPOUTFILE @RM_SWITCHES @RM_PARAMS
            @OTHER_SWITCHES %OK_FIELD);

use strict;
use Bio::PrimarySeqI;
use Bio::SeqFeature::FeaturePair;
use Bio::Root::Root;
use Bio::Root::IO;


# Let the code begin...

@ISA = qw(Bio::Root::Root Bio::Root::IO );

BEGIN {

    if (defined $ENV{REPEATMASKERDIR}) {
        $PROGRAMDIR = $ENV{REPEATMASKERDIR} || '';
        $PROGRAM = Bio::Root::IO->
	    catfile($PROGRAMDIR."/src/bin/",
		    'RepeatMasker'.($^O =~ /mswin/i ?'.exe':''));
    }
    else {
        $PROGRAM = 'RepeatMasker';
    }
    @RM_PARAMS = qw(DIV LIB CUTOFF PARALLEL GC FRAG );

    @RM_SWITCHES = qw(NOLOW LOW L NOINT INT NORNA ALU M MUS ROD RODENT MAM MAMMAL COW AR 
                       ARABIDOPSIS DR DROSOPHILA EL ELEGANS IS_ONLY IS_CLIP NO_IS RODSPEC
		       PRIMSPEC W WUBLAST S Q QQ GCCALC NOCUT NOISY QUIET); 

    # Authorize attribute fields
    foreach my $attr ( @RM_PARAMS, @RM_SWITCHES,
                       @OTHER_SWITCHES) { $OK_FIELD{$attr}++; }
}

sub AUTOLOAD {
    my $self = shift;
    my $attr = $AUTOLOAD;
    $attr =~ s/.*:://;
    $attr = uc $attr;
    $self->throw("Unallowed parameter: $attr !") unless $OK_FIELD{$attr};
    $self->{$attr} = shift if @_;
    return $self->{$attr};
}

=head2 new

 Title   : new
 Usage   : $rm->new($seq)
 Function: creates a new wrapper
 Returns:  Bio::Tools::Run::RepeatMasker
 Args    : self

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  # to facilitiate tempfile cleanup
  $self->_initialize_io();

  my ($attr, $value);
  (undef,$TMPDIR) = $self->tempdir(CLEANUP=>1);
  (undef,$TMPOUTFILE) = $self->tempfile(-dir => $TMPDIR);
  while (@args) { 
    $attr =   shift @args;
    $value =  shift @args;
    next if( $attr =~ /^-/ ); # don't want named parameters
    if ($attr =~/'PROGRAM'/i) {
      $self->program($value);
      next;
    }
    $self->$attr($value);
  }
  if (! defined $self->program) {
    $self->program($PROGRAM);
  }
  unless ($self->exists_rm()) {
    if( $self->verbose >= 0 ) {
      warn "RepeatMasker program not found as ".$self->program.
	  " or not executable. \n"; 
    }
  }

  return $self;
}

=head2  exists_rm()

 Title   : exists_rm
 Usage   : $rmfound= Bio::Tools::Run::RepeatMasker->exists_rm()
 Function: Determine whether RepeatMasker program can be found on 
           current host
 Example :
 Returns : 1 if repeatmasker program found at expected location, 
           0 otherwise.
 Args    : none

=cut

sub exists_rm{
    my $self = shift;
    if( my $f = Bio::Root::IO->exists_exe($PROGRAM) ) {
  $PROGRAM = $f if( -e $f );
  return 1;
    }
}

=head2 program

 Title   : program
 Usage   : $obj->program($newval)
 Function:
 Returns : value of program
 Args    : newvalue (optional)

=cut

sub program{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'program'} = $value;
    }
    return $self->{'program'};

}

=head2  version

 Title   : version
 Usage   : 
 Function: Determine the version number of the program
 Example :
 Returns : float or undef
 Args    : none

=cut

sub version {
    my ($self) = @_;

    return undef unless $self->exists_rm;
    my $string = `RepeatMasker -- ` ;
    $string =~ /\(([\d.]+)\)/;
    return $1 || undef;

}

=head2  mask

 Title   : mask
 Usage   : $rm->mask($seq)
 Function: carry out repeat mask
 Example :
 Returns : an array of repeat features that are
           Bio::SeqFeature::FeaturePairs
 Args    : Bio::PrimarySeqI compliant object

=cut

sub mask {
  my ($self,$seq) = @_;
  my ($infile);
  $infile = $self->_setinput($seq);

  my $param_string = $self->_setparams();
  my @repeat_feats = $self->_run($infile,$param_string);

  return @repeat_feats;

}

=head2  _run

 Title   : _run
 Usage   : $rm->_run ($filename,$param_string)
 Function: internal function that runs the repeat masker
 Example :
 Returns : an array of repeat features
 Args    : the filename to the input sequence and the parameter string

=cut

sub _run {
  my ($self,$infile,$param_string) = @_;
  my $instring;
  $self->debug( "Program ".$self->program."\n");

  my $outfile = $infile.".out";
  my $cmd_str = $self->program." $param_string ". $infile;
  $self->debug("repeat masker command = $cmd_str");
  my $status = system($cmd_str);

  $self->throw("Repeat Masker Call($cmd_str) crashed: $?\n") 
      unless $status == 0;
  my $rpt_feat = $self->_parse_results($outfile);
  $self->repeat_features($rpt_feat);

  #get masked sequence
  my $masked = $infile.".masked";
  my $seqio = Bio::SeqIO->new(-file=>$masked,-format=>'FASTA');
  $self->masked_seq($seqio->next_seq);

  return @{$rpt_feat};
}

=head2  masked_seq

 Title   : masked_seq
 Usage   : $rm->masked_seq($seq)
 Function: get/set for masked sequence
 Example :
 Returns : the masked sequence
 Args    : Bio::Seq object

=cut

sub masked_seq {
  my ($self,$seq) = @_;
  if($seq){
    $self->{'_masked_seq'} = $seq;
  }
  return $self->{'_masked_seq'};
}

=head2  repeat_features

 Title   : repeat_features
 Usage   : $rm->repeat_features(\@rf)
 Function: get/set for repeat features array
 Example :
 Returns : the array of repeat features
 Args    : 

=cut

sub repeat_features {
  my ($self,$rf) = @_;
  if($rf) {
    $self->{'_rf'} = $rf;
  }
  return @{$self->{'_rf'}};
}

=head2  _parse_results()

 Title   : _parse_results
 Usage   : Internal function, not to be called directly
 Function: parses the results from RepeatMasker output 
           (largely copied from Ensembl RepeatMasker Runnable)
 Example :
 Returns : array of repeat features 
 Args    : the name of the output file

=cut

sub _parse_results {
  my ($self,$outfile) = @_;
  $outfile || $self->throw("No outfile specified");
  open(REPOUT,"<$outfile") || $self->throw("Error opening $outfile\n");
  my $filehandle = \*REPOUT;
  my @repeat_features;
  #check if no repeats found
    if (<$filehandle> =~ /no repetitive sequences detected/)
    {
        print STDERR "RepeatMasker didn't find any repetitive sequences\n";
        close $filehandle;
        return;
    }
    #extract values
     while (<$filehandle>) {  
       if (/\d+/) { #ignore introductory lines
            my @element = split;
            # ignore features with negatives 
            next if ($element[11-13] =~ /-/); 
            my (%feat1, %feat2);
            my ($score, $query_name, $query_start, $query_end, $strand,
		$repeat_name, $repeat_class ) = (split)[0, 4, 5, 6, 8, 9, 10];
            my ($hit_start,$hit_end);
	    if ($strand eq '+') {
            ($hit_start, $hit_end) = (split)[11, 12];
            $strand = 1;
	    }
       elsif ($strand eq 'C') {
	        ($hit_start, $hit_end) = (split)[12, 13];
          $strand = -1;
	    }
     #my $rc = $self->_get_consensus($repeat_name, $repeat_class);

	    my $rf = Bio::SeqFeature::FeaturePair->new;
	    $rf->score            ($score);
	    $rf->start            ($query_start);
	    $rf->end              ($query_end);
	    $rf->strand           ($strand);
	    $rf->hstart           ($hit_start);
	    $rf->hend             ($hit_end);
	    $rf->add_tag_value("repeat_name", $repeat_name);
	    $rf->add_tag_value("repeat_class",$repeat_class);

	    #$rf->repeat_consensus ($rc);

            push @repeat_features, $rf;
        }
    }
    close $filehandle;

    return \@repeat_features;
}


=head2  _setparams()

 Title   : _setparams
 Usage   : Internal function, not to be called directly
 Function:  Create parameter inputs for repeatmasker program
 Example :
 Returns : parameter string to be passed to repeatmasker
 Args    : name of calling object

=cut

sub _setparams {
    my ($attr, $value, $self);

    $self = shift;

    my $param_string = "";
    for  $attr ( @RM_PARAMS ) {
      $value = $self->$attr();
      next unless (defined $value);

      my $attr_key = lc $attr; #put params in format expected by dba

      $attr_key = ' -'.$attr_key;
      $param_string .= $attr_key.' '.$value;
    }

    for  $attr ( @RM_SWITCHES) {
      $value = $self->$attr();
      next unless ($value);
      my $attr_key = lc $attr; #put switches in format expected by dba
      $attr_key = ' -'.$attr_key;
      $param_string .= $attr_key ;
    }


    if ($self->quiet() || $self->verbose() < 0) {
      $param_string .= '  >/dev/null ';
    }
    return $param_string;
}


=head2  _setinput()

 Title   : _setinput
 Usage   : Internal function, not to be called directly
 Function: writes input sequence to file and return the file name
 Example :
 Returns : string 
 Args    : a Bio::PrimarySeqI compliant object

=cut

sub _setinput {
  my ($self,$seq) = @_;
  $seq->isa("Bio::PrimarySeqI") || 
      $self->throw("Need a Bio::PrimarySeq compliant object for RepeatMasker");
#  my  $in  = Bio::SeqIO->new(-file => $infilename , '-format' => 'Fasta');
  my ($tfh1,$outfile1) = $self->tempfile(-dir=>$TMPDIR);
  my $out1 = Bio::SeqIO->new(-fh=> $tfh1 , '-format' => 'Fasta');
  $out1->write_seq($seq);

  return ($outfile1);
}


1;


