# BioPerl module for Bio::Tools::Run::Pseudowise
#
# Cared for by
#
# Copyright Kiran 
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::Pseudowise - Object for prediting pseudogenes in a
given sequence given a protein and a cdna sequence 

=head1 SYNOPSIS

  #  Build a pseudowise alignment factory
  my  $factory = Bio::Tools::Run::Pseudowise->new();

  #  Pass the factory 3 Bio:SeqI objects (in the order of query peptide and cdna and
     target_genomic) 
  #@genes is an array of GenericSeqFeature objects
  my @genes = $factory->predict_genes($seq1, $seq2, $seq3); 

=head1 DESCRIPTION

Pseudowise is a pseudogene predition program developed by Ewan Birney
http://www.sanger.ac.uk/software/wise2.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Kiran 

Email kiran@fugu-sg.org 

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


package Bio::Tools::Run::Pseudowise;
use vars qw($AUTOLOAD @ISA $PROGRAM $PROGRAMDIR
            $TMPDIR $TMPOUTFILE @DBA_SWITCHES 
            @OTHER_SWITCHES %OK_FIELD);
use strict;
use Bio::SeqIO;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::Root::Root;
use Bio::Root::IO;
use Bio::Factory::ApplicationFactoryI;
use Bio::Search::HSP::GenericHSP;

@ISA = qw(Bio::Root::Root Bio::Root::IO Bio::Factory::ApplicationFactoryI);

# You will need to enable pseudowise to find the pseudowise program. This
# can be done in (at least) two ways:
#
# 1. define an environmental variable WISEDIR 
# export WISEDIR =/usr/local/share/wise2.2.0
# where the wise2.2.20 package is installed
#
# 2. include a definition of an environmental variable WISEDIR in
# every script that will use DBA.pm
# $ENV{WISEDIR} = '/usr/local/share/wise2.2.20';

BEGIN {

    if (defined $ENV{WISEDIR}) {
        $PROGRAMDIR = $ENV{WISEDIR} || '';
        $PROGRAM = Bio::Root::IO->catfile($PROGRAMDIR."/src/bin/",
                                          'dba'.($^O =~ /mswin/i ?'.exe':''));
    }
    else {
        $PROGRAM = 'pseudowise';
    }
}

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
  unless ($self->exists_pseudowise()) {
    if( $self->verbose >= 0 ) {
      warn "Pseudowise program not found as ".$self->program." or not executable. \n Pseudowise can be obtained from http://www.sanger.ac.uk/software/wise2\n"; 
    }
  }

  return $self;
}


=head2  exists_pseudowise()

 Title   : exists_pseudowis_
 Usage   : $pseudowisefound= Bio::Tools::Run::Pseudowise->exists_pseudowise()
 Function: Determine whether pseudowise program can be found on current host
 Example :
 Returns : 1 if pseudowise program found at expected location, 0 otherwise.
 Args    :  none

=cut


sub exists_pseudowise{
    #my $self = shift;
    my ($self) = @_;
    my $rt = Bio::Root::IO->new();
    #if( my $f = Bio::Root::IO->exists_exe($PROGRAM) ) {
    if( my $f = $rt->exists_exe($PROGRAM) ) {
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
 Usage   : exit if $prog->version() < 1.8
 Function: Determine the version number of the program
 Example :
 Returns : float or undef
 Args    : none

=cut

sub version {
    my ($self) = @_;

    return undef unless $self->exists_pseudowise;
    my $string = `pseudowise -- ` ;
    $string =~ /\(([\d.]+)\)/;
    return $1 || undef;

}

=head2 predict_genes 

 Title   : predict_genes 
 Usage   :
            3 sequence objects 
            @feats = $factory->predict_genes($seq1, $seq2, $seq3);
 
Function: Predict pseudogenes


 Returns : An array of Bio::Seqfeature::Generic objects 
 Args    : Name of a file containing a set of 3 fasta sequences in the order of 
           peptide, cdna and genomic sequences
           or else 3  Bio::Seq objects.

 Throws an exception if argument is not either a string (eg a
 filename) or 3 Bio::Seq objects.  If
 arguments are strings, throws exception if file corresponding to string
 name can not be found. 

=cut

sub predict_genes {

    my ($self,$seq1, $seq2, $seq3)=@_; 
    my ($attr, $value, $switch);

# Create input file pointer
    my ($infile1,$infile2,$infile3)= $self->_setinput($seq1, $seq2, $seq3);
    if (!($infile1 && $infile2 && $infile3)) {$self->throw("Bad input data (sequences need an id ) ");}


# run pseudowise 
    my @feats = $self->_run($infile1,$infile2,$infile3);
    return @feats;
}


=head2  _run

 Title   :  _run
 Usage   :  Internal function, not to be called directly
 Function:   makes actual system call to a pseudowise program
 Example :
 Returns : nothing; pseudowise  output is written to a
           temporary file $TMPOUTFILE
 Args    : Name of a files containing 3 sequences in the order of peptide, cdna and genomic 

=cut

sub _run {
    my ($self,$infile1,$infile2,$infile3) = @_;
    my $instring;
    $self->debug( "Program ".$self->program."\n");
    #my $outfile = $self->outfile() || $TMPOUTFILE ;
    my $outfile =  $TMPOUTFILE ;
    my $commandstring = $self->program." $infile1 $infile2 $infile3> $outfile";
    $self->debug( "pseudowise command = $commandstring");
    my $status = system($commandstring);
    $self->throw( "Pseudowise call ($commandstring) crashed: $? \n") unless $status==0;

    #parse the outpur and return a Bio::Seqfeature array
    my $genes   = $self->_parse_results($outfile);

    return @{$genes};
}

=head2  _parse_results

 Title   :  __parse_results
 Usage   :  Internal function, not to be called directly
 Function:  Parses pseudowise output 
 Example :
 Returns : an reference to an array of Seqfeatures 
 Args    : the name of the output file 

=cut

sub _parse_results {
    my ($self,$outfile) = @_;
    $outfile||$self->throw("No outfile specified");
    my ($self) = @_;

    print "Parsing the file\n";

    my $filehandle;
    if (ref ($outfile) !~ /GLOB/)
    {
        open (PSEUDOWISE, "<".$outfile)
            or $self->throw ("Couldn't open file ".$outfile.": $!\n");
        $filehandle = \*PSEUDOWISE;
    }
    else
    {
        $filehandle = $outfile;
    }

    my @genes;
    #The big parsing loop - parses exons and predicted peptides

            while (<$filehandle>)
            {
                my $gene;
                if (/Gene/i)
                {
                  $gene = new Bio::SeqFeature::Generic (
                                -primary => 'pseudogene',
                                -source => 'pseudowise');
                  push @genes, $gene;

                  while(<$filehandle>) {
                    my @gene_elements = split;
                    my $no = scalar(@gene_elements);
                    if ((/Gene/i) && $no == 3) {
                      my @element = split;
                      my $no = scalar(@element);
                      my $gene_start = $element[1];
                      my $gene_end = $element[2];
                      print "No. of elements in Gene-- $no\n";
                      print "Gene start - $gene_start\n";
                      print "Gene End - $gene_end\n";
                      $gene->start($gene_start);
                      $gene->end($gene_end);
                    }
                    elsif (/Exon/i) {
                      my @element = split;
                      my $no = scalar(@element);
                      my $exon_start = $element[1];
                      my $exon_end = $element[2];
                      my $exon_phase = $element[4];
                      print "No. of elements in Exon-- $no\n";
                      print "Exon start - $exon_start\n";
                      print "Exon End - $exon_end\n";
                      print "Exon Phase - $exon_phase\n";
                      my $exon = new Bio::SeqFeature::Generic (
                                -start => $exon_start,
                                -end => $exon_end,
                                -primary => 'exon',
                                -source => 'pseudowise',
                                -frame  => $exon_phase);
                      $gene->add_sub_SeqFeature($exon);
                    }
                    elsif ((/Gene/i) && $no != 3) {
                       $gene = new Bio::SeqFeature::Generic (
                                -primary => 'pseudogene',
                                -source => 'pseudowise');
                       push @genes, $gene;

                    }
                  }
                }
            }
    return \@genes;
}
   
=head2  _setinput()

 Title   :  _setinput
 Usage   :  Internal function, not to be called directly
 Function:   Create input files for pseudowise program
 Example :
 Returns : name of file containing dba data input
 Args    : Seq objects or input file name

=cut

sub _setinput {
  my ($self, $seq1, $seq2, $seq3) = @_;
  my ($tfh1,$tfh2,$tfh3,$outfile1,$outfile2,$outfile3);

    if(!($seq1->isa("Bio::PrimarySeqI") && $seq2->isa("Bio::PrimarySeqI")&& $seq2->isa("Bio::PrimarySeqI"))) 
      { $self->throw("One or more of the sequences are nor Bio::PrimarySeqI objects\n"); }

    ($tfh1,$outfile1) = $self->tempfile(-dir=>$TMPDIR);
    ($tfh2,$outfile2) = $self->tempfile(-dir=>$TMPDIR);
    ($tfh3,$outfile3) = $self->tempfile(-dir=>$TMPDIR);

    my $out1 = Bio::SeqIO->new(-fh=> $tfh1 , '-format' => 'Fasta');
    my $out2 = Bio::SeqIO->new(-fh=> $tfh2 , '-format' => 'Fasta');
    my $out3 = Bio::SeqIO->new(-fh=> $tfh3 , '-format' => 'Fasta');

    $out1->write_seq($seq1);
    $out2->write_seq($seq2);
    $out3->write_seq($seq3);
    $self->_query_pep_seq($seq1);
    $self->_query_cdna_seq($seq2);
    $self->_subject_dna_seq($seq3);
    return $outfile1,$outfile2,$outfile3;
  
}


=head2  _query_pep_seq()

 Title   :  _query_pep_seq
 Usage   :  Internal function, not to be called directly
 Function:  get/set for the query sequence 
 Example :
 Returns : 
 Args    : 

=cut

sub _query_pep_seq {
  my ($self,$seq) = @_;
  if(defined $seq){
    $self->{'_query_pep_seq'} = $seq;
  }
  return $self->{'_query_pep_seq'};
}

=head2  _query_cdna_seq()

 Title   :  _query_cdna_seq
 Usage   :  Internal function, not to be called directly
 Function:  get/set for the query sequence
 Example :
 Returns :
 Args    :

=cut

sub _query_cdna_seq {
  my ($self,$seq) = @_;
  if(defined $seq){
    $self->{'_query_cdna_seq'} = $seq;
  }
  return $self->{'_query_cdna_seq'};
}

=head2  _subject_dna_seq()

 Title   :  _subject_dna_seq
 Usage   :  Internal function, not to be called directly
 Function:  get/set for the subject sequence
 Example :
 Returns :

 Args    :

=cut

sub _subject_dna_seq {
  my ($self,$seq) = @_;
  if(defined $seq){
    $self->{'_subject_dna_seq'} = $seq;
  }
  return $self->{'_subject_dna_seq'};
}
1; # Needed to keep compiler happy
  
