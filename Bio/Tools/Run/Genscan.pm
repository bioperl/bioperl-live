# BioPerl module for Bio::Tools::Run::Genscan
#
# Cared for by
#
# Copyright Balamurugan Kumarasamy
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::Genscan - Object for identifying genes in a
given sequence given a matrix(for appropriate organisms).

=head1 SYNOPSIS

  # Build a Genscan factory
  my $param = ('MATRIX'=>HumanIso.smat);
  my  $factory = Bio::Tools::Run::Genscan->new($param);

  # Pass the factory a Bio::Seq object
  #@genes is an array of Bio::Tools::Predictions::Gene objects
  my @genes = $factory->predict_genes($seq);

=head1 DESCRIPTION

 Genscan is a gene identifying program developed by Christopher Burge
http://genes.mit.edu/burgelab/

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

=head1 AUTHOR - Bala

Email savikalpa@fugu-sg.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Tools::Run::Genscan;
 
use vars qw($AUTOLOAD @ISA $PROGRAM  $PROGRAMDIR                
            $TMPDIR $TMPOUTFILE  @GENSCAN_PARAMS
             %OK_FIELD);
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::Root;
use Bio::Root::IO;
use Bio::Factory::ApplicationFactoryI;
use Bio::Tools::Genscan;

@ISA = qw(Bio::Root::Root Bio::Root::IO Bio::Factory::ApplicationFactoryI);



BEGIN {

   if (defined $ENV{GENSCANDIR}) {
        $PROGRAMDIR = $ENV{GENSCANDIR} || '';
        $PROGRAM = Bio::Root::IO->catfile($PROGRAMDIR,
                                          'genscan'.($^O =~ /mswin/i ?'.exe':''));
    }
    else {                                                                 
        $PROGRAM = 'genscan';
    }
              
     @GENSCAN_PARAMS=qw(MATRIX);
      foreach my $attr ( @GENSCAN_PARAMS)
                        { $OK_FIELD{$attr}++; }
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

sub new {                                                       
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    # to facilitiate tempfile cleanup
    $self->_initialize_io();

    my ($attr, $value);
    (undef,$TMPDIR) = $self->tempdir(CLEANUP=>1);
    (undef,$TMPOUTFILE) = $self->tempfile(-dir => $TMPDIR);
    while (@args)  {
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
    unless ($self->exists_genscan()) {
        if( $self->verbose >= 0 ) {
            warn "Genscan program not found as ".$self->program." or not executable. \n ";
      }
    }

    return $self;
}


=head2  exists_genscan()

 Title   : exists_genscan
 Usage   : $genscanfound = Bio::Tools::Run::Genscan->exists_genscan()
 Function: Determine whether genscan program can be found on current host
 Returns : 1 if genscan program found at expected location, 0 otherwise.
 Args    : none

=cut


sub exists_genscan {
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


=head2 predict_genes()
    Title   :   predict_genes()
    Usage   :   $obj->predictgenes($matrixFile,$seqFile)
    Function:   Runs genscan and creates an array of Genes
    Returns :   An array of Bio::Tools::Prediction::Gene objects
    Args    :   A Bio::PrimarySeqI    

=cut
sub predict_genes() {
    my ($self,$seq) = @_;
    my $infile1 = $self->_writeSeqFile($seq);  
    $self->_set_input($infile1);  
    print"matrix". $self->MATRIX;
    my @feat = $self->_run();
    return @feat;
}

=head2 _run

    Title   :   _run
    Usage   :   $obj->_run()
    Function:   Internal(not to be used directly)
    Returns :   An array of Bio::Tools::Prediction::Gene objects
    Args    :   

=cut

sub _run {                                                           

    my ($self) = @_;
    my @genes;
    my $gene;
    my $result = $TMPOUTFILE;

    my $str = $self->program.' '.$self->MATRIX.' '.$self->{'input'}." > $result";
    system($str);
    $self->throw($result." not created by Genscan\n") unless (-e $result);
    my $genScanParser = Bio::Tools::Genscan->new(-file => $result);
    
    while( $gene = $genScanParser->next_prediction()){
       push(@genes, $gene);
    }     
   return @genes;
}

=head2 _set_input()

    Title   :   _set_input
    Usage   :   obj->_set_input($matrixFile,$seqFile)
    Function:   Internal(not to be used directly)
    Returns :  
    Args    :

=cut
sub _set_input() {
   my ($self,$infile1) = @_;
   $self->{'input'}=$infile1;
}

=head2 _writeSeqFile()

    Title   :   _writeSeqFile
    Usage   :   obj->_writeSeqFile($seq)
    Function:   Internal(not to be used directly)
    Returns :   
    Args    : 

=cut


sub _writeSeqFile(){
  my ($self,$seq) = @_;
  my($tfh,$inputfile);
 ($tfh,$inputfile) = $self->tempfile(-dir=>$TMPDIR);
  my $in  = Bio::SeqIO->new(-fh => $tfh , '-format' => 'Fasta');
  $in->write_seq($seq);
 
  return $inputfile;
}

1;
