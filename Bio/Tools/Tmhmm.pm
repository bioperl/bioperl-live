
# BioPerl module for Bio::Tools::Tmhmm
#
# Copyright Balamurugan Kumarasamy
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#
# Copyright 
#
# You may distribute this module under the same terms as perl itself

=head1 NAME

Bio::Tools::Tmhmm - parse TmHMM output (transmembrane HMM)

=head1 SYNOPSIS

  use Bio::Tools::Tmhmm;
  my $parser = new Bio::Tools::Tmhmm(-fh =>$filehandle );
  while( my $tmhmm_feat = $parser->next_result ) {
     #do something
     #eg
     push @tmhmm_feat, $tmhmm_feat;
  }

=head1 DESCRIPTION

Parser for Tmhmm output

=head1 FEEDBACK

=head2 Mailing Lists

user feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Bala

Email savikalpa@fugu-sg.org


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::Tmhmm;
use vars qw(@ISA);
use strict;

use Bio::Tools::AnalysisResult;
use Bio::Root::Root;
use Bio::SeqFeature::FeaturePair;
use Bio::Root::IO;
use Bio::SeqFeature::Generic;
@ISA = qw(Bio::Root::Root Bio::Root::IO Bio::Tools::AnalysisResult);



=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Tmhmm();
 Function: Builds a new Bio::Tools::Tmhmm object
 Returns : Bio::Tools::Tmhmm
 Args    : -fh/-file => $val, # for initing input, see Bio::Root::IO


=cut

sub new {
      my($class,@args) = @_;

      my $self = $class->SUPER::new(@args);
      $self->_initialize_io(@args);

      return $self;
}


=head2 next_result

 Title   : next_result
 Usage   : my $feat = $Tmhmm->next_result
 Function: Get the next result set from parser data
 Returns : Bio::SeqFeature::Generic
 Args    : none


=cut

sub next_result {
        my ($self) = @_;

        my $line;

        # parse
        my $id;
        while ($_=$self->_readline()) { 
           $line = $_;
           chomp $line;


           next if /^$/;
           if ($line=~/^#\s+(\S+)/) { 
                   #if the line starts with a '#' for example in # 13 Length: 522 
                   #assign 13 as the id.

                    $id = $1;
                    my ($junk, $values) = split /:/;
                   $self->seqname($id);
                    next;
           }

           elsif ($line=~/^(\S+)\s+(\S+)\s+(\w+)\s+(\d+)\s+(\d+)/) {

                    # Example :-  13      TMHMM2.0        inside       1   120
                    # assign $orien(inside) $start(1) and $end(120)


                    my $orien = $3;
                    my $start = $4;
                    my $end = $5;
                    $orien = uc ($orien);

                    if ($orien eq "TMHELIX") {
                         my (%feature);
                         $feature{name} = $self->seqname;
                         $feature{start} = $start;
                         $feature{end} = $end;
                         $feature{source} ='tmhmm';
                         $feature{primary}= 'transmembrane';
                         $feature{program} ='tmhmm';
                         $feature{logic_name} = 'TMHelix';
                         my $new_feat= $self->create_feature(\%feature);
                         return $new_feat;
                    }
                    next;
           }
           next;
        }
}

=head2 create_feature

 Title   : create_feature
 Usage   : obj->create_feature(\%feature)
 Function: Internal(not to be used directly)
 Returns : A Bio::SeqFeature::Generic object
 Args    :

=cut

sub create_feature {
       my ($self, $feat) = @_;


       # create feature object
       my $feature = Bio::SeqFeature::Generic->new(-seqname     => $feat->{name},
                                                 -start       => $feat->{start},
                                                 -end         => $feat->{end},
                                                 -score       => $feat->{score},
                                                 -source      => $feat->{source},
                                                 -primary     => $feat->{primary},
                                                 -logic_name  => $feat->{logic_name}, 
                                               );
       return $feature;
   }

=head2 seqname

 Title   :   seqname
 Usage   :   obj->seqname($seqname)
 Function:   Internal(not to be used directly)
 Returns :
 Args    :   seqname

=cut

sub seqname{
    my ($self,$seqname)=@_;

    if (defined$seqname){

        $self->{'seqname'}=$seqname;
    }

    return $self->{'seqname'};

}


1;


