
# BioPerl module for Bio::Tools::Seg
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
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Seg

=head1 SYNOPSIS

  use Bio::Tools::Seg;
  my $parser = new Bio::Tools::Seg(-fh =>$filehandle );
  while( my $seg_feat = $parser->next_result ) {
        #do something
        #eg
        push @seg_feat, $seg_feat;
  }

=head1 DESCRIPTION

Parser for Seg output

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
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


=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::Seg;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SeqFeature::FeaturePair;
use Bio::Root::IO;
use Bio::SeqFeature::Generic;
@ISA = qw(Bio::Root::Root Bio::Root::IO);





=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Seg();
 Function: Builds a new Bio::Tools::Seg object
 Returns : Bio::Tools::Seg
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
 Usage   : my $feat = $seg->next_result
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
           if ($line=~/^\>/) { #if it is a line starting with a ">"
               $line=~/^\>\s*(\S+)\s*\((\d+)\-(\d+)\)\s*complexity=(\S+)/;
               my $id = $1;
               my $start = $2;
               my $end = $3;
               my $score = $4;

               #for example in this line test_prot(214-226) complexity=2.26 (12/2.20/2.50)
               #$1 is test_prot  $2 is 214 $3 is 226 and $4 is 2.26

               my (%feature);
               $feature{name} = $id;
               $feature{score} = $score;
               $feature{start} = $start;
               $feature{end} = $end;
               $feature{source} = "Seg";
               $feature{primary} = 'low_complexity';
               $feature{program} = "Seg";
               $feature{logic_name} = 'low_complexity';
               my $new_feat =  $self->create_feature (\%feature);
               return $new_feat;
            }
          next;
        }

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

1;


