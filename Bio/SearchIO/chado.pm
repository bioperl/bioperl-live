# $Id$
#
# BioPerl module for Bio::SearchIO::chado
#
# Chris Mungall <cjm@fruitfly.org>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::chado - chado sequence input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the SearchIO handler system. Go:

    $stream = Bio::SearchIO->new(-file => $filename, -format => 'chado');

    while ( my $seq = $stream->next_seq() ) {
	# do something with $seq
    }

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from chado flat
file databases. CURRENTLY ONLY TO


=head2 Optional functions

=over 3

=item _show_dna()

(output only) shows the dna or not

=item _post_sort()

(output only) provides a sorting func which is applied to the FTHelpers
before printing


=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://www.bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Chris Mungall

Email cjm@fruitfly.org

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SearchIO::chado;
use vars qw(@ISA);
use strict;

use Bio::SearchIO;
use Bio::SeqFeature::Generic;
use Bio::Seq::SeqFactory;
use Bio::Annotation::Collection;
use Bio::Annotation::Comment;
use Bio::Annotation::Reference;
use Bio::Annotation::DBLink;


use Bio::SeqIO::chado;

use Data::Stag qw(:all);

# should really inherit off of a chado helper...
@ISA = qw(Bio::SearchIO Bio::SeqIO::chado);
 
sub _initialize {
    my($self,@args) = @_;
    
    $self->SUPER::_initialize(@args); 
    my $wclass = $self->default_handler_class;
    $self->handler($wclass->new);
    $self->{_end_of_data} = 0;
    $self->handler->S("chado");
    return;
}

sub DESTROY {
    my $self = shift;
    $self->end_of_data();
    $self->SUPER::DESTROY();
}

sub end_of_data {
    my $self = shift;
    $self->{_end_of_data} = 1;
    $self->handler->E("chado");
}

sub default_handler_class {
    return "Data::Stag::BaseHandler";
} 

=head2 write_result

 Title   : write_result
 Usage   : $stream->write_result($result)
 Function: writes the $result object (must be result) to the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Result


=cut

sub write_result {
    my ($self,$result) = @_;
    
    if( !defined $result ) {
	$self->throw("Attempting to write with no result!");
    }
    
    my $w = $self->handler;
    $w->S("result");
#    my $result_temp_uid = $self->get_temp_uid($result);

    my @stats =
      (map {
          [analysisprop=>[
                          [pkey=>$_],
                          [pval=>$result->get_statistic($_)]]]
      } $result->available_statistics);
    my @params =
      (map {
          [analysisprop=>[
                          [pkey=>$_],
                          [pval=>$result->get_parameter($_)]]]
      } $result->available_parameters);

    my $cid = $self->get_temp_uid($result);
    $w->ev(companalysis=>[
                          [companalysis_id=>$cid],
                          [datasource=>$result->database_name],
                          @stats,
                          @params,
                         ]
          );
    while( my $hit = $result->next_hit ) {
	# process the Bio::Search::Hit::HitI object
        $self->write_hit($hit, $cid);
    }
    $w->E("result");
    return 1;
}

sub write_hit {
    my $self = shift;
    my $hit = shift;
    my $cid = shift;

    my $w = $self->handler;
    my $hit_id = $self->get_temp_uid($hit);

    # we should determine the type by the type of blast;
    # eg blastx gives CDS for hit and CDS_exon for HSP
    my $fnode =
      [feature=> [
                  [feature_id=>$hit_id],
                  [name=>$hit->name],
                  [typename=>"hit"],
                  [analysisfeature=>[
                                     [rawscore=>$hit->raw_score],
                                     [significance=>$hit->significance],
                                     [analysis_id=>$cid]]]]];
    $w->ev(@$fnode);
    foreach my $hsp ( $hit->hsps) {
        $self->write_hsp($hsp, $hit_id);
    }
    return 1;
}

sub write_hsp {
    my $self = shift;
    my $hsp = shift;
    my $hid = shift;

    my $w = $self->handler;
    my $hsp_id = $self->get_temp_uid($hsp);
    my $order = 0;
    my @lnodes =
      map {
          my ($nbeg, $nend, $strand) = 
            $self->bp2ib([$hsp->start($_),
                          $hsp->end($_),
                          $hsp->strand($_)
                         ]);
          my $src = $_ eq 'query' ? $hsp->query->seq_id : $hsp->hit->seq_id;
          [featureloc=>[
                        [nbeg=>$nbeg],
                        [nend=>$nend],
                        [strand=>$strand],
                        [srcfeature=>$src],
                        [group=>0],
                        [order=>$order++],
                       ]
          ]
      } qw(query subject);
    my $fnode =
      [feature => [
       
                   [feature_id=>$hsp_id],
                   [typename=>"hsp"],
                   [analysisfeature=>[
                                      [rawscore=>$hsp->score],
                                      [significance=>$hsp->significance],
                                     ]
                   ],
                   @lnodes,
                  ]
      ];
    $w->ev(@$fnode);
    $w->ev(feature_relationship=>[
                                  [subjfeature_id=>$hsp_id],
                                  [objfeature_id=>$hid]
                                 ]
          );
    return 1;
}


1;
