# $Id$
#
# BioPerl module for Bio::SeqIO::chado
#
# Chris Mungall <cjm@fruitfly.org>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::chado - chado sequence input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the SeqIO handler system. Go:

    $stream = Bio::SeqIO->new(-file => $filename, -format => 'chado');

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

package Bio::SeqIO::chado;
use vars qw(@ISA);
use strict;

use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::Species;
use Bio::Seq::SeqFactory;
use Bio::Annotation::Collection;
use Bio::Annotation::Comment;
use Bio::Annotation::Reference;
use Bio::Annotation::DBLink;

use Data::Stag qw(:all);

@ISA = qw(Bio::SeqIO);
 
sub _initialize {
    my($self,@args) = @_;
    
    $self->SUPER::_initialize(@args); 
    if( ! defined $self->sequence_factory ) {
	$self->sequence_factory(new Bio::Seq::SeqFactory
				(-verbose => $self->verbose(), 
				 -type => 'Bio::Seq::RichSeq'));
    }
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

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    :

=cut

sub next_seq {
    my ($self,@args) = @_;
    my $seq = $self->sequence_factory->create
	(
         #         '-verbose' =>$self->verbose(), 
         #	 %params,
         #	 -seq => $seqc,
         #	 -annotation => $annotation,
         #	 -features => \@features
        );
    return $seq;
}

sub handler {
    my $self = shift;
    $self->{_handler} = shift if @_;
    return $self->{_handler};
}


=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq($seq)
 Function: writes the $seq object (must be seq) to the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq


=cut

sub write_seq {
    my ($self,$seq) = @_;
    
    if( !defined $seq ) {
	$self->throw("Attempting to write with no seq!");
    }
    
    if( ! ref $seq || ! $seq->isa('Bio::SeqI') ) {
	$self->warn(" $seq is not a SeqI compliant module. Attempting to dump, but may fail!");
    }

    # get a handler - must inherit from Data::Stag::BaseHandler;
    my $w = $self->handler;

    # start of data
    $w->S("seqset");

    #    my $seq_temp_uid = $self->get_temp_uid($seq);

    my $seq_temp_uid = $seq->accession . '.' . ($seq->can('seq_version') ? $seq->seq_version : $seq->version);

    # data structure representing the core sequence for this record
    my $seqnode =
      Data::Stag->new(feature=>[
                                    [feature_id=>$seq_temp_uid],
                                    [dbxrefstr=>$seq->accession_number],
                                    [name=>$seq->display_name],
                                    [residues=>$seq->seq],
                                   ]);
    
    # soft properties
    my %prop = ();

    my ($div, $mol);
    my $len = $seq->length();

    if ( $seq->can('division') ) {
	$div=$seq->division;
    } 
    if( !defined $div || ! $div ) { $div = 'UNK'; }

    if( !$seq->can('molecule') || ! defined ($mol = $seq->molecule()) ) {
	$mol = $seq->alphabet || 'DNA';
    }
    
    
    my $circular = 'linear  ';
    $circular = 'circular' if $seq->is_circular;

    # cheeky hack - access symbol table
    no strict 'refs';
    map {
        $prop{$_} = 
          $ {*$_};
    } qw(mol div circular);
    use strict 'refs';

    map {
        $prop{$_} = $seq->$_() if $seq->can($_);
    } qw(desc keywords);

    local($^W) = 0;   # supressing warnings about uninitialized fields.
    
    # Organism lines
    if (my $spec = $seq->species) {
        my ($species, $genus, @class) = $spec->classification();
	my $OS;
        if( $spec->common_name ) {
	    $OS = $spec->common_name;
	} else { 
	    $OS = "$genus $species";
	}
        if (my $ssp = $spec->sub_species) {
            $OS .= " $ssp";
        }
    }
    
    # Reference lines
    my $count = 1;
    foreach my $ref ( $seq->annotation->get_Annotations('reference') ) {
        # TODO
    }
    # Comment lines
    
    foreach my $comment ( $seq->annotation->get_Annotations('comment') ) {
        $seqnode->add_featureprop([[pkey=>'comment'],[pval=>$comment->text]]);
    }

    # throw the writer an event
    $w->ev(@$seqnode);

    $seqnode = undef;      # free memory

    # make events for all the features within the record
    foreach my $sf ( $seq->top_SeqFeatures ) {
        $self->write_sf($sf, $seq_temp_uid);
    }

    # data end
    $w->E("seqset");
    return 1;
}

# ----
# writes a seq feature
# ----

sub write_sf {
    my $self = shift;
    my $sf = shift;
    my $seq_temp_uid = shift;

    my $w = $self->handler;

    my %props =
      map {
          $_=>[$sf->each_tag_value($_)]
      } $sf->all_tags;

    my $loc = $sf->location;
    my $name = $sf->display_name;
    my $type = $sf->primary_tag;
    my @subsfs = $sf->sub_SeqFeature;
    my @locnodes = ();
    my $sid = $loc->is_remote ? $loc->seq_id : $seq_temp_uid;
    if( $loc->isa("Bio::Location::SplitLocationI") ) {
        # turn splitlocs into subfeatures
        my $n = 1;
        push(@subsfs,
             map {
                 my $ssf =
                   Bio::SeqFeature::Generic->new(

                                                 -start=>$_->start,
                                                 -end=>$_->end,
                                                 -strand=>$_->strand,
                                                 -primary=>$self->subpartof($type),
                                              );
                 if ($_->is_remote) {
                     $ssf->location->is_remote(1);
                     $ssf->location->seq_id($_->seq_id);
                 }
                 $ssf;
               } $loc->each_Location);
    }
    elsif( $loc->isa("Bio::Location::RemoteLocationI") ) {
        # turn splitlocs into subfeatures
        my $n = 1;
        push(@subsfs,
             map {
                 Bio::SeqFeature::Generic->new(
#                                               -name=>$name.'.'.$n++,
                                               -start=>$_->start,
                                               -end=>$_->end,
                                               -strand=>$_->strand,
                                               -primary=>$self->subpartof($type),
                                              )
               } $loc->each_Location);
    }
    else {
        my ($beg, $end, $strand) = $self->bp2ib($loc);
        @locnodes = (
                     [featureloc=>[
                                   [nbeg=>$beg],
                                   [nend=>$end],
                                   [strand=>$strand],
                                   [srcfeature_id=>$sid],
                                   [group=>0],
                                   [rank=>0],
                                  ]
                     ]
                    );
    }
    my $feature_id = $self->get_temp_uid($sf);
    
    my $fnode =
      [feature=>[
                 [feature_id=>$feature_id],
                 [name=>$name],
                 [typename=>$type],
                 @locnodes,
                 (map {
                     my $k = $_;
                     map { [featureprop=>[[pkey=>$k],[pval=>$_]]] } @{$props{$k}}
                 } keys %props),
                ]];
    $w->ev(@$fnode);

    foreach my $ssf (@subsfs) {
        my $ssfid = $self->write_sf($ssf, $sid);
        $w->ev(feature_relationship=>[
                                      [subjfeature_id=>$ssfid],
                                      [objfeature_id=>$feature_id]
                                     ]
              );
    }
    return $feature_id;
}

# private;
# an ID for this session that should be
# unique... hmm
sub session_id {
    my $self = shift;
    $self->{_session_id} = shift if @_;
    if (!$self->{_session_id}) {
        $self->{_session_id} = $$.time;
    }
    return $self->{_session_id};
}


our $next_id = 1;
our %obj2id_hash = ();
sub get_temp_uid {
    my $self = shift;
    my $ob = shift;
    my $id = $obj2id_hash{$ob};
    if (!$id) {
        $id = $next_id++;
        $obj2id_hash{$ob} = $id;
    }
    return $self->session_id.'.'.$id;
}

# interbase and directional semantics
sub bp2ib {
    my $self = shift;
    my $loc = shift;
    my ($s, $e, $str) = 
      ref($loc) eq "ARRAY" ? (@$loc) : ($loc->start, $loc->end, $loc->strand);
    if ($str < 0) {
        ($s, $e) = ($e, $s);
    }
    $s--;
    return ($s, $e, $str);
}

sub subpartof {
    my $self = shift;
    my $type = 'partof_'.shift;
    $type =~ s/partof_CDS/CDS_exon/;
    $type =~ s/partof_\wRNA/exon/;
    return $type;
}

1;
