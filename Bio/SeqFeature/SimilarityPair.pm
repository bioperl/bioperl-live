# $Id$
#
# BioPerl module for Bio::SeqFeature::SimilarityPair
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::SimilarityPair - Sequence feature based on the similarity
                  of two sequences.

=head1 SYNOPSIS

  $sim_pair = Bio::SeqFeature::SimilarityPair->from_searchResult($blastHit);

  $sim = $sim_pair->query(); # a Bio::SeqFeature::Similarity object - the query
  $sim = $sim_pair->hit();   # dto - the hit.

  # some properties for the similarity pair
  $expect = $sim_pair->significance();
  $score = $sim_pair->score();
  $bitscore = $sim_pair->bits();

  # this will not write the description for the sequence (only its name)
  print $sim_pair->query()->gff_string(), "\n";

=head1 DESCRIPTION

Lightweight similarity search result as a pair of Similarity
features. This class inherits off Bio::SeqFeature::FeaturePair and
therefore implements Bio::SeqFeatureI, whereas the two features of the
pair are descendants of Bio::SeqFeature::Generic, with better support
for representing similarity search results in a cleaner way.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Hilmar Lapp

Email hlapp@gmx.net or hilmar.lapp@pharma.novartis.com

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::SimilarityPair;
use vars qw(@ISA);
use strict;

use Bio::SeqFeature::FeaturePair;
use Bio::SeqFeature::Similarity;
use Bio::Tools::Blast::Sbjct;

@ISA = qw(Bio::SeqFeature::FeaturePair);

=head2 new

 Title   : new
 Usage   : my $similarityPair = new Bio::SeqFeature::SimilarityPair
                                 (-hit   => $hit,
				  -query => $query,
				  -source => 'blastp');
 Function: Initializes a new SimilarityPair object
 Returns : Bio::SeqFeature::SimilarityPair
 Args    : -query => The query in a Feature pair 
           -hit   => (formerly '-subject') the subject/hit in a Feature pair


=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);

    # Hack to deal with the fact that SimilarityPair calls strand()
    # which will lead to an error in Bio::Search::HSP::BlastHSP 
    # because parsing hasn't yet occurred.
    # TODO: Remove this when BlastHSP doesn't do lazy parsing.
    $self->{'_initializing'} = 1;

    my ($hit, $query, $fea1, $source,$sbjct) =
	$self->_rearrange([qw(HIT
			      QUERY
			      FEATURE1
                              SOURCE
			      SUBJECT
			      )],@args);
    
    if( $sbjct ) { 
	$self->deprecated("use of -subject deprecated: SimilarityPair now uses 'hit'");
	if(! $hit) { $hit = $sbjct } 
	else { 
	    $self->warn("-hit and -subject were specified, using -hit and ignoring -subject");
	}
    }

    # make sure at least the query feature exists -- this refers to feature1
    if($query && ! $fea1) { $self->query( $query);  } 
    else { $self->query('null'); } # call with no args sets a default value for query
    
    $hit && $self->hit($hit);
    # the following refer to feature1, which has been ensured to exist
    $self->primary_tag('similarity') unless( defined $self->primary_tag() );
    $source && $self->source_tag($source);
    $self->strand(0) unless( defined $self->strand() );

    $self->{'_initializing'} = 0;  # See "Hack" note above
    return $self;
}

=head2 from_searchResult

 Title   : from_searchResult
 Usage   : $sim_pair = Bio::SeqFeature::SimilarityPair->from_searchResult($blast_obj);
           $sim_pair->from_searchResult($blast_obj);
 Function: This method creates or fills SimilarityPair objects from objects
           representing similarity search results.

           Since there is no public interface search result objects are
           required to implement, this method basically checks for the type
           of the object and dispatches the actual SimilarityPair creation
           to a method capable of this.

           At present, the following classes are recognized:
           Bio::Tools::Blast::Sbjct
           Bio::Tools::Blast::HSP
           An exception will be thrown if an object of an unrecognized class
           is passed.

           Note that this is probably the point where you will want to add
           your class if you have a method for creating SimilarityPair
           objects from it.

           Note that passing an object that has already previously been
           filled is potentially error-prone, because undefined fields
           will not be (re-)set to an undef value.

 Returns : The object created or filled.
 Args    : 

=cut

sub from_searchResult {
    my ($obj, $blastobj) = @_;

    if(! ref($obj)) {
	$obj = $obj->new();
    }
    if($blastobj->isa('Bio::Tools::Blast::Sbjct')) {
	return $obj->_from_BlastObj($blastobj);
    } elsif($blastobj->isa('Bio::Tools::Blast::HSP')) {
	return $obj->_from_BlastObj($blastobj);
    } else {
	$obj->throw("don't know how to handle object of " . ref($blastobj));
    }
}

=head2 _from_blastObject

 Title   : from_blastObject
 Usage   : $sim_pair = Bio::SeqFeature::SimilarityPair->_from_blastObj($blast_obj);
           $sim_pair->_from_blastObj($blast_obj);
 Function: See documentation for from_searchResult(). This one handles
           Bio::Tools::Blast::Sbjct and Bio::Tools::Blast::HSP objects.
 Returns : A SimilarityPair object.
 Args    : 

=cut

sub _from_BlastObj {
    my ($obj, $blastObj) = @_;

    if(! ref($obj)) {
	$obj = $obj->new();
    }
    my $simqu = $obj->query();
    my $simsb = $obj->hit();
    my $report = $blastObj->parent();
    my $blastSbjct = $blastObj;
    if($blastObj->isa('Bio::Tools::Blast::HSP')) {
	$report = $report->parent();
	$blastSbjct = $blastSbjct->parent();
    }
    #
    # set the overall properties: score, bits, E value
    #
    $obj->score($blastObj->score());
    $obj->bits($blastObj->bits());
    $obj->significance($blastObj->signif());
    $obj->source_tag($report->program());
    #
    # set the query and hit specific properties
    #
    # seq names
    $simqu->seqname($report->query());
    $simsb->seqname($blastSbjct->name());
    # seq descriptions
    $simqu->seqdesc($report->query_desc());
    $simsb->seqdesc($blastSbjct->desc());
    # start and end points
    $simqu->start($blastObj->start('query'));
    $simsb->start($blastObj->start('sbjct'));
    $simqu->end($blastObj->end('query'));
    $simsb->end($blastObj->end('sbjct'));
    # frame and strand
    if($blastObj->isa('Bio::Tools::Blast::HSP')) {
	$simqu->strand($blastObj->strand('query'));
	$simsb->strand($blastObj->strand('sbjct'));
	my $frm = $blastObj->frame();
	if($frm) {
	    # convert to 0,1,2 format
	    $frm = ($frm * $simqu->strand($blastObj->strand('query'))) -1;
	    $simqu->frame($frm);
	}
    } else {
	# scan through all HSPs, and if they agree set the respective field
	my @hsps = $blastObj->hsps();
	my $strandq = $hsps[0]->strand('query');
	my $strands = $hsps[0]->strand('sbjct');
	my $frm = $hsps[0]->frame();
	foreach my $hsp (@hsps) {
	    if(!$hsp->strand('query') || $hsp->strand('query') != $strandq) {
		$strandq = undef;
		last; # frames can hardly agree in this case
	    }
	    if(defined($strands) && ($hsp->strand('sbjct') != $strands)) {
		$strands = undef;
	    }
	    # we can assume that either all hsps have a frame or none
	    if($frm && ($frm != $hsp->frame())) {
		$frm = undef;
	    }
	}
	$simqu->strand($strandq) if defined($strandq);
	$simsb->strand($strands) if defined($strands);
	$simqu->strand($frm) if $frm;
    }
    # 'percent' identity
    $simqu->frac_identical($blastObj->frac_identical('query'));
    $simsb->frac_identical($blastObj->frac_identical('sbjct'));
    #
    # return the created or fille object
    #
    return $obj;
}

#
# Everything else is just inherited from SeqFeature::FeaturePair.
#

=head2 query

 Title   : query
 Usage   : $query_feature = $obj->query();
           $obj->query($query_feature);
 Function: 
 Returns : 
 Args    : 


=cut

sub query {
    my ($self, @args) = @_;
    my $f = $self->feature1();
    if( ! @args || ( !ref($args[0]) && $args[0] eq 'null') ) {
	if( ! defined( $f) ) {
	    @args = Bio::SeqFeature::Similarity->new();	    
	} elsif( ! $f->isa('Bio::SeqFeature::Similarity') && 
		 $f->isa('Bio::SeqFeatureI') ) {
	    # a Bio::SeqFeature::Generic was placeholder for feature1	    
	    my $newf = new 
	      Bio::SeqFeature::Similarity( -start   => $f->start(),
					   -end     => $f->end(),
					   -strand  => $f->strand(),
					   -primary => $f->primary_tag(),
					   -source  => $f->source_tag(),
					   -seqname => $f->seqname(),
					   -score   => $f->score(),
					   -frame   => $f->frame(),
					   );
	    foreach my $tag ( $newf->all_tags ) {
		$tag->add_tag($tag, $newf->each_tag($tag));
	    }
	    @args = $newf;	   
	} else {
	    @args = ();
	}
    }
    return $self->feature1(@args);
}




=head2 subject

 Title   : subject
 Usage   : $sbjct_feature = $obj->subject();
           $obj->subject($sbjct_feature);
 Function: Get/Set Subject for a SimilarityPair 
 Returns : Bio::SeqFeature::Similarity
 Args    : [optional] Bio::SeqFeature::Similarity
 Notes   : Deprecated.  Use the method 'hit' instead

=cut

sub subject { 
    my $self = shift;
    $self->deprecated("Method subject deprecated: use hit() instead");
    $self->hit(@_); 
}

=head2 hit

 Title   : hit
 Usage   : $sbjct_feature = $obj->hit();
           $obj->hit($sbjct_feature);
 Function: Get/Set Hit for a SimilarityPair 
 Returns : Bio::SeqFeature::Similarity
 Args    : [optional] Bio::SeqFeature::Similarity


=cut

sub hit {
    my ($self, @args) = @_;
    my $f = $self->feature2();
    if(! @args || (!ref($args[0]) && $args[0] eq 'null') ) {
	if( ! defined( $f) ) {
	    @args = Bio::SeqFeature::Similarity->new();
	} elsif( ! $f->isa('Bio::SeqFeature::Similarity') && 
		 $f->isa('Bio::SeqFeatureI')) {
	    # a Bio::SeqFeature::Generic was placeholder for feature2	    
	    my $newf = new 
	      Bio::SeqFeature::Similarity( -start   => $f->start(),
					   -end     => $f->end(),
					   -strand  => $f->strand(),
					   -primary => $f->primary_tag(),
					   -source  => $f->source_tag(),
					   -seqname => $f->seqname(),
					   -score   => $f->score(),
					   -frame   => $f->frame(),
					   );
	    foreach my $tag ( $newf->all_tags ) {
		$tag->add_tag($tag, $newf->each_tag($tag));
	    }
	    @args = $newf;
	}
    }
    return $self->feature2(@args);
}

=head2 source_tag

 Title   : source_tag
 Usage   : $source = $obj->source_tag(); # i.e., program
           $obj->source_tag($evalue);
 Function: 
 Returns : 
 Args    : 


=cut

sub source_tag {
    my ($self, @args) = @_;

    if(@args) {
	$self->hit()->source_tag(@args);
    }
    return $self->query()->source_tag(@args);
}

=head2 significance

 Title   : significance
 Usage   : $evalue = $obj->significance();
           $obj->significance($evalue);
 Function: 
 Returns : 
 Args    : 


=cut

sub significance {
    my ($self, @args) = @_;

    if(@args) {
	$self->hit()->significance(@args);
    }
    return $self->query()->significance(@args);
}

=head2 score

 Title   : score
 Usage   : $score = $obj->score();
           $obj->score($value);
 Function: 
 Returns : 
 Args    : 


=cut

sub score {
    my ($self, @args) = @_;

    if(@args) {
	$self->hit()->score(@args);
    }
    return $self->query()->score(@args);
}

=head2 bits

 Title   : bits
 Usage   : $bits = $obj->bits();
           $obj->bits($value);
 Function: 
 Returns : 
 Args    : 


=cut

sub bits {
    my ($self, @args) = @_;

    if(@args) {
	$self->hit()->bits(@args);
    }
    return $self->query()->bits(@args);
}

1;
