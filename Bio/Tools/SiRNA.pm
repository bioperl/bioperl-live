# $Id$
#
# BioPerl module for Bio::Tools::SiRNA
#
# Cared for by Donald Jackson, donald.jackson@bms.com
#
# Copyright Donald Jackson
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

SiRNA - Perl object for designing small inhibitory RNAs.

=head1 SYNOPSIS

  use Bio::Tools::SiRNA;

  my $sirna_designer = Bio::Tools::SiRNA->new( -target => $bio_rich_seq );
  my @pairs = $sirna_designer->design;

  foreach $pair (@pairs) {
      my $sense_oligo_sequence = $pair->sense->seq;
      my $antisense_oligo_sequence = $pair->antisense->seq;

      # print out results
      print join ("\t", $pair->start, $pair->end, $pair->rank,
                  $sense_oligo_sequence, $antisense_oligo_sequence), "\n";
  }

=head1 DESCRIPTION

Package for designing SiRNA reagents.

Input is a Bio::RichSeq object (the target).

Output is a list of Bio::SeqFeature::SiRNA::Pair objects, which are
added to the feature table of the target sequence.  Each
Bio::SeqFeature::SiRNA::Pair contains two subfeatures
(Bio::SeqFeature::Oligo objects) which correspond to the individual
oligos.  These objects provide accessors for the information on the
individual reagent pairs.

This module implements the design rules described by Tuschl and
colleagues (see
http://www.mpibpc.gwdg.de/abteilungen/100/105/sirna.html).  They
describe three rules for RNAi oligos, which I label as rank 1 (best),
2, and 3 (worst).

I added two modifications: SiRNAs that overlap known SNPs (identified
as SeqFeatures with primary tag = variation) are avoided, and other
regions (with primary tag = 'Excluded') can also be skipped.  I use
this with Bio::Tools::Run::Mdust to avoid low-complexity regions (must
be run separately), but other programs could also be used.

=head2 EXPORT

None.

=head1 SEE ALSO

L<Bio::Tools::Run::Mdust>, L<Bio::SeqFeature::SiRNA::Pair>,
L<Bio::SeqFeature::SiRNA::Oligo>, L<perl>.

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
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR

Donald Jackson (donald.jackson@bms.com)

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::SiRNA;

require 5.005_62;
use strict;
use warnings;

use vars qw($AUTOLOAD);

use Bio::Seq::RichSeq;
use Bio::SeqFeature::Generic;
use Bio::Root::Root;
use Bio::SeqFeature::SiRNA::Oligo;
use Bio::SeqFeature::SiRNA::Pair;


our @ISA = qw(Bio::Root::Root);


our $VERSION = '1.0';

our %PATTERNS = ( 1 	=> '(AA.{19}TT)',
		  2 	=> '(AA.{19}[ACG][ACG])',
		  3 	=> '([CGT]A.{21})',
		  Pol3	=> '(.A[AG].{17}[CT]..)'
		  );

our %COMP = ( A => 'T',
	      T => 'A',
	      C => 'G',
	      G => 'C',
	      N => 'N',
	      );

our @ARGNAMES = qw(START_PAD END_PAD MIN_GC CUTOFF OLIGOS AVOID_SNPS
		   GSTRING TMPDIR TARGET DEBUG);


=head2 new

 Title		: new
 Usage		: my $sirna_designer = Bio::Tools::SiRNA->new();
 Function	: Constructor for designer object
 Returns	: Bio::Tools::SiRNA object
 Args		: target - the target sequence for the SiRNAs as a Bio::Seq::RichSeq
                  start_pad - distance from the CDS start to skip (default 75)
                  end_pad - distance from the CDS end to skip (default 50)
                  min_gc - minimum GC fraction (NOT percent) (default 0.4)
                  max_gc - maximum GC fraction (NOT percent) (default 0.6)
                  cutoff - worst 'rank' accepted(default 3)
                  avoid_snps - boolean - reject oligos that overlap a variation
                     SeqFeature in the target (default true)
                  gstring - maximum allowed consecutive Gs.
                     Too many can cause problems in synthesis (default 4)
  Note		: All arguments can also be changed/accessed using autoloaded 
                 methods such as:

    my $start_pad = $sirna_designer->start_pad().

=cut

sub new {
    my ($proto, @args) = @_;
    my $pkg = ref($proto) || $proto;

    my %args;
    my $self = {};
    bless ($self, $pkg);

    @args{@ARGNAMES} = $self->_rearrange(\@ARGNAMES, @args); 

    $self->{'start_pad'} = $args{'START_PAD'} || 75; # nt from start to mask
    $self->{'end_pad'} = $args{'END_PAD'} || 50; # nt from end to mask
    $self->{'min_gc'} = $args{'MIN_GC'} || 0.40;
    $self->{'max_gc'} = $args{'MAX_GC'} || 0.60;
    $self->{'cutoff'} = $args{'CUTOFF'} || 3; # highest (worst) rank wanted
    $self->{'oligos'} = [];
    defined($args{'AVOID_SNPS'}) ? $self->{'avoid_snps'} = $args{'AVOID_SNPS'} :  
	$self->{'avoid_snps'} = 1; # (t/f to avoid or include reagents that cover SNPs)
    $self->{'gstring'} = $args{'GSTRING'} || 4; # maximum allowed consecutive Gs - too many can cause problems in oligo synthesis
    $self->{'tmpdir'} = $args{'TMPDIR'}  || $ENV{'TMPDIR'} || $ENV{'TMP'} || '';
    $self->{'debug'} = $args{'DEBUG'} || 0;

    $self->target($args{'TARGET'}) if ($args{'TARGET'});


    return $self;
}


=head2 target 

  Title		: target
  Usage		: my $target_seq = $sirna_designer->target(); # get the current target
                  OR 
                  $sirna_designer->target($new_target_seq); # set a new target 
  Function	: Set/get the target as a Bio::Seq::RichSeq object
  Returns	: a Bio::Seq::RichSeq object
  Args		: a Bio::Seq::RichSeq object (optional)

=cut

sub target {
    my ($self, $target) = @_;

    if ($target) {
	unless ($target->isa('Bio::Seq::RichSeq')) {
	    $self->throw(  -class => 'Bio::Root::BadParameter',
			   -text  => "Target must be passed as a Bio::Seq::RichSeq object" );
	}
	unless ( grep { uc($target->molecule) eq $_ } qw(DNA MRNA CDNA)) {
	    $self->throw(  -class => 'Bio::Root::BadParameter',
			   -text  =>  "Sequences of type ". $target->molecule. " are not supported"
			   );
        }

	$self->{'target'} = $target;
	return 1;

    }
    elsif ($self->{'target'}) {
	return $self->{'target'};
    }
    else {
	$self->throw("Target sequence not defined");
    }
}

=head2 design

  Title		: design
  Usage		: my @pairs = $sirna_designer->design();
  Purpose	: Design SiRNA oligo pairs.  
  Returns	: A list of SiRNA pairs as Bio::SeqFeature::SiRNA::Pair objects
  Args		: none

=cut

sub design {	
    my ($self) = @_;

    unless ( grep { $_->primary_tag eq 'Target' } $self->target->top_SeqFeatures ) {
	$self->_define_target();
    }

    foreach ( 1 .. $self->cutoff ) {
	$self->_get_oligos($_);
    }
       
    return ( grep { $_->isa('Bio::SeqFeature::SiRNA::Pair') } $self->target->top_SeqFeatures );
}
    
sub _define_target {
    my ($self) = @_;
    my ($feat, $cds, $left, $right);

    my $target = $self->target or 
	$self->throw("Unable to design oligos - no target provided");

    ($cds) = grep { $_->primary_tag eq 'CDS' } $target->top_SeqFeatures;
    
    if ($cds) {
	$left = $cds->start + $self->end_pad;
	$right = $cds->end - $self->start_pad;	
    }
    else {
	$left = $target->start + $self->end_pad;
	$right = $target->end - $self->start_pad;
    }
    # define target region 
    my $targregion = Bio::SeqFeature::Generic->new( -start 		=> $left,
						    -end 		=> $right,
						    -primary		=> 'Target' );
    $self->target->add_SeqFeature($targregion);
}

sub _regex {
    my ($self, $rank) = @_;
    return $PATTERNS{$rank};
}

sub _get_oligos {
    # use regular expressions to pull out oligos

    my ($self, $rank) = @_;
    my $regex = $self->_regex($rank);
    my @exclude;


    my ($targregion) = grep { $_->primary_tag eq 'Target' } $self->target->top_SeqFeatures;
    my $seq = $targregion->seq->seq;
    # but this way I loose start info
    my $targstart = $targregion->start;
    
    # exclude masked region
    push(@exclude, grep { $_->primary_tag eq 'Excluded' } $self->target->top_SeqFeatures);

    # add SNP checking
    if ($self->avoid_snps) {
	my @snps =  grep { $_->primary_tag eq 'variation' } $self->target->top_SeqFeatures;
	push(@exclude, @snps);
    }

    while ( $seq =~ /$regex/gi ) {
	my $target = $1;

	# check for too many Gs (or Cs on the other strand)
	next if ( $target =~ /G{ $self->gstring,}/io );
	next if ( $target =~ /C{ $self->gstring,}/io );
	# skip Ns (for filtering)
	next if ( $target =~ /N/i);

	my $start = length($`) + $targstart;
	my $stop = $start + length($target) -1;

	my @gc = ( $target =~ /G|C/gi);
	my $fxGC = sprintf("%2.2f", (scalar(@gc) / length($target)));
	next if ($fxGC < $self->min_gc);
	next if ($fxGC > $self->max_gc);

	my $sense = Bio::SeqFeature::SiRNA::Oligo->new( -start 		=> $start,
							-end 		=> $stop,
							-strand 	=> 1,
							-seq 		=> _get_sense($target),
							-source_tag	=> ref($self),
						       );	

	my $asense = Bio::SeqFeature::SiRNA::Oligo->new( -start 	=> $start,
							 -end		=> $stop,
							 -strand	=> -1,
							 -seq 		=> _get_anti($target), 
							 -source_tag	=> ref($self),
							 );

  	my $sirna = Bio::SeqFeature::SiRNA::Pair->new( -rank 		=> $rank,
						       -fxGC		=> $fxGC,
						       -sense 		=> $sense,
						       -antisense 	=> $asense,     
						       -source_tag	=> ref($self),
						       );

	unless ($self->_has_overlap($sirna, \@exclude)) {
	    $self->target->add_SeqFeature($sirna);
	}
    }
}    

sub _has_overlap {
    # flag any pairs that overlap an UNDESIRED feature (eg SNP)
    # return true if there is overlap, false if not

    my ($self, $test, $flist) = @_;
    print STDERR "Checking oligo at ", $test->start, " to ",$test->end, "\n" 
	if ($self->debug);
    
    foreach my $feat (@$flist) {
	if (($test->start <= $feat->end) and ($test->end >= $feat->start)) {
	    print STDERR "Overlaps ", $feat->primary_tag, " at ",
	    $feat->start, " to ", $feat->end, "\n" if ($self->debug);
	    return 1;
	}
    }
    return 0; # default - no overlap
}
    
	 
sub _get_sense {
    my ($target) = @_;
    # trim off 1st 2 nt to get overhang
    $target =~ s/^..//;
    # convert T's to U's (transcribe)
    $target =~ s/T/U/g;
    # force last 2 nt to be T's
    $target =~ s/..$/TT/;

    return $target;
}

sub _get_anti {
    my ($target) = @_;
    my @target = split(//, $target);
    my ($nt,@antitarget);

    while ($nt = pop @target) {
	push(@antitarget, $COMP{$nt});
    }
    my $anti = join('', @antitarget);
    # trim off 1st 2 nt to get overhang
    $anti =~ s/^..//;
    # convert T's to U's
    $anti =~ s/T/U/g;
    # convert last 2 NT's to T
    $anti =~ s/..$/TT/;

    return $anti;
}


sub AUTOLOAD {
    my ($self, $value) = @_;
    my $name = $AUTOLOAD;
    $name =~ s/.+:://;

    return if ($name eq 'DESTROY');


    if (defined $value) {
	$self->{$name} = $value;
    }

    unless (exists $self->{$name}) {
	$self->throw("Attribute $name not defined for ". ref($self));
  	#return undef;
    }

    return $self->{$name};
}

1;
