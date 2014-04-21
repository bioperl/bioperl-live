#
# BioPerl module for Bio::Tools::SiRNA
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Donald Jackson, donald.jackson@bms.com
#
# Copyright Bristol-Myers Squibb
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

SiRNA - Perl object for designing small inhibitory RNAs.

=head1 SYNOPSIS

  use Bio::Tools::SiRNA;

  my $sirna_designer = Bio::Tools::SiRNA->new( -target => $bio_seq,
                                               -rules  => 'saigo'
    );
  my @pairs = $sirna_designer->design;

  foreach $pair (@pairs) {
      my $sense_oligo_sequence = $pair->sense->seq;
      my $antisense_oligo_sequence = $pair->antisense->seq;

      # print out results
      print join ("\t", $pair->start, $pair->end, $pair->rank,
                  $sense_oligo_sequence, $antisense_oligo_sequence), "\n";
  }

=head1 DESCRIPTION

Package for designing siRNA reagents.

Input is a L<Bio::SeqI>-compliant object (the target).

Output is a list of Bio::SeqFeature::SiRNA::Pair objects, which are
added to the feature table of the target sequence.  Each
Bio::SeqFeature::SiRNA::Pair contains two subfeatures
(Bio::SeqFeature::Oligo objects) which correspond to the individual
oligos.  These objects provide accessors for the information on the
individual reagent pairs.

This verion of Bio::Tools::SiRNA represents a major change in architecture.
Specific 'rulesets' for siRNA selection as developed by various groups are
implemented as Bio::Tools::SiRNA::Ruleset objects, which inherit from
Bio::Tools::SiRNA.  This will make it easier to add new rule sets or modify
existing approaches. Currently the Tuschl and Ui-Tei (2004) rules are 
implemented. For consistency, the Tuschl rules are implemented by default.

In addition, this module provides three 'extra' rules which can be added
above and beyond any ruleset.

=over 3

=item 1.

SiRNAs that overlap known SNPs (identified as SeqFeatures with 
primary tag = variation) can be avoided.

=item 2.

Other regions (with primary tag = 'Excluded') can also be skipped.  I
use this with Bio::Tools::Run::Mdust to avoid low-complexity regions
(must be run separately), but other programs could also be used.

=item 3.

SiRNAs may also be selected in the 3 prime UTR of a gene by setting
$sirna_designer-E<gt>include_3pr() to true.

=back

=head2 EXPORT

None.

=head1 SEE ALSO

L<Bio::Tools::Run::Mdust>, L<Bio::SeqFeature::SiRNA::Pair>,
L<Bio::SeqFeature::SiRNA::Oligo>..

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
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Donald Jackson (donald.jackson@bms.com)

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::SiRNA;

use strict;
use warnings;

use vars qw($AUTOLOAD);

use Bio::Seq::RichSeq;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::SiRNA::Oligo;
use Bio::SeqFeature::SiRNA::Pair;


use base qw(Bio::Root::Root);


our %COMP = ( A => 'T',
	      T => 'A',
	      C => 'G',
	      G => 'C',
	      N => 'N',
	      );

our @ARGNAMES = qw(RULES START_PAD END_PAD MIN_GC CUTOFF OLIGOS AVOID_SNPS
		   GSTRING TMPDIR TARGET DEBUG);


=head2 new

 Title		: new
 Usage		: my $sirna_designer = Bio::Tools::SiRNA->new();
 Function	: Constructor for designer object
 Returns	: Bio::Tools::SiRNA object
 Args		: target - the target sequence for the SiRNAs as a Bio::Seq::RichSeq
                  start_pad - distance from the CDS start to skip (default 75)
                  end_pad - distance from the CDS end to skip (default 50)
                  include_3pr - set to true to include SiRNAs in the 3prime UTR (default false)
                  rules - rules for selecting siRNAs, currently supporting saigo and tuschl
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

    my $self = {};
    bless ($self, $pkg);

    my %args;

    @args{@ARGNAMES} = $self->_rearrange(\@ARGNAMES, @args); 
    
    if ($args{'RULES'}) {
	$self->rules($args{'RULES'});
    }

    $self->{'start_pad'} = $args{'START_PAD'} || 75; # nt from start to mask
    $self->{'end_pad'} = $args{'END_PAD'} || 50; # nt from end to mask
    $self->{'include_3pr'} = $args{'INCLUDE_3PR'} || 0; # look for oligos in 3prime UTR
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
  Function	: Set/get the target as a Bio::SeqI-compliant object
  Returns	: a Bio::SeqI-compliant object
  Args		: a Bio::SeqI-compliant object (optional)

=cut

sub target {
    my ($self, $target) = @_;

    if ($target) {
	unless ($target->isa('Bio::SeqI')) {
	    $self->throw(  -class => 'Bio::Root::BadParameter',
			   -text  => "Target must be passed as a Bio::Seq object" );
	}
	if ($target->can('molecule')) {
	    ( grep { uc($target->molecule) eq $_ } qw(DNA MRNA CDNA)) or
		$self->throw(  -class => 'Bio::Root::BadParameter',
			       -text  =>  "Sequences of type ". $target->molecule. " are not supported"
			       );
	}
	else {
	    ($target->alphabet eq 'dna') or 
		$self->throw(  -class => 'Bio::Root::BadParameter',
			       -text  =>  "Sequences of alphabet ". $target->alphabet. " are not supported"
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

=head2 rules

    Title	: rules
    Usage	: $sirna->rules('ruleset')
    Purpose	: set/get ruleset to use for selecting SiRNA oligo pairs.
    Returns	: not sure yet
    Args	: a ruleset name (currently supported: Tuschl, Saigo)
                  or a Bio::Tools::SiRNA::RulesetI compliant object

=cut

sub rules {
    my ($self, $rules) = @_;

    if ($rules) {
	$self->_load_ruleset($rules);
    }
    # default: use tuschl rules
    unless ($self->{_rules}) {
	$self->_load_ruleset('tuschl');
    }
    return $self->{_rules};
}

sub _load_ruleset {
    my ($self, $ruleset) = @_;

    my $rule_module = join('::', ref($self), 'Ruleset', lc($ruleset));

    eval "require $rule_module";
    
    if ($@) {
	#warn join("\n", '@INC contains:', @INC, undef);
	$self->throw("Unable to load $rule_module: $@");
	return;
    }

    else {
	$self->{_rules} = $rule_module;
	bless($self, $rule_module); # recast as subclass
    }
	
    return 1;
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

    ($self->rules) or $self->throw('Unable to design siRNAs: no rule set specified');

#     unless ( grep { $_->primary_tag eq 'Target' } $self->target->top_SeqFeatures ) {
# 	$self->_define_target();
#     }

    my @oligos = $self->_get_oligos();
       
    return ( grep { $_->isa('Bio::SeqFeature::SiRNA::Pair') } $self->target->top_SeqFeatures );
}
    
sub _define_target {
    my ($self) = @_;
    my ($feat, $cds, $left, $right);

    my $target = $self->target or 
	$self->throw("Unable to design oligos - no target provided");

    ($cds) = grep { $_->primary_tag eq 'CDS' } $target->top_SeqFeatures if ($target->can('top_SeqFeatures'));
    
    if ($cds) {
	$left = $cds->start + $self->start_pad;
	if (!$self->include_3pr) {
	    $right = $cds->end - $self->end_pad;
	}
	else {
	    $right = $target->length - $self->end_pad;
	}
    }
    else {
	$left = 0 + $self->start_pad;
	$right = $target->length - $self->end_pad;
    }


    # is there anything left?
    if (($right - $left) < 20) {
	$self->throw("There isn't enough sequence to design oligos.  Please reduce start_pad and end_pad or supply more sequence");
    }
    # define target region 
    my $targregion = Bio::SeqFeature::Generic->new( -start 		=> $left,
						    -end 		=> $right,
						    -primary		=> 'Target' );
    $self->target->add_SeqFeature($targregion);

    # locate excluded regions
    my @excluded = grep { $_->primary_tag eq 'Excluded' } $self->target->top_SeqFeatures;

    if ($self->avoid_snps) {
	my @snps =  grep { $_->primary_tag eq 'variation' } $self->target->top_SeqFeatures;
	push(@excluded, @snps);
    }
    
    $self->excluded(\@excluded);

    return $targregion;
}

sub _get_targetregion {
    my ($self) = @_;
    
    my ($targregion) = grep { $_->primary_tag eq 'Target' } $self->target->top_SeqFeatures;
    $targregion ||= $self->_define_target;

    $self->throw("Target region for SiRNA design not defined") unless ($targregion);

    my $seq = $targregion->seq->seq;
    # but this way I loose start info
     my $targstart = $targregion->start;

    return ($seq, $targstart);
}   

# MOVE to SiRNA::Ruleset::tuschl
# sub _regex {
#     my ($self, $rank) = @_;
#     return $PATTERNS{$rank};
# }

# sub _get_oligos {
#     # use regular expressions to pull out oligos

#     my ($self, $rank) = @_;
#     my $regex = $self->_regex($rank);
#     my @exclude;


#     my ($targregion) = grep { $_->primary_tag eq 'Target' } $self->target->top_SeqFeatures;
#     my $seq = $targregion->seq->seq;
#     # but this way I loose start info
#     my $targstart = $targregion->start;
    
#     # exclude masked region
#     push(@exclude, grep { $_->primary_tag eq 'Excluded' } $self->target->top_SeqFeatures);

#     # add SNP checking
#     if ($self->avoid_snps) {
# 	my @snps =  grep { $_->primary_tag eq 'variation' } $self->target->top_SeqFeatures;
# 	push(@exclude, @snps);
#     }

#     while ( $seq =~ /$regex/gi ) {
# 	my $target = $1;

# 	# check for too many Gs (or Cs on the other strand)
# 	next if ( $target =~ /G{ $self->gstring,}/io );
# 	next if ( $target =~ /C{ $self->gstring,}/io );
# 	# skip Ns (for filtering)
# 	next if ( $target =~ /N/i);

# 	my $start = length($`) + $targstart;
# 	my $stop = $start + length($target) -1;

# 	my @gc = ( $target =~ /G|C/gi);
# 	my $fxGC = sprintf("%2.2f", (scalar(@gc) / length($target)));
# 	next if ($fxGC < $self->min_gc);
# 	next if ($fxGC > $self->max_gc);

# 	my $sense = Bio::SeqFeature::SiRNA::Oligo->new( -start 		=> $start,
# 							-end 		=> $stop,
# 							-strand 	=> 1,
# 							-seq 		=> _get_sense($target),
# 							-source_tag	=> ref($self),
# 						       );	

# 	my $asense = Bio::SeqFeature::SiRNA::Oligo->new( -start 	=> $start,
# 							 -end		=> $stop,
# 							 -strand	=> -1,
# 							 -seq 		=> _get_anti($target), 
# 							 -source_tag	=> ref($self),
# 							 );

#   	my $sirna = Bio::SeqFeature::SiRNA::Pair->new( -rank 		=> $rank,
# 						       -fxGC		=> $fxGC,
# 						       -sense 		=> $sense,
# 						       -antisense 	=> $asense,     
# 						       -source_tag	=> ref($self),
# 						       );

# 	unless ($self->_has_overlap($sirna, \@exclude)) {
# 	    $self->target->add_SeqFeature($sirna);
# 	}
#     }
# }    

=head2 add_oligos

  Title		: add_oligos
  Usage	 	: $sirna_designer->add_oligos($sequence, $start, $rank);
  Purpose	: Add SiRNA olgos to target Bio::Seq as Bio::SeqFeature::SiRNA::Pair objects
  Args		: Oligo sequence and start position (required), rank/score (optional)

=cut

sub add_oligos {
    my ($self, $seq, $start, $rank) = @_;

    ($seq) or throw ('No sequence supplied for add_oligos');
    (defined $start) or throw ('No start position specified for  add_oligos');
    
    my ($end) = $start + length($seq);

    my ($sseq) = $self->_get_sense($seq);
    my $sense = Bio::SeqFeature::SiRNA::Oligo->new( -start 		=> $start,
						    -end 		=> ($start + length($sseq)),
						    -strand 	=> 1,
						    -seq 		=> $sseq,
						    -source_tag	=> ref($self),
						    );	

    my $aseq = $self->_get_anti($seq);
    my $asense = Bio::SeqFeature::SiRNA::Oligo->new( -start 		=> $end,
						     -end		=> ($end - length($aseq)),
						     -strand		=> -1,
						     -seq 		=> $aseq, 
						     -source_tag	=> ref($self),
						     );

    my $sirna = Bio::SeqFeature::SiRNA::Pair->new( -rank 		=> $rank,
						  # -fxGC		=> $fxGC,
						   -sense 		=> $sense,
						   -antisense 	=> $asense,     
						   -source_tag	=> ref($self),
						   );

    unless ($self->_has_overlap($sirna, $self->excluded)) {
	$self->target->add_SeqFeature($sirna);
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
    
# MOVE to SiRNA::Ruleset::tuschl
	 
# sub _get_sense {
#     my ($target) = @_;
#     # trim off 1st 2 nt to get overhang
#     $target =~ s/^..//;
#     # convert T's to U's (transcribe)
#     $target =~ s/T/U/gi;
#     # force last 2 nt to be T's
#     $target =~ s/..$/TT/;

#     return $target;
# }

# sub _get_anti {
#     my ($target) = @_;
#     my @target = split(//, $target);
#     my ($nt,@antitarget);

#     while ($nt = pop @target) {
# 	push(@antitarget, $COMP{$nt});
#     }
#     my $anti = join('', @antitarget);
#     # trim off 1st 2 nt to get overhang
#     $anti =~ s/^..//;
#     # convert T's to U's
#     $anti =~ s/T/U/gi;
#     # convert last 2 NT's to T
#     $anti =~ s/..$/TT/;

#     return $anti;
# }


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
    }

    return $self->{$name};
}

sub _comp {
    my ($self, $char) = @_;

    return unless ($char);
    $char = uc($char);
    return $COMP{ $char };
}
1;
