#
# BioPerl module for Bio::Tools::EMBOSS::Palindrome
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::EMBOSS::Palindrome - parse EMBOSS palindrome output

=head1 SYNOPSIS

  # a simple script to turn palindrome output into GFF3
  use Bio::Tools::EMBOSS::Palindrome;
  use Bio::Tools::GFF;

  my $parser = Bio::Tools::EMBOSS::Palindrome->new(-file => $filename);
  my $out    = Bio::Tools::GFF->new(-gff_version => 3,
                                   -file => ">$filename.gff");
  while( my $seq = $parser->next_seq ) {
     for my $feat ( $seq->get_SeqFeatures ) {
        $out->write_feature($feat);
     }
  }

=head1 DESCRIPTION

This is a parser for the EMBOSS tool 'palindrome'.  It will produce a
L<Bio::Seq> object for each sequence analyzed.  The sequence will be
empty (but will be of the correct length) and will have attached to it
L<Bio::SeqFeature::FeaturePair> objects which wil


=head2 FUTURE WORK

It may be consolidated into another framework at a later time, but for
the time being it will stay a separate modules.

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
email or the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::EMBOSS::Palindrome;
use vars qw($DEFAULT_SOURCETAG);
use strict;

use Bio::SeqFeature::FeaturePair;
use Bio::SeqFeature::Generic;

use base qw(Bio::Root::IO);
$DEFAULT_SOURCETAG = 'palindrome';

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::EMBOSS::Palindrome->new();
 Function: Builds a new Bio::Tools::EMBOSS::Palindrome object 
 Returns : an instance of Bio::Tools::EMBOSS::Palindrome
 Args    : -file/-fh  => a filename or filehandle for
                         initializing the parser

=cut

=head2 next_seq

 Title   : next_seq
 Usage   : my $seq = $parser->next_seq;
 Function: Get the next feature set from the 
 Returns : L<Bio::SeqI> object
 Args    : none


=cut

sub next_seq {
    my ($self) = @_;
    my (%searching, $seq,$state);
    my $source = $self->source_tag;
    $state = 0;
    while(defined($_ = $self->_readline)) {
	if( /^\s+$/ ) {
	    next;
	} elsif( /^Palindromes\s+of\s*:\s+(\S+)/o ) {
	    $state = 0;
	    if( $seq )  {
		$self->_pushback($_);
		return $seq;
	    } 
	    $seq = Bio::Seq->new(-display_id => $1);
	    # now get ready to store for the next record
	    $searching{'-seq_id'} = $1;
	} elsif( /^Sequence\s+length\s+is\s*:\s+(\d+)/o ) {
	    $seq->length($1);
	    $searching{'-tag'}->{'seqlength'} = $1;
	} elsif( /^(Start|End)\s+at\s+position\s*:\s+(\d+)/ ) {
	    $searching{'-tag'}->{lc($1)} = $2;
	} elsif( m/^(Maximum|Minimum)\s+length\s+of\s+Palindromes\s+
		 is\s*:\s+(\d+)/ox) {
	    $searching{'-tag'}->{lc($1).'_length'} = $2;
	} elsif( /^(Maximum\s+gap)\s+between\s+elements\s+is\s*:\s+(\d+)/o ) {
	    $searching{'-tag'}->{lc($1)} = $2;
	} elsif( m/^Number\s+of\s+mismatches\s+allowed\s+
		 in\s+Palindrome\s*:\s+(\d+)/ox ) {
	    $searching{'-tag'}->{'allowed_mismatches'} = $1;
	} elsif( /^Palindromes:/o ) {
	    $state = 1;
	} elsif( $state == 1 ) {
	    my $feature = Bio::SeqFeature::FeaturePair->new
		(-primary_tag  => 'similarity',
		 -source_tag   => $source);
	    for(my $i = 0; $i < 3; $i++ ) {
		if ($i != 1) {
		    if( /^(\d+)\s+(\S+)\s+(\d+)/o ) {
			my ($start,$match,$end) = ($1,$2,$3);
			my $type = $i == 0 ? 'feature1' : 'feature2';
			($start,$end) = sort { $a <=> $b } ($start,$end);
			$feature->$type(
					Bio::SeqFeature::Generic->new
					(%searching,
					 -start       => $start,
					 -end         => $end,
					 -strand      => $i == 0 ? 1 : -1,
					 -primary_tag => 'similarity',
					 -source_tag  => $source)
					);
		    } else { 
			chomp;
			warn("Out of sync, line did not match:'$_'\n");
		    }

		}
		$_ = $self->_readline;
	    }
	    $seq->add_SeqFeature($feature);
	}
    }
    return $seq;
}

=head2 source_tag

 Title   : source_tag
 Usage   : $obj->source_tag($newval)
 Function: Get/Set Source Tag ('palindrome') by default
 Returns : value of source_tag (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub source_tag{
    my $self = shift;

    return $self->{'source_tag'} = shift if @_;
    return $self->{'source_tag'} || $DEFAULT_SOURCETAG;
}

1;
