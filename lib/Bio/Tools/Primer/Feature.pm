#
# BioPerl module for Bio::Tools::Primer::Feature
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Primer::Feature - position of a single primer

=head1 SYNOPSIS

    use Bio::Tools::Primer::Feature;

    my $pf = Bio::Tools::Primer::Feature->new( -start => $start, -end => $end, -strand => $strand);
    $pf->attach_seq($seq);

    # is a SeqFeatureI

    print "primer starts at ",$pf->start," with sequence ",$pf->seq->seq(),"\n";

    # helper functions

    print "GC percentage ",$pf->gc(),"\n";
    print "has inversion of size 4 at ",$pf->inversion(4),"\n";



=head1 DESCRIPTION

Primer Features represents one primer in a primer pair. This object is
mainly for designing primers, and probably principly used in the
primer design system

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ewan Birney

Email birney-at-ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...



package Bio::Tools::Primer::Feature;

use base qw(Bio::SeqFeature::Generic);



sub new {
    my ( $caller, @args) = @_;   
    my ($self) = $caller->SUPER::new(@args); 

    # done - we hope
    return $self;
}

sub gc_percent {
    my $self = shift;

    my $seq = $self->seq();

    if( !defined $seq ) {
	$self->throw("Primer feature has no attached sequence, can't calculate GC");
    }

    my $str = $seq->seq();

    my $count = $str =~ tr/GCgc/GCgc/;

    return $count*100.0 / $seq->length;
}

sub inversion {
    my $self = shift;
    my $size = shift;

    if( !defined $size ) {
	$self->throw("Must have size paramter in inversion");
    }

    my $seq = $self->seq();

    if( !defined $seq ) {
	$self->throw("Primer feature has no attached sequence, can't calculate inversion");
    }

    my $len = $seq->length - $size;

    my $str = $seq->seq();

    foreach my $i ( 0 .. $len ) {
	my $revstr = substr($str,$i,$size);
	my $orig = $revstr;
	$revstr = reverse $revstr;
	$revstr = s/[^ATGCNatgcn]/N/g;

	$revstr =~ tr/ATGCNatgcn/TACGNtacgn/;

	if( $str =~ /$revstr/ ) {
	    return $orig;
	}
    }

    return;
}

1;
