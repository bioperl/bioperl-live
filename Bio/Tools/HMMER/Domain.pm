
#
# Perl Module for HMMUnit
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
#Copyright Genome Research Limited (1997). Please see information on licensing in LICENSE

package Bio::Tools::HMMER::Domain;

use vars qw($AUTOLOAD @ISA @EXPORT_OK);
use Exporter;
use Carp;
use strict;

#
# Place functions/variables you want to *export*, ie be visible from the caller package into @EXPORT_OK
#

@EXPORT_OK = qw();

#
# @ISA has our inheritance.
#

@ISA = ( 'Exporter' );



my %fields = (
    #Insert field names here as field => undef,
	      seqname => undef,
	      seq_range => undef,
	      hmmname => undef,
	      hmmacc  => undef,
	      hmm_range => undef,
	      bits => undef,
	      evalue => undef,
	      prob => undef,
	      seqbits => undef,
	      alignlines => undef, #raw HMMer2 alignment alignment lines for printing out
);


sub new {
    my $ref = shift;
    my $class = ref($ref) || $ref;
    my $self = {
	'_permitted' => \%fields,
	%fields, };

    $self->{'seq_range'} = new Range;
    $self->{'hmm_range'} = new Range;
    $self->{'alignlines'} = [];
    bless $self, $class;
    return $self;
}


sub AUTOLOAD {
    my $self = shift;
    my $type = ref($self) || carp "$self is not an object - can't therefore find a member!";
    my $name = $AUTOLOAD;
    $name =~ /::DESTROY/ && return;
    $name =~ s/.*://;
    unless (exists $self->{'_permitted'}->{$name} ) {
	carp "In type $type, can't access $name - probably passed a wrong variable into HMMUnit";
    }
    if (@_) {
	return $self->{$name} = shift;
    } else {
	return $self->{$name};
    }
}

sub add_alignment_line {
    my $self = shift;
    my $line = shift;
    push(@{$self->{'alignlines'}},$line);
}

sub each_alignment_line {
    my $self = shift;
    return @{$self->{'alignlines'}};
}

sub get_nse {
    my $self = shift;
    my $sep1 = shift;
    my $sep2 = shift;

    if( !defined $sep2 ) {
	$sep2 = "-";
    }
    if( !defined $sep1 ) {
	$sep1 = "/";
    }

    return sprintf("%s%s%d%s%d",$self->seqname,$sep1,$self->start_seq,$sep2,$self->end_seq);
}


sub start_seq {
    my $self = shift;
    my $start = shift;

    if( !defined $start ) {
	return $self->{'seq_range'}->start();
    }
    $self->{'seq_range'}->start($start);
    return $start;
}

sub end_seq {
    my $self = shift;
    my $end = shift;

    if( !defined $end ) {
	return $self->{'seq_range'}->end();
    }
    $self->{'seq_range'}->end($end);
    return $end;

}


sub start_hmm {
    my $self = shift;
    my $start = shift;

    if( !defined $start ) {
	return $self->{'hmm_range'}->start();
    }
    $self->{'hmm_range'}->start($start);
    return $start;
}

sub end_hmm {
    my $self = shift;
    my $end = shift;

    if( !defined $end ) {
	return $self->{'hmm_range'}->end();
    }
    $self->{'hmm_range'}->end($end);
    return $end;

}



1;  # says use was ok
__END__

=head1 NAME

HMMUnit

=head1 DESCRIPTION

Description for B<HMMUnit>

=head1 AUTHOR

B<Ewan Birney> Email birney@sanger.ac.uk

=over

=item get_nse

No current documentation

=item start_seq

No current documentation

=item end_seq

No current documentation

=item add_alignment_line

No current documentation

=item each_alignment_line

No current documentation
