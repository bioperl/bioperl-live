# $Id$
#
# BioPerl module for Bio::Tools::Analysis::SimpleAnalysisBase
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Richard Adams
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Analysis::SimpleAnalysisBase - abstact superclass for
SimpleAnalysis implementations

=head1 SYNOPSIS

  # no to be run directly

=head1 DESCRIPTION

This class is a generic implementation of SimpleAnalysisI and should
be used as a base class for specific implementations.

SimpleAnalysis predictions just need to provide a specific run()
result() and _init() methods.

=head1 SEE ALSO

L<Bio::SimpleAnalysisI>, 
L<Bio::WebAgent>

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                       - General discussion
  http://bio.perl.org/MailList.html           - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHORS

Richard Adams, Richard.Adams@ed.ac.uk, 
Heikki Lehvaslaiho, heikki@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Analysis::SimpleAnalysisBase;

use vars qw(@ISA);
use strict;
use Data::Dumper;
use Bio::SimpleAnalysisI;
use Bio::WebAgent;

my $FLOAT = '[+-]?\d*\.\d*';

my %STATUS =  map { $_ => 1 } qw(CREATED COMPLETED TERMINATED_BY_ERROR);

@ISA = qw(Bio::WebAgent Bio::SimpleAnalysisI );

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(); #WebAGent new
    $self->_init;      #this line has to be before the attributes are filled in
    while ( @_ ) {
        my $key = lc shift;
        $key =~ s/^-//;
        $self->$key(shift);
    }
    return $self;
}


sub seq {
    my ($self,$value) = @_;
    if ( defined $value) {
        $self->throw("I need a Bio::PrimarySeqI, not  [". $value. "]")
            unless $value->isa('Bio::PrimarySeqI');

        my $mol_type = $self->analysis_spec->{'type'};
        $self->throw("I need a [" . $mol_type . "]  seq, not a  [". $value->alphabet. "]")
            unless $value->alphabet =~/$mol_type/i;
        $self->{'_seq'} = $value;
        return $self;
    }
    return $self->{'_seq'} ;
}



# the next 4 subs are just 'getters'
sub analysis_name {
    my $self = shift;
    return $self->{'_ANALYSIS_NAME'};
}

sub analysis_spec {
    my $self = shift;
    return $self->{'_ANALYSIS_SPEC'};
}

sub input_spec {
    my $self = shift;
    return $self->{'_INPUT_SPEC'};
}

sub result_spec {
    my $self = shift;
    return $self->{'_RESULT_SPEC'};
}

sub run {
    my ($self, $args) = @_;
    $self->_process_arguments ($args) if $args;

    # check input
    $self->throw("Need a sequence object as an input") unless $self->seq;
    $self->debug(Data::Dumper->Dump([$self],[$self]));

    # internal run()
    $self->_run;
    return $self;
}

sub wait_for {
    my ($self, $args) = @_;
    $self->run($args);
}

sub status {
    my ($self,$value) = @_;

    if( defined $value) {
        no strict 'refs';
        my $class = ref($self);
        $self->throw("Not a valid status value [$value]\n".
                     "Valid values are ". join(", ", keys %STATUS ))
            unless defined $STATUS{$value};
        $self->{'_status'} = $value;
        use strict;
    }
    return $self->{'_status'} || 'CREATED' ;
}

sub _process_arguments {
    my ($self, $args) = @_;

    my %spec;
    map {$spec{ $_->{'name'} } = $_ } @{$self->input_spec};

    $self->debug(Data::Dumper->Dump([\%spec, $args],[\%spec, $args]));
    foreach my $key (keys %$args) {
        my $value = $args->{$key};

        $self->throw("Unknown argument [$key]")
            unless $spec{$key};
        $self->$key($value);
    }

    foreach my $key (keys %spec) {
        $self->throw("Mandatory argument [$key] is not set")
            if $spec{$key}{'mandatory'} eq 'true' and not defined $self->$key;
    }
}


sub _run { shift->throw_not_implemented();}
	


