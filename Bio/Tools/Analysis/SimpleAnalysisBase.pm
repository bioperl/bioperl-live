#
# BioPerl module for Bio::Tools::Analysis::SimpleAnalysisBase
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Richard Adams
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Analysis::SimpleAnalysisBase - abstract superclass for
SimpleAnalysis implementations

=head1 SYNOPSIS

# not to be run directly

=head1 DESCRIPTION

This class is a generic implementation of SimpleAnalysisI and should
be used as a base class for specific implementations.

Modules implementing SimpleAnalysisBase only need to provide specific 
_init(), _run() and result() methods, plus any get/set methods for 
parameters to the analysis program.

=head1 SEE ALSO

L<Bio::SimpleAnalysisI>, 
L<Bio::WebAgent>

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHORS

Richard Adams, Richard.Adams@ed.ac.uk, 
Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Analysis::SimpleAnalysisBase;

use strict;
use Data::Dumper;

my $FLOAT = '[+-]?\d*\.\d*';

my %STATUS =  map { $_ => 1 } qw(CREATED COMPLETED TERMINATED_BY_ERROR);

use base qw(Bio::WebAgent Bio::SimpleAnalysisI);

=head2 new

 Usage   : $job->new(...)
 Returns : a new analysis object, 
 Args    : none (but an implementation may choose
           to add arguments representing parameters for the analysis
           program. Each key value of must have a method implemented
           for it in a subclass. A seq () method is provided here as
           this will probably be needed by all sequence analysis programs

=cut

sub new {
    my $class = shift;

    my $self = $class->SUPER::new(); #WebAgent new
    $self->_init; #this line has to be before the attributes are filled in
    while ( @_ ) {
        my $key = lc shift;
        $key =~ s/^-//;
        $self->$key(shift);
    }
    return $self;
}

=head2 seq

 Usage   : $job->seq()
 Returns : a Bio::PrimarySeqI implementing sequence object, or void
 Args    : None, or a Bio::PrimarySeqI implementing object 

=cut

sub seq {
    my ($self,$value) = @_;
    if ( defined $value) {
        $self->throw("I need a Bio::PrimarySeqI, not  [". $value. "]")
            unless $value->isa('Bio::PrimarySeqI');
		$self->throw(" I need a PrimarySeq object, not a BioSeq object ")
			if $value->isa('Bio::SeqI');

        my $mol_type = $self->analysis_spec->{'type'};
        $self->throw("I need a [" . $mol_type . "]  seq, not a  [". $value->alphabet. "]")
            unless $value->alphabet =~/$mol_type/i;
        $self->{'_seq'} = $value;
        return $self;
    }
    return $self->{'_seq'} ;
}

=head2  analysis_name

    Usage     : $analysis->analysis_name();
    Returns   : The analysis name
    Arguments : none

=cut

sub analysis_name {
    my $self = shift;
    return $self->{'_ANALYSIS_NAME'};
}

=head2  analysis_spec

    Usage    :  $analysis->analysis_spec();
    Returns  :  a hash reference to  a hash of analysis parameters. See
                Bio::SimpleAnalysisI for a list of recommended key values.
    Arguments:  none

=cut

sub analysis_spec {
    my $self = shift;
    return $self->{'_ANALYSIS_SPEC'};
}

=head2 clear

    Usage     : $analysis->clear();
    Returns   : true value on success
    Arguments : none
    Purpose   : to remove raw results from a previous analysis so that
                an analysis can be repeated with different parameters.

=cut

sub clear {
	my $self= shift;
	if (defined($self->{'_result'})) {
		delete $self->{'_result'};
		}
	if (defined ($self->{'_parsed'})) {
		delete $self->{'_parsed'};
		}
	return 1;
}
		 


=head2  input_spec

    Usage     : $analysis->input_spec();
    Returns   : a  reference to  an array of  hashes of analysis parameters. See
                Bio::SimpleAnalysisI for a list of recommended key values.
    Arguments : none

=cut

sub input_spec {
    my $self = shift;
    return $self->{'_INPUT_SPEC'};
}

=head2  result_spec

    Usage     : $analysis->result_spec();
    Returns   : a  reference to  a   hashes of resultformats. See
                Bio::SimpleAnalysisI for a list of recommended key values. 
                The key values can be used as parameters to the result() 
                method, the values provide descriptions.
    Arguments : none

=cut

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
	
1;
