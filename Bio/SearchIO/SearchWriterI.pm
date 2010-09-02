#-----------------------------------------------------------------
#
# BioPerl module Bio::SearchIO::SearchWriterI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

=head1 NAME

Bio::SearchIO::SearchWriterI - Interface for outputting parsed Search results

=head1 SYNOPSIS

Bio::SearchIO::SearchWriterI objects cannot be instantiated since this
module defines a pure interface.

Given an object that implements the Bio::SearchIO::SearchWriterI interface,
you can do the following things with it:

    print $writer->to_string( $result_obj, @args );

=head1 DESCRIPTION

This module defines abstract methods that all subclasses must implement
to be used for outputting results from L<Bio::Search::Result::ResultI>
objects.

=head1 AUTHOR

Steve Chervitz E<lt>sac-at-bioperl.orgE<gt>

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

package Bio::SearchIO::SearchWriterI;


use base qw(Bio::Root::RootI);

=head2 to_string

 Purpose   : Produces data for each Search::Result::ResultI in a string.
           : This is an abstract method. For some useful implementations,
           : see ResultTableWriter.pm, HitTableWriter.pm, 
           : and HSPTableWriter.pm.
 Usage     : print $writer->to_string( $result_obj, @args );
 Argument  : $result_obj = A Bio::Search::Result::ResultI object
           : @args = any additional arguments used by your implementation.
 Returns   : String containing data for each search Result or any of its
           : sub-objects (Hits and HSPs).
 Throws    : n/a

=cut

sub to_string {
    my ($self, $result, @args) = @_;
    $self->throw_not_implemented;
}

=head2 start_report

 Title   : start_report
 Usage   : $self->start_report()
 Function: The method to call when starting a report. You can override it
           to make a custom header
 Returns : string
 Args    : none

=cut

sub start_report { return '' }

=head2 end_report

 Title   : end_report
 Usage   : $self->end_report()
 Function: The method to call when ending a report, this is
           mostly for cleanup for formats which require you to 
           have something at the end of the document (</BODY></HTML>)
           for HTML
 Returns : string
 Args    : none


=cut

sub end_report {  return '' }

=head2 filter

 Title   : filter
 Usage   : $writer->filter('hsp', \&hsp_filter);
 Function: Filter out either at HSP,Hit,or Result level
 Returns : none
 Args    : string => data type,
           CODE reference


=cut

# yes this is an implementation in the interface, 
# yes it assumes that the underlying class is hash-based
# yes that might not be a good idea, but until people
# start extending the SearchWriterI interface I think
# this is an okay way to go

sub filter {
    my ($self,$method,$code) = @_;    
    return unless $method;
    $method = uc($method);
    if( $method ne 'HSP' &&
	$method ne 'HIT' &&
	$method ne 'RESULT' ) {
	$self->warn("Unknown method $method");
	return;
    }
    if( $code )  {
	$self->throw("Must provide a valid code reference") unless ref($code) =~ /CODE/;
	$self->{$method} = $code;
    }
    return $self->{$method};
}

1;


