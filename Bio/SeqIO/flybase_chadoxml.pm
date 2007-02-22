package Bio::SeqIO::flybase_chadoxml;
use strict;

=pod

This is a simple subclass of L<Bio::SeqIO::chadoxml>; please see
its documentation for details.

=cut

use base 'Bio::SeqIO::chadoxml';

sub _initialize {
    my($self,%args) = @_;
    $self->SUPER::_initialize(%args);

    #default for standard chado is polypeptide
    $Bio::SeqIO::chadoxml::feattype_args2so{'CDS'} = 'protein';
    $Bio::SeqIO::chadoxml::cv_name{'sequence'} = 'SO';
    $Bio::SeqIO::chadoxml::cv_name{'relationship'} = 'relationship type';
    $Bio::SeqIO::chadoxml::cv_name{'feature_property'} = 'property type';

    return;
} 


=head2 return_ftype_hash

=over

=item Usage

  $obj->return_ftype_hash()

=item Function

A simple hash where returning it has be factored out of the main
code to allow subclasses to override it.

=item Returns

A hash that indicates what the name of the SO term is and what
the name of the Sequence Ontology is in the cv table.

=item Arguments

The string that represents the SO term.

=back

=cut

sub return_ftype_hash {
    my $self  = shift;
    my $ftype = shift;
    my %ftype_hash = ( "name" => $ftype,
                       "cv_id" => {"name" => $Bio::SeqIO::chadoxml::cv_name{'sequence'} });
    return %ftype_hash;
}

=head2 return_reltypename

=over

=item Usage

  $obj->return_reltypename()

=item Function

Return the appropriate relationship type name depending on the 
feature type (typically partof, but producedby for proteins).

=item Returns

A relationship type name.

=item Arguments

A SO type name.

=back

=cut

sub return_reltypename {
    my $self   = shift;
    my $sftype = shift;

    my $reltypename;
    if ($sftype eq 'protein' || $sftype eq 'polypeptide') {
        $reltypename = 'producedby';
    } else {
        $reltypename = 'partof';
    }

    return $reltypename;
}


1;
