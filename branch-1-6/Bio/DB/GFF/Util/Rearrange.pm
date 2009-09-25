=head1 NAME

Bio::DB::GFF::Util::Rearrange - rearrange utility

=head1 SYNOPSIS

 use Bio::DB::GFF::Util::Rearrange 'rearrange';

 my ($arg1,$arg2,$arg3,$others) = rearrange(['ARG1','ARG2','ARG3'],@args);

=head1 DESCRIPTION

This is a different version of the _rearrange() method from
Bio::Root::Root.  It runs as a function call, rather than as a method
call, and it handles unidentified parameters slightly differently.

It exports a single function call:

=over 4

=item @rearranged_args = rearrange(\@parameter_names,@parameters);

The first argument is an array reference containing list of parameter
names in the desired order.  The second and subsequent arguments are a
list of parameters in the format:

  (-arg1=>$arg1,-arg2=>$arg2,-arg3=>$arg3...)

The function calls returns the parameter values in the order in which
they were specified in @parameter_names.  Any parameters that were not
found in @parameter_names are returned in the form of a hash reference
in which the keys are the uppercased forms of the parameter names, and
the values are the parameter values.

=back

=head1 BUGS

None known yet.

=head1 SEE ALSO

L<Bio::DB::GFF>,

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

package Bio::DB::GFF::Util::Rearrange;

use strict;
require Exporter;
use vars qw(@EXPORT @EXPORT_OK);
use base qw(Exporter);
@EXPORT_OK = qw(rearrange);
@EXPORT = qw(rearrange);
use Bio::Root::Version;

# default export
sub rearrange {
    my($order,@param) = @_;
    return unless @param;
    my %param;

    if (ref $param[0] eq 'HASH') {
      %param = %{$param[0]};
    } else {
      return @param unless (defined($param[0]) && substr($param[0],0,1) eq '-');

      my $i;
      for ($i=0;$i<@param;$i+=2) {
        $param[$i]=~s/^\-//;     # get rid of initial - if present
        $param[$i]=~tr/a-z/A-Z/; # parameters are upper case
      }

      %param = @param;                # convert into associative array
    }

    my(@return_array);

    local($^W) = 0;
    my($key)='';
    foreach $key (@$order) {
        my($value);
        if (ref($key) eq 'ARRAY') {
            foreach (@$key) {
                last if defined($value);
                $value = $param{$_};
                delete $param{$_};
            }
        } else {
            $value = $param{$key};
            delete $param{$key};
        }
        push(@return_array,$value);
    }
    push (@return_array,\%param) if %param;
    return @return_array;
}

1;
