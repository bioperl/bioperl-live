#!/usr/bin/env perl
# Copyright (c) 2016 by Weigang Qiu Lab

package Bio::BPWrapper;

our $VERSION = '1.13';
use strict; use warnings;
use 5.010;

use constant PROGRAM => 'Bio::BPWrapper';

=pod

=head2 Common Subroutines

=head3 print_version $program_name

Show program name and version and exit

=cut


sub print_version($)
{
    my $program = shift;
    say "${program}, version $Bio::BPWrapper::VERSION";
    exit;
}

use Pod::Usage;

=head3 common_opts $opts

Handle some common options:

=over 4

=item C<--help>, C<-h>

Show usage help

=item C<--man>

Show a manual page via L<Pod::Usage>

=back

=cut

sub common_opts($)
{
    my $opts = shift;
    pod2usage(1) if $opts->{"help"};
    pod2usage(-exitstatus => 0, -verbose => 2) if $opts->{"man"};
}


unless (caller) {
    print "Pssst... this is a module. Invoke via bioaln, bioseq, biopop, or biotree.\n";
}
1;
