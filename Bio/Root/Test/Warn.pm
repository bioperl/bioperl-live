#
# BioPerl module for Bio::Root::Test::Warn
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Root::Test::Warn - Perl extension to test Bioperl methods for warnings

=head1 SYNOPSIS

  use Bio::Root::Test::Warn;

  warning_is {$bio_object->method()} 'Must supply a parameter', "a missing parameter test";
  warning_like {$bio_object->method()} qr/Must supply a parameter/i, "a missing parameter test";

=head1 DESCRIPTION

This module provides a few convenience methods for testing warning based code.

See Test::Warn for details.

You will normally not use this module directly, but have it auto-loaded for you
by Bio::Root::Test.

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

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Root::Test::Warn;

use strict;
use warnings;
use Exporter qw(import);

use Test::Builder;
use Test::Warn;

our @EXPORT = qw(warning_is
                 warnings_are
                 warning_like
                 warnings_like);

{
    my $Tester = Test::Builder->new;
    
    no warnings 'redefine';
    sub Test::Warn::_canonical_got_warning {
        my ($called_from, $msg) = @_;
        my $warn_kind = $called_from eq 'Carp' ? 'carped' : ($called_from =~ /Bio::/ ? 'Bioperl' : 'warn');
        
        my $warning;
        if ($warn_kind eq 'Bioperl') {
            ($warning) = $msg =~ /\n--------------------- WARNING ---------------------\nMSG: (.+)\n---------------------------------------------------\n$/m;
            $warning ||= $msg; # shouldn't ever happen
        }
        else {
            my @warning_stack = split /\n/, $msg;   # some stuff of uplevel is included
            $warning = $warning_stack[0];
        }
        
        return {$warn_kind => $warning}; # return only the real message
    }
    
    sub Test::Warn::_diag_found_warning {
        foreach (@_) {
            if (ref($_) eq 'HASH') {
                ${$_}{carped} ? $Tester->diag("found carped warning: ${$_}{carped}")
                              : (${$_}{Bioperl} ? $Tester->diag("found Bioperl warning: ${$_}{Bioperl}")
                                 : $Tester->diag("found warning: ${$_}{warn}"));
            } else {
                $Tester->diag( "found warning: $_" );
            }
        }
        $Tester->diag( "didn't find a warning" ) unless @_;
    }
    
    sub Test::Warn::_cmp_got_to_exp_warning {
        my ($got_kind, $got_msg) = %{ shift() };
        my ($exp_kind, $exp_msg) = %{ shift() };
        return 0 if ($got_kind eq 'warn') && ($exp_kind eq 'carped');
        
        my $cmp;
        if ($got_kind eq 'Bioperl') {
            $cmp = $got_msg =~ /^\Q$exp_msg\E$/;
        }
        else {
            $cmp = $got_msg =~ /^\Q$exp_msg\E at \S+ line \d+\.?$/;
        }
        
        return $cmp;
    }
}
 
1;
