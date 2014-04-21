#
# BioPerl module for Bio::Restriction::IO::prototype
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chris Fields 
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Restriction::IO::prototype - prototype enzyme set

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::Restriction::IO class.

=head1 DESCRIPTION

This is a parser for the proto/neo file REBASE format, which contains
prototype information as well as (in the neo file) neoschizomer data.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists Your participation is much appreciated.

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

=head1 AUTHOR

Rob Edwards, redwards@utmem.edu

=head1 CONTRIBUTORS

Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Restriction::IO::prototype;

use vars qw(%WITH_REFM_FIELD);
use strict;

#use Bio::Restriction::IO;
use Bio::Restriction::Enzyme;
use Bio::Restriction::EnzymeCollection;

use Data::Dumper;

use base qw(Bio::Restriction::IO::base);

=head2 read

 Title   : read
 Usage   : $renzs = $stream->read
 Function: reads all the restrction enzymes from the stream
 Returns : a Bio::Restriction::Restriction object
 Args    : none

=cut

sub read {
    my $self = shift;
    my $coll = Bio::Restriction::EnzymeCollection->new(-empty => 1);
    my ($seentop, $last_type);
    while (defined (my $line = $self->_readline)) {
        chomp $line;
        next unless $line;
        if ($line =~ /TYPE\s+(I)+/) {
            $last_type = $1;
            $seentop ||= 1;
            next;
        }
        next unless $seentop;
        my @data = split /\s+/,$line,2;
        next if $data[0] =~ /^[-\s]*$/;
        # neo
        my ($enzyme, $is_neo, $is_proto, $site);
        if ($data[0] =~ /^\s+(\S+)\s+(\S+)/) {
            ($enzyme, $site, $is_proto, $is_neo) = ($1, $2, 0, 1);
        } else {
            ($enzyme, $site, $is_proto, $is_neo) = ($data[0], $data[1], 1, 0);
        }
        $site =~ s/\s+//g;
        
        my $precut;
        if ($site =~ m/^\((\d+\/\d+)\)[RYATGCN]+/) {
            $precut=$1;
            $site =~ s/\($precut\)//;
        }
        
        my ($cut, $comp_cut);
        ($site, $cut, $comp_cut) = $self->_cuts_from_site($site);
        
        my $re = Bio::Restriction::Enzyme->new(
            -type => $last_type,
            -site => $site,
            -name => $enzyme,
            -is_prototype => $is_proto,
            -is_neoschizomer => $is_neo);
        
        if ($cut) {
            $re->cut($self->_coordinate_shift_to_cut(length($site), $cut));
            $re->complementary_cut($self->_coordinate_shift_to_cut(length($site), $comp_cut));
        }
        $coll->enzymes($re);
    }
    return $coll->enzymes;
}

=head2 write

 Title   : write
 Usage   : $stream->write($renzs)
 Function: writes restriction enzymes into the stream
 Returns : 1 for success and 0 for error
 Args    : a Bio::Restriction::Enzyme
           or a Bio::Restriction::EnzymeCollection object

=cut

sub write {
    my ($self,@h) = @_;
    $self->throw_not_implemented;
}

1;
