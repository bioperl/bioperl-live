# $Id$
#
# BioPerl module for Bio::SearchIO::chadosxpr
#
# Chris Mungall <cjm@fruitfly.org>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::chadosxpr - chadosxpr sequence input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the SearchIO handler system. Go:

    $stream = Bio::SearchIO->new(-file => $filename, -format => 'chadosxpr');

    while ( my $seq = $stream->next_seq() ) {
	# do something with $seq
    }

=head1 DESCRIPTION

This object can transform Bio::Search objects to and from chadosxpr flat
file databases. CURRENTLY ONLY TO


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://www.bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Chris Mungall

Email cjm@fruitfly.org

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SearchIO::chadosxpr;
use Bio::SearchIO::chado;
use vars qw(@ISA);
use strict;

use Data::Stag::SxprWriter;

@ISA = qw(Bio::SearchIO::chado);

sub default_handler_class {
    return "Data::Stag::SxprWriter";
} 

1;
