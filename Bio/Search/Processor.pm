
#
# BioPerl module for Bio::Search::Processor
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Processor - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

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

=head1 AUTHOR - Aaron Mackey

Email amackey@virginia.edu

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::Processor;

use strict;

use Bio::Root::Version;

=head2 new

 Title   : new
 Usage   : $proc = Bio::Search::Processor->new -file      => $filename,
                                              -algorithm => 'Algorithm' ;
 Function: Used to specify and initialize a data processor of search
           algorithm results.
 Returns : A processor specific to the algorithm type, if it exists.
 Args    : -file => filename
           -algorithm => algorithm specifier
           -fh => filehandle to attach to (file or fh required)

=cut

sub new {

    my $type = shift;
    my $proc;
    my ($module, $load, $algorithm);

    my %args = @_;

    exists $args{'-algorithm'} or do { 
	print STDERR "Must supply an algorithm!";
	return;
    };

    $algorithm = $args{'-algorithm'} || $args{'-ALGORITHM'};

    $module = "_<Bio/Search/Processor/$algorithm.pm";
    $load = "Bio/Search/Processor/$algorithm.pm";

    unless ( $main::{$module} ) {
	eval { require $load; };
	if ( $@ ) {
	    print STDERR <<"EOF";
$load: $algorithm cannot be found
Exception $@
For more information about the Search/Processor system please see the
Processor docs.  This includes ways of checking for processors at 
compile time, not run time
EOF
	    return;
	}
    }

    $proc = "Bio::Search::Processor::$algorithm"->new(@_);
    return $proc;
}

1;
