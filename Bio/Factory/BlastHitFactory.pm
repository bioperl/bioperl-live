#-----------------------------------------------------------------
# $Id$
#
# BioPerl module for Bio::Factory::BlastHitFactory
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

# POD documentation - main docs before the code

=head1 NAME

Bio::Factory::BlastHitFactory - Factory for Bio::Search::Hit::BlastHit objects

=head1 SYNOPSIS

    use Bio::Factory::BlastHitFactory;

    my $hit_fact = Bio::Factory::BlastHitFactory->new();

    my $hit = $hit_fact->create_hit( %parameters );

See documentation for create_hit() for information about C<%parameters>.

=head1 DESCRIPTION

This module encapsulates code for creating Bio::Search::Hit::BlastHit
and Bio::Search::HSP::BlastHSP objects from traditional BLAST report
data (i.e., non-XML formatted).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                - General discussion
  http://bio.perl.org/MailList.html    - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR 

Steve Chervitz E<lt>sac@bioperl.orgE<gt>

See L<the FEEDBACK section | FEEDBACK> for where to send bug reports and comments.

=head1 COPYRIGHT

Copyright (c) 2001 Steve Chervitz. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

#'

package Bio::Factory::BlastHitFactory;

use strict;
use Bio::Root::Root;
use Bio::Factory::HitFactoryI;
use Bio::Search::Hit::BlastHit;

use vars qw(@ISA);

@ISA = qw(Bio::Root::Root Bio::Factory::HitFactoryI); 

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    return $self;
}

=head2 create_hit

 Title   : create_hit
 Usage   : $hit = $factory->create_hit( %params );
 Function: Creates a new Bio::Search::Hit::BlastHit object given 
           raw BLAST report data, formatted in traditional BLAST report format.
 Returns : A single Bio::Search::Hit::BlastHit object
 Args    : Named parameters to be passed to the BlastHit object.
           Parameter keys are case-insensitive.
           See Bio::Search::Hit::BlastHit::new() documentation for 
           details about these parameters.
           The only additional parameter required is:
              -RESULT    => a Bio::Search::Result::BlastResult object.
           From this result object, the program, query length, 
           and iteration are obtained and passed on to the BlastHit.

=cut

sub create_hit {
    my ($self, @args) = @_;

    my ($blast, $raw_data, $shallow_parse) =
      $self->_rearrange( [qw(RESULT
			     RAW_DATA
			     SHALLOW_PARSE)], @args);

    my %args = @args;
    $args{'-PROGRAM'}   = $blast->analysis_method;
    $args{'-QUERY_LEN'} = $blast->query_length;
    $args{'-ITERATION'} = $blast->iterations;

    my $hit = Bio::Search::Hit::BlastHit->new( %args );
    
    unless( $shallow_parse ) {
      $self->_add_hsps( $hit, 
			$args{'-PROGRAM'}, 
			$args{'-QUERY_LEN'}, 
			$blast->query_name, 
			@{$raw_data} );
    }

    return $hit;
}

#=head2 _add_hsps
#
# Usage     : Private method; called automatically by create_hit().
# Purpose   : Creates BlastHSP.pm objects for each HSP in a BLAST hit alignment.
#           : Also collects the full description of the hit from the
#           : HSP alignment section.
# Returns   : n/a
# Argument  : (<$BlastHit_object>, <$program_name>, <$query_length>, <$query_name>, <@raw_data>
#             'raw data list' consists of traditional BLAST report 
#             format for a single HSP, supplied as a list of strings.
# Throws    : Warnings for each BlastHSP.pm object that fails to be constructed.
#           : Exception if no BlastHSP.pm objects can be constructed.
#           : Exception if can't parse length data for hit sequence.
# Comments  : Requires Bio::Search::HSP::BlastHSP.pm.
#           : Sets the description using the full string present in 
#           : the alignment data.
#=cut

#--------------
sub _add_hsps { 
#--------------
    my( $self, $hit, $prog, $qlen, $qname, @data ) = @_;
    my $start     = 0;
    my $hspCount  = 0;

    require Bio::Search::HSP::BlastHSP;

#    printf STDERR "\nBlastHit \"$hit\" _process_hsps(). \nDATA (%d lines) =\n@data\n", scalar(@data);

    my( @hspData, @hspList, @errs, @bad_names );
    my($line, $set_desc, @desc);
    $set_desc = 0;
    my $hname = $hit->name;
    my $hlen;

    hit_loop:
   foreach $line( @data ) {

       if( $line =~ /^\s*Length = ([\d,]+)/ ) {
	   $hit->_set_description(@desc);
	   $set_desc = 1;
	   $hit->_set_length($1);
           $hlen = $hit->length;
	   next hit_loop;
       } elsif( !$set_desc) {
	   $line =~ s/^\s+|\s+$//g;
	   push @desc, $line;
	   next hit_loop;
       } elsif( $line =~ /^\s*Score/ ) {
	   ## This block is for setting multiple HSPs.

	   if( not scalar @hspData ) {
	       $start = 1; 
	       push @hspData, $line; 
	       next hit_loop;

	    } elsif( scalar @hspData) {  
		$hspCount++;
		$self->verbose and do{ print STDERR +( $hspCount % 10 ? "+" : "+\n" ); };

#		print STDERR "\nBlastHit: setting HSP #$hspCount \n@hspData\n";
		my $hspObj =  Bio::Search::HSP::BlastHSP->new
				      (-RAW_DATA   => \@hspData, 
				       -RANK       => $hspCount,
				       -PROGRAM    => $prog,
				       -QUERY_NAME => $qname,
				       -HIT_NAME   => $hname,
				      ); 
		push @hspList, $hspObj;
		@hspData = ();
		push @hspData, $line;
		next;
	   } else {
	       push @hspData, $line;
	   }
       } elsif( $start ) {
	   ## This block is for setting the last HSP (which may be the first as well!).
	   if( $line =~ /^(end|>|Parameters|CPU|Database:)/ ) {
	       $hspCount++;
	       $self->verbose and do{ print STDERR +( $hspCount % 10 ? "+" : "+\n" ); };

#	       print STDERR "\nBlastHit: setting HSP #$hspCount \n@hspData"; 

	       my $hspObj = Bio::Search::HSP::BlastHSP->new
				     (-RAW_DATA   => \@hspData, 
				      -RANK       => $hspCount,
				      -PROGRAM    => $prog,
				      -QUERY_NAME => $qname,
				      -HIT_NAME   => $hname,
				     );
	       push @hspList, $hspObj;
	   } else {
	       push @hspData, $line;
	   }
       }
   }		

    $hit->{'_length'} or $self->throw( "Can't determine hit sequence length.");

    # Adjust logical length based on BLAST flavor.
    if($prog =~ /TBLAST[NX]/) {
	$hit->{'_logical_length'} = $hit->{'_length'} / 3;
    }

    $hit->{'_hsps'} = [ @hspList ];

#    print STDERR "\n--------> Done building HSPs for $hit (total HSPS: ${\$hit->num_hsps})\n";

}



1;
