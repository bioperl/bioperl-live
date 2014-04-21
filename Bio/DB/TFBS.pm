# $Id: TFBS.pm,v 1.11 2006/08/12 11:00:03 sendu Exp $
#
# BioPerl module for Bio::DB::TFBS
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

Bio::DB::TFBS - Access to a Transcription Factor Binding Site database

=head1 SYNOPSIS

  use Bio::DB::TFBS;

  my $db = Bio::DB::TFBS->new(-source => 'transfac');
  my ($factor_id) = $db->get_factor_ids('PPAR-gamma1');
  my ($matrix_id) = $db->get_matrix_ids('PPAR-gamma1');

  # get a Bio::Map::TranscriptionFactor with all the positions of a given factor
  my $factor = $db->get_factor(-factor_id => $factor_id);

  # get a Bio::Map::GeneMap containing all the factors that bind near a given gene
  my $gene_map = $db->get_gene_map(-gene_name => 'AQP 7');

  # get a PSM (Bio::Matrix::PSM) of a given matrix
  my $psm = $db->get_matrix(-matrix_id => $matrix_id);

  # get the aligned sequences (Bio::SimpleAlign) that were used to build a given
  # matrix
  my $align = $db->get_alignment(-matrix_id => $matrix_id);

  # get a specific instance sequence (Bio::LocatableSeq)
  my $seq = $db->get_seq($id);

=head1 DESCRIPTION

This is a front end module for access to a Transcription Factor Binding Site
database.

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

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 CONTRIBUTORS

Based on Bio::DB::Taxonomy by Jason Stajich

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::TFBS;
use strict;

use Bio::Root::Root;

use base qw(Bio::Root::Root);

our $DefaultSource = 'transfac';

=head2 new

 Title   : new
 Usage   : my $obj = Bio::DB::TFBS->new(-source => 'transfac');
 Function: Builds a new Bio::DB::TFBS object.
 Returns : an instance of Bio::DB::TFBS
 Args    : -source => which database source: currently only 'transfac_pro'

=cut

sub new {
    my ($class, @args) = @_;
  
    if ($class =~ /Bio::DB::TFBS::(\S+)/) {
        my ($self) = $class->SUPER::new(@args);
        $self->_initialize(@args);
        return $self;
    }
    else { 
        my %param = @args;
        @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
        my $source = $param{'-source'} || $DefaultSource;
        
        $source = "\L$source";	# normalize capitalization to lower case
        
        # normalize capitalization
        return unless( $class->_load_tax_module($source) );
        return "Bio::DB::TFBS::$source"->new(@args);
    }
}

# empty for now
sub _initialize { }

=head2 _load_tax_module

 Title   : _load_tax_module
 Usage   : *INTERNAL Bio::DB::TFBS stuff*
 Function: Loads up (like use) a module at run time on demand

=cut

sub _load_tax_module {
    my ($self, $source) = @_;
    my $module = "Bio::DB::TFBS::" . $source;
    my $ok;

    eval { $ok = $self->_load_module($module) };
    if ( $@ ) {
	print STDERR $@;
	print STDERR <<END;
$self: $source cannot be found
Exception $@
For more information about the Bio::DB::TFBS system please see
the Bio::DB::TFBS docs.  This includes ways of checking for 
formats at compile time, not run time.
END
  ;
    }
    return $ok;
}

1;
