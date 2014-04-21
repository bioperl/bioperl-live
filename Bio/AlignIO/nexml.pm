#
# BioPerl module for Bio::AlignIO::nexml
#
# Copyright Chase Miller
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::nexml - NeXML format sequence alignment input/output stream driver

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::AlignIO> class.

=head1 DESCRIPTION

This object can transform L<Bio::SimpleAlign> objects to and from
NeXML format. For more information on NeXML, visit L<http://www.nexml.org>.

=head1 FEEDBACK

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

=head1 AUTHORS

Chase Miller

=head1 CONTRIBUTORS

Mark Jensen, maj@fortinbras.us
Rutger Vos, rutgeraldo@gmail.com

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::nexml;

use strict;
use lib '../..';
use Bio::Nexml::Factory;
use Bio::Phylo::IO qw(parse unparse);

use base qw(Bio::AlignIO);


sub _initialize {
  my($self,@args) = @_;
  $self->SUPER::_initialize(@args);  
  $self->{_doc} = undef;
}

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln
 Function: returns the next alignment in the stream.
 Returns : Bio::Align::AlignI object - returns 0 on end of file
	        or on error
 Args    : 

See L<Bio::Align::AlignI>

=cut

sub next_aln {
	my ($self) = @_;
    unless ( $self->{'_parsed'} ) {
    	#use a parse function to load all the alignment objects found in the nexml file at once
        $self->_parse;
    }
    return $self->{'_alns'}->[ $self->{'_alnsiter'}++ ];
}

=head2 rewind

 Title   : rewind
 Usage   : $alnio->rewind
 Function: Resets the stream
 Returns : none
 Args    : none


=cut

sub rewind {
    my $self = shift;
    $self->{'_alniter'} = 0;
}

=head2 doc

 Title   : doc
 Usage   : $treeio->doc
 Function: Returns the biophylo nexml document object
 Returns : Bio::Phylo::Project
 Args    : none or Bio::Phylo::Project object

=cut

sub doc {
	my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_doc'} = $value;
	}
	return $obj->{'_doc'};
}

sub _parse {
	my ($self) = @_;

    $self->{'_parsed'}   = 1;
    $self->{'_alnsiter'} = 0;
    my $fac = Bio::Nexml::Factory->new();
	
	$self->doc(parse(
 	'-file'       => $self->{'_file'},
 	'-format'     => 'nexml',
 	'-as_project' => '1'
 	));

	$self->{'_alns'} = $fac->create_bperl_aln($self);
 	if(@{ $self->{'_alns'} } == 0)
 	{
 		self->debug("no seqs in $self->{_file}");
 	}
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in nexml format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object

See L<Bio::Align::AlignI>

=cut

sub write_aln {
	my ($self, $aln) = @_;
	
    my $fac = Bio::Nexml::Factory->new();
    my $taxa = $fac->create_bphylo_taxa($aln);
	my ($matrix) = $fac->create_bphylo_aln($aln, $taxa);
	$matrix->set_taxa($taxa);
	
	$self->doc(Bio::Phylo::Factory->create_project());
	$self->doc->insert($matrix);
	my $ret = $self->_print($self->doc->to_xml());
	$self->flush;
	return $ret;
}








1;
