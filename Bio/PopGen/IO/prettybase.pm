#
# BioPerl module for Bio::PopGen::IO::prettybase
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PopGen::IO::prettybase - Extract individual allele data from PrettyBase format

=head1 SYNOPSIS

Do not use directly, use through the Bio::PopGen::IO driver

=head1 DESCRIPTION

This object will parse comma delimited PrettyBase output.  PrettyBase
is defined by the SeattleSNPs http://pga.gs.washington.edu/

This is expected to be tab delimited (you can vary with the
field_delimiter flag SITE SAMPLE ALLELE1 ALLELE2

There are 2 initialization parameters, the delimiter
(-field_delimiter) [default 'tab'] and a boolean -no_header which
specifies if there is no header line to read in.  All lines starting
with '#' will be skipped

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

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 CONTRIBUTORS

Matthew Hahn, matthew.hahn-at-duke.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PopGen::IO::prettybase;
use vars qw($FieldDelim $Header);
use strict;

($FieldDelim,$Header) =( '\t',0);


use Bio::PopGen::Individual;
use Bio::PopGen::Population;
use Bio::PopGen::Genotype;

use base qw(Bio::PopGen::IO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::PopGen::IO::prettybase->new();
 Function: Builds a new Bio::PopGen::IO::prettybase object 
 Returns : an instance of Bio::PopGen::IO::prettybase
 Args    : -field_delimiter      => a field delimiter character or regexp (default is /\t/ ) 
           -header               => boolean if the file will have a header and parser should
                                    skip first line in the file (default is false)
           -convert_indel_states => convert alleles which are longer than one character
                                    to an 'I' meaning insert state, and alleles which are
                                    '-' to a delete state.
                                    (default is false)

=cut

sub _initialize {
    my($self, @args) = @_;
    my ($fieldsep,
	$conv_indels,
	$header) = $self->_rearrange([qw(FIELD_DELIMITER
					 CONVERT_INDEL_STATES
					 HEADER)],@args);

    $self->flag('header', defined $header ? $header : $Header);
    $self->flag('field_delimiter',defined $fieldsep ? $fieldsep : $FieldDelim);
    $self->{'_header'} = undef;
    $self->{'_parsed_individiuals'} = [];
    $self->{'_parsed'} = 0;
    $self->flag('convert_indel',$conv_indels || 0);
    return 1;
}

=head2 flag

 Title   : flag
 Usage   : $obj->flag($flagname,$newval)
 Function: Get/Set the flag value
 Returns : value of a flag (a boolean)
 Args    : A flag name, currently we expect 
           'header', 'field_delimiter', or 'allele_delimiter' 
           on set, new value (a boolean or undef, optional)


=cut

sub flag{
    my $self = shift;
    my $fieldname = shift;
    return unless defined $fieldname;
    
    return $self->{'_flag'}->{$fieldname} = shift if @_;
    return $self->{'_flag'}->{$fieldname};
}


=head2 next_individual

 Title   : next_individual
 Usage   : my $ind = $popgenio->next_individual;
 Function: Retrieve the next individual from a dataset
 Returns : Bio::PopGen::IndividualI object
 Args    : none


=cut

sub next_individual {
    my ($self) = @_;
    unless( $self->{'_parsed'} ) {
	$self->_parse_prettybase;
    }
    return $self->{'_parsed_individiuals'}->[$self->{'_iterator'}++];
}



=head2 next_population

 Title   : next_population
 Usage   : my $ind = $popgenio->next_population;
 Function: Retrieve the next population from a dataset
 Returns : Bio::PopGen::PopulationI object
 Args    : none
 Note    : Many implementation will not implement this

=cut

# Plan is to just return the whole dataset as a single population by 
# default I think - people would then have each population in a separate
# file.

sub next_population{
    my ($self) = @_;
    my @inds;
    while( my $ind = $self->next_individual ) {
	push @inds, $ind;
    }
    return unless @inds;
    Bio::PopGen::Population->new(-individuals => \@inds);
}


sub _parse_prettybase {
    my $self = shift;
    my %inds;
    my $convert_indels = $self->flag('convert_indel');
    while( defined( $_ = $self->_readline) ) {
	next if( /^\s*\#/ || /^\s+$/ || ! length($_) );
	
	my ($site,$sample,@alleles) = split($self->flag('field_delimiter'),$_);
	if( ! defined $sample ) { 
	    warn("sample id is undefined for $_");
	    next;
	}
	for my $allele ( @alleles ) {
	    $allele =~ s/^\s+//;
	    $allele =~ s/\s+$//;
	    if( $convert_indels ) {
		if( length($allele) > 1 ) {
		    # we have an insert state
		    $allele = 'I';
		} elsif( $allele eq '-' ) {
		    # have a delete state
		    $allele = 'D';
		}
	    }
	}
	
	my $g = Bio::PopGen::Genotype->new(-alleles      => \@alleles,
					  -marker_name  => $site,
					  -individual_id=> $sample); 
	

	if( ! defined $inds{$sample} ) {
	    $inds{$sample} = Bio::PopGen::Individual->new(-unique_id => $sample);
	}
	$inds{$sample}->add_Genotype($g);
    }
    $self->{'_parsed_individiuals'} = [ values %inds ];
    $self->{'_parsed'} = 1;
    return;
}


=head2 write_individual

 Title   : write_individual
 Usage   : $popgenio->write_individual($ind);
 Function: Write an individual out in the file format
 Returns : none
 Args    : L<Bio::PopGen::PopulationI> object(s)

=cut

sub write_individual{
    my ($self,@inds) = @_;
    foreach my $ind ( @inds ) {
	if (! ref($ind) || ! $ind->isa('Bio::PopGen::IndividualI') ) {
	    $self->warn("Cannot write an object that is not a Bio::PopGen::IndividualI object");
	    next;
	}
	foreach my $marker ( $ind->get_marker_names ) { 
	    my $g = $ind->get_Genotypes(-marker=> $marker);
	    next unless defined $g;
	    $self->_print( join("\t", $marker, $ind->unique_id, 
				$g->get_Alleles), "\n");	    
	}
    }
    
}



=head2 write_population

 Title   : write_population
 Usage   : $popgenio->write_population($pop);
 Function: Write a population out in the file format
 Returns : none
 Args    : L<Bio::PopGen::PopulationI> object(s)
 Note    : Many implementation will not implement this

=cut

sub write_population{
    my ($self,@pops) = @_;
    foreach my $pop ( @pops ) {
	if (! ref($pop) || ! $pop->isa('Bio::PopGen::PopulationI') ) {
	    $self->warn("Cannot write an object that is not a Bio::PopGen::PopulationI object");
	    next;
	}
	my @mnames = $pop->get_marker_names;
	foreach my $ind ( $pop->get_Individuals ) {
	    if (! ref($ind) || ! $ind->isa('Bio::PopGen::IndividualI') ) {
		$self->warn("Cannot write an object that is not a Bio::PopGen::IndividualI object");
		next;
	    }
	    foreach my $marker ( @mnames ) { 
		my $g = $ind->get_Genotypes(-marker=> $marker);
		next unless defined $g;
		$self->_print( join("\t", $marker, $ind->unique_id, 
				    $g->get_Alleles), "\n");
			   
	    }
	}
    }
}

1;
