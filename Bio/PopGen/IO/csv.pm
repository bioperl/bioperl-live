# $Id$
#
# BioPerl module for Bio::PopGen::IO::csv
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PopGen::IO::csv -Extract individual allele data from a CSV parser 

=head1 SYNOPSIS

Do not use directly, use through the Bio::PopGen::IO driver

=head1 DESCRIPTION

This object will parse comma delimited format (CSV) or whatever
delimiter you specify. It currently doesn't handle the more complex
quote escaped CSV format.  There are 3 initialization parameters, 
the delimiter (-field_delimiter) [default ','], (-allele_delimiter) 
[default ' '].    The third initialization parameter is a boolean 
-no_header which specifies if there is no header line to read in.  All lines starting with '#' will be skipped

When no_header is not specific the data is assumed to be of the following form.
Having a header line this
SAMPLE,MARKERNAME1,MARKERNAME2,...

and each data line having the form (diploid data)
SAMP1,101 102,100 90,a b
or for haploid data
SAMP1,101,100,a

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 CONTRIBUTORS

Matthew Hahn, matthew.hahn-at-duke.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PopGen::IO::csv;
use vars qw(@ISA $FieldDelim $AlleleDelim $NoHeader);
use strict;

($FieldDelim,$AlleleDelim,$NoHeader) =( ',', '\s+',0);

# Object preamble - inherits from Bio::Root::Root

use Bio::PopGen::IO;

use Bio::PopGen::Individual;
use Bio::PopGen::Population;
use Bio::PopGen::Genotype;

@ISA = qw(Bio::PopGen::IO );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::PopGen::IO::csv();
 Function: Builds a new Bio::PopGen::IO::csv object 
 Returns : an instance of Bio::PopGen::IO::csv
 Args    :


=cut

sub _initialize {
    my($self, @args) = @_;
    my ($fieldsep,$all_sep, 
	$noheader) = $self->_rearrange([qw(FIELD_DELIMITER
					   ALLELE_DELIMITER
					   NO_HEADER)],@args);

    $self->flag('no_header', defined $noheader ? $noheader : $NoHeader);
    $self->flag('field_delimiter',defined $fieldsep ? $fieldsep : $FieldDelim);
    $self->flag('allele_delimiter',defined $all_sep ? $all_sep : $AlleleDelim);

    $self->{'_header'} = undef;
    return 1;
}

=head2 flag

 Title   : flag
 Usage   : $obj->flag($flagname,$newval)
 Function: Get/Set the flag value
 Returns : value of a flag (a boolean)
 Args    : A flag name, currently we expect 
           'no_header', 'field_delimiter', or 'allele_delimiter' 
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

sub next_individual{
    my ($self) = @_;
    while( defined( $_ = $self->_readline) ) {
	next if( /^\s*\#/ || /^\s+$/ || ! length($_) );
	last;
    }
    return if ! defined $_; 
    if( $self->flag('no_header') || 
	defined $self->{'_header'} ) {
	my ($samp,@marker_results) = split($self->flag('field_delimiter'),$_);
	my $i = 1;
	foreach my $m ( @marker_results ) {
	    $m =~ s/^\s+//;
	    $m =~ s/\s+$//;
	    my $markername;
	    if( defined $self->{'_header'} ) {
		$markername = $self->{'_header'}->[$i];
	    } else { 
		$markername = "Marker$i";
	    }
	    $self->debug( "markername is $markername alleles are $m\n");
	    my @alleles = split($self->flag('allele_delimiter'), $m);
	    $m = new Bio::PopGen::Genotype(-alleles      => \@alleles,
					   -marker_name  => $markername,
					   -individual_id=> $samp); 
	    $i++; 
	}
	return new Bio::PopGen::Individual(-unique_id => $samp,
					   -genotypes => \@marker_results);
    } else {
	chomp;
	$self->{'_header'} = [split($self->flag('field_delimiter'),$_)];
	return $self->next_individual; # rerun loop again
    }
    return undef;
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
    $self->warn("Sorry current this implementation doesn't know how to extract populations from CSV data");
    return undef;
}



1;
