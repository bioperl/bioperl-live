#-----------------------------------------------------------------
#
# BioPerl module Bio::AnalysisResultI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# Derived from Bio::Tools::AnalysisResult by Hilmar Lapp <hlapp@gmx.net>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

# POD documentation - main docs before the code

=head1 NAME

Bio::AnalysisResultI - Interface for analysis result objects

=head1 SYNOPSIS

Bio::AnalysisResultI defines an interface that must be implemented by
a subclass. So you cannot create Bio::AnalysisResultI objects,
only objects that inherit from Bio::AnalysisResultI. 

=head1 DESCRIPTION

The AnalysisResultI module provides an interface for modules
encapsulating the result of an analysis that was carried out with a
query sequence and an optional subject dataset.

The notion of an analysis represented by this base class is that of a unary or
binary operator, taking either one query or a query and a subject and producing
a result. The query is e.g. a sequence, and a subject is either a sequence,
too, or a database of sequences. 

This interface defines methods to access analysis result data and does
not impose any constraints on how the analysis result data is acquired.

Note that this module does not provide support for B<running> an analysis.
Rather, it is positioned in the subsequent parsing step (concerned with
turning raw results into BioPerl objects).

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
the bugs and their resolution.  Bug reports can be submitted via the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Steve Chervitz, Hilmar Lapp

Email sac@bioperl.org
Email hlapp@gmx.net (author of Bio::Tools::AnalysisResult on which this module is based)

=head1 COPYRIGHT

Copyright (c) 2001 Steve Chervitz. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::AnalysisResultI;
use strict;


use base qw(Bio::Root::RootI);


=head2 analysis_query

 Usage     : $query_obj = $result->analysis_query();
 Purpose   : Get a Bio::PrimarySeqI-compatible object representing the entity 
             on which the analysis was performed. Lacks sequence information.
 Argument  : n/a
 Returns   : A Bio::PrimarySeqI-compatible object without sequence information.
             The sequence will have display_id, description, moltype, and length data.

=cut

#---------------------
sub analysis_query {
#---------------------
    my ($self) = @_;
    $self->throw_not_implemented;
}


=head2 analysis_subject

 Usage     : $obj = $result->analyis_subject();
 Purpose   : Get the subject of the analysis against which it was
             performed. For similarity searches it will probably be a database,
             and for sequence feature predictions (exons, promoters, etc) it
             may be a collection of models or homologous sequences that were
             used, or undefined.
 Returns   : An object of a type the depends on the implementation
             May also return undef for analyses that don\'t involve subjects.
 Argument  : n/a
 Comments  : Implementation of this method is optional.
             AnalysisResultI provides a default behavior of returning undef.

=cut

#---------------
sub analysis_subject { 
#---------------
    my ($self) = @_; 
    return;
}

=head2 analysis_subject_version

 Usage     : $vers = $result->analyis_subject_version();
 Purpose   : Get the version string of the subject of the analysis.
 Returns   : String or undef for analyses that don\'t involve subjects.
 Argument  : n/a
 Comments  : Implementation of this method is optional.
             AnalysisResultI provides a default behavior of returning undef.

=cut

#---------------
sub analysis_subject_version { 
#---------------
    my ($self) = @_; 
    return;
}


=head2 analysis_date

 Usage     : $date = $result->analysis_date();
 Purpose   : Get the date on which the analysis was performed.
 Returns   : String
 Argument  : n/a

=cut

#---------------------
sub analysis_date {
#---------------------
    my ($self) = @_;
    $self->throw_not_implemented;
}

=head2 analysis_method

 Usage     : $meth = $result->analysis_method();
 Purpose   : Get the name of the sequence analysis method that was used
             to produce this result (BLASTP, FASTA, etc.). May also be the
             actual name of a program.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method { 
#-------------
    my ($self) = @_;  
    $self->throw_not_implemented;
}

=head2 analysis_method_version

 Usage     : $vers = $result->analysis_method_version();
 Purpose   : Get the version string of the analysis program.
           : (e.g., 1.4.9MP, 2.0a19MP-WashU).
 Returns   : String
 Argument  : n/a

=cut

#---------------------
sub analysis_method_version {
#---------------------
    my ($self) = @_; 
    $self->throw_not_implemented;
}

=head2 next_feature

 Title   : next_feature
 Usage   : $seqfeature = $obj->next_feature();
 Function: Returns the next feature available in the analysis result, or
           undef if there are no more features.
 Example :
 Returns : A Bio::SeqFeatureI implementing object, or undef if there are no
           more features.
 Args    : none

=cut

#---------------------
sub next_feature {
#---------------------
    my ($self);
    $self->throw_not_implemented;
}


1;
