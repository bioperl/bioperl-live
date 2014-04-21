#
# BioPerl module for Bio::Matrix::Scoring
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Matrix::Scoring - Object which can hold scoring matrix information

=head1 SYNOPSIS

  use Bio::Matrix::Scoring;

=head1 DESCRIPTION

An object which can handle AA or NT scoring matrix information.  Some
transformation properties are available too.

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

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Matrix::Scoring;
use strict;


use base qw(Bio::Matrix::Generic);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Matrix::Scoring->new();
 Function: Builds a new Bio::Matrix::Scoring object 
 Returns : an instance of Bio::Matrix::Scoring
 Args    :


=cut


sub new { 
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($entropy,$expected,$scale,$scaleval,$database,
		  $lowestscore,$highestscore,$lambda,$H) = 
			 $self->_rearrange([qw(
        ENTROPY EXPECTED SCALE SCALE_VALUE DATABASE
		  LOWEST_SCORE HIGHEST_SCORE LAMBDA H)], @args);

    $self->entropy  ($entropy);
    $self->expected_score($expected);
    $self->scale    ($scale);
    $self->scale_value($scaleval);
    $self->database ($database);
    $self->lowest_score($lowestscore);
    $self->highest_score($highestscore);
    $self->lambda($lambda);
    $self->H($H);
				    
    return $self;
}

=head2 entropy

 Title   : entropy
 Usage   : $obj->entropy($newval)
 Function: 
 Example : 
 Returns : value of entropy (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub entropy{
    my $self = shift;

    return $self->{'entropy'} = shift if @_;
    return $self->{'entropy'};
}

=head2 expected_score

 Title   : expected_score
 Usage   : $obj->expected_score($newval)
 Function: 
 Example : 
 Returns : value of expected (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub expected_score{
    my $self = shift;

    return $self->{'expected'} = shift if @_;
    return $self->{'expected'};
}

=head2 scale

 Title   : scale
 Usage   : $obj->scale($newval)
 Function: 
 Example : 
 Returns : value of scale (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub scale{
    my $self = shift;

    return $self->{'scale'} = shift if @_;
    return $self->{'scale'};
}

=head2 scale_value

 Title   : scale_value
 Usage   : $obj->scale_value($newval)
 Function: 
 Example : 
 Returns : value of scale_value (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub scale_value{
    my $self = shift;

    return $self->{'scale_value'} = shift if @_;
    return $self->{'scale_value'};
}

=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: 
 Example : 
 Returns : value of description (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub description{
    my $self = shift;

    return $self->{'description'} = shift if @_;
    return $self->{'description'};
}

=head2 database

 Title   : database
 Usage   : $obj->database($newval)
 Function: 
 Example : 
 Returns : value of database (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub database{
    my $self = shift;

    return $self->{'database'} = shift if @_;
    return $self->{'database'};
}

=head2 lowest_score

 Title   : lowest_score
 Usage   : $obj->lowest_score($newval)
 Function: 
 Example : 
 Returns : value of lowest_score (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub lowest_score{
    my $self = shift;

    return $self->{'lowest_score'} = shift if @_;
    return $self->{'lowest_score'};
}

=head2 highest_score

 Title   : highest_score
 Usage   : $obj->highest_score($newval)
 Function: 
 Example : 
 Returns : value of highest_score (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub highest_score{
    my $self = shift;

    return $self->{'highest_score'} = shift if @_;
    return $self->{'highest_score'};
}

=head2 lambda

 Title   : lambda
 Usage   : $obj->lambda($newval)
 Function: 
 Example : 
 Returns : value of lambda (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub lambda{
    my $self = shift;

    return $self->{'lambda'} = shift if @_;
    return $self->{'lambda'};
}

=head2 H

 Title   : H
 Usage   : $obj->H($newval)
 Function: 
 Example : 
 Returns : value of H (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub H{
    my $self = shift;
    return $self->{'H'} = shift if @_;
    return $self->{'H'};
}

1;
