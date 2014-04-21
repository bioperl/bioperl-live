#
# BioPerl module for Bio::LiveSeq::Mutation
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::LiveSeq::Mutation - Mutation event descriptor class

=head1 SYNOPSIS

  # full descrition of a point mutation
  $mutation1a = Bio::LiveSeq::Mutation->new ( -seq => 'A',
					      -seqori => 'T',
					      -pos  => 100,
					      -len => 1 # optional, defaults to length(seq)
					     );

  # minimal information for a point mutation
  $mutation1b = Bio::LiveSeq::Mutation->new ( -seq => 'A',
					      -pos  => 100
					      );
  # insertion
  $mutation2 = Bio::LiveSeq::Mutation->new ( -seq => 'ATT',
					     -pos  => 100,
					     -len => 0
					     );
  # deletion
  $mutation3 = Bio::LiveSeq::Mutation->new ( -seq => '',  # optional
					     -seqori => 'TTG',  # optional
					     -pos  => 100
					     -len => 3
					     );
  # complex
  $mutation4 = Bio::LiveSeq::Mutation->new ( -seq => 'CC', 
					     -seqori => 'TTG',  # optional
					     -pos  => 100
					     -len => 3
					     );


=head1 DESCRIPTION

This class describes a local mutation event using minimalistic
description.  It is not necessary to know anything about the original
sequence. You need to give the changed sequence, the position of the
mutation in the (unidentified) reference sequence, and the length of
the affected subsequence in the reference sequence. If the original
allele sequence is given, the objects applying the mutation into the
reference sequence (e.g. L<Bio::LiveSeq::Mutator>) might check for its
validity.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the 
Bioperl mailing lists  Your participation is much appreciated.

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

report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::LiveSeq::Mutation;
use strict;

# Object preamble - inheritance


use base qw(Bio::Root::Root);

sub new {
    my($class,@args) = @_;
    my $self;
    $self = {};
    bless $self, $class;

    my ($seq, $seqori, $pos, $len, $label) =
	    $self->_rearrange([qw(SEQ
				  SEQORI
				  POS
				  LEN
				  )],
			      @args);

    $seq && $self->seq($seq);
    $seqori && $self->seqori($seqori);
    $pos && $self->pos($pos);
    defined($len) && $self->len($len); # defined() added otherwise won't work for len==0
    
    return $self; # success - we hope!
}


=head2 seq

 Title   : seq
 Usage   : $obj->seq();
 Function: 

            Sets and returns the mutated sequence. No checking is done
            to validate the symbols.

 Example : 
 Returns : string
 Args    : integer

=cut


sub seq {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'seq'} = $value;
    }
    return $self->{'seq'} || '';
}


=head2 seqori

 Title   : seqori
 Usage   : $obj->seqori();
 Function: 

            Sets and returns the original subsequence in the reference
            sequence. No checking is done to validate the symbols.
            Optional value.

 Example : 
 Returns : string
 Args    : string

=cut


sub seqori {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'seqori'} = $value;
    }
    return $self->{'seqori'} || '';
}


=head2 pos

 Title   : pos
 Usage   : $obj->pos();
 Function: 

            Sets and returns the position of the first element in the
            sequence.

 Example : 
 Returns : string
 Args    : integer

=cut


sub pos {
    my ($self,$value) = @_;
    if( defined $value) {
	if ( $value !~ /^([+-])?\d+$/ ) {
	    $self->throw("[$value] for pos has to be an integer\n");
	} else {
	    $self->{'pos'} = $value;
	}
    }
    return $self->{'pos'};
}

=head2 len

 Title   : len
 Usage   : $obj->len();
 Function: 

            Sets and returns the len of the affected original allele
            sequence.  If value is not set, defaults to the length of
            the mutated sequence (seq).

 Example : 
 Returns : string
 Args    : string

=cut

sub len {
   my ($self,$value) = @_;
   if ( defined $value) {
       $self->{'len'} = $value;
  }
   if ( ! exists $self->{'len'} ) {
       return length $self->{'seq'};
   }
   return $self->{'len'};
}

=head2 label

 Title   : label
 Usage   : $obj->label();
 Function: 

            Sets and returns the label of the affected original allele
            location. Label is a stable identifier whereas location
            can be changed by mutations. Label comes from
            l<Bio::LiveSeq::Gene>.

 Example : 
 Returns : string
 Args    : string

=cut

sub label {
   my ($self,$value) = @_;
   if ( defined $value) {
       $self->{'label'} = $value;
   }
   if ( ! exists $self->{'label'} ) {
       return;
   } 
   return $self->{'label'};
}


=head2 transpos

 Title   : transpos
 Usage   : $obj->transpos();
 Function: 

            Sets and returns the transcript position of the mutation.
            Set when associated with a reference sequence. Value
            depends on reference molecule and the co-ordinate system
            used.

 Example : 
 Returns : string
 Args    : integer

=cut


sub transpos {
    my ($self,$value) = @_;
    if( defined $value) {
	if ( $value !~ /^([+-])?\d+$/ ) {
	    $self->throw("[$value] for transpos has to be an integer\n");
	} else {
	    $self->{'transpos'} = $value;
	}
    }
    return $self->{'transpos'};
}


=head2 issue

 Title   : issue
 Usage   : $obj->issue();
 Function: 

            Sets and returns the position of the mutation in an array
            of mutations to be issued. Set after the validity of the
            mutation has been confirmed.

 Example : 
 Returns : string
 Args    : integer

=cut


sub issue {
    my ($self,$value) = @_;
    if( defined $value) {
	if ( $value !~ /^([+-])?\d+$/ ) {
	    $self->throw("[$value] for issue has to be an integer\n");
	} else {
	    $self->{'issue'} = $value;
	}
    }
    return $self->{'issue'};
}


=head2 prelabel

 Title   : prelabel
 Usage   : $obj->prelabel();
 Function: 

            Sets and returns the prelabel of the affected original allele
            location. Prelabel is a stable identifier whereas location
            can be changed by mutations. Prelabel comes from
            l<Bio::LiveSeq::Gene>.

 Example : 
 Returns : string
 Args    : string

=cut

sub prelabel {
   my ($self,$value) = @_;
   if ( defined $value) {
       $self->{'prelabel'} = $value;
   }
   if ( ! exists $self->{'prelabel'} ) {
       return;
   } 
   return $self->{'prelabel'};
}


=head2 postlabel

 Title   : postlabel
 Usage   : $obj->postlabel();
 Function: 

            Sets and returns the postlabel of the affected original allele
            location. Postlabel is a stable identifier whereas location
            can be changed by mutations. Postlabel comes from
            l<Bio::LiveSeq::Gene>.

 Example : 
 Returns : string
 Args    : string

=cut

sub postlabel {
   my ($self,$value) = @_;
   if ( defined $value) {
       $self->{'postlabel'} = $value;
   }
   if ( ! exists $self->{'postlabel'} ) {
       return;
   } 
   return $self->{'postlabel'};
}


=head2 lastlabel

 Title   : lastlabel
 Usage   : $obj->lastlabel();
 Function: 

            Sets and returns the lastlabel of the affected original allele
            location. Lastlabel is a stable identifier whereas location
            can be changed by mutations. Lastlabel comes from
            l<Bio::LiveSeq::Gene>.

 Example : 
 Returns : string
 Args    : string

=cut

sub lastlabel {
   my ($self,$value) = @_;
   if ( defined $value) {
       $self->{'lastlabel'} = $value;
   }
   if ( ! exists $self->{'lastlabel'} ) {
       return;
   } 
   return $self->{'lastlabel'};
}

1;
