package Bio::Seq::SimulatedRead;


=head1 NAME

Bio::Seq::SimulatedRead - Read with sequencing errors taken from a reference sequence

=head1 SYNOPSIS

    use Bio::Seq::SimulatedRead;
    use Bio::PrimarySeq;
  
    # Create a reference sequence
    my $genome = Bio::PrimarySeq->new( -id    => 'human_chr2',
                                       -seq   => 'TAAAAAAACCCCTG',
                                       -desc  => 'The human genome' );

    # A 10-bp error-free read taken from a genome
    my $read = Bio::Seq::SimulatedRead->new(
       -reference => $genome  , # sequence to generate the read from
       -id        => 'read001', # read ID
       -start     => 3        , # start of the read on the genome forward strand
       -end       => 12       , # end of the read on the genome forward strand
       -strand    => 1        , # genome strand that the read is on
    );

    # Display the sequence of the read
    print $read->seq."\n";

    # Add a tag or MID to the beginning of the read
    $read->mid('ACGT');

    # Add sequencing errors (error positions are 1-based and relative to the
    # error-free MID-containing read)
    my $errors = {};
    $errors->{'8'}->{'+'} = 'AAA'; # insertion of AAA after residue 8
    $errors->{'1'}->{'%'} = 'G';   # substitution of residue 1 by a G
    $errors->{'4'}->{'-'} = undef; # deletion of residue 4
    $read->errors($errors);

    # Display the sequence of the read with errors
    print $read->seq."\n";

    # String representation of where the read came from and its errors
    print $read->desc."\n";

=head1 DESCRIPTION

This object is a simulated read with sequencing errors. The user can provide a
reference sequence to take a read from, the position and orientation of the 
read on the reference sequence, and the sequencing errors to generate.

The sequence of the read is automatically calculated based on this information.
By default, the description of the reads contain tracking information and will
look like this (Bioperl-style):

  reference=human_chr2 start=3 end=12 strand=-1 mid=ACGT errors=1%G,4-,8+AAA description="The human genome"

or Genbank-style:

  reference=human_chr2 position=complement(3..12) mid=ACGT errors=1%G,4-,8+AAA description="The human genome"

Creating a simulated read follows these steps:
  1/ Define the read start(), end(), strand() and qual_levels() if you want
     quality scores to be generated. Do not change these values once set because
     the read will not be updated.
  2/ Specify the reference sequence that the read should be taken from. Once
     this is done, you have a fully functional read. Do not use the reference()
     method again after you have gone to the next step.
  3/ Use mid() to input a MID (or tag or barcode) to add to the beginning of the
     read. You can change the MID until you go to next step.
  4/ Give sequencing error specifications using errors() as the last step. You
     can do that as many times as you like, and the read will be updated.

=head1 AUTHOR

Florent E Angly E<lt>florent . angly @ gmail-dot-comE<gt>.

Copyright (c) 2011 Florent E Angly.

This library is free software; you can redistribute it under the GNU General
Public License version 3.

=cut


use strict;
use warnings;
use Bio::LocatableSeq;
use base qw( Bio::Seq::Quality Bio::LocatableSeq );


=head2 new

 Title    : new
 Function : Create a new simulated read object
 Usage    : my $read = Bio::Seq::SimulatedRead->new(
               -id        => 'read001',
               -reference => $seq_obj ,
               -errors    => $errors  ,
               -start     => 10       ,
               -end       => 135      ,
               -strand    => 1        ,
            );
 Arguments: -reference   => Bio::SeqI, Bio::PrimarySeqI object representing the
                            reference sequence to take the read from. See
                            reference().
            -errors      => Hashref representing the position of errors in the read
                            See errors().
            -mid         => String of a MID to prepend to the read. See mid().
            -track       => Track where the read came from in the read description?
                            See track().
            -coord_style => Define what coordinate system to use. See coord_style().
            All other methods from Bio::LocatableSeq are available.
 Returns  : new Bio::Seq::SimulatedRead object

=cut

sub new {
   my ($class, @args) = @_;
   my $self = $class->SUPER::new(@args);
   my ($qual_levels, $reference, $mid, $errors, $track, $coord_style) =
      $self->_rearrange([qw(QUAL_LEVELS REFERENCE MID ERRORS TRACK COORD_STYLE)], @args);
   $coord_style = defined $coord_style ? $coord_style : 'bioperl';
   $self->coord_style($coord_style);
   $track = defined $track ? $track : 1;
   $self->track($track);
   $qual_levels = defined $qual_levels ? $qual_levels : [];
   $self->qual_levels($qual_levels) if defined $qual_levels;
   $self->reference($reference) if defined $reference;
   $self->mid($mid) if defined $mid;
   $self->{_mutated} = 0;
   $self->errors($errors) if defined $errors;
   return $self;
}


=head2 qual_levels

 Title    : qual_levels
 Function : Get or set the quality scores to give to the read. By default, if your
            reference sequence does not have quality scores, no quality scores
            are generated for the simulated read. The generated quality scores
            are very basic. If a residue is error-free, it gets the quality score
            defined for good residues. If the residue has an error (is an
            addition or a mutation), the residue gets the quality score specified
            for bad residues. Call the qual_levels() method before using the
            reference() method.
 Usage    : my $qual_levels = $read->qual_levels( );
 Arguments: Array reference containing the quality scores to use for:
            1/ good residues (e.g. 30)
            2/ bad residues (e.g. 10)
 Returns  : Array reference containing the quality scores to use.

=cut

sub qual_levels {
   my ($self, $qual_levels) = @_;
   if (defined $qual_levels) {
      if ( (scalar @$qual_levels != 0) && (scalar @$qual_levels != 2) ) {
         $self->throw("The quality score specification must define the score".
            " to use for good and for bad residues\n");
      }
      $self->{qual_levels} = $qual_levels;
   }
   return $self->{qual_levels};
}


=head2 reference

 Title    : reference
 Function : Get or set the reference sequence that the read comes from. Once the
            reference has been set, you have a functional simulated read which
            supports all the Bio::LocatableSeq methods. This method must be
            called after qual_levels() but before mid() or errors().
 Usage    : my $seq_obj = $read->reference();
 Arguments: Bio::SeqI or Bio::PrimarySeqI object
 Returns  : Bio::SeqI or Bio::PrimarySeqI object

=cut

sub reference {
   my ($self, $reference) = @_;
   if (defined $reference) {
      # Sanity check 1
      if ( (not $reference->isa('Bio::SeqI')) && (not $reference->isa('Bio::PrimarySeqI')) ) {
         $self->throw("Expected a Bio::SeqI object as reference, but got: $reference\n");
      }
      # Sanity check 2
      if ($self->{mid} || $self->{errors}) {
         $self->throw("Cannot change the reference sequence after an MID or ".
            "sequencing errors have been added to the read\n");
      }
      # Use beginning of reference sequence as start default
      if (not defined $self->start) {
         $self->start(1);
      }
      # Use end of reference sequence as end default
      if (not defined $self->end) {
         $self->end($reference->length);
      }
      # Use strand 1 as strand default
      if (not defined $self->strand) {
         $self->strand(1);
      }
      # Set the reference sequence object
      $self->{reference} = $reference;
      # Create a sequence, quality scores and description from the reference
      $self->_create_seq;
      $self->_create_qual if scalar @{$self->qual_levels};
      $self->_create_desc if $self->track;
   }
   return $self->{reference};
}


sub _create_seq {
   my $self = shift;
   # Get a truncation of the reference sequence
   my $reference = $self->reference;
   my $read_obj = $reference->trunc( $self->start, $self->end );
   # Reverse complement the read if needed
   if ($self->strand == -1) {
      $read_obj = $read_obj->revcom();
   }
   $self->seq($read_obj->seq);
   return 1;
}


sub _create_qual {
   my $self = shift;
   $self->qual([ ($self->qual_levels->[0]) x ($self->end - $self->start + 1) ]);
   return 1;
}


sub _create_desc {
   # Create the read description of the error-free read
   my $self = shift;
   # Reference sequence ID
   my $desc_str = '';
   my $ref_id = $self->reference->id;
   if (defined $ref_id) {
      $desc_str .= 'reference='.$ref_id.' ';
   }
   # Position of read on reference sequence: start, end and strand
   my $strand = $self->strand;
   if ($self->coord_style eq 'bioperl') {
     $desc_str .= 'start='.$self->start.' end='.$self->end.' ';
     if (defined $strand) {
        # Strand of the reference sequence that the read is on
        $strand = '+1' if $strand == 1;
        $desc_str .= 'strand='.$strand.' ';
     }
   } else {
     if ( (defined $strand) && ($strand == -1) ) {
       # Reverse complemented
       $desc_str .= 'position=complement('.$self->start.'..'.$self->end.') ';
     } else {
       # Regular (forward) orientation
       $desc_str .= 'position='.$self->start.'..'.$self->end.' ';
     }
   }
   # Description of the original sequence
   my $ref_desc = $self->reference->desc;
   if ( (defined $self->reference->desc) && ($self->reference->desc !~ m/^\s*$/) ) {
      $ref_desc =~ s/"/\\"/g; # escape double-quotes to \"
      $desc_str .= 'description="'.$ref_desc.'" ';
   }
   $desc_str =~ s/\s$//g;
   # Record new description
   $self->desc($desc_str);
   return 1;
}


=head2 mid

 Title    : mid
 Function : Get or set a multiplex identifier (or MID, or tag, or barcode) to
            add to the read. By default, no MID is used. This method must be
            called after reference() but before errors().
 Usage    : my $mid = read->mid();
 Arguments: MID sequence string (e.g. 'ACGT')
 Returns  : MID sequence string

=cut

sub mid {
   my ($self, $mid) = @_;
   if (defined $mid) {
      # Sanity check 1
      if (not defined $self->reference) {
         $self->throw("Cannot add MID because the reference sequence was not ".
            "set\n");
      }
      # Sanity check 2
      if ($self->{errors}) {
         $self->throw("Cannot add an MID after sequencing errors have been ".
            "introduced in the read\n");
      }
      # Sanity check 3
      if (not $self->validate_seq($mid)) {
         $self->throw("MID is not a valid DNA sequence\n");
      }
      # Update sequence, quality scores and description with the MID
      $self->_update_seq_mid($mid);
      $self->_update_qual_mid($mid) if scalar @{$self->qual_levels};
      $self->_update_desc_mid($mid) if $self->track;
      # Set the MID value
      $self->{mid} = $mid;
   }
   return $self->{mid}
}


sub _update_seq_mid {
   # Update the MID of a sequence
   my ($self, $mid) = @_;
   # Remove old MID
   my $seq = $self->seq;
   my $old_mid = $self->{mid};
   if (defined $old_mid) {
      $seq =~ s/^$old_mid//;
   }
   # Add new MID
   $seq = $mid . $seq;
   $self->seq( $seq );
   return 1;
}


sub _update_qual_mid {
   # Update the MID of a quality scores
   my ($self, $mid) = @_;
   # Remove old MID
   my $qual = $self->qual;
   my $old_mid = $self->{mid};
   if (defined $old_mid) {
      splice @$qual, 0, length($old_mid);
   }
   $qual =  [($self->qual_levels->[0]) x length($mid), @$qual];
   $self->qual( $qual );
   return 1; 
}


sub _update_desc_mid {
   # Update MID specifications in the read description
   my ($self, $mid) = @_;
   if ($mid) {
      # Sequencing errors introduced in the read
      my $mid_str = "mid=".$mid;
      my $desc_str = $self->desc;
      $desc_str =~ s/((position|strand)=\S+)( mid=\S+)?/$1 $mid_str/g;
      $self->desc( $desc_str );
   }
   return 1;
}


=head2 errors

 Title    : errors
 Function : Get or set the sequencing errors and update the read. By default, no
            errors are made. This method must be called after the mid() method.
 Usage    : my $errors = $read->errors();
 Arguments: Reference to a hash of the position and nature of sequencing errors.
            The positions are 1-based relative to the error-free MID-containing
            read (not relative to the reference sequence). For example:
               $errors->{34}->{'%'} = 'T'  ; # substitution of residue 34 by a T
               $errors->{23}->{'+'} = 'GG' ; # insertion of GG after residue 23
               $errors->{45}->{'-'} = undef; # deletion of residue 45
            Substitutions and deletions are for a single residue, but additions
            can be additions of several residues.
            An alternative way to specify errors is by using array references
            instead of scalar for the hash values. This allows to specify
            redundant mutations. For example, the case presented above would
            result in the same read sequence as the example below:
               $errors->{34}->{'%'} = ['C', 'T']    ; # substitution by a C and then a T
               $errors->{23}->{'+'} = ['G', 'G']    ; # insertion of G and then a G
               $errors->{45}->{'-'} = [undef, undef]; # deletion of residue, and again
 Returns  : Reference to a hash of the position and nature of sequencing errors.

=cut

sub errors {
   my ($self, $errors) = @_;
   if (defined $errors) {
      # Verify that we have a hashref
      if ( (not defined ref $errors) || (not ref $errors eq 'HASH') ) {
         $self->throw("Error specification has to be a hashref. Got: $errors\n");
      }
      # Verify that we have a reference sequence
      if (not defined $self->reference) {
         $self->throw("Cannot add errors because the reference sequence was not set\n");
      }
      # Convert scalar error specs to arrayref specs
      $errors = $self->_scalar_to_arrayref($errors);
      # Check validity of error specifications
      $errors = $self->_validate_error_specs($errors);
      # Set the error specifications
      $self->{errors} = $errors;
      # Need to recalculate the read from the reference if previously mutated
      if ($self->{_mutated}) {
         $self->_create_seq;
         $self->_create_qual if scalar @{$self->qual_levels};
         $self->_create_desc if $self->track;
      }
      # Now mutate the read, quality score and description
      $self->_update_seq_errors;
      $self->_update_qual_errors if scalar @{$self->qual_levels};
      $self->_update_desc_errors if $self->track;

   }
   return $self->{errors};
}


sub _scalar_to_arrayref {
   # Replace the scalar values in the error specs by more versatile arrayrefs
   my ($self, $errors) = @_;
   while ( my ($pos, $ops) = each %$errors ) {
      while ( my ($op, $res) = each %$ops ) {
         if (ref $res eq '') {
            my $arr = [ split //, ($res || '') ];
            $arr = [undef] if scalar @$arr == 0;
            $$errors{$pos}{$op} = $arr;
         }
      }
   }
   return $errors;
}


sub _validate_error_specs {
   # Clean error specifications and warn of any issues encountered
   my ($self, $errors) = @_;
   my %valid_ops = ('%' => undef, '-' => undef, '+' => undef); # valid operations

   # Calculate read length
   my $read_length = $self->length;
   while ( my ($pos, $ops) = each %$errors ) {

      # Position cannot be no longer than the read length
      if ( (defined $read_length) && ($pos > $read_length) ) {
         $self->warn("Position $pos is beyond end of read ($read_length ".
            "residues). Skipping errors specified at this position.\n");
         delete $errors->{$pos};
      }

      # Position has to be 0+ for addition, 1+ for substitution and deletion
      if ( $pos < 1 && (exists $ops->{'%'} || exists $ops->{'-'}) ) {
         $self->warn("Positions of substitutions and deletions have to be ".
            "strictly positive but got $pos. Skipping substitution or deletion".
            " at this position\n");
         delete $ops->{'%'};
         delete $ops->{'-'};
      }
      if ( $pos < 0 && exists $ops->{'+'}) {
         $self->warn("Positions of additions have to be zero or more. ".
            "Skipping addition at position $pos.\n");
         delete $ops->{'+'};
      }

      # Valid operations are '%', '+' and '-'
      while ( my ($op, $res) = each %$ops ) {
         if (not exists $valid_ops{$op}) {
            $self->warn("Skipping unknown error operation '$op' at position".
               " $pos\n");
            delete $ops->{$op};
         } else {
            # Substitutions: have to have at least one residue to substitute
            if ( ($op eq '%') && (scalar @$res < 1) ) {
               $self->warn("At least one residue must be provided for substitutions,".
                  "but got ".scalar(@$res)." at position $pos.\n");
            }
            # Additions: have to have at least one residue to add
            if ( ($op eq '+') && (scalar @$res < 1) ) {
               $self->warn("At least one residue must be provided for additions,".
                  "but got ".scalar(@$res)." at position $pos.\n");
            }
            # Deletions
            if ( ($op eq '-') && (scalar @$res < 1) ) {
               $self->warn("At least one 'undef' must be provided for deletions,".
                  "but got ".scalar(@$res)." at position $pos.\n");
            }
         }
      }

      delete $errors->{$pos} unless scalar keys %$ops;
   }

   return $errors;
}


sub _update_seq_errors {
   my $self = shift;
   my $seq_str = $self->seq;
   my $errors = $self->errors;
   if (scalar keys %$errors > 0) {
      my $off = 0;
      for my $pos ( sort {$a <=> $b} (keys %$errors) ) {
         # Process sequencing errors at that position
         for my $type ( '%', '-', '+' ) {
            next if not exists $$errors{$pos}{$type};
            my $arr = $$errors{$pos}{$type};
            if ($type eq '%') {
               # Substitution at residue position. If there are multiple
               # substitutions to do, directly skip to the last one.
               substr $seq_str, $pos - 1 + $off, 1, $$arr[-1];
            } elsif ($type eq '-') {
               # Deletion at residue position
               substr $seq_str, $pos - 1 + $off, 1, '';
               $off--;
            } elsif ($type eq '+') {
               # Insertion after residue position
               substr $seq_str, $pos + $off, 0, join('', @$arr);
               $off += scalar @$arr;
            }
         }
      }
      $self->{_mutated} = 1;
   } else {
      $self->{_mutated} = 0;
   }
   $self->seq($seq_str);
   return 1;
}


sub _update_qual_errors {
   my $self = shift;
   my $qual = $self->qual;
   my $errors = $self->errors;
   my $bad_qual = $self->qual_levels->[1];
   if (scalar keys %$errors > 0) {
      my $off = 0;
      for my $pos ( sort {$a <=> $b} (keys %$errors) ) {
         # Process sequencing errors at that position
         for my $type ( '%', '-', '+' ) {
            next if not exists $$errors{$pos}{$type};
            my $arr = $$errors{$pos}{$type};
            if ($type eq '%') {
               # Substitution at residue position
               splice @$qual, $pos - 1 + $off, 1, $bad_qual;
            } elsif ($type eq '-') {
               # Deletion at residue position
               splice @$qual, $pos - 1 + $off, 1;
               $off--;
            } elsif ($type eq '+') {
               # Insertion after residue position
               splice @$qual, $pos + $off, 0, ($bad_qual) x scalar(@$arr);
               $off += scalar @$arr;
            }
         }
      }
   }
   $self->qual($qual);
   return 1;
}


sub _update_desc_errors {
   # Add or update error specifications in the read description
   my $self = shift;
   my $errors = $self->errors;
   if (defined $errors and scalar keys %$errors > 0) {
      # Sequencing errors introduced in the read
      my $err_str = 'errors=';
      for my $pos ( sort {$a <=> $b} (keys %$errors) ) {
         # Process sequencing errors at that position
         for my $type ( '%', '-', '+' ) {
            next if not exists $$errors{$pos}{$type};
            for my $val ( @{$$errors{$pos}{$type}} ) {
               $val = '' if not defined $val;
               $err_str .= $pos . $type . $val . ',';
            }
         }
      }
      $err_str =~ s/,$//;
      my $desc_str = $self->desc;
      $desc_str =~ s/((position|strand)=\S+( mid=\S+)?)( errors=\S+)?/$1 $err_str/g;
      $self->desc( $desc_str );
   }
   return 1;
}


=head2 track

 Title    : track
 Function : Get or set the tracking status in the read description. By default,
            tracking is on. This method can be called at any time.
 Usage    : my $track = $read->track();
 Arguments: 1 for tracking, 0 otherwise
 Returns  : 1 for tracking, 0 otherwise

=cut

sub track {
   my ($self, $track) = @_;
   if (defined $track) {
      if (defined $self->reference) {
         if ($track == 1) {
            $self->_create_desc;
            $self->_update_desc_mid($self->mid);
            $self->_update_desc_errors;
         } else {
            $self->desc(undef);
         }
      }
      $self->{track} = $track;
   }
   return $self->{track};
}


=head2 coord_style

 Title    : coord_style
 Function : When tracking is on, define which 1-based coordinate system to use
            in the read description:
              * 'bioperl' uses the start, end and strand keywords (default),
                similarly to the GFF3 format. Example:
                  start=1 end=10 strand=+1
                  start=1 end=10 strand=-1
              * 'genbank' does only provide the position keyword. Example:
                  position=1..10
                  position=complement(1..10)
 Usage    : my $coord_style = $read->track();
 Arguments: 'bioperl' or 'genbank'
 Returns  : 'bioperl' or 'genbank'

=cut

sub coord_style {
   my ($self, $coord_style) = @_;
   my %styles = ( 'bioperl' => undef, 'genbank' => undef );
   if (defined $coord_style) {
      if (not exists $styles{$coord_style}) {
         die "Error: Invalid coordinate style '$coord_style'\n";
      }
      $self->{coord_style} = $coord_style;
   }
   return $self->{coord_style};
}


1;
