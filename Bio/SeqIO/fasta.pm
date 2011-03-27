# BioPerl module for Bio::SeqIO::fasta
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#          and Lincoln Stein <lstein@cshl.org>
#
# Copyright Ewan Birney & Lincoln Stein
#
# You may distribute this module under the same terms as perl itself
# _history
# October 18, 1999  Largely rewritten by Lincoln Stein

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::fasta - fasta sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::SeqIO class.

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from fasta flat
file databases.

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

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHORS - Ewan Birney & Lincoln Stein

Email: birney@ebi.ac.uk
       lstein@cshl.org

=head1 CONTRIBUTORS

Jason Stajich, jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::fasta;
use vars qw($WIDTH @SEQ_ID_TYPES $DEFAULT_SEQ_ID_TYPE);
use strict;

use Bio::Seq::SeqFactory;
use Bio::Seq::SeqFastaSpeedFactory;

use base qw(Bio::SeqIO);

@SEQ_ID_TYPES = qw(accession accession.version display primary);
$DEFAULT_SEQ_ID_TYPE = 'display';

BEGIN { $WIDTH = 60}

sub _initialize {
  my($self,@args) = @_;
  $self->SUPER::_initialize(@args);  
  my ($width) = $self->_rearrange([qw(WIDTH)], @args);
  $width && $self->width($width);
  unless ( defined $self->sequence_factory ) {
      $self->sequence_factory(Bio::Seq::SeqFastaSpeedFactory->new());
  }
}

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object, or nothing if no more available
 Args    : NONE

=cut

sub next_seq {
    my( $self ) = @_;
    my $seq;
    my $alphabet;
    local $/ = "\n>";
    return unless my $entry = $self->_readline;

    chomp($entry);
    if ($entry =~ m/\A\s*\Z/s)  { # very first one
	return unless $entry = $self->_readline;
	chomp($entry);
    }
    
    # this just checks the initial input; beyond that, due to setting $/ above,
    # the > is part of the record separator and is removed
    $self->throw("The sequence does not appear to be FASTA format ".
        "(lacks a descriptor line '>')") if $. == 1 && $entry !~ /^>/;
    
    $entry =~ s/^>//;
    
    my ($top,$sequence) = split(/\n/,$entry,2);
    defined $sequence && $sequence =~ s/>//g;
#    my ($top,$sequence) = $entry =~ /^>?(.+?)\n+([^>]*)/s
#	or $self->throw("Can't parse fasta entry");

    my ($id,$fulldesc);
    if( $top =~ /^\s*(\S+)\s*(.*)/ ) {
	($id,$fulldesc) = ($1,$2);
    }
    
    if (defined $id && $id eq '') {$id=$fulldesc;} # FIX incase no space 
                                                   # between > and name \AE
    defined $sequence && $sequence =~ tr/ \t\n\r//d;	# Remove whitespace

    # for empty sequences we need to know the mol.type
    $alphabet = $self->alphabet();
    if(defined $sequence && length($sequence) == 0) {
	if(! defined($alphabet)) {
	    # let's default to dna
	    $alphabet = "dna";
	}
    } else {
	# we don't need it really, so disable
	# we want to keep this if SeqIO alphabet was set by user
	# not sure if this could break something
	#$alphabet = undef;
    }

    $seq = $self->sequence_factory->create(
					   -seq         => $sequence,
					   -id          => $id,
					   # Ewan's note - I don't think this healthy
					   # but obviously to taste.
					   #-primary_id  => $id,
					   -desc        => $fulldesc,
					   -alphabet    => $alphabet,
					   -direct      => 1,
					   );




    # if there wasn't one before, set the guessed type
    #unless ( defined $alphabet ) {
	# don't assume that all our seqs are the same as the first one found
	#$self->alphabet($seq->alphabet());
    #}
    return $seq;

}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq(@seq)
 Function: Writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Array of 1 or more Bio::PrimarySeqI objects

=cut

sub write_seq {
   my ($self,@seq) = @_;
   my $width = $self->width;
   foreach my $seq (@seq) {
		$self->throw("Did not provide a valid Bio::PrimarySeqI object") 
		  unless defined $seq && ref($seq) && $seq->isa('Bio::PrimarySeqI');

		my $top;

		# Allow for different ids 
		my $id_type = $self->preferred_id_type;
		if( $id_type =~ /^acc/i ) {
			$top = $seq->accession_number();
			if( $id_type =~ /vers/i ) {
				$top .= "." . $seq->version();
			}
		} elsif($id_type =~ /^displ/i ) { 
			$self->warn("No whitespace allowed in FASTA ID [". $seq->display_id. "]")
			  if defined $seq->display_id && $seq->display_id =~ /\s/;
			$top = $seq->display_id();
			$top = '' unless defined $top;
			$self->warn("No whitespace allowed in FASTA ID [". $top. "]")
			  if defined $top && $top =~ /\s/;
		} elsif($id_type =~ /^pri/i ) {
			$top = $seq->primary_id();
		}

		if ($seq->can('desc') and my $desc = $seq->desc()) {
			$desc =~ s/\n//g;
			$top .= " $desc";
		}
		
		if( $seq->isa('Bio::Seq::LargeSeqI') ) {
		  $self->_print(">$top\n");
		  # for large seqs, don't call seq(), it defeats the
		  # purpose of the largeseq functionality.  instead get
		  # chunks of the seq, $width at a time
		  my $buff_max = 2000;
		  my $buff_size = int($buff_max/$width)*$width; #< buffer is even multiple of widths
		  my $seq_length = $seq->length;
		  my $num_chunks = int($seq_length/$buff_size+1);
		  for( my $c = 0; $c < $num_chunks; $c++ ) {
		    my $buff_end = $buff_size*($c+1);
		    $buff_end = $seq_length if $buff_end > $seq_length;
		    my $buff = $seq->subseq($buff_size*$c+1,$buff_end);
		    if($buff) {
		      $buff =~ s/(.{1,$width})/$1\n/g;
		      $self->_print($buff);
		    } else {
		      $self->_print("\n");
		    }
		  }
		} else {
		  my $str = $seq->seq;
		  if(defined $str && length($str) > 0) {
		    $str =~ s/(.{1,$width})/$1\n/g;
		  } else {
		    $str = "\n";
		  }
		  $self->_print (">",$top,"\n",$str) or return;
		}
   }

   $self->flush if $self->_flush_on_write && defined $self->_fh;
   return 1;
}


=head2 width

 Title   : width
 Usage   : $obj->width($newval)
 Function: Get/Set the line width for FASTA output
 Returns : value of width
 Args    : newvalue (optional)


=cut

sub width{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'width'} = $value;
    }
    return $self->{'width'} || $WIDTH;
}

=head2 preferred_id_type

 Title   : preferred_id_type
 Usage   : $obj->preferred_id_type('accession')
 Function: Get/Set the preferred type of identifier to use in the ">ID" position
           for FASTA output.
 Returns : string, one of values defined in @Bio::SeqIO::fasta::SEQ_ID_TYPES.
           Default = $Bio::SeqIO::fasta::DEFAULT_SEQ_ID_TYPE ('display').
 Args    : string when setting. This must be one of values defined in 
           @Bio::SeqIO::fasta::SEQ_ID_TYPES. Allowable values:
           accession, accession.version, display, primary
 Throws  : fatal exception if the supplied id type is not in @SEQ_ID_TYPES.

=cut

sub preferred_id_type {
    my ($self,$type) = @_;
    if( defined $type ) {
	if( ! grep lc($type) eq $_, @SEQ_ID_TYPES) {
	    $self->throw(-class=>'Bio::Root::BadParameter',
			 -text=>"Invalid ID type \"$type\". Must be one of: @SEQ_ID_TYPES");
	}
	$self->{'_seq_id_type'} = lc($type);
#	print STDERR "Setting preferred_id_type=$type\n";
    }
    $self->{'_seq_id_type'} || $DEFAULT_SEQ_ID_TYPE;
}

1;
