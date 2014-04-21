#
# BioPerl module for Bio::Tools::QRNA
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

Bio::Tools::QRNA - A Parser for qrna output

=head1 SYNOPSIS

  use Bio::Tools::QRNA;
  my $parser = Bio::Tools::QRNA->new(-file => $qrnaoutput);
  while( my $feature = $parser->next_feature ) {
    # do something here
  }

=head1 DESCRIPTION

Parses QRNA output (E.Rivas:
http://selab.janelia.org/software.html
ftp://selab.janelia.org/pub/software/qrna/).

This module is not complete, but currently it packs information from
each QRNA alignment into a single Bio::SeqFeature::Generic object.

Not all options for QRNA output have been tested or tried.  It has
been tested on sliding window output (-w -x) and shuffled output (-b
or -B).

See t/QRNA.t for example usage.

At some point we may have more complicated feature object which will
support this data rather than forcing most of the information into
tag/value pairs in a SeqFeature::Generic.

Running with -verbose =E<gt> 1 will store extra data in the feature.  The
entire unparsed entry for a particular feature will be stored as a
string in the tag 'entry' it is accessible via:

  my ($entry) = $f->each_tag_value('entry');

The winning model for any given alignment test will be the name stored
in the primary_tag field of feature.  The bit score will stored in the
score field.  The logoddpost is available via the a tag/value pair.
This example code will show how to print out the score and log odds
post for each model.

  # assuming you got a feature already
  print "model score logoddspost\n";
  foreach my $model ( qw(OTH COD RNA) ) {
    my ($score)       = $f->get_tag_values("$model\_score");
    my ($logoddspost) = $f->get_tag_values("$model\_logoddspost");
    print "$model $score $logoddspost\n";
  }

The start and end of the alignment for both the query and hit sequence
are available through the L<Bio::SeqFeature::FeaturePair> interface,
specifically L<Bio::SeqFeature::FeaturePair::feature1> and
L<Bio::SeqFeature::FeaturePair::feature2>.  Additionally if you have
run QRNA with an input file which has the location of the alignment
stored in the FASTA filename as in (ID/START-END) which is the default
output format from L<Bio::AlignIO::fasta> produced alignment output,
this module will re-number start/end for the two sequences so they are
in the actual coordinates of the sequence rather than the relative
coordinates of the alignment.  You may find the bioperl utillity
script search2alnblocks useful in creating your input files for QRNA.

Some other words of warning, QRNA uses a 0 based numbering system for
sequence locations, Bioperl uses a 1 based system.  You'll notice that
locations will be +1 they are reported in the raw QRNA output.

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


package Bio::Tools::QRNA;
use vars qw(@Models);
use strict;

use Bio::SeqFeature::Generic;
use Bio::SeqFeature::FeaturePair;

use base qw(Bio::Root::IO Bio::SeqAnalysisParserI);
@Models = qw(OTH COD RNA);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::QRNA->new();
 Function: Builds a new Bio::Tools::QRNA object 
 Returns : an instance of Bio::Tools::QRNA
 Args    : -fh/-file filehandle/filename standard input for 
                     Bio::Root:IO objects

=cut

=head2 next_feature

 Title   : next_feature
 Usage   : my $feature = $parser->next_feature
 Function: Get the next QRNA feature
 Returns : 
 Args    :


=cut

sub next_feature {
    my ($self) = @_;
    my $f = shift @{$self->{'_parsed_features'} || []};
    if( ! defined $f && $self->_parse_pair ) {
	$f = shift @{$self->{'_parsed_features'} || []};
    }
    return $f;
}

sub _parse_pair {
   my ($self,@args) = @_;
   my (@features,%data);
   my $seenstart = 0;
   while( defined( $_ = $self->_readline) ) {
       next if( /^\#\-\-/o );
       if( /^\#\s+(qrna)\s+(\S+)\s+\(([^\)]+)\)/o ) {
	   $self->program_name($1);
	   $self->program_version($2);
	   $self->program_date($3);
       } elsif( /^\#\s+(PAM model)\s+\=\s+(.+)\s+$/o ) {
	   $self->PAM_model($2);
       } elsif( /^\#\s+(RNA model)\s+\=\s+(\S+)/o ) {
	   $self->RNA_model($2);
       } elsif( /^\#\s+(seq file)\s+\=\s+(.+)\s+$/o ) {
	   $self->seq_file($2);	   
       } elsif( /^\#\s+(\d+)\s+\[([\-+])\s+strand\]/o ) {
	   if( $seenstart ) { 
	       if( $data{'alignment_len'} ) {
		   push @features, $self->_make_feature(\%data);
	       }
	       $self->_pushback($_);
	       last;
	   }
	   $seenstart = 1;
       } elsif( /^\#/ ) {
	   next;
       } elsif( />(\S+)\s+\((\d+)\)/ ) {
	   if( @{$data{'seqs'} || []} == 2 ) { 
	       $self->warn( "already seen seqs ".join(' ', ,map { $_->[0] } 
						      @{$data{'seqs'}}). "\n");
	   } else { 
	       push @{$data{'seqs'}}, [$1,$2];
	   }
       } elsif( /^length alignment:\s+(\d+)\s+\(id\=(\d+(\.\d+)?)\)/o ) {
	   
	   if( $data{'alignment_len'} ) {
	       push @features, $self->_make_feature(\%data);	
	       # reset all the data but the 'seqs' field
	       %data  = ( 'seqs' => $data{'seqs'} );
	   }
	   
	   if( /\(((sre_)?shuffled)\)/ ) { 
	       $data{'shuffled'} = $1;
	   }
	   $data{'alignment_len'} = $1;
	   $data{'alignment_pid'} = $2;
       } elsif ( /^pos([XY]):\s+(\d+)\-(\d+)\s+\[(\d+)\-(\d+)\]\((\d+)\)\s+
		 \-\-\s+\((\S+\s+\S+\s+\S+\s+\S+)\)/ox ) {
	   $data{"seq\_$1"}->{'aln'} = [ $2,$3, $4,$5, $6];
	   @{$data{"seq\_$1"}->{'base_comp'}} = split(/\s+/,$7);
       } elsif( /^winner\s+\=\s+(\S{3})/ ) {
	   $data{'winning_model'} = $1;
       } elsif( /^(\S{3})\s+ends\s+\=\s+(\-?\d+)\s+(\-?\d+)/ ) {
	   # QRNA is 0-based
	   # Bioperl is 1 based
	   $data{'model_location'}->{$1} = [ $2,$3 ];
       }  elsif( /^\s+(logoddspost)?OTH\s+\=\s+/ox ) {
	   while( /(\S+)\s+\=\s+(\-?\d+(\.\d+))/g ) {
	       my ($model,$score)= ($1,$2);
	       if( $model =~ s/^logoddspost// ) {
		   $data{'model_scores'}->{'logoddspost'}->{$model} = $score;
	       } else {
		   $data{'model_scores'}->{'bits'}->{$model} = $score;
	       }
	   }
       }
       $data{'entry'} .= $_;
   }
   if( @features ) {
       push @{$self->{'_parsed_features'}}, @features;
       return scalar @features;
   }
   return 0;
}

=head2 PAM_model

 Title   : PAM_model
 Usage   : $obj->PAM_model($newval)
 Function: 
 Example : 
 Returns : value of PAM_model (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub PAM_model{
    my $self = shift;
    return $self->{'PAM_model'} = shift if @_;
    return $self->{'PAM_model'};
}

=head2 RNA_model

 Title   : RNA_model
 Usage   : $obj->RNA_model($newval)
 Function: 
 Example : 
 Returns : value of RNA_model (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub RNA_model{
    my $self = shift;

    return $self->{'RNA_model'} = shift if @_;
    return $self->{'RNA_model'};
}

=head2 seq_file

 Title   : seq_file
 Usage   : $obj->seq_file($newval)
 Function: 
 Example : 
 Returns : value of seq_file (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub seq_file{
    my $self = shift;

    return $self->{'seq_file'} = shift if @_;
    return $self->{'seq_file'};
}


=head2 program_name

 Title   : program_name
 Usage   : $obj->program_name($newval)
 Function: 
 Example : 
 Returns : value of program_name (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub program_name{
    my $self = shift;

    return $self->{'program_name'} = shift if @_;
    return $self->{'program_name'} || 'qrna';
}

=head2 program_version

 Title   : program_version
 Usage   : $obj->program_version($newval)
 Function: 
 Example : 
 Returns : value of program_version (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub program_version{
    my $self = shift;

    return $self->{'program_version'} = shift if @_;
    return $self->{'program_version'};
}

=head2 program_date

 Title   : program_date
 Usage   : $obj->program_date($newval)
 Function: 
 Example : 
 Returns : value of program_date (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub program_date{
    my $self = shift;
    return $self->{'program_date'} = shift if @_;
    return $self->{'program_date'};
}

sub _make_feature { 
    my ($self,$data) = @_; 
    my ($qoffset,$hoffset) = (1,1);
    # when you run qrna and have produced ID/START-END
    # formatted input strings we can remap the location
    # to the original

    # name is stored as the first entry in the seq array ref
    my ($qid,$hid) = ( $data->{'seqs'}->[0]->[0],
		       $data->{'seqs'}->[1]->[0]);
    if( $qid =~ /(\S+)\/(\d+)\-(\d+)/ ) {
	($qid,$qoffset) = ($1,$2);
    }
    if( $hid =~ /(\S+)\/(\d+)\-(\d+)/ ) {
	($hid,$hoffset) = ($1,$2);
    }

    my $f = Bio::SeqFeature::FeaturePair->new();

    my ($s,$e) = @{$data->{'model_location'}->{$data->{'winning_model'}}};
    my $qf = Bio::SeqFeature::Generic->new
	( -primary_tag => $data->{'winning_model'},
	  -source_tag  => $self->program_name,
	  -score       => $data->{'model_scores'}->{'bits'}->{$data->{'winning_model'}},
	  -start       => $s+$qoffset,
	  -end         => $e+$qoffset,
	  -seq_id      => $qid,
	  -strand      => ($s < $e ) ? 1 : -1,
	  );

    my $hf = Bio::SeqFeature::Generic->new
	( -primary_tag => $qf->primary_tag,
	  -source_tag  => $qf->source_tag,
	  -score       => $qf->score,
	  -seq_id      => $hid,
	  -start       => $s + $hoffset,
	  -end         => $e + $hoffset,
	  -strand      => $qf->strand,
	  );
    $f->feature1($qf);
    $f->feature2($hf);
    $f->add_tag_value('alignment_len', $data->{'alignment_len'});
    $f->add_tag_value('alignment_pid', $data->{'alignment_pid'});
    # store the other model scores and data
    foreach my $model ( @Models ) {
	$f->add_tag_value("$model\_score", $data->{'model_scores'}->{'bits'}->{$model});
	$f->add_tag_value("$model\_logoddspost", $data->{'model_scores'}->{'logoddspost'}->{$model});
	if( ! $data->{'model_location'}->{$model} ) {
	    if( $self->verbose > 0 ) {
		$self->debug( $data->{'entry'} );
	    }
	    $self->throw("no location parsed for $model in ",
	    (map { @$_ } @{$data->{'seqs'}}), " ", $f->start, " ", $f->end);
	} else { 
	    $f->add_tag_value("$model\_positions", 
			      join("..",@{$data->{'model_location'}->{$model} }));
	}
    }
    # probably a better way to store this - as 
    # a seq object perhaps
    $f->add_tag_value('seq1', @{$data->{'seqs'}->[0]});
    $f->add_tag_value('seq2', @{$data->{'seqs'}->[1]});
    $f->add_tag_value('entry', $data->{'entry'}) if $self->verbose > 0;
    if( $data->{'shuffled'} ) {
	$f->add_tag_value('shuffled', $data->{'shuffled'});
    }
    return $f;					       
}
1;
