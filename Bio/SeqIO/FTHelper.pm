#
# BioPerl module for Bio::SeqIO::FTHelper
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::FTHelper - Helper class for Embl/Genbank feature tables

=head1 SYNOPSIS

Used by Bio::SeqIO::EMBL to help process the Feature Table

=head1 DESCRIPTION

Represents one particular Feature with the following fields

      key - the key of the feature
      loc - the location string of the feature
      <other fields> - other fields

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

   bioperl-l@bioperl.org             - General discussion
   bioperl-guts-l@bioperl.org        - Automated bug and CVS messages
   http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::FTHelper;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my ($self, @args) = @_;

  my $make = $self->SUPER::_initialize;
  $self->{'_field'} = {};
# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 from_SeqFeature

 Title   : from_SeqFeature
 Usage   : @fthelperlist = Bio::SeqIO::FTHelper::from_SeqFeature($sf,$context_annseq);
 Function: constructor of fthelpers from SeqFeatures
         :
         : The additional annseq argument is to allow the building of FTHelper
         : lines relevant to particular sequences (ie, when features are spread over
         : enteries, knowing how to build this)
 Returns : an array of FThelpers
 Args    : seq features

=cut

=head2 _generic_seqfeature

 Title   : _generic_seqfeature
 Usage   : $fthelper->_generic_seqfeature($annseq, "GenBank")
 Function: processes fthelper into a generic seqfeature
 Returns : TRUE on success and otherwise FALSE
 Args    : Bio::Seq, string indicating the source (GenBank/EMBL/SwissProt)

=cut

sub _generic_seqfeature {
    my ($fth, $annseq, $source) = @_;
    my ($sf);

    # print "Converting from", $fth->key, "\n";

    # set a default if not specified
    if(! defined($source)) {
	$source = "EMBL/GenBank/SwissProt";
    }
    
    $sf = new Bio::SeqFeature::Generic;
    
    # Trap 23.46 type fuzzy locations
    if ($fth->loc =~ /\d+\.\d+/) {
        $fth->warn("unable to parse fuzzy location ('" .
		   $fth->loc . "') ignoring feature (seqid=" .
		   $annseq->id() . ")");
        $sf = undef;
    }
    # Parse compound features
    elsif ( $fth->loc =~ /join/ ) {
	my $strand;
	if ( $fth->loc =~ /complement/ ) {
	    $strand = -1;
	} else {
	    $strand = 1;
	}

	$sf->strand($strand);
	$sf->primary_tag($fth->key . "_span");
	$sf->source_tag($source);
	$sf->has_tag("parent", 1);
	$sf->_parse->{'parent_homogenous'} = 1;

	# we need to make sub features
	my $loc = $fth->loc;
	while ( $loc =~ /(\<?\d+[ \W]{1,3}\>?\d+)/g ) {
	    my $next_loc = $1;
	    my $sub = new Bio::SeqFeature::Generic;
	    $sub->primary_tag($fth->key);
	    $sub->strand($strand);
	    if($fth->_parse_loc($sub, $next_loc)) {
		$sub->source_tag($source);
		$sf->add_sub_SeqFeature($sub,'EXPAND');
	    } else {
		$fth->warn("unable to parse location successfully out of " .
			   $next_loc . ", ignoring feature (seqid=" .
			   $annseq->id() . ")");
                $sf = undef;
	    }
	}

    }
    # Parse simple locations
    else {
	$sf->source_tag($source);
	$sf->primary_tag($fth->key);
	if ( $fth->loc =~ /complement/ ) {
	    $sf->strand(-1);
	} else {
	    $sf->strand(1);
	}
	if(! $fth->_parse_loc($sf, $fth->loc())) {
	    $annseq->warn("unexpected location line [" . $fth->loc() .
			  "] in reading $source, ignoring feature " .
			  $fth->key() . " (seqid=" . $annseq->id() . ")");
	    $sf = undef;
	}
    }
    #print "Adding B4 ", $sf->primary_tag , "\n";

    if(defined($sf)) {
	foreach my $key ( keys %{$fth->field} ){
	    foreach my $value ( @{$fth->field->{$key}} ) {
		$sf->add_tag_value($key,$value);
	    }
	}
	$annseq->add_SeqFeature($sf);
	return 1;
    } else {
	$fth->warn("unable to parse feature " . $fth->key() .
		   " in $source sequence entry (id=" .
		   $annseq->id() . "), ignoring");
	return 0;
    }
}

=head2 _parse_loc

 Title   : _parse_loc
 Usage   : $fthelper->_parse_loc($feature, $loc_string)
 Function: Parses the given location string and sets start() and end() in the
           given feature object appropriately. As side effects, tag values
           (private) for features where only partial locations are given
           may be set.  These tags are prefixed with an underscore, and are
           not included in the FTHelper object generated by the from_SeqFeature
           method.

           Note that this method is private.
 Returns : TRUE on success and otherwise FALSE
 Args    : Bio::SeqFeatureI, string

=cut

sub _parse_loc {
    my ($self, $sf, $loc) = @_;
    my ($start,$end,$fea_type,$tagval);
    my %compl_of = ("5" => "3", "3" => "5");
    my $fwd = ($sf->strand() > -1) ? "5" : "3";

    #print "Processing $loc\n";

    # Two numbers separated by anything of '.', '^', and spaces (SRS puts a
    # space between the two dots), optionally surrounded by parentheses and a
    # qualifier.  (Qualifiers in locations are theoretically illegal in the
    # EMBL/GenBank feature table, but have been observed in a file from
    # GenBank.)
    # The numbers may optionally be preceded by '<' (first) or '>' second), or
    # replaced by '?' (sometimes in SwissProt), in which case you'll find a
    # tag called '_part_feature' with the value of which end is missing.
    # After the numbers there may be text, separated from the numbers by any
    # of [,;" ], which will be recorded as the value of the tag with the same
    # as the qualifier, or "note" if there is no qualifier.
    # Examples: 10..70
    #           10^11     # you'll find a '_zero_width_feature' tag
    #           <10..>70
    #           ?..70
    #           10. .70   # I've seen SRS doing such a thing
    #           replace(10..12, "ACT")
    #           ^^^^^^^ ^^  ^^   ^^^
    #         qualifier from/to  note (Example from the evil GenBank file)
    #
    # Fuzzy locations like 200.202 or (200.202)..220 are neither covered
    # correctly by this method nor by the Feature object itself. So, it should
    # *not* be passed to this method.
    #
    # HL 05/16/2000
    if($loc =~ /^\s*(\w+[A-Za-z])?\(?([\<\?\d]\d*)[.\^\s]{1,3}([\>\?\d]\d*)[,;\" ]*([A-Za-z]\w*)?\"?\)?\s*$/) {
	#print "1 = \"$1\", 2 = \"$2\", 3 = \"$3\", 4 = \"$4\"\n";
	$fea_type = $1 if $1;
	$start = $2;
	$end   = $3;
	$tagval = $4 if $4;
    } 
    # like before, but only one number
    elsif($loc =~ /^\s*(\w+[A-Za-z])?\(?([\<\>\?\d]\d*)[,;\" ]*([A-Za-z]\w*)?\"?\)?\s*$/) {
	#print "1 = \"$1\", 2 = \"$2\", 3 = \"$3\"\n";
	$fea_type = $1 if $1;
	$start = $end = $2;
	$tagval = $3 if $3;
    } else {
	#print "didn't match\n";
	return 0;
    }
    if ( $start =~ s/^[\<\>\?]// ) {
	#print "partial feature 5' (loc=$loc; start=$start)\n";
	$sf->add_tag_value('_part_feature', $fwd . '_prime_missing');
    }
    if ( $end =~ s/^[\<\>\?]// ) {
	#print "partial feature 3' (loc=$loc; end=$end)\n";
	$sf->add_tag_value('_part_feature',
			    $compl_of{$fwd} . '_prime_missing');
    }
    $sf->start($start) if(length($start) > 0);
    $sf->end($end) if(length($end) > 0);
    
    # Commented out as this is taken care of by the next block - JGRG
    #if(defined($fea_type) && ($fea_type ne "complement")) {
    #    #print "featype: $fea_type\n";
    #    $sf->add_tag_value('_feature_type', $fea_type);
    #}
    if(defined($tagval)) {
	if(! $fea_type) {
	    $fea_type = "note";
	}
	$sf->add_tag_value($fea_type, $tagval);
    }
    
    if ( $loc =~ /\d+\^\d+/ ) {
	$sf->add_tag_value('_zero_width_feature', 1);
	$sf->strand(0);
	if ($sf->start + 1 != $sf->end ) {
	    # FIXME this is a uncertainty condition, which is not yet retained
	    # in the feature object
	}
    }
    return 1;
}

sub from_SeqFeature {
    my ($sf, $context_annseq) = @_;
    my @ret;
    my $key;

    #
    # If this object knows how to make FThelpers, then let it
    # - this allows us to store *really* weird objects that can write
    # themselves to the EMBL/GenBank...
    #

    if ( $sf->can("to_FTHelper") ) {
	return $sf->to_FTHelper($context_annseq);
    }

    # build something sensible...
    # if the parent homogenous flag is set, build things from the
    # sub level
    my $loc;
    my ($start_mod, $end_mod, $delim) = ( "", "", ".." );
    my $ph_flag;

    if ( $sf->can('_parse') ) {
	$ph_flag = $sf->_parse->{'parent_homogenous'} || 0;
    } else {
	$ph_flag = 0;
    }

    if ( $ph_flag == 1 ) {
	$key = $sf->primary_tag();
	$key =~ s/_span//g;
	$loc = "join(";
	my $touch = 0;
	foreach my $sub ( $sf->sub_SeqFeature() ) {
	    if ( $touch == 1 ) {
		$loc .= ",";
	    } else {
		$touch = 1;
	    }

	    # decide which symbols to use to describe locations
	    #     Supported: .. ^ < >
	    # Not supported: . (fuzzy ranges)

	    if ( $sf->has_tag('_part_feature') ) {

		if ( grep /5_prime_missing/, $sub->each_tag_value('_part_feature') ) {

		    if ( $sub->strand() == 1 ) {
			$start_mod = '<';
		    } elsif ( $sub->strand() == -1 ) {
			$end_mod = '>';
		    }
		}

		if ( grep /3_prime_missing/, $sub->each_tag_value('_part_feature') ) {

		    if ( $sub->strand() == 1 ) {
			$end_mod = '>';
		    } elsif ( $sub->strand() == -1 ) {
			$start_mod = '<';
		    }
		}
	    }

	    if ( $sub->has_tag('_zero_width_feature')) {
		$delim = "^";
	    }

	    $loc .= $start_mod . $sub->start() . $delim . $end_mod . $sub->end();
	}

	$loc .= ")";
	if ( $sf->strand() == -1 ) {
	    $loc = "complement($loc)";
	}
    } else {

	# decide which symbols to use to describe locations
	#     Supported: .. ^ < >
	# Not supported: . (fuzzy ranges)

	if ( $sf->has_tag('_part_feature')) {

	    if ( grep /5_prime_missing/, $sf->each_tag_value('_part_feature') ) {

		if ( $sf->strand() == 1 ) {
		    $start_mod = '<';
		}
		elsif ( $sf->strand() == -1 ) {
		    $end_mod = '>';
		}
	    }

	    if ( grep /3_prime_missing/, $sf->each_tag_value('_part_feature') ) {

		if ( $sf->strand() == 1 ) {
		    $end_mod = '>';
		}
		elsif ( $sf->strand() == -1 ) {
		    $start_mod = '<';
		}
	    }
	}

	if ( $sf->has_tag('_zero_width_feature')) {
	    $delim = "^";
	}

	$loc = $start_mod . $sf->start() . $delim . $end_mod . $sf->end();

	$key = $sf->primary_tag();

	if ( $sf->strand() == -1 ) {
	    $loc = "complement($loc)";
	}
	# going into sub features
	foreach my $sub ( $sf->sub_SeqFeature() ) {
	    my @subfth = &Bio::SeqIO::FTHelper::from_SeqFeature($sub);
	    push(@ret, @subfth);
	}
    }

    my $fth = Bio::SeqIO::FTHelper->new();

    $fth->loc($loc);
    $fth->key($key);
    $fth->field->{'note'} = [];
    #$sf->source_tag && do { push(@{$fth->field->{'note'}},"source=" . $sf->source_tag ); };
    $sf->score && do { push(@{$fth->field->{'note'}},"score=" . $sf->score ); };
    $sf->frame && do { push(@{$fth->field->{'note'}},"frame=" . $sf->frame ); };
    #$sf->strand && do { push(@{$fth->field->{'note'}},"strand=" . $sf->strand ); };

    foreach my $tag ( $sf->all_tags ) {
        # Tags which begin with underscores are considered
        # private, and are therefore not printed
        next if $tag =~ /^_/;
	if ( !defined $fth->field->{$tag} ) {
	    $fth->field->{$tag} = [];
	}
	foreach my $val ( $sf->each_tag_value($tag) ) {
	    push(@{$fth->field->{$tag}},$val);
	}
    }

    push(@ret, $fth);

    unless (@ret) {
	$context_annseq->throw("Problem in processing seqfeature $sf - no fthelpers. Error!");
    }
    foreach my $ft (@ret) {
	if ( !$ft->isa('Bio::SeqIO::FTHelper') ) {
	    $sf->throw("Problem in processing seqfeature $sf - made a $fth!");
	}
    }

    return @ret;

}

=head2 key

 Title   : key
 Usage   : $obj->key($newval)
 Function:
 Example :
 Returns : value of key
 Args    : newvalue (optional)

=cut

sub key {
   my ($obj, $value) = @_;
   if ( defined $value ) {
      $obj->{'key'} = $value;
    }
    return $obj->{'key'};

}

=head2 loc

 Title   : loc
 Usage   : $obj->loc($newval)
 Function:
 Example :
 Returns : value of loc
 Args    : newvalue (optional)

=cut

sub loc {
   my ($obj, $value) = @_;
   if ( defined $value ) {
      $obj->{'loc'} = $value;
    }
    return $obj->{'loc'};

}

=head2 field

 Title   : field
 Usage   :
 Function:
 Example :
 Returns :
 Args    :

=cut

sub field {
   my ($self) = @_;

   return $self->{'_field'};
}

=head2 add_field

 Title   : add_field
 Usage   :
 Function:
 Example :
 Returns :
 Args    :

=cut

sub add_field {
   my ($self, $key, $val) = @_;

   if ( !exists $self->field->{$key} ) {
       $self->field->{$key} = [];
   }
   push( @{$self->field->{$key}} , $val);

}

