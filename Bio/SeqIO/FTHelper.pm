# $Id$
#
# BioPerl module for Bio::SeqIO::FTHelper
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
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

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://www.bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 CONTRIBUTORS

Jason Stajich jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqIO::FTHelper;
use vars qw(@ISA);
use strict;

use Bio::SeqFeature::Generic;
use Bio::Location::Simple;
use Bio::Location::Fuzzy;
use Bio::Location::Split;


use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);

sub new {
    my ($class, @args) = @_;

    # no chained new because we make lots and lots of these. 
    my $self = {};
    bless $self,$class;
    $self->{'_field'} = {};
    return $self; 
}

=head2 _generic_seqfeature

 Title   : _generic_seqfeature
 Usage   : $fthelper->_generic_seqfeature($annseq, "GenBank")
 Function: processes fthelper into a generic seqfeature
 Returns : TRUE on success and otherwise FALSE
 Args    : Bio::Seq, string indicating the source (GenBank/EMBL/SwissProt)


=cut

sub _generic_seqfeature {
    my ($fth, $seqid, $source) = @_;
    my ($sf);


    # print "Converting from", $fth->key, "\n";

    # set a default if not specified
    if(! defined($source)) {
	$source = "EMBL/GenBank/SwissProt";
    }

    $sf = Bio::SeqFeature::Generic->direct_new();

    my $strand = ( $fth->loc =~ /complement/ ) ? -1 : 1;    
    $sf->strand($strand);

    # Parse compound features
    if ( $fth->loc =~ /(join)/i || $fth->loc =~ /(order)/i  || 
	 $fth->loc =~ /(bond)/i ) {

	my $combotype=$1;
	$sf->primary_tag($fth->key);
	$sf->source_tag($source);
	my $splitlocation = new Bio::Location::Split(-strand=>$strand, 
						     -seq_id => $seqid,
						     -splittype => $combotype);
	# we need to make sub features
	my $loc = $fth->loc;
	$loc =~ s/^.*$combotype\((\S+)\)/$1/;
	foreach my $next_loc ( split(/\s*,\s*/, $loc) ) {
	    my $remote=0;
	    if ( $next_loc =~ s/\(?\s*([A-Za-z\d\_]+(\.\d+)?):// ) {
		$seqid = $1;
		$remote=1;
	    }
	    if( my $location = $fth->_parse_loc($sf,$next_loc)) {
		$fth->debug( "I got remote: ". 
			      join(",", ($location->start(), 
					 $location->end(), 
					 $location->strand()))
			      . " for $next_loc\n");
		$location->seq_id($seqid);
		if ($remote) {
		    $location->is_remote(1);
		}
		$splitlocation->add_sub_Location($location);
	    } else {
		$fth->warn("unable to parse location successfully out of " .
			   $next_loc . ", ignoring feature (seqid=" .
			   $seqid . ")");		
                $sf = undef;
		last;
	    }
	}
	# see bug #930
        # we'll skip this SeqFeature if we can't parse the location 
	$sf->location($splitlocation) if( defined $sf);


    }     
    # Parse simple locations and fuzzy locations
    else {

      $sf->source_tag($source);
      $sf->primary_tag($fth->key);	
      my $loc = $fth->loc();
      my $seqid;
      if( $loc =~ /^(\d+)\.\.(\d+)$/ ) {
	my $start = $1;
	my $end   = $2;

	# hard core object building to accelerate time
	my $location = {};
	bless $location,'Bio::Location::Simple';
	$location->{'_start'}  = $start;
	$location->{'_end'}    = $end;
	$location->{'_strand'} = 1;

	$sf->location($location);
      } elsif ( $loc =~ /^complement\((\d+)\.\.(\d+)\)$/ ) {
	my $start = $1;
	my $end   = $2;

	# hard core object building to accelerate time
	my $location = {};
	bless $location,'Bio::Location::Simple';
	$location->{'_start'}  = $start;
	$location->{'_end'}    = $end;
	$location->{'_strand'} = -1;
	$sf->location($location);

      } else {
	
	if ( $loc =~ s/\(?\s*([A-Za-z\d\_]+(\.\d+)?):// ) {
	  ($seqid) = $1;
	}
	
	
	if( my $location = $fth->_parse_loc($sf,$loc) ) {
	  $location->seq_id($seqid) if ( $seqid);
	  $sf->location($location);
	} else {
	  $fth->warn("unexpected location line [" . $loc .
		      "] in reading $source, ignoring feature " .
		      $fth->key() . " (seqid=" . $seqid . ")");
	  $sf = undef;
	}

      }

    }
    #print "Adding B4 ", $sf->primary_tag , "\n";

    if(defined($sf)) {

	#print "dogfood location is ", $sf->location->to_FTstring(), "\n";
	foreach my $key ( keys %{$fth->field} ){
	    foreach my $value ( @{$fth->field->{$key}} ) {
		$sf->add_tag_value($key,$value);
	    }
	}
	return $sf;
    } else {
	$fth->warn("unable to parse feature " . $fth->key() .
		   " in $source sequence entry (id=" .
		   $seqid . "), ignoring");
	return ();
    }
}

=head2 _parse_loc

 Title   : _parse_loc
 Usage   : $fthelper->_parse_loc( $loc_string)

 Function: Parses the given location string and returns a location object 
           with start() and end() and strand() set appropriately.
           Note that this method is private.
 Returns : location object or 0 on fail
 Args    : location string

=cut

sub _parse_loc {
    my ($self, $sf,$locstr) = @_;
#my ($start,$end,$fea_type,$tagval)
#    my %compl_of = ("5" => "3", "3" => "5");
    my ($fea_type, $tagval) = ('','');
    my ($strand,$start,$end) = (1);


    $self->debug( "Location parse, processing $locstr\n");

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
    #
    # FuzzyLocation management works now, 
    # however the location strings
    # (5.12)..17 
    # (5.18)..(300.305)
    # will not be parsed by the regex below.  Something to work on
    if( $locstr =~ /complement\((.+)/ ) {
	$locstr = $1;
	$strand = -1;
    }
    my ($delim) = '';
    if($locstr =~ /^\s*(\w+[A-Za-z])?\({0,2}([\<\>\?]?\d*[\<\>\?]?([\.\^]\d+)?)\)?([\.\^\s]{1,3})\(?([\<\>\?]?\d*[\<\>\?]?([\.\^]\d+)?)\){0,2}[,;\" ]*([A-Za-z]\w*)?\"?\)?\s*$/) {
#	print "1 = \"$1\", 2 = \"$2\", 3 = \"$3\", 4 = \"$4\", 5 = \"$5\", 6 = \"$6\", 7 = \"$7\"\n";
	$fea_type = $1 if $1;
	$start = $2;
	$delim = $4;
	$end   = $5;
	$tagval = $7 if $7;
    } 
    # like before, but only one number
    elsif($locstr =~ /^\s*(\w+[A-Za-z])?\(?([\<\>\?]?\d*[\<\>\?]?([\.\^]\d+)?)\)?[,;\" ]*([A-Za-z]\w*)?\"?\)?\s*$/) {
#	print "1 = \"$1\", 2 = \"$2\", 3 = \"$3\"\n";	
	$fea_type = $1 if $1;
	$start = $end = $2;
	$tagval = $4 if $4;
    } else  {
	$self->warn( "$locstr didn't match\n") if( $self->verbose > 0);
	return 0;
    }
    
    my $type = 'Bio::Location::Simple';
    my @args = ('-start'=> $start, '-end' => $end,
		'-strand' => $strand);
    if ( $start =~ /[\>\<\?]/ || 
	 $end    =~ /[\>\<\?]/ || 
	 $delim =~ /^[\.^]$/ )
    {
	unless ($start =~ /^\d+$/  && $end =~ /^\d+$/  && --$end == $start) {
	    $type = 'Bio::Location::Fuzzy';
	}
	push @args, ('-location_type' => $delim); 
    } 
    my $location = $type->new(@args);
    if(defined($tagval) && $tagval ne '') {
	if(! $fea_type) {
	    $fea_type = "note";
	}
	$sf->add_tag_value($fea_type, $tagval);
    }


    return $location;
}

=head2 from_SeqFeature

 Title   : from_SeqFeature
 Usage   : @fthelperlist = Bio::SeqIO::FTHelper::from_SeqFeature($sf,
						     $context_annseq);
 Function: constructor of fthelpers from SeqFeatures
         :
         : The additional annseq argument is to allow the building of FTHelper
         : lines relevant to particular sequences (ie, when features are spread over
         : enteries, knowing how to build this)
 Returns : an array of FThelpers
 Args    : seq features


=cut

sub from_SeqFeature {
    my ($sf, $context_annseq) = @_;
    my @ret;

    #
    # If this object knows how to make FThelpers, then let it
    # - this allows us to store *really* weird objects that can write
    # themselves to the EMBL/GenBank...
    #

    if ( $sf->can("to_FTHelper") ) {
	return $sf->to_FTHelper($context_annseq);
    }

    my $fth = Bio::SeqIO::FTHelper->new();
    my $key = $sf->primary_tag();
    my $locstr = $sf->location->to_FTstring;
    
    # ES 25/06/01 Commented out this code, Jason to double check
    #The location FT string for all simple subseqfeatures is already 
    #in the Split location FT string

    # going into sub features
    #foreach my $sub ( $sf->sub_SeqFeature() ) {
	#my @subfth = &Bio::SeqIO::FTHelper::from_SeqFeature($sub);
	#push(@ret, @subfth);    
    #}

    $fth->loc($locstr);
    $fth->key($key);
    $fth->field->{'note'} = [];
    #$sf->source_tag && do { push(@{$fth->field->{'note'}},"source=" . $sf->source_tag ); };
    
    ($sf->can('score') && $sf->score) && do { push(@{$fth->field->{'note'}},
						   "score=" . $sf->score ); };
    ($sf->can('frame') && $sf->frame) && do { push(@{$fth->field->{'note'}},
						   "frame=" . $sf->frame ); };
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

1;
