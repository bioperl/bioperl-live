# You may distribute this module under the same terms as perl itself

=head1 NAME

Bio::Tools::SwissProtParser - Lightweight SwissProt Parser, loosely adapted
from FASTAParser.

=head1 SYNOPSIS

 my $report = Bio::Tools::SwissProtParser->new(\*STDIN) or die;

 while(my $sbjct = $report->nextSbjct) {
     #same thing
     $sbjct->AC;
     $sbjct->accession;

     #same again
     $sbjct->ID;
     $sbjct->identification;

     #you also can call in scalar or list context
     my @gene_names = $sbjct->GN; # the list itself
     my $gene_names = $sbjct->GN; # space joined list
 }

=head1 DESCRIPTION

Provides OO access to SwissProt records.

=head2 Constructor

To create a new report, you pass a filehandle reference to the
SwissProtParser constructor.

 my $report = new Bio::Tools::SwissProtParser(\*STDIN);

=head Accesors

See the Swissprot Manual for the lowdown:
 reference: http://us.expasy.org/sprot/userman.html

Here I provide a mapping of tags to verbose accessor names.  
Useful if you like to type.  All methods are read/write.

 SwissProt Tag		Method	Verbose Method
 ---------------------------------------------
 ID			ID	identification
 AC			AC	accession
 DT			DT	date
 DE			DE	description
 GN			GN	gene_name
 OS			OS	organism_species
 OG			OG	organelle
 OC			OC	organism_classification
 OX			OX	organism_taxonomy_cross_reference
 KW			KW	keyword
 DR			DR	database_cross_reference
 FT			FT	feature_table
 SQ			SQ	sequence_header
 			SS	sequence
 RN			RR	reference
 RP			RR	reference
 RC			RR	reference
 RX			RR	reference
 RA			RR	reference
 RT			RR	reference
 RL			RR	reference
 CC			CC	comment

=head1 AUTHOR

Allen Day <allenday@ucla.edu>

=head1 SEE ALSO

L<Bio::Tools::FASTAParser>

=head1 COPYRIGHT

Copyright (c) 2002 Allen Day

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut

package Bio::Tools::SwissProtParser;
use strict;
use Bio::Root::Root;
use vars qw(@ISA);
@ISA = qw(Bio::Root::Root);

=head2 new

 Title   : new
 Function: Create a new Bio::Tools::SwissParser object
 Returns : Bio::Tools::SwissParser object
 Args    : Filehandle (ie, \*STDIN)

=cut

sub new {
        my ($class, $fh) = @_;
        if (ref $fh !~ /GLOB/)
                {die "SwissProtParser error: new expects a GLOB reference not $fh\n"}
        my $this = bless {};
        $this->{FH} = $fh;
        return $this;
}

sub FH { return shift->{FH} }

=head2 nextSbjct

 Title    : nextSbjct
 Usage    : $sbjct = $obj->nextSbjct();
 Function : Method of iterating through all the Sbjct retrieved
            from parsing the report
 Example  : while ( my $sbjct = $obj->nextSbjct ) {}
 Returns  : next Sbjct object or 0 if finished
 Args     :

=cut

sub nextSbjct {
	my $self = shift;
	local $/ = "//\n";
	my $fh = $self->FH;
	my $record = <$fh>;
	return 0 unless $record;
        return Bio::Tools::SwissProtParser::Sbjct::new ($self,$record);
}

###############################################################################
# SwissProtParser::Subjct
###############################################################################

package Bio::Tools::SwissProtParser::Sbjct;
use strict;
use Bio::Tools::SwissProtParser;
use vars qw(@ISA);
@ISA = qw (Bio::Tools::SwissProtParser);

sub new {
        my $this = bless {};
	my $parent = shift;
	my $record = shift;
	return 0 unless $record;

	$this->parse($record);
	return $this;
}

sub parse {
	my ($self,$record) = @_;
	$self->throw("record undefined: $!") unless defined $record;

	chomp $record;
	my @lines = split "\n", $record;
	foreach my $line (@lines) {
		my($tag,$content) = $line =~ /^(..)   (.+)$/;

		#reference: http://us.expasy.org/sprot/userman.html
		my $func = 
		  $tag eq 'ID' ? 'ID'		#one
		: $tag eq 'AC' ? 'AC'		#many delimited by ';'
		: $tag eq 'DT' ? 'DT'		#one
		: $tag eq 'DE' ? 'DE'		#many
		: $tag eq 'GN' ? 'GN'		#many line delited by 'AND|OR'
		: $tag eq 'OS' ? 'OS'		#one
		: $tag eq 'OG' ? 'OG'		#many
		: $tag eq 'OC' ? 'OC'		#many delimited by ';'
		: $tag eq 'OX' ? 'OX'		#many delimited by ';'
		: $tag eq 'KW' ? 'KW'		#many delimited by ';'
		: $tag eq 'DR' ? 'DR'		#many; these are special
		: $tag eq 'FT' ? 'FT'		#these are special
		: $tag eq 'SQ' ? 'SQ'		#one
		: $tag eq '  ' ? 'SS'		#i made SS up for Sequence Sequence
						#i made RR up for Reference Reference
		: $tag eq 'RN' ? 'RR'		#reference_number
		: $tag eq 'RP' ? 'RR'		#reference_position
		: $tag eq 'RC' ? 'RR'		#reference_comment
		: $tag eq 'RX' ? 'RR'		#reference_cross_reference
		: $tag eq 'RA' ? 'RR'		#reference_author
		: $tag eq 'RT' ? 'RR'		#reference_title
		: $tag eq 'RL' ? 'RR'		#reference_location
		: $tag eq 'CC' ? 'CC'		#many, multiline
		: $tag eq '//' ? '//'
		: undef;

		next if $func eq '__';
		warn "$line\nunknown tag: '$tag'" and next unless defined $func;

		$self->$func($tag,$content);
	}
}

#this is the piece that figures out what to return
#based on calling context
sub omatic {
	my($self,$tag) = @_;
	$tag = lc $tag;

	return wantarray 
		? defined $self->{$tag}
			? @{$self->{$tag}}
			: []
		: defined $self->{$tag}
			? join " ",@{$self->{$tag}}
			: undef;
}

sub identification { return shift->ID(@_) }
sub ID {
	my($self,$tag,$datum)=@_;

	return $self->{id} unless defined $tag;

	$datum =~ s!^(\S+).+!$1!;
	$self->{id} = $datum;
}

sub accession { return shift->AC(@_) }
sub AC {
	my($self,$tag,$datum)=@_;
	return $self->omatic($tag) unless defined $tag ;
	my @data = split ';',$datum;
	push @{$self->{ac}}, @data;
}

sub date { return shift->DT(@_) }
sub DT {
	my($self,$tag,$datum)=@_;
	return $self->{dt} unless defined $tag;
	$self->{dt} = $datum;
}

sub description { return shift->DE(@_) }
sub DE {
	my($self,$tag,$datum)=@_;
	return $self->omatic($tag) unless defined $tag ;
	push @{$self->{de}}, $datum;
}

sub gene_name { return shift->GN(@_) }
sub GN {
	my($self,$tag,$datum)=@_;
	return $self->omatic($tag) unless defined $tag ;
	my @data = split /AND|OR/ ,$datum;
	push @{$self->{gn}}, @data;
}

sub organism_species { return shift->OS(@_) }
sub OS {
	my($self,$tag,$datum)=@_;
	return $self->omatic($tag) unless defined $tag ;
	$self->{os} = $datum;
}

sub organelle { return shift->OG(@_) }
sub OG {
	my($self,$tag,$datum)=@_;
	return $self->omatic($tag) unless defined $tag ;
	push @{$self->{og}}, $datum;
}

sub organism_classification { return shift->OC(@_) }
sub OC {
	my($self,$tag,$datum)=@_;
	return $self->omatic($tag) unless defined $tag ;
	my @data = split ';',$datum;
	push @{$self->{oc}}, @data;
}

sub organism_taxonomy_cross_reference { return shift->OX(@_) }
sub OX {
	my($self,$tag,$datum)=@_;
	return $self->omatic($tag) unless defined $tag ;
	my @data = split ';',$datum;
	push @{$self->{ox}}, @data;
}

sub keyword { return shift->KW(@_) }
sub KW {
	my($self,$tag,$datum)=@_;
	return $self->omatic($tag) unless defined $tag ;
	my @data = split ';',$datum;
	push @{$self->{kw}}, @data;
}

sub database_cross_reference { return shift->DR(@_) }
sub DR {
	my($self,$tag,$datum)=@_;
	return $self->omatic($tag) unless defined $tag ;
	my @data = split ';',$datum;
	push @{$self->{dr}}, [@data];
}

sub feature_table { return shift->FT(@_) }
sub FT {
	my($self,$tag,$datum)=@_;
	return $self->omatic($tag) unless defined $tag ;
	my @data = $datum =~ /^(.{8}) (.{6}) (.{6})       (.+)$/;
	chomp $_ foreach @data;
	push @{$self->{ft}}, [@data];
}

sub sequence_header { return shift->SQ(@_) }
sub SQ {
	my($self,$tag,$datum)=@_;
	return $self->omatic($tag) unless defined $tag ;
	my @data = split ';',$datum;
	push @{$self->{sq}}, [@data];
}

sub sequence { return shift->SS(@_) }
sub SS {
	my($self,$tag,$datum)=@_;
	return $self->{ss} unless defined $tag;
	$datum =~ s/\s//gs;
	$self->{ss} .= $datum;
}

sub reference { return shift->RR(@_) }
sub RR {
	my($self,$tag,$datum)=@_;
	$self->{rr} ||= ();

	unless(defined $tag){
		return wantarray
			? $self->{rr}
			: $self->prettyRR($self->{rr});
	}

	$self->refcount(1) if $tag eq 'RN';
	%{$self->{rr}[$self->refcount]}->{$tag} = $datum;
}

sub comment { return shift->CC(@_) }
sub CC {
	my($self,$tag,$datum)=@_;
	$self->{rr} ||= ();

	unless(defined $tag){
		return wantarray
			? $self->{cc}
			: $self->prettyCC($self->{cc});
	}


	if($datum =~ /-\!- (.+?):(.+)/){
		push @{$self->{cc}}, "$1\t$2";
	}
}

sub prettyCC {
	my $self = shift;
	my $cc   = shift;

	my $return = '';

	foreach my $item (sort @$cc) {
		$return .= $item."\n";
	}

	return $return;
}

sub prettyRR {
	my $self = shift;
	my @refs = @_;

	my $return = '';

	foreach my $ref (@refs){
		foreach my $element (@$ref){  
			foreach my $line (keys %$element){
				$return .= $line."\t".$element->{$line}."\n";
			}
		}
	}

	return $return;
}

sub refcount {
	my($self,$arg) = @_;
	return $self->{refcount} unless defined $arg;
	if(defined $self->{refcount}){
		$self->{refcount}++;
		return $self->{refcount};
	} else {
		$self->{refcount} = 0;
		return $self->{refcount};
	}
}

sub comholder {
	my($self,$arg) = @_;
	return $self->{comholder} unless defined $arg;
	$self->{comholder} = $arg;
}

1;
