# $Id$
#
# BioPerl module for Bio::DB::CUTG
#
# Cared for by Richard Adams (richard.adams@ed.ac.uk)
#
# Copyright Richard Adams
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::CUTG - for access to the Codon usage Database
at http://www.kazusa.or.jp/codon.

=head1 SYNOPSIS

       use Bio::CodonUsage::Table; 
       use Bio::DB::CUTG;

       my $db = Bio::DB::CUTG->new(-sp =>'Pan troglodytes');
       $db->get_web_request();
       $db->write_data(-file => '>OUT');
       $db->get_local_request(-file => 'OUT');

       my $cdtable = $db->next_data;

=head1 DESCRIPTION


This class retrieves and objectifies codon usage tables either from
the web database or from a local file. The idea is that you can 
initially retrieve a CUT from the web database, and write it to file
in a way that can be read in later. The next_data method returns a 
Bio::CodonUsage::Table object. 

For  a web query, two parameters need to be specified: species(sp) and 
genetic code id (gc). The database is searched 
using regular expressions therefore the full latin name 
must be given to specify the organism. If the species name
is ambiguous the first CUT in the list is retrieved. 
Defaults are Homo sapiens and 1(standard genetic code).
If you are retrieving CUTs from organisms using other genetic codes this needs
to be  put in as a parameter. Parameters can be entered in the constructor
or in the get_web_request ()method. Allowable parameters are listed in the
$QUERY_KEYS hash reference variable.
 
=head1 SEE ALSO

L<Bio::Tools::CodonTable>, 
L<Bio::WebAgent>,
L<Bio::CodonUsage::Table>

=head1 FEEDBACK

=head2 Mailing Lists


User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                       - General discussion
  http://bio.perl.org/MailList.html           - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHORS

Richard Adams, Richard.Adams@ed.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...



package Bio::DB::CUTG;
use Bio::WebAgent;
use Bio::Root::IO;
use vars qw($URL @ISA $QUERY_KEYS);

@ISA = qw(Bio::WebAgent);

$QUERY_KEYS = { 
				sp => 'full Latin species name',	
				gc => 'genetic code id'
			 };

BEGIN {
		 $URL = "http://www.kazusa.or.jp"
	}


=head2 new

 Title   : new
 Usage   : my $db = Bio::DB::CUTG->new()
 Returns : a reference to a new Bio::DB::CUTG 
 Args    : hash of optional values for db query

=cut

sub new {
	my ($class, @args ) =@_;
	_check_args(@args);
	my $self = $class->SUPER::new(@args);
	return $self;
}
=head2 query_keys

 Title   : query_keys
 Usage   : $db->query_keys()
 Purpose : To determine valid keys for parameters for db query.
 Returns : a reference to a hash describing valid query keys
 Args    : none

=cut

sub query_keys {
	return $QUERY_KEYS;
	}

=head2  sp

 Title  : sp
 Usage  : my $sp = $db->sp();
 Purpose: Get/set method for species name
 Returns: void or species name string
 Args   : None or species name string

=cut

sub sp {
	my $self = shift;
	if (@_) {
		my $name = shift;
		if ($name =~ /[^\w\s]/) {
			$self->warn (" contains non-word characters, setting to default
							of Homo sapiens");
			$self->{'_sp'} = "Homo sapiens";
				}
		else{
			$self->{'_sp'} = $name;
			}
		}
	return $self->{'_sp'}|| "Homo sapiens";
	
}

=head2  gc

 Title  : gc
 Usage  : my $gc = $db->gc();
 Purpose: Get/set method for genetic code id
 Returns: void or genetic code  integer
 Args   : None or genetic code integer

=cut

sub gc {
	#### genetic code id for translations ####
	my $self = shift;
	if (@_) {
		if($_[0] =~ /^\d+$/ && $_[0] >= 1 && $_[0] <=15 && $_[0] != 7 
				&& $_[0] != 8) {
			$self->{'_gc'} = shift;
			}
		else {
			$self->warn("invalid genetic code index - setting to default");
			$self->{'_gc'} = 1;
			}
		}
	return $self->{'_gc'} || 1; #return 1 if not defined

	}

=head2  get_local_request

 Title  : get_local_request
 Usage  : $db->get_local_request();
 Purpose: To read in a formatted file of a CUT
 Returns: void 
 Args   : File name(mandatory), genetic code id(optional)

=cut

sub get_local_request {

	my ($self, @args) = @_;
	$self->{'_source'} = 'file';
	$self->{'_raw_cud'} = ''; #reinitialise
	my $ct = Bio::Root::IO->new(@args);
	while (my $l = $ct->_readline() ) {
		$self->{'_raw_cud'}.=  $l;
	}
	
}

=head2  write_data

 Title       : write_data
 Usage       : $db->write_data(-file =>">filename");
 Purpose     : To write  a formatted file of a CUT
 Returns     : void 
 Args        : File name(mandatory)
 Description : At present just writes raw DB entry to file. NO functionality
               for directly writing Bio::CodonUsage::Table objects yet.

=cut

sub write_data {

	##currently just writes text from web retrieval to a text file.
	my ($self, @args)= @_;
	my $ct = Bio::Root::IO->new(@args);
	if (defined ($self->{'_raw_cud'} )) {
		my $title = "#Codon usage table for " . $self->sp . "\n";
		my ($content) = $self->{'_raw_cud'}=~/.*<PRE>(.+)<\/PRE>/s ;
		$content 			||= $self->{'_raw_cud'};
		$ct->_print("$title" . "$content");
		return;
		}
	else {
		$self->warn(" cannot yet write object to file, need raw text file");
		return;
		}
}

=head2  get_web_request

 Title  : get_web_request
 Usage  : $db->get_web_request();
 Purpose: To query remote CUT with a species name
 Returns: void 
 Args   : species  name(mandatory), genetic code id(optional)

=cut
	

sub get_web_request {
	my ($self, @args) = @_;
	_check_args(@args);
	shift;
	### can out in parameters here as well
	while( @_ ) {
	my $key = shift;
        $key =~ s/^-//;
        $self->$key(shift);
    }	
	$self->url($URL);
	$self->{'_source'} = 'www';

	###1st of all search DB to check species exists and is unique
	my $nameparts =  join "+", $self->sp =~ /(\w+)/g;
	my $search_url = $self->url . "/codon/cgi-bin/spsearch.cgi?species=" 
					. $nameparts . "&c=s";
	my $rq = HTTP::Request->new(GET=>$search_url);
	my $reply = $self->request($rq);

	my $content = $reply->content;
	#####  if no matches, assign defaults - or can throw here?  ######
	if ($content =~ /not found/i) {
		$self->warn ("organism not found -selecting human as default");
		$self->sp("Homo sapiens");
		$self->_db("gbpri");
	
	}

	
	else {
		 my @names = $content =~ /(species)/g;
		### get 1st species data from report ####
		my ($sp, $db)  = $content =~ /species=(.*)\+\[(\w+)\]"/;
		
		$sp =~ s/\+/ /g;
		## warn if  more than 1 matching species ##
		## if multiple species retrieved, choose first one by default ##
		if (@names >1 ){
			$self->warn ("too many species - not a unique species id - selecting $sp  ");
			}
		### now assign species and database value
		$self->sp($sp);
		$self->_db($db);
		}


	######## now get codon table , all defaults established now
	 $nameparts =  join "+", $self->sp =~ /(\w+)/g;
	my $CT_url = $self->url . "/codon/cgi-bin/showcodon.cgi?species="
				. $nameparts . "+%5B" . $self->_db . "%5D&aa=" . $self->gc . "&style=GCG";
	my $rq2 = HTTP::Request->new(GET=>$CT_url);
	my $reply2 = $self->request($rq2);

	my $content2 = $reply2->content;
	$self->{'_raw_cud'} = ''; #reinitialise
	($self->{'_raw_cud'}) = $content2 =~ /\[.+\]:(.*)Genetic/s;
	$self->{'_raw_cud'}=~ s/End/Ter/g;
	  
	}

=head2  next_data

 Title  : next_data
 Usage  : my $cut = $db->next_data();
 Purpose: To  obtain a Bio:CodonUsage::Table object from a returned query
 Returns: A  Bio:CodonUsage::Table object
 Args   : none

=cut

sub next_data {
	my $self = shift;
	my $cut = $self->_parse;
	return $cut;
	}


sub _check_args {

	###checks parameters for matching $QUERYKEYS
	my @args = @_;
	while (my $key = lc(shift @args)) {
		$key =~ s/\-//;
		
		if (!exists ($QUERY_KEYS->{$key})) {
			Bio::Root::Root->throw("invalid parameter - must be one of [" .
						(join "] [", keys %$QUERY_KEYS) . "]");
		}
		shift @args;
	}
}

#### internal URL parameter not specifiable ######
sub _db {
	my $self = shift;
	if (@_) {
		$self->{'_db'} = shift;
		}
	return $self->{'_db'};
}

sub _parse {
	my $self = shift;
	my $cdtableobj = Bio::CodonUsage::Table->new();
	#### get stats from outside main table, just for web retrieved version
	if ($self->{'_source'} eq 'www') {
		my ($gcstring) = $self->{'_raw_cud'} =~ /(Coding.+)Genetic/s;
		($cdtableobj->{'_coding_gc'}{'all'},
		$cdtableobj->{'_coding_gc'}{'1'},
		$cdtableobj->{'_coding_gc'}{'2'},
		$cdtableobj->{'_coding_gc'}{'3'}) = $gcstring =~ /(\d\d\.\d\d)/sg;

		($cdtableobj->{'_cds_count'}) = $self->{'_raw_cud'} =~ /\s(\d+)\sCDS/;

	### now just keep main table ########
		$self->{'_raw_cud'} =~ s/.*<PRE>(.+)<\/PRE>.*/$1/s;
		}
	elsif ($self->{'_source'} eq 'file') {
		my $sp = $self->{'_raw_cud'} =~ /^#Codon.+for\s(.+)$/;
		$self->sp($sp);
		}
															
		
	## now parse main table #####
	my @lines = split "\n", $self->{'_raw_cud'};
	for my $line (@lines) {
		## ignore  spurious lines ##
		next if $line =~ /^$/ || $line !~ /^\w\w\w\s/;
		 $line =~ /(\w\w\w)\s+(\w+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/;
		if (defined ($1)) {
			$cdtableobj->{'_table'}{$1}{$2} = {
						'abs_count'=>$3,
						 'per1000'=> $4,
						 'rel_freq'=> $5,
						};
			}
		else {
			print "line is $line";
			}
		}
		## check has been parsed ok ##
		if (scalar keys %{$cdtableobj->{'_table'}} != 21) {
			$cdtableobj->warn("probable parsing error - should be 21 entries for 20aa + stop codon");
		}
		return $cdtableobj;
		
}


	return 1;

