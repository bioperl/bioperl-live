#!perl


# command-line script for HIV sequence queries
# using HIV.pm and HIVQuery.pm

=head1 NAME

bp_hivq.PL - an interactive command-line interface to L<Bio::DB::HIV> and L<Bio::DB::Query::HIVQuery>

=head1 SYNOPSIS

 $ perl bp_hivq.PL
 hivq> query C[subtype] SI[phenotype]
 hivq> prerun
 80 sequences returned
 Query: C[subtype] SI[phenotype]
 hivq> outfile csi.fas
 hivq> run
 Download complete.
 hivq> outfile dsi.fas
 hivq> run D[subtype] SI[phenotype]
 Download complete.
 hivq> count
 25 sequences returned
 Query: D[subtype] SI[phenotype]
 hivq> exit
 $

=head1 DESCRIPTION

The BioPerl modules L<Bio::DB::HIV> and L<Bio::DB::Query::HIVQuery>
together allow batch queries against the Los Alamos National
Laboratories' HIV Sequence Database using a simple query
language. C<bp_hivq.PL> provides both an example of the use of these
modules, and a standalone interactive command-line interface to the
LANL HIV DB. Simple commands allow the user to retrieve HIV sequences
and annotations using the query language implemented in
L<Bio::DB::Query::HIVQuery>. Visit the man pages for those modules for
more details.

=head1 USAGE

Run the script using C<perl bp_hivq.PL> or, in Unix, C<./bp_hivq.PL>. You will see
the

 hivq>

prompt. Type commands with queries to retrieve sequence and annotation
data.  See the SYNOPSIS for a sample session. Available commands
are described below.

=head2 TIPS

The LANL database is pretty complex and extensive. Use the C<find> facility to
explore the available database tables and fields. To identify aliases for a particular field, use C<find alias [fieldname]>. For example, to find a short alias to the 
weirdly named field C<seq_sample.ssam_second_receptor>, do

 hivq> find alias seq_sample.ssam_second_receptor

which returns

 coreceptor             second_receptor 

Now, instead of the following query

 hivq> run C[subtype] CCR5[seq_sample.ssam_second_receptor]

you know you can do

 hivq> run C[subtype] CCR5[coreceptor]

Use the C<outfile> command to set the file that receives the retrieved
sequences. You can change the current output file simply by issuing a
new C<outfile> command during the session. The output file defaults to 
standard output.

Use the C<query> command to validate a query without hitting the
database. Use the C<prerun> or C<count> commands to get a count of
sequence hits for a query without retrieving the data. Use C<run> or
C<do> to perform a complete query, retrieving sequence data into the
currently set output files.

To process C<bp_hivq.PL> commands in batch, create a text file
(C<bp_hivq.cmd>, for example) containing desired commands one per
line. Then execute the following from the shell:

 $ cat bp_hivq.cmd | perl bp_hivq.PL

=head1 COMMANDS

Here is a complete list of commands. Options in single brackets (C<[req_option]>) are 
required; options in double brackets (C<[[opt_option]]>) are optional.

 confirm            : Toggle interactive confirmation before 
                      executing queries
 exit               : Quit script
 find               : Explore database schema
  find tables                 Display all database tables
  find fields                 Display all database fields (columns)
  find fields [table]         Display all fields in [table]
  find alias [field]          Display valid aliases for [field]
 help [[command]]   : Show command help
                      if [[command]] not specified, list all 
                      available commands
 id                 : Display current session id
 outfile [filename] : Set file for collecting retrieved data
 ping               : Check if LANL DB is available
 prerun [[query]]   : Execute query but retreive hit count only
                      if [[query]] not specified, use current query
 query [query]      : Validate and set current query
 run [[query]]      : Execute query and retrieve data
                      if [[query]] not specified, use current query
 state              : Display current state of the script

 bye                : Alias for 'exit'
 config             : Alias for 'state'
 count              : Alias for 'prerun'
 do                 : Alias for 'run'
 out                : Alias for 'outfile'
 quit               : Alias for 'exit'

=head1 OPTIONS

 -v : verbose; turns on the internal debug() function

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Mark A. Jensen

Mark A. Jensen  E<lt>maj@fortinbras.usE<gt>

=cut

$|=1;
use strict;
use warnings;
use Term::ReadLine;
use Bio::DB::HIV;
use Bio::DB::Query::HIVQuery;
use Bio::SeqIO;
use Error qw(:try);

use Getopt::Std;
our ($opt_v);
getopts('v');
my $db = new Bio::DB::HIV;

my $q = new Bio::DB::Query::HIVQuery();
my $schema = $q->_schema;
my $seqio;
my @cmds = qw(
bye
config
confirm
count
do
exit
find
help
id
out
outfile
ping
prerun
query
quit
run
state
);

my %state = (
    'confirm'    => 0,
    'outfile'    => {'text'=>'[stdout]','handle'=>\*STDOUT},
    'query'      => '',
    'session_id' => '',
    'timeout'    => '',
    'curct'      => 0
);

my %help = (
    'bye' => 'bye: End script',
    'config' => 'config: Print current configuration state of hivq',
    'confirm' => 'confirm: Toggle user confirmation of query execution',
    'count' => 'count [[query string]]: Get number of sequences available for query',
    'do' => 'do [[query string]]: Run query to completion',
    'exit' => 'exit: End script',
    'find' => join("\n","find tables: Print all LANL table names",
               "find fields: Print all field names",
               "find fields [tablename]: Print all fields in specified table",
               "find aliases [fieldname]: Print all valid aliases for specified field",
		   "find options [fieldname/alias]: Print valid match options for specified field"),
    'help' => 'help [[command]]: Get description of command',
    'id' => 'id: Print current session id',
    'out' => 'out [output filename]: Set filename to hold returned sequences',
    'outfile' => 'outfile [output filename]: Set filename to hold returned sequences',
    'ping' => 'ping: Ping LANL database',
    'prerun' => 'prerun [[query string]]: Get number of sequences available',
    'query' => 'query [[query string]]: Validate and set current query string',
    'quit' => 'quit: End script',
    'run' => 'run [[query string]]: Run query to completion',
    'state' => 'state: Print current configuration state of hivq',
    );
my $done = 0;
debug("We have arrived.");
#### main loop
my ($cmd, @tok);
my $prompt = "hivq> ";
my $term = new Term::ReadLine 'hivq';
CMD:
while ( !$done ) {
    @tok=();$cmd='';
    for ( $cmd = $term->readline($prompt) ) {
	chomp;
	s/^\s+//; #trim whsp
	$term->addhistory($_);
	@tok = split(/\s+/);
	$tok[0] = lc $tok[0]; # normalize
	for ($tok[0]) {
	    (!$_ || /^\#/) && do {
		next CMD;
	    };
	    !(grep(/^$tok[0]$/, @cmds)) && do {
		error("Unrecognized command \'$tok[0]\'");
		next CMD;
	    };
	    (m/^quit$/ || m/^bye$/ || m/^exit$/) && do {
		my $fh = $state{outfile}->{handle};
		unless ($fh == \*STDOUT) {
		    close $fh;
		}
		$done=1;
		next CMD;
	    };
	    m/^confirm$/ && do {
		$state{confirm} = !$state{confirm};
		print "User confirmation before running query: ". ($state{confirm} ? "yes" : "no")."\n";
		next CMD;
	    };
	    (m/^outfile$/ || m/^out$/) && do {
		unless ($tok[1]) {
		    error('Output filename not specified');
		    next CMD;
		}
		my $fh;
		eval {
		    open $fh, ">$tok[1]" or die $!;
		};
		if ($@) {
		    error($@);
		}
		else {
		    my $oldfh = $state{outfile}->{handle};
		    if ($state{outfile}->{text} ne '[stdout]') {
			close($oldfh);
		    }
		    $state{outfile}->{handle} = $fh;
		    $state{outfile}->{text} = $tok[1];
		}
		next CMD;
	    };
	    m/^query$/ && do {
		my $query = join(' ', @tok[1..$#tok]);
		if ($query) {
		    $q->_reset;
		    $q->query($query);
		}
		else {
		    if ($state{query}) {
			print "Current query: $state{query}\n";
		    }
		    else {
			print "No query set\n";
		    }
		    next CMD;
		}
		try { #catch errors
		    $q->_do_query(0);
		}
		catch Error with {
		    my $E = shift;
		    error($E->text);
		    next CMD;
		};
		$state{session_id} = $q->_session_id;
		$state{query} = $query;
		next CMD;
	    };
	    (m/^prerun$/ || m/^count$/) && do {
		my $query = join(' ', @tok[1..$#tok]);
		if (!$query) {
		    if (!$q->query()) {
			error('No query currently set');
			next CMD;
		    }
		    elsif ($q->{'_RUN_LEVEL'} >= 1 &&
			($q->query() eq $state{query})) {
			outputPrerun();
			next CMD;
		    }
		    # else, fallthrough to prerun current query
		}
		else { # new query
		    $q->_reset;
		    $q->query($query);
		    try {
			$q->_do_query(0);
		    }
		    catch Error with {
			my $E = shift;
			error($E->text);
			next CMD;
		    };
		}
		# on success:
		$state{confirm} && do { next CMD unless getConfirm(); };
		$state{session_id} = $q->_session_id;
		$state{query} = $query;
		try {
		    $q->_do_query(1);
		}
		catch Error with {
		    my $E = shift;
		    error($E->text);
		    next CMD;
		};
		# on success:
		$state{curct} = $q->count;
		$state{query} = $q->query;
		outputPrerun();
		next CMD;

	    };
	    (m/^run$/ || m/^do$/) && do {
		my $query = join(' ', @tok[1..$#tok]);
		if (!$query) {
		    unless (defined $q->{'_RUN_LEVEL'} && $q->{'_RUN_LEVEL'} < 2) {
			error('No query currently set');
			next CMD;
		    }
		}
		else { # new query
		    $q->_reset;
		    $q->query($query);
		    try {
			$q->_do_query(0);
		    }
		    catch Error with {
			my $E = shift;
			error($E->text);
			next CMD;
		    };
		    $state{query} = $query;
		    $state{curct} = 0;
		    $state{session_id} = $q->_session_id;
		}
		# on success:
		$state{confirm} && do { next CMD unless getConfirm(); };
		try {
		    $q->_do_query(2);
		    $seqio = $db->get_Stream_by_query($q);
		}
		catch Error with {
		    my $E = shift;
		    error($E->text);
		    next CMD;
		};
		# on success:
		$state{curct} = $q->count;
		# output the seqs
		if ($q->count) {
		    outputSeqs($seqio);
		}
		else {
		    error('No sequences returned');
		}
		next CMD;
	    };
	    m/^find$/ && do {
		outputFind(@tok);
		next CMD;
	    };
	    m/^help$/ && do {
		unless ($tok[1]) {
		    outputUsage();
		    next CMD;
		}
		if (grep(/$tok[1]/, @cmds)) {
		    print $help{$tok[1]}, "\n";
		}
		else {
		    error("Command \'$tok[1]\' unrecognized\n");
		}
		next CMD;
	    };
	    (m/^state$/ || m/^config$/) && do {
		foreach my $k (sort keys %state) {
		    if (ref($state{$k}) eq 'HASH') {
			print "$k: ".$state{$k}->{text}."\n";
		    }
		    else {
			print "$k: ";
			if ($k eq 'confirm') {
			    print $state{$k} ? "yes" : "no";
			    print "\n";
			}
			else {
			    print "$state{$k}\n";
			}
		    }
		}
		next CMD;
	    };
	    do {
		error("Command currently unimplemented\n");
		next CMD;
	    };
	}
    }
}

### end main

sub error {
    my $msg = shift;
    if ($msg =~ /MSG:/) {
	($msg) = grep (/^MSG:/, split(/\n|\r/,$msg));	
	$msg =~ s/^MSG: *//;
	$msg =~ s/\sat\s.*$//;
    }
    print STDERR "hivq: $msg\n";
    return 1;
}

sub outputPrerun {
    print (($state{curct} ? $state{curct} : "No")
	. " sequence"
	. ($state{curct}>1 ? "s" : "") 
	. " returned\n");
    print "Query: ".$state{query}."\n";
    return 1;
}

sub outputSeqs {
    my ($seqio) = @_;
    my @ret;
    while (my $seq = $seqio->next_seq) {
	my $nameline = ">";
	$nameline .= $seq->annotation->get_value('Special', 'accession').
	    '('.$seq->id.') ';
	# loop through categories:
	foreach my $cat ($seq->annotation->get_keys()) {
	    foreach my $an ($seq->annotation->get_keys($cat)) {
		next if ($an eq 'accession');
		my $value = $seq->annotation->get_value($cat, $an);
		# next line: kludge to skip if there's an annotation 
		# object instead of a value (I believe this is a bug)
		next if ref($value);
		$nameline .= "\t".join('=', "'$an'", "'".$value."'");
	    }
	}
	push @ret, $nameline."\n";
	push @ret, $seq->seq()."\n";
    }
    doOutputSeqs(@ret);
}

sub doOutputSeqs {
    my @lines = @_;
    if (!@lines) {
	error("No sequences to output");
    }
    my $fh = $state{outfile}->{handle};
    print $fh @lines;
    print "Download complete\n";
    return 1;
}

sub getConfirm {
    print "Run query? [y/n]";
    my $ans = <STDIN>;
    ($ans =~ /^[yY]/) && do {return 1;};
    return 0;
}

sub outputUsage {
    print "Available commands:\n";
    outputInColumns(\@cmds);
    return 1;
}

sub outputFind {
    my @tok = @_;
    for my $arg ($tok[1]) {
	!$arg && do {
	    print $help{find},"\n";
	    return;
	};
	($arg =~ m/^tables$/) && do {
	    outputInColumns([$schema->tables]);
	    return;
	};
	($arg =~ m/^fields$/) && do {
	    my $tbl = $tok[2];
	    if (!$tbl) {
		outputInColumns([$schema->fields]);
	    }
	    else {
		unless (grep /^$tbl$/, $schema->tables) {
		    error("Table \'$tbl\' not valid");
		    return;
		}
		outputInColumns([grep /^$tbl/, $schema->fields()]);
	    }
	    return;
	};
	($arg =~ /^options$/) && do {
	    my $fld = $tok[2];
	    my @aliases = $schema->aliases($fld);
	    if (!@aliases) {
		unless (grep /^$fld$/, $schema->aliases) {
		    error("Field \'$tok[2]\' not valid");
		    return;
		}
		foreach ($schema->fields) {
		    if ( grep /$fld/, $schema->aliases($_) ) {
			$fld = $_;
			last;
		    }
		}
		# on success: $fld is set to valid field
	    }
	    my @options = sort {$a cmp $b}  $schema->options($fld);
	    @options = (@options ? @options : ('Free text'));
	    outputInColumns(\@options);
	    return;
	};
	do {
	    ($arg =~ /^alias/) && ($arg = $tok[2]);
	    if (grep /$arg/, $schema->fields) {
		my @aliases = sort $schema->aliases($arg);
		if (@aliases) {
		    outputInColumns(\@aliases);
		}
		else {
		    error("No aliases to field \'$arg\'");
		}
	    }
	    else {
		error("Field \'$arg\' not valid");
	    }
	    return;
	};
    }
}




sub outputInColumns {
    my ($items, $n, $w) = @_;
    my $i = 0;
    my @items = @$items;
    my $maxl = length([sort {length($b)<=>length($a)} @items]->[0]);
    $n ||= int(60/($maxl+3)) || 1;
    $w ||= $maxl+3;
    my $coll = int(@items/$n);
    $coll == @items/$n ? $coll : ++$coll;
    my @t;
    for my $j (0..$n-1) {
	if ($j<$n-1) {
	    $t[$j] = [@items[$j*$coll..$j*$coll+$coll-1]];
	}
	else {
	    $t[$j] = [@items[$j*$coll..$#items]];
	}
    }
    @items = map { my $j = $_; map { $t[$_]->[$j] || () } (0..$#t) } (0..scalar(@{$t[0]})-1);
    foreach (@items) {
	$_ .= (' ' x ($w-length($_)));
    }
    while ($i < @items) {
	print join('', map { $_ || ()  } @items[$i..$i+$n-1] ),"\n";
	$i+=$n;
    }
    return 1;
}
    
sub debug {
    print STDERR shift()."\n" if $opt_v;
    return;
}
