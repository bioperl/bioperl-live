#!/usr/bin/perl

# bioperl.pl
# cjm@fruitfly.org

use strict;
use lib '.';
no strict "vars";
use Data::Dumper;
use Bio::Perl;
use Bio::SeqIO;
use Getopt::Long;
my $h = {};

GetOptions($h,
           "file|f=s",
           );
my @cmds = get_default_cmds();
shell($h, \@cmds, @ARGV);

# prepare for some seriously hacky code....
sub shell {
    my $h = shift;
    my @cmds = @{shift || []};
    my @args = @_;
    my $prompt = $ENV{BIOPERL_PROMPT} || "BioPerl> ";
    my $quit = 0;
    my @lines = ();
    my $r;
    my $rv;
    my $seq;
    my @pseqs = ();
    my $seqio;
    my $wseqio;
    my $fastadb;
    my $options =
      {echo=>0, chatty=>10};

    my $loadfn = $h->{'file'};
    if ($loadfn) {
        @lines = ("load '$loadfn'");
    }

    sub hr {
	print "\n===============================\n";
    }

    sub nl {
	print "\n";
    }

    sub demo {
        if (! -d 't/data') {
            print "To run the demo, you must be in the bioperl directory\n";
        }
        @lines = 
          split(/\n/,
                q[
                  %keep = %$options;
                  +format ''
                  +outformat ''
                  +echo 1
                  # BioPerl shell utility - Demo
                  #
                  # We're now going to take a tour
                  # through some of the features of
                  # this tool.
                  #
                  #
                  # This demo will go through some of
                  # the major commands, feeding you
                  # the commands as you go. all you have
                  # to do is hit <ENTER> every time
                  # you see the prompt $prompt
                  # you will then see the output of
                  # the command on your terminal window.
                  # type 'q' to end the tour
                  # at any time.
                  #
                  waitenter
                  # PARSING GENBANK RECORDS
                  # -----------------------
                  # to parse genbank files, use
                  # the read_seq() method, or
                  # simply use the '<' command.
                  #
                  # First of all we're going to
                  # take a look at the file
                  # 't/data/test.genbank'
                  # Let's examine the file itself
                  # using the unix command "cat"
                  # (you can use any unix command
                  #  using the ! at the beginning
                  #  of a line)
                  ^!cat t/data/test.genbank
                  waitenter
                  # Ok, you can see this is a
                  # typical file of genbank records.
                  # Let's get the first sequence
                  # from the file
                  ^<t/data/test.genbank
                  waitenter
                  # we have parsed the first
                  # record of the file, and placed
                  # the sequence object into
                  # the variable $seq
                  #
                  # if you are familiar with perl
                  # objects and the bioperl object
                  # model, you can interact with
                  # the object; for instance, to
                  # display the residues we use the
                  # seq() method like this:
                  ^print $seq->seq()
                  waitenter
                  #
                  # we can cycle through all the
                  # sequences in the file using
                  # the ',' command.
                  ^,
                  waitenter
                  # this fetched the second sequence
                  # and placed it in the $seq variable
                  #
                  # we can change the output format
                  # by setting the 'outformat' parameter
                  # like this:
                  ^+outformat fasta
                  ^,
                  waitenter
                  # now the sequences are output in
                  # fasta format
                  # to change to embl format:
                  ^+outformat embl
                  ^,
                  waitenter
                  # we can also fetch _all_ seqs from
                  # a file; for this example we will
                  # use t/data/swiss.dat, which is in
                  # swiss format. usually bioperl can guess
                  # the file format from the file extension
                  # but this isn't possible here, so we
                  # must help by setting the input format:
                  ^+format swiss
                  # now lets get all the sequences, like this:
                  ^<*t/data/swiss.dat
                  waitenter
                  # typing <* is equivalent to
                  # using the read_seqs() function,
                  # like this:
                  ^read_seqs('t/data/swiss.dat')
                  waitenter
                  # we now have all the sequences in
                  # the array @seqs
                  # we can write these all out as fasta
                  ^+outformat fasta
                  ^>*
                  # we can also write these out to a file:
                  ^>*myfile.tmp
                  ^!cat myfile.tmp
                  #
                  # RANDOM ACCESS OF FASTA FILES
                  # END
                  +echo 0
                  %$options = %keep
                 ]);

        @lines = 
          map {
              s/^ *//;
              $_;
          } @lines;
    }

    sub error {
        if ($error) {
            print "Error:\n$error";
        }
        else {
            print "No errors have been reported\n";
        }
    }

    sub fmt {
        $options->{format} = shift if @_;
        print "format=$options->{format}\n";
    }

    # should this move to Bio::Perl ?
    sub seqio {
        my $filename = shift;
        $options->{format} = shift if @_;

        if( !defined $filename ) {
            warn "read_sequence($filename) - usage incorrect";
        }

        if( defined $options->{format} ) {
            $seqio = Bio::SeqIO->new( '-file' => $filename, '-format' => $options->{format});
        } else {
            $seqio = Bio::SeqIO->new( '-file' => $filename);
        }

        $seqio;
    }

    sub wseqio {
        my $filename = shift;
        $options->{format} = shift if @_;

        my @args = ();
        if ($filename && $filename !~ /^\>/) {
            $filename = ">$filename";
        }
        push(@args, -file => "$filename") if $filename;
        push(@args, -fh => \*STDOUT) unless $filename;
        push(@args, -format => $options->{outformat}) if $options->{outformat};
        $wseqio = Bio::SeqIO->new( @args );

        $wseqio;
    }

    sub show_seq {
        return unless $seq;
        if ($wseqio) {
            $wseqio->write_seq($seq);
        }
        else {
            printf "seq display id: %s\n", $seq->display_id;
        }
    }

    sub addseq {
        push(@pseqs, @_);
        while (scalar(@pseqs) > 50) {
            # todo - history variable
            shift @pseqs;
        }
    }

    sub next_seq {
        if ($seqio) {
            eval {
                $seq = $seqio->next_seq;
            };
            if ($@) {
                $error = $@;
                print "There was an error getting the seq. Type 'error'\n";
                print "for full details\n";
                print "(Maybe you have to explicitly set the format?)";
            }
            addseq($seq);
        }
        else {
            print STDERR "use read_seq first\n";
        }
        show_seq;
        $seq;
    }

    sub next_seqs {
        @seqs = ();
        if ($seqio) {
            while ($seq = $seqio->next_seq) {
                printf "%s\n", $seq->display_id;
                push(@seqs, $seq);
            }
        }
        $seq = $seqs[$#seqs] if @seqs;
        @seqs
    }

    sub read_seq {
        seqio(@_);
        next_seq();
    }

    sub read_seqs {
        seqio(@_);
        next_seqs();
    }

    sub write_seq {
        wseqio(@_);
        $wseqio->write_seq($seq) if $seq;
    }

    sub write_seqs {
        wseqio(@_);
        map {
            $wseqio->write_seq($_)
        } @seqs;
    }

    sub pod {
        if (!-d "Bio") {
            print "You need to be in the bioperl directory!\n";
        }
        else {
            my $mod = shift;
            unix("pod2text", "Bio/$mod.pm");
        }
    }
    sub fastadb {
        require "Bio/DB/Fasta.pm";
        my $f = shift;
        $fastadb = Bio::DB::Fasta->new($f);
        print "Set \$fastadb object\n";
        $fastadb;
    }

    sub subseq {
        if (!$fastadb) {
            fastadb(shift);
        }
        $seq = $fastadb->get_Seq_by_id(shift);
        if (@_) {
            printf "%s\n",
              $seq->subseq(@_);
        }
        $seq;
    }

    sub load {
        open(F, shift);
        @lines = map {chomp;$_} <F>;
        close(F);
    }

    sub waitenter {
        print "<hit ENTER to continue>";
        <STDIN>;
    }

    sub showintro {
        hr;
        print "This is a text-based commandline interface to BioPerl;\n";
        print "\n";
    }

    sub checkoptions {
    }

    sub showoptions {
        my $k = shift;
        my @k = defined $k ? ($k) : keys %$options;
        foreach my $ok ($k) {
            my $v = sprintf("%s", $options->{$k});
            if ($v =~ /HASH/) {
                # hide perl internal details
                # from user; if they are experienced
                # perlhackers they can just
                # type "x $options" to see the
                # gory details
                $v = "undisplayable";
            }
            printf("%20s:%s\n",
                   $ok,
                   $b);
        }
    }

    sub set {
        my ($k,$v) = @_;
        if (defined($v)) {
            $options->{$k} = $v;
            checkoptions;
        }
        else {
            showoptions($k);
        }
#        if ($k eq "format") {
#            seqio;
#        }
        if ($k eq "outformat") {
            wseqio;
        }
    }

    sub echo {
        my $e = shift;
        if (defined($e)) {
            set("echo", $e);
        }
        else {
            set("echo", !$options->{echo});
        }
    }

    sub options {
        map {print "$_ = $options->{$_}\n"} keys%$options;
    }

    sub showcommands {
        hr;
        print "BioPerl Shell Commands:\n";
        my $layout = "%5s : %-20s - %s\n";
        printf $layout, "cmd", "function", "summary";
        printf "%s\n", ("-" x 40);
        foreach my $c (@cmds) {
            my $sc = $c->{shortcut};
            $sc =~ s/\\//g;
            printf($layout,
                   $sc,
                   $c->{'func'} . "()",
                   $c->{'summary'}
                   );
        }
        
    }

    sub showexamples {
        print "\nExamples:\n-------\n";
    }

    sub showvariables {
        hr;
        print "Shell variables:\n";
        print q[
                $seq     : Bio::SeqI object
                $seqio   : Bio::SeqIO object
                @pseqs   : array of previous Bio::SeqI objects
               ];
        nl;
    }

    sub welcome {
	print "Welcome to the BioPerl shell interface!\n\n";
        print "\n\nType 'help' for instructions\n";
        print "\n\nType 'demo' for demonstration\n";
        print "\n\nThis is ALPHA software - commands may change\n";
        print "-lots more commands need to be added to take full\n";
        print "advantage of the bioperl functionality\n\n";
    }

    sub help {
        my $topic = shift;
        my $c;
        if ($topic) {
            ($c) = grep {$_->{func} eq $topic} @cmds;
        }
        if ($c) {
            print "Function:   $c->{func}\n";
            print "Shortcut:   $c->{shortcut}\n" if $c->{shortcut};
            print "Summary:    $c->{summary}\n" if $c->{summary};
            print "=======\n";
            print "$c->{docs}\n" if $c->{docs};
        }
        elsif ($topic eq "advanced") {
            hr;
            nl;
        }
        else {
            hr;
            print "\nBioPerl Shell Help\n\n";
            showintro;
            waitenter;
            showcommands;
            waitenter;
            showvariables;
            waitenter;
            showexamples;
            nl;
            nl;
            nl;
            print "Type \"demo\" for an interactive demo of commands\n\n";
            print "Type \"help advanced\" for advanced options\n\n";
            hr;
            nl;
        }
    }

    sub p {
	print shift;
	print "\n";
    }

    sub x {
	print Dumper shift;
	print "\n";
    }

    # trick to allow barewords as keywords...
    sub advanced {"advanced"}

    sub unix {
        my @cmds = @_;
        my $c = join(" ", @cmds);
        print `$c`;
    }

    welcome;
    require Term::ReadLine;
    require Shell;

    checkoptions;
    print "\n";
    my $termline = shift || Term::ReadLine->new($prompt);

    my $rcfile = "$ENV{HOME}/.goshellrc";
    if (-f $rcfile) {
	open(F, $rcfile);
	@lines = <F>;
	close(F);
	
    }
    my $end_signal = "";
    while (!$quit) {
	if ($end_signal) {
	    @lines = ($lines);
	    while ($end_signal && ($line = $termline->readline("? "))) {
		if($line !~ /$end_signal/) {
		    $lines[0].= "\n$line";
		}
		else {
		    $end_signal = "";
		}
	    }
	    next;
	}
	my $line = 
	  @lines ? shift @lines : $termline->readline($prompt);
        if ($line =~ /^\^/) {
            $line =~ s/^\^//;
            print "$prompt$line";
            my $e = <STDIN>;
            if ($e =~ /^q/) {
                $line = "";
                @lines = ();
            }
        }
        if ($options->{echo} && $line !~ /\+?wait/) {
            if ($line =~ /^\#/) {
                print "$line\n";
            }
            else {
                print "$prompt$line\n";
            }
            if ($options->{sleep}) {
                sleep $options->{sleep};
            }
            if ($options->{wait}) {
                sleep $options->{wait};
            }
        }
	my ($cmd, @w) = split(' ',$line);

	$_ = $cmd;
	if (/^\<\<(.*)/) {
	    $end_signal = $1;
	}

        # check for shortcuts
        my $selected;
        foreach my $c (@cmds) {
            my $shortcut = $c->{'shortcut'};
            next unless $shortcut;
            if ($line =~ /^$shortcut(.*)/) {
                if (!defined($selected) ||
                    length($shortcut) > length($selected->{shortcut} || "")) {
                    # get the most specific match
                    $selected = $c;
                }
            }
        }
        if ($selected) {
            my $shortcut = $selected->{'shortcut'};
            if ($line =~ /^$shortcut(.*)/) {
                my @w = map {"'".$_."'" } split(' ', $1);
                $line = $selected->{'func'}." ".join(", ", @w);
            }
        }

	$rv = eval $line;
#        print "\n";
#        print "RV=$rv;;;\n";
	if ($@) {
	    print STDERR $@;
	}
        if ($options->{sleep}) {
            sleep $options->{sleep};
        }
        if ($options->{wait}) {
            sleep $options->{wait};
            $options->{wait} = 0;
        }

    }
}

sub get_default_cmds {
    my @cmds =
      (
       {
        func         =>  'read_seq',
        shortcut     =>  '\<',
        summary      =>  'read a Seq from a file', 
       },

       {
        func         =>  'next_seq',
        shortcut     =>  ',',
        summary      =>  'get the next Seq', 
       },

       {
        func         =>  'read_seqs',
        shortcut     =>  '\<\*',
        summary      =>  'read all Seqs from a file', 
       },

       {
        func         =>  'write_seq',
        shortcut     =>  '\>',
        summary      =>  'write a Seq to screen/file', 
       },

       {
        func         =>  'write_seqs',
        shortcut     =>  '\>\*',
        summary      =>  'write a Seq to screen/file', 
       },

       {
        func         =>  'fastadb',
        shortcut     =>  'fa',
        summary      =>  'fast fasta access', 
       },

       {
        func         =>  'subseq',
        summary      =>  'get a subseq from a fastadb', 
       },

       {
        func         =>  'set',
        shortcut     =>  '\+',
        summary      =>  'set a shell parameter', 
       },

       {
        func         =>  'unix',
        shortcut     =>  '\!',
        summary      =>  'run a unix command', 
       },

       {
        func         =>  'x',
        summary      =>  'display variable (and internals) using dumper', 
       },

      );
    return @cmds;
}
