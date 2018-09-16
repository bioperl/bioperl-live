#!/usr/bin/perl
# 
# Counts up all the used and inherited modules in a directory of modules to
# help indicate which the most important modules are, graphs it also
#
# Written by Sendu Bala, using much code directly from
# http://www.perlmonks.org/?displaytype=displaycode;node_id=87329
# and also
# http://search.cpan.org/src/NEILB/pmusage-1.2/pmusage

use strict;
use warnings;

use IO::File;
use File::Find;
use Getopt::Std;
use GraphViz;

sub usage
{
	print <<EOF;
$0 -- Shows which modules are used the most in a directory of modules

usage: $0 [-f outfile] [-l listfile] [-h] [-v] [dir]
-f outfile   specify output file (default=module_usage)
-h           get this help message
-i fmt       set image format to fmt (default=jpeg)
             also available: canon,text,ps,hpgl,pcl,mif,pic,gd,gd2,gif,jpeg,
             png,wbmp,vrml,vtx,mp,fig,svg,plain *nb* only jpeg is tested as
			 working correctly
-l listfile  get filenames/options from listfile
-v           list filenames to STDERR

If directory names are given, all the *.p[lm] files in the directory will
be processed. The default is to do all the Perl files in the current directory.


The graphical output is a *simplification with information loss*, since BioPerl
is too complex to graph 'raw'. The simplification proceeds as follows:
# It is determined which modules each BioPerl package (aka class, module)
  'uses' (what modules does it load via 'use', 'require' or inherit from
  via 'use base', excluding external (non-BioPerl) modules).
# Packages with identical usage (ignoring the type of usage) are grouped
  together.
# The graph shows all the groups with more than one member as nodes, with edges
  from them pointing to the individual packages that they use.
# The set of those individual packages pointed to by groups also have edges
  showing their use-relationship to other members of the set (only).
# Members of that set are also shaded in red. The saturation of the shade
  indicates how many packages use that package (so dark red packages are used a
  lot).
# Edges are coloured green to show inheritance-type usage, or blue to show
  'use'/'require'-type usage. It is possible that some members of a group
  inherit from a particular package, whilst other members only 'use' it. In that
  case, the colour is based on which is most common for the group.
# Groups with members all from the same subdirectory, and individual packages
  pointed to by groups, are 'clustered' if they are from the same subdirectory
EOF
	exit shift;
}

# process cmdline options
my $opts = 'f:l:hvi:';
my %opts;
getopts($opts, \%opts) || usage(1);
usage(0) if defined($opts{h});
while (defined($opts{l})) {
	my $lFile = IO::File->new($opts{l}) or die "can't open -l file $opts{l} : $!\n";
	my @largs = <$lFile>;
	chomp(@largs);
	splice(@ARGV, 0, 0, @largs);
	delete($opts{l});
	getopts($opts, \%opts) || usage(1);
	$lFile->close();
}

my $outfile = defined($opts{f}) ? $opts{f} : "module_usage";
my $format = defined($opts{i}) ? $opts{i} : 'jpeg';

# now filenames are in @ARGV
push(@ARGV, '.') if !@ARGV;

my @files;
my %sections;

sub findPerlFiles {
	-f $_ && /^.*\.p[ml]\z/si && push(@files, $File::Find::name);
}

# process directories
foreach my $top (@ARGV) {
	File::Find::find({wanted => \&findPerlFiles}, $top);
}

my %usage;
my %users;
my %inheritance;
my %packages;

sub store_package_usage {
    my ($package, $used) = @_;
    my %used = %{$used};
    
    STDERR->print("package $package used (".join(' ', keys %used).")\n") if $opts{v};
    
    $packages{$package} = \%used;
    
    foreach my $module (keys %used) {
        $usage{$module}++;
        push (@{$users{$module}}, $package);
    }
}

foreach my $file (@files) {
	$file =~ s#^./##;
	STDERR->print("processing $file\n") if $opts{v};
	my $f = IO::File->new($file) or warn "can't open $file: $!\n", next;
    
	my ($package, %used);
    
	my $pod = 0;
	while (<$f>) {
		if (/^=cut/) {
			$pod=0;
			next;
		}
		if (/^=[a-zA-Z]+/) {
			$pod=1;
			next;
		}
		next if $pod;
        
		if (/^\s*package\s+([[:word:]:]+)\s*;/) {
            if ($package) {
                store_package_usage($package, \%used);
                %used = ();
            }
            
			$package = $1;
			next;
		}
		if (/use base\s*(.*)/) {
			my $tmp = $1;
			while (!/;/)	# accumulate ISA value for multiple lines
			{
				$_ = <$f>;
				$tmp .= $_;
			}
			my @use_base = eval $tmp;
			if ($@) { warn "Unparseable 'use base' line for $package: $tmp"; next }
			
            foreach my $module (@use_base) {
                $used{$module} = 1;
				$inheritance{$package}->{$module} = 1;
            }
		}
		elsif (/^\s*use\s+([^\s;()]+)/ || /^\s*require\s+([^\s;()'"]+)/) {
            $used{$1} = 1;
		}
	}
	$f->close();
    
    if ($package) {
        store_package_usage($package, \%used);
    }
}

# simplify so we can view a graph of usage: we group all packages that have
# identical usage. NB: this doesn't look at external modules at all
my %groups;
while (my ($package, $used_hash) = each %packages) {
    my @used_packages;
    foreach my $used_module (sort keys %{$used_hash}) {
        next unless defined $packages{$used_module};
        push(@used_packages, $used_module);
    }
    @used_packages || next;
    
    push(@{$groups{join('|', @used_packages)}}, $package);
}

# we're going to shade boxes based on usage later, figure out an appropriate
# shade range by ranking
my %counts;
while (my ($group, $pack_list) = each %groups) {
    my @children = @{$pack_list};
    
    @children > 1 || next;
    
    my @parents = split(/\|/, $group);
    foreach my $parent (@parents) {
        my $count = $usage{$parent};
        $counts{$parent} = $count;
    }
}
my %ranks;
my $rank = 0;
my $prev_count;
foreach my $parent (sort { $counts{$a} <=> $counts{$b} } keys %counts) {
    my $this_count = $counts{$parent};
    $ranks{$parent} = $prev_count && $prev_count != $this_count ? ++$rank : $rank;
    $prev_count = $this_count;
}

sub class_to_subdir {
    my $class = shift;
    $class =~ s/::[^:]+$//;
    return $class;
}

my $g = GraphViz->new(concentrate => 1,
                      node => {shape => 'box'},
                      $format eq 'ps' ? (pagewidth => 46.81, pageheight => 33.11) : ()); # A0 for ps output
my $inherited_edge_colour = 'green';
my $used_edge_colour = 'blue';
my $cluster_colour = 'black'; #*** darkgray, 0,0,0.31 don't work, why?!
my $child_id = 0;
my $group_definitions = '';
my %parents;
while (my ($group, $pack_list) = each %groups) {
    my @children = @{$pack_list};
    
    # ignore single child groups (required or graph gets too wide to jpeg)
    @children > 1 || next;
    
    # we'll cluster if all children belong to the same subdirectory
    my %subdirs;
    foreach my $child (@children) {
        $subdirs{class_to_subdir($child)} = 1;
    }
    my $subdir;
    if (keys %subdirs == 1) {
        ($subdir) = keys %subdirs;
        undef $subdir if $subdir eq 'Bio';
    }
    
    my $this_child = 'group'.++$child_id;
    $g->add_node($this_child,
                 style => 'dashed',
                 label => "$this_child:\n".join("\n", @children),
                 $subdir ? (cluster => {name => $subdir, style => 'dotted', color => $cluster_colour}) : ());
    
    my @parents = split(/\|/, $group);
    
    $group_definitions .= "  $this_child consists of ".scalar(@children)." packages: ".join(', ', @children)."\n  $this_child members use ".scalar(@parents)." other packages: ".join(', ', @parents)."\n\n";
    
    foreach my $parent (@parents) {
        # we'll shade the parent box based on how many packages use it
        my $this_rank = $ranks{$parent};
        my $shade = (1 / $rank) * $this_rank;
        
		# we'll colour the edge based on if we inherited this parent or just
		# used it, going by the most common for the group
		my ($inherited, $used) = (0, 0);
		foreach my $child (@children) {
			if (defined $inheritance{$child}->{$parent}) {
				$inherited++;
			}
			else {
				$used++;
			}
		}
		my $edge_colour = $inherited > $used ? $inherited_edge_colour : $used_edge_colour;
		
        # we'll cluster if this isn't a base Bio::x class
        my $subdir = class_to_subdir($parent);
        undef $subdir if $subdir eq 'Bio';
        
        $g->add_node($parent,
                     style => 'filled',
                     fillcolor => "0,$shade,1",
                     $subdir ? (cluster => {name => $subdir, style => 'dotted', color => $cluster_colour}) : ());
        $parents{$parent} = 1;
        $g->add_edge($this_child => $parent, color => $edge_colour);
    }
}

# show links between parents
foreach my $parent (keys %parents) {
    my %used = %{$packages{$parent}};
    
    foreach my $used (keys %used) {
        next unless defined $parents{$used};
        $g->add_edge($parent => $used, color => defined $inheritance{$parent}->{$used} ? $inherited_edge_colour : $used_edge_colour);
    }
}

# write out graph
my $output = IO::File->new($outfile.".$format", 'w') or die "can't open $outfile.$format: $!\n";
$output->print(eval "\$g->as_$format()");
$output->close();


my $package_count = keys %packages;
my $total_used = keys %usage;

my $results_str = "Packages investigated: $package_count\nTotal modules used: $total_used\n\n";


# descriptive text output
# list by popularity
my @internal;
my @external;
foreach my $module (sort { $usage{$b} <=> $usage{$a} || $a cmp $b } keys %usage) {
    my $count = $usage{$module};
    
    if (defined $packages{$module}) {
        push(@internal, "  $module => used $count times");
    }
    else {
        my $by = '';
        if ($count <= 5) {
            $by = " by ".join(", ", @{$users{$module}});
        }
        push(@external, "  $module => used $count times$by");
    }
}

$results_str .= "External module usage:\n".join("\n", @external);
$results_str .= "\n\nPackage usage:\n".join("\n", @internal);

# list the packages that aren't used by any other package
$results_str .= "\n\nPackages not used by any other:\n";
foreach my $package (sort keys %packages) {
    next if $usage{$package};
    $results_str .= "  $package\n";
}

# define the groups referenced in the graph
$results_str .= "\nGroup definitions:\n$group_definitions";

# write out descriptive text file
$output = IO::File->new($outfile.'.txt', 'w') or die "can't open $outfile.txt: $!\n";
$output->print($results_str);
$output->close();

exit;
