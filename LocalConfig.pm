package LocalConfig;

use strict;
use Exporter;
use vars qw(@ISA %Local @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw( %Local );

open(CONF,'bioperl.conf');

while(<CONF>){
	chomp;
	s/\s*#.*//g; next unless $_; ## clean out comments;
	my($k,$v) = /^.*?(\S+)\s+(.+)$/;

	if($v =~ m!^file://(.+)$!){ #value is an external file reference
		open(F,$1) or die "couldn't open $v";
		local $/;
		$Local{$k} = <F>;
		close(F);
	} else { # a normal value
		$Local{$k} = $v;
	}
}
close(CONF);

1;
