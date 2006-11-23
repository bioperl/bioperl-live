#!/usr/bin/perl -w

# This is a stub that simply tells you to use Build.PL instead

warn <<'WARN';

The installation process for Bioperl no longer uses Makefile.PL

Please use Build.PL instead:
perl Build.PL
./Build test
./Build install

For futher details, see:
http://search.cpan.org/~kwilliams/Module-Build/lib/Module/Build.pm

If you want to stick with using Makefile.PL and make syntax, you can
generate a working Makefile.PL by:
perl Build.PL
./Build manifest
./Build distmeta

Now you can:
perl Makefile.PL
make test
make install

WARN

exit;
