package Bio::DB::GFF::Adaptor::dbi::pg_fts;


=head1 NAME

Bio::DB::GFF::Adaptor::dbi::pg_fts -- Database adaptor for a specific postgres schema with a TSearch2 implementation

=head1 SYNOPSIS

    #create new GFF database connection
    my $db      = Bio::DB::GFF->new( -adaptor => 'dbi::pg_fts',
                                     -dsn     => 'dbi:Pg:dbname=worm');

    #add full text indexing 'stuff'
    #assumes that TSearch2 is available to PostgreSQL
    #this will take a VERY long time for a reasonably large database
    $db->install_TSearch2();

    ...some time later...
    #we don't like full text searching...
    $db->remove_TSearch2();

=head1 DESCRIPTION

This adaptor is based on Bio::DB::GFF::Adaptor::dbi::pg but it implements
the TSearch2 PostgreSQL contrib module for fast full text searching.  To
use this module with your PostgreSQL GFF database, you need to make
TSearch2 available in the database. 

To use this adaptor, follow these steps:

=over

=item Install TSearch2 contrib module for Pg

Can be as easy as `sudo yum install postgresql-contrib`, or you may
need to recompile PostgreSQL to include it.  See
L<http://www.sai.msu.su/~megera/postgres/gist/tsearch/V2/docs/tsearch-V2-intro.html>
for more details

=item Load the TSearch2 functions to you database

  % cat tsearch2.sql | psql <your database>

=item Load your data using the pg adaptor:

 % bp_pg_bulk_load_gff.pl -c -d yeast saccharomyces_cerevisiae.gff

or

 % bp_load_gff.pl -c -d yeast -a dbi::pg saccharomyces_cerevisiae.gff

=item Add GFF/TSearch2 specific modifications

Execute a perl script like this one:

  #!/usr/bin/perl -w
  use strict;

  use Bio::DB::GFF;

  my $db = Bio::DB::GFF->new(
      -adaptor   => 'dbi::pg_fts',
      -dsn       => 'dbi:Pg:dbname=yeast',
      -user      => 'scott',
    );

  print "Installing TSearch2 columns...\n";

  $db->install_TSearch2();

  print "Done\n";

=back

Note that this last step will take a long time.  For a S. cerevisiae
database with 15K rows, it took over an hour on my laptop, and
with a C. elegans database (~10 million rows) it took well over a day.

If at some point you add more data you your database, you need to run
a similar script to the one above, only executing the update_TSearch2()
method.  Finally, if you want to remove the TSearch2 columns from your 
database and go back to using the pg adaptor, you can execute a script
like the one above, only executing the remove_TSearch2() method.

=head1 NOTES ABOUT TSearch2 SEARCHING

You should know a few things about how searching with TSearch2 works in
the GBrowse environment:

=over

=item 1

TSearch2 does not do wild cards, so you should encourage your users not
to use them.  If wild cards are used, the adaptor will fall back on 
an ILIKE search, which will be much slower.

=item 2

However, TSearch2 does do 'word stemming'.  That is, if you search
for 'copy', it will find 'copy', 'copies', and 'copied'.

=item 3

TSearch2 does not do phrase searching; all of the terms in the
search string are ANDed together.

=back

=head1 ACKNOWLEDGEMENTS

Special thanks to Russell Smithies and Paul Smale at AgResearch in
New Zealand for giving me their recipe for doing full text indexing
in a GFF database.

=head1 BUGS

Please report bugs to the BioPerl and/or GBrowse mailing lists
(L<mailto:bioperl-l@lists.open-bio.org> and L<mailto:gmod-gbrowse@lists.sourceforge.net>
respectively).

=head1 SEE ALSO

Please see L<Bio::DB::GFF::Adaptor::dbi::pg> for more information
about tuning your PostgreSQL server for GFF data, and for general
information about GFF database access, see L<Bio::DB::GFF>.

=head1 AUTHOR

Scott Cain, cain@cshl.edu

=head1 APPENDIX

=cut

# a simple postgres adaptor
use strict;
use Bio::DB::GFF::Adaptor::dbi;
use base qw(Bio::DB::GFF::Adaptor::dbi::pg);

use constant FULLTEXTSEARCH => <<END;
SELECT distinct gclass,gname,fattribute_value,fmethod,fsource
    FROM fgroup,fattribute_to_feature,fdata,ftype
     WHERE fgroup.gid=fdata.gid
       AND fdata.fid=fattribute_to_feature.fid
       AND fdata.ftypeid=ftype.ftypeid
       AND (fattribute_to_feature.idxfti @@ to_tsquery('default', ?))
END
;

use constant FULLTEXTWILDCARD => <<END;
SELECT distinct gclass,gname,fattribute_value,fmethod,fsource
    FROM fgroup,fattribute_to_feature,fdata,ftype
     WHERE fgroup.gid=fdata.gid
       AND fdata.fid=fattribute_to_feature.fid
       AND fdata.ftypeid=ftype.ftypeid
       AND lower(fattribute_to_feature.fattribute_value) LIKE lower(?)
END
;

sub new {
  my $class = shift;
  my $self  = $class->SUPER::new(@_);
  return $self;
}

=head2 search_notes

 Title   : search_notes
 Usage   : @search_results = $db->search_notes("full text string",$limit)
 Function: Search the notes for a text string, using PostgreSQL TSearch2
 Returns : array of results
 Args    : full text search string, and an optional row limit
 Status  : public

This is based on the mysql-specific method that makes use of the TSearch2
functionality in PosgreSQL's contrib directory. Given a search string,
it performs a full-text search of the notes table and returns an array
of results.  Each row of the returned array is a arrayref containing
the following fields:

  column 1   A Bio::DB::GFF::Featname object, for passing to segment()
  column 2   The text of the note
  column 3   A relevance score.

=cut

sub search_notes {
  my $self = shift;
  my ($search_string,$limit) = @_;

  my @terms = split /\s+/, $search_string;

  my $sth;
  if ($search_string =~ /\*/) {
      $search_string =~ tr/*/%/s;
      my $query = FULLTEXTWILDCARD;
      $query   .= " limit $limit" if defined $limit;
      $sth      = $self->dbh->do_query($query,$search_string);
  }
  elsif (@terms == 1) {
      my $query = FULLTEXTSEARCH;
      $query   .= " limit $limit" if defined $limit;
      $sth      = $self->dbh->do_query($query,$search_string);
  }
  else {
      my $query = FULLTEXTSEARCH;
      my $andstring = join (' & ', @terms);
#      $query   .= qq{ AND (fattribute_to_feature.fattribute_value ILIKE '\%$search_string%')};
      $query   .= " LIMIT $limit" if defined $limit;
      $sth      = $self->dbh->do_query($query,$andstring);
  } 
  
  my @results;
  while (my ($class,$name,$note,$method,$source) = $sth->fetchrow_array) {

     next unless $class && $name;    # sorry, ignore NULL objects
     my $featname = Bio::DB::GFF::Featname->new($class=>$name);
     my $type     = Bio::DB::GFF::Typename->new($method,$source);
     push @results,[$featname,$note,0,$type]; #gbrowse expects a score, but
                                              #pg doesn't give one, thus the 0
  }

  return @results;
}

=head2 make_features_by_name_where_part

 Title   : make_features_by_name_where_part
 Function: constructs a TSearch2-compliant WHERE clause for a name search
 Status  : protected

=cut

#need a make_features_by_name_where_part method to override pg
sub make_features_by_name_where_part {
  my $self = shift;
  my ($class,$name) = @_;

  my @terms = split /\s+/, $name; 

  if ($name =~ /\*/) {
    $name =~ tr/*/%/s;
    return ("fgroup.gclass=? AND lower(fgroup.gname) LIKE lower(?)",$class,$name);
  }
  else {
    my $where_str = "fgroup.gclass=? AND (fgroup.idxfti @@ to_tsquery('default', ?)) ";
    if (@terms == 1) {
      return ($where_str,$class,$name);
    }
    else {
      my $andstring = join (' & ', @terms);
#      $where_str .= qq{ AND (fgroup.gname ILIKE '\%$name%')};
      return ($where_str,$class,$andstring); 
    }
  }
}

=head2 install_TSearch2

 Title   : install_TSearch2
 Function: installs schema modifications for use with TSearch2
 Usage   : $db->install_TSearch2
 Status  : public

=cut


#needs method for installing TSearch2 (does that mean that the SQL for
#creating the tables and functions should go in here?  That would be
#the safest and easiest thing to do
sub install_TSearch2 {
  my $self = shift;

  my $dbh = $self->features_db;

  $dbh->do('ALTER TABLE fattribute_to_feature ADD COLUMN idxFTI tsvector') 
     or $self->throw('adding FTI column to f_to_f failed');

  $dbh->do('ALTER TABLE fgroup ADD COLUMN idxFTI tsvector')
     or $self->throw('adding FTI column to fgroup failed');

  $self->update_TSearch2();

  return;
}

=head2 update_TSearch2

 Title   : update_TSearch2
 Function: Updates TSearch2 columns
 Usage   : $db->update_TSearch2
 Status  : public

=cut


sub update_TSearch2 {
  my $self = shift;

  my $dbh = $self->features_db;

  $self->warn('updating full text column; this may take a very long time...');
  $dbh->do("UPDATE fattribute_to_feature "
          ."SET idxFTI= to_tsvector('default', fattribute_value) "
          ."WHERE idxFTI IS NULL") 
       or $self->throw('updating fti column failed');
  $dbh->do("UPDATE fgroup "
          ."SET idxFTI= to_tsvector('default', gname) "
          ."WHERE idxFTI IS NULL")
       or $self->throw('updating fgroup fti column failed');

  $self->warn('Preliminary optimization of database; this may also take a long time...');
  $dbh->do('VACUUM FULL ANALYZE')
       or $self->throw('vacuum failed');

  $self->warn('Updating full text index; again, this may take a long time');
  $dbh->do('CREATE INDEX idxFTI_idx ON fattribute_to_feature '
          .'USING gist(idxFTI)')
       or $self->warn('creating full text index failed');
  $dbh->do('CREATE INDEX fgroup_idxFTI_idx ON fgroup '
          .'USING gist(idxFTI)')
       or $self->warn('creating fgroup full text index failed');

  $self->warn('Optimizing database; hopefully, this will not take as long as other steps');
  $dbh->do('VACUUM FULL ANALYZE');
  $dbh->do("SELECT set_curcfg('default')");

  return;
}

=head2 remove_TSearch2

 Title   : remove_TSearch2
 Function: Removes TSearch2 columns
 Usage   : $db->remove_TSearch2
 Status  : public

=cut

sub remove_TSearch2 {
  my $self = shift;

  my $dbh = $self->features_db;

  $self->warn('Removing full text search capabilities');
  $dbh->do('DROP INDEX idxFTI_idx')
     or $self->throw('dropping full text index failed');
  $dbh->do('DROP INDEX fgroup_idxFTI_idx')
     or $self->throw('dropping full text index failed');

  $dbh->do('ALTER TABLE fattribute_to_feature DROP COLUMN idxFTI')
     or $self->throw('dropping full text column failed');
  $dbh->do('ALTER TABLE fgroup DROP COLUMN idxFTI')
     or $self->throw('dropping full text column failed');


  return;
}


1;
