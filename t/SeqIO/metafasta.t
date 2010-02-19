<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">     <html>
  <head>
    <title>view: /bioperl/bioperl-live/trunk/t/SeqIO/metafasta.t (Rev: 15112, via SVN::Web)</title>  <link rel="stylesheet" type="text/css" href="/svnweb/css/trac/svnweb.css" />   </head>
  <body>
    <div>  <div id="navpath">
	<h1><a class="first" href="/svnweb/index.cgi/">Repository List</a> /   <a href="/svnweb/index.cgi/bioperl/">bioperl</a>    <span class="sep">/</span>    <a href="/svnweb/index.cgi/bioperl/browse/bioperl-live/?rev=15112">bioperl-live</a>     <span class="sep">/</span>    <a href="/svnweb/index.cgi/bioperl/browse/bioperl-live/trunk/?rev=15112">trunk</a>     <span class="sep">/</span>    <a href="/svnweb/index.cgi/bioperl/browse/bioperl-live/trunk/t/?rev=15112">t</a>     <span class="sep">/</span>    <a href="/svnweb/index.cgi/bioperl/browse/bioperl-live/trunk/t/SeqIO/?rev=15112">SeqIO</a>      <span class="sep">/</span> metafasta.t    @ r15112</h1>
      </div>  <div id="language-selection">
        <form action="http://code.open-bio.org/svnweb/index.cgi/bioperl/view/bioperl-live/trunk/t/SeqIO/metafasta.t?rev=16733">
          <select name="lang" onchange="this.form.submit();">             <option value="en" selected="yes">English</option>             <option value="fr">Fran&ccedil;ais</option>             <option value="zh_cn">Chinese (Simplified)</option>             <option value="zh_tw">Chinese (Traditional)</option>  </select>   <input type="hidden" name="rev" value="16733" />  <noscript>
            <input type="submit" value="Go" />
          </noscript>
        </form>
      </div>
    </div>

    <div id="content">
                            <div class="actions">
  <ul>  <li><a href="/svnweb/index.cgi/bioperl/blame/bioperl-live/trunk/t/SeqIO/metafasta.t?rev=15112">Blame/Annotate</a></li>  <li><a href="/svnweb/index.cgi/bioperl/checkout/bioperl-live/trunk/t/SeqIO/metafasta.t?rev=15112">Checkout</a></li>  <li><a href="/svnweb/index.cgi/bioperl/log/bioperl-live/trunk/t/SeqIO/metafasta.t?rev=15112">View Revision Log</a></li>  </ul>
</div>         <table id="info" summary="Revision Log">
  <tr>
    <th scope="row"> Revision <a href="/svnweb/index.cgi/bioperl/revision?rev=15112">15112</a> (by sendu, 2008/12/08 18:12:38)</th>
    <td class="message">BioperlTest -&gt; Bio::Root::Test</td>
  </tr>
</table>

<div id="preview">  <pre class="code-block"># -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests =&gt; 6);
	
	use_ok('Bio::SeqIO::metafasta');
}

my $verbose = test_debug();

my $io = Bio::SeqIO-&gt;new(-format =&gt; 'metafasta',
								 -verbose =&gt; $verbose,
								 -file =&gt; test_input_file('test.metafasta'));

isa_ok($io, 'Bio::SeqIO');
ok(my $seq = $io-&gt;next_seq);
isa_ok($seq, 'Bio::Seq::Meta');
is($seq-&gt;seq, &quot;ABCDEFHIJKLMNOPQRSTUVWXYZ&quot;);
is($seq-&gt;display_id,'test');
</pre>  </div>     </div>
    <div id="footer">
      <hr />
      <p class="right"><em>Back to <a href="http://code.open-bio.org">code.open-bio.org</a> | Hosting donated by <a href="http://www.bioteam.net">the BioTeam</a> | <a href="http://search.cpan.org/dist/SVN-Web/">Powered by SVN::Web</a></em></p>
    </div>
  </body>
</html>
 