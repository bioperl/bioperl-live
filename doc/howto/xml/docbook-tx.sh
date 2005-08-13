#!/bin/tcsh
# Docbook XML to HTML, PDF, and text using e-novative
# stylesheets, Saxon version 6, Apache FOP, and lynx.

setenv HTML_STYLE '~/Documents/programming/docbook/stylesheet/e-novative_article_html.xsl'
setenv PDF_STYLE '~/Documents/programming/docbook/stylesheet/e-novative_article_fo.xsl'
setenv XML '~/bioperl-live/doc/howto/sgml/Beginners.xml'
setenv HTML 'test.html'
setenv TXT 'test.txt'
setenv PDF 'test.pdf'
setenv FOP_JAR '/usr/local/fop-0.20.5/build/fop.jar:/usr/local/fop-0.20.5/lib/avalon-framework-cvs-20020806.jar:/usr/local/fop-0.20.5/lib/batik.jar'
setenv SAXON_JAR '/usr/local/saxon6-5-4/saxon.jar'
setenv FOP 'org.apache.fop.apps.Fop'
setenv JAVA '/usr/bin/java'
setenv LYNX '/sw/bin/lynx'


# Create HTML file
$JAVA -jar $SAXON_JAR -o $HTML $XML $HTML_STYLE use.extensions=1 

# Create text from HTML file
$LYNX -dump $HTML > $TXT

# Create PDF
$JAVA -classpath $FOP_JAR $FOP -q -xml $XML -xsl $PDF_STYLE -pdf $PDF
