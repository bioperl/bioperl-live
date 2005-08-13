<?xml version="1.0" encoding="iso-8859-1"?>

<!-- >e-novative> DocBook Environment (eDE)                                  -->
<!-- (c) 2002 e-novative GmbH, Munich, Germany                               -->
<!-- http://www.e-novative.de                                                -->

<!-- Single HTML File Generation                                             -->

<!-- This file is part of eDE                                                -->

<!-- eDE is free software; you can redistribute it and/or modify             -->
<!-- it under the terms of the GNU General Public License as published by    -->
<!-- the Free Software Foundation; either version 2 of the License, or       -->
<!-- (at your option) any later version.                                     -->

<!-- eDE is distributed in the hope that it will be useful,                  -->
<!-- but WITHOUT ANY WARRANTY; without even the implied warranty of          -->
<!-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           -->
<!-- GNU General Public License for more details.                            -->

<!-- You should have received a copy of the GNU General Public License       -->
<!-- along with eDe; if not, write to the Free Software                      -->
<!-- Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA -->

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                version="1.0"
                xmlns="http://www.w3.org/TR/xhtml1/transitional"
>

<xsl:import href="../xsl/htmlhelp/profile-htmlhelp.xsl" />
<xsl:import href="e-novative.xsl" />
<xsl:import href="e-novative_article.xsl" />
<xsl:import href="custom.xsl" />
<xsl:import href="custom_article_htmlhelp.xsl" />

<xsl:output method="html" encoding="iso-8859-1" indent="no" />

<!-- make all gui elements appear bold                                       -->
<!-- This adds html bold tags to each gui element (guibutton, guiicon,       -->
<!-- guilabel, guimenu, guimenuitem, guisubmenu)                             -->
<!-- applies to: html output only                                            -->
<!-- Comment the whole block out to make gui elements appear as regular      -->
<!-- text                                                                    -->

<xsl:template match="guibutton">
<b><xsl:call-template name="inline.charseq"/></b>
</xsl:template>

<xsl:template match="guiicon">
<b><xsl:call-template name="inline.charseq"/></b>
</xsl:template>

<xsl:template match="guilabel">
<b><xsl:call-template name="inline.charseq"/></b>
</xsl:template>

<xsl:template match="guimenu">
<b><xsl:call-template name="inline.charseq"/></b>
</xsl:template>

<xsl:template match="guimenuitem">
<b><xsl:call-template name="inline.charseq"/></b>
</xsl:template>

<xsl:template match="guisubmenu">
<b><xsl:call-template name="inline.charseq"/></b>
</xsl:template>

</xsl:stylesheet>
