<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: EGG_NOM_1b_data_block
Version: 1.0
Date: 23 Dec 2010
-->

<xsl:stylesheet id="EGG_NOM_1b_data_block" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:strip-space elements="*"/>
  <xsl:output name="EGG_NOM_1b" method="text" encoding="US-ASCII" media-type="text/plain"/>

  <xsl:include href="eef2txt_defaults.xsl"/>

  <xsl:variable name="Symbol_Field_Length"         select="20"/>
  <xsl:variable name="Empty_Symbol_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Symbol_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>


  <xsl:include href="eef2txt_data_records_EGG_NOM_1b.xsl"/>

</xsl:stylesheet>
