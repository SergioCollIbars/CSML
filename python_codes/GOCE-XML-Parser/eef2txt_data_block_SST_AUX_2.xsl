<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: SST_AUX_2_data_block
Version: 1.0
Date: 20 Jun 2008
-->

<xsl:stylesheet id="SST_AUX_2_data_block" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:strip-space elements="*"/>
  <xsl:output name="SST_AUX_2" method="text" encoding="US-ASCII" media-type="text/plain"/>

  <xsl:include href="eef2txt_defaults.xsl"/>

  <xsl:variable name="Value_Prefix" select="' '"/>

  <xsl:variable name="Header_Keyword_Field_Length" select="30"/>

  <xsl:variable name="Empty_Header_Keyword_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Header_Keyword_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
      
  <xsl:include href="eef2txt_data_header_SST_AUX_2.xsl"/>
  <xsl:include href="eef2txt_data_records_SST_AUX_2.xsl"/>

</xsl:stylesheet>
