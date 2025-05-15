<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: header_L1b
Version: 1.0
Date: 23 Dec 2010
-->

<xsl:stylesheet id="header_L1b" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:strip-space elements="*"/>
  <xsl:output name="header_L1b" method="text" encoding="US-ASCII" media-type="text/plain"/>

  <xsl:include href="eef2txt_defaults.xsl"/>

  <xsl:include href="eef2txt_header_product_record.xsl"/>
  <xsl:include href="eef2txt_header_dsd_records_L1b.xsl"/>

</xsl:stylesheet>
