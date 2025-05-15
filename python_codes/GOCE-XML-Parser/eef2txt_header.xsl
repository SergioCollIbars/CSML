<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: header
Version: 1.0
Date: 20 Jun 2008
-->

<xsl:stylesheet id="header" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:strip-space elements="*"/>
  <xsl:output name="header" method="text" encoding="US-ASCII" media-type="text/plain"/>

  <xsl:include href="eef2txt_defaults.xsl"/>

  <xsl:include href="eef2txt_header_product_record.xsl"/>
  <xsl:include href="eef2txt_header_dsd_records.xsl"/>

</xsl:stylesheet>
