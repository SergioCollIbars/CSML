<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: SST_PSO_2_data_block
Version: 1.0
Date: 20 Jun 2008
-->

<xsl:stylesheet id="SST_PSO_2_data_block" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:strip-space elements="*"/>
  <xsl:output name="SST_PSO_2" method="text" encoding="US-ASCII" media-type="text/plain"/>

  <xsl:include href="eef2txt_defaults.xsl"/>

  <xsl:variable name="Comment_Indicator" select="'#'"/>
  <xsl:variable name="Value_Prefix"      select="' '"/>

  <xsl:variable name="Header_Keyword_Field_Length" select="33"/>
  <xsl:variable name="Symbol_Field_Length"         select="2"/>

  <xsl:variable name="Year_Start_Field_Length"         select="4"/>
  <xsl:variable name="Month_Start_Field_Length"        select="$Symbol_Field_Length"/>
  <xsl:variable name="Day_Of_Month_Start_Field_Length" select="$Symbol_Field_Length"/>
  <xsl:variable name="Hour_Start_Field_Length"         select="$Symbol_Field_Length"/>
  <xsl:variable name="Minute_Start_Field_Length"       select="$Symbol_Field_Length"/>
  <xsl:variable name="Second_Start_Field_Length"       select="11"/>

  <xsl:variable name="Empty_Header_Keyword_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Header_Keyword_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Symbol_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Symbol_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="Empty_Year_Start_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Year_Start_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Month_Start_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Month_Start_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Day_Of_Month_Start_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Day_Of_Month_Start_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Hour_Start_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Hour_Start_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Minute_Start_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Minute_Start_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Second_Start_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Second_Start_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:include href="eef2txt_data_header_SST_PSO_2.xsl"/>
  <xsl:include href="eef2txt_data_records_SST_PSO_2.xsl"/>

</xsl:stylesheet>
