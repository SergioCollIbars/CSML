<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: EGM_GOC_2_data_block
Version: 1.0
Date: 20 Jun 2008
-->

<xsl:stylesheet id="EGM_GOC_2_data_block" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:strip-space elements="*"/>
  <xsl:output name="EGM_GOC_2" method="text" encoding="US-ASCII" media-type="text/plain"/>

  <xsl:include href="eef2txt_defaults.xsl"/>

  <xsl:variable name="Value_Prefix" select="' '"/>

  <xsl:variable name="Header_Keyword_Field_Length"   select="30"/>
  <xsl:variable name="Number_Of_Values_Field_Length" select="7"/>

  <xsl:variable name="Empty_Header_Keyword_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Header_Keyword_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Number_Of_Values_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Number_Of_Values_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="EGM_GEO_2_Data_Format" select="*/Earth_Explorer_Header/Variable_Header/SPH/EGM_GOC_2/EGM_GEO_2/Original_Source/Format/Fortran_Notation"/>
  <xsl:variable name="EGM_GAN_2_Data_Format" select="*/Earth_Explorer_Header/Variable_Header/SPH/EGM_GOC_2/EGM_GAN_2/Original_Source/Format/Fortran_Notation"/>
  <xsl:variable name="EGM_GVE_2_Data_Format" select="*/Earth_Explorer_Header/Variable_Header/SPH/EGM_GOC_2/EGM_GVE_2/Original_Source/Format/Fortran_Notation"/>
  <xsl:variable name="EGM_GVN_2_Data_Format" select="*/Earth_Explorer_Header/Variable_Header/SPH/EGM_GOC_2/EGM_GVN_2/Original_Source/Format/Fortran_Notation"/>
  <xsl:variable name="EGM_GER_2_Data_Format" select="*/Earth_Explorer_Header/Variable_Header/SPH/EGM_GOC_2/EGM_GER_2/Original_Source/Format/Fortran_Notation"/>
      
  <xsl:include href="eef2txt_data_header_EGM_GOC_2.xsl"/>
  <xsl:include href="eef2txt_data_records_EGM_GOC_2.xsl"/>

</xsl:stylesheet>
