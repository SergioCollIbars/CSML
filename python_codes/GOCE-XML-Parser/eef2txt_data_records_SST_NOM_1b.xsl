<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: SST_NOM_1b_data_records
Version: 1.0
Date: 23 Dec 2010
-->

<xsl:stylesheet id="SST_NOM_1b_data_records" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:variable name="GPS_Time_Field_Length"      select="20"/>
  <xsl:variable name="Cov_Pt_Field_Length"        select="71"/>
  <xsl:variable name="Cov_V_Field_Length"         select="53"/>
  <xsl:variable name="SST_Pos_Field_Length"       select="53"/>
  <xsl:variable name="Rec_Clock_Err_Field_Length" select="15"/>   

  <xsl:variable name="Empty_GPS_Time_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$GPS_Time_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="Empty_Cov_Pt_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Cov_Pt_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="Empty_Cov_V_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Cov_V_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="Empty_SST_Pos_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$SST_Pos_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="Empty_Rec_Clock_Err_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Rec_Clock_Err_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>



  <xsl:template match="*/Data_Block/SST_COV_DS">
  
    <xsl:for-each select="SST_COV_1i">

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Tt_GPS"/>
        <xsl:with-param name="Empty_Field" select="$Empty_GPS_Time_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="SST_Cov_Pt/Row1"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Cov_Pt_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="SST_Cov_Pt/Row2"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Cov_Pt_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="SST_Cov_Pt/Row3"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Cov_Pt_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="SST_Cov_Pt/Row4"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Cov_Pt_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="SST_Cov_V/Row1"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Cov_V_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="SST_Cov_V/Row2"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Cov_V_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>
  
      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="SST_Cov_V/Row3"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Cov_V_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

  </xsl:template>



  <xsl:template match="*/Data_Block/SST_PVT_DS">
  
    <xsl:for-each select="SST_PVT_1i">

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Tt_GPS"/>
        <xsl:with-param name="Empty_Field" select="$Empty_GPS_Time_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="SST_Pos/Position"/>
        <xsl:with-param name="Empty_Field" select="$Empty_SST_Pos_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="SST_Pos/Rec_Clock_Err"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Rec_Clock_Err_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="SST_Vel"/>
        <xsl:with-param name="Empty_Field" select="$Empty_SST_Pos_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

  </xsl:template>




</xsl:stylesheet>
