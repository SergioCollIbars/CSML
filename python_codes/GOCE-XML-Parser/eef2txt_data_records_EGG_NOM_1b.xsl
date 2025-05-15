<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: EGG_NOM_1b_data_records
Version: 1.0
Date: 23 Dec 2010
-->

<xsl:stylesheet id="EGG_NOM_1b_data_records" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:variable name="GPS_Time_Field_Length"         select="20"/>
  <xsl:variable name="Gravity_Gradient_Field_Length" select="95"/>
  <xsl:variable name="Flag_Field_Length"             select="17"/>
  <xsl:variable name="Common_Mode_Field_Length"      select="140"/>
  <xsl:variable name="Diff_Mode_Field_Length"        select="140"/>
  <xsl:variable name="Q_Grad_Length"                 select="63"/>

  <xsl:variable name="Empty_GPS_Time_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$GPS_Time_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="Empty_Gravity_Gradient_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Gravity_Gradient_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="Empty_Flag_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Flag_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="Empty_Common_Mode_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Common_Mode_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="Empty_Diff_Mode_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Diff_Mode_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="Empty_Q_Grad_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Q_Grad_Length"/>
    </xsl:call-template>
  </xsl:variable>



  <xsl:template match="*/Data_Block/EGG_GGT_DS">
  
    <xsl:for-each select="EGG_GGT_1i">

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Tt_GPS"/>
        <xsl:with-param name="Empty_Field" select="$Empty_GPS_Time_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Gravity_Grad_Tensor"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Gravity_Gradient_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Qual_Flag"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

  </xsl:template>

  <xsl:template match="*/Data_Block/EGG_CCD_DS">
  
    <xsl:for-each select="EGG_CCD_1i">

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Tt_GPS"/>
        <xsl:with-param name="Empty_Field" select="$Empty_GPS_Time_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>
          
      <xsl:call-template name="Recombine_3_Gradients">
         <xsl:with-param name="XX" select="Acc_Ccm/X"/>
         <xsl:with-param name="YY" select="Acc_Ccm/Y"/>
         <xsl:with-param name="ZZ" select="Acc_Ccm/Z"/>        
      </xsl:call-template>

      <xsl:call-template name="Recombine_3_Gradients">
         <xsl:with-param name="XX" select="Acc_Cdm/X"/>
         <xsl:with-param name="YY" select="Acc_Cdm/Y"/>
         <xsl:with-param name="ZZ" select="Acc_Cdm/Z"/>        
      </xsl:call-template>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

  </xsl:template>


  <xsl:template match="*/Data_Block/EGG_IAQ_DS">
  
    <xsl:for-each select="EGG_IAQ_1i">

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Tt_GPS"/>
        <xsl:with-param name="Empty_Field" select="$Empty_GPS_Time_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corr_Quat/Q_Grad"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Q_Grad_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

  </xsl:template>


</xsl:stylesheet>
