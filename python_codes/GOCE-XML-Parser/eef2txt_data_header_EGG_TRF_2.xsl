<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: EGG_TRF_2_data_header
Version: 1.1
Date: 14 Jul 2008
-->

<xsl:stylesheet id="EGG_TRF_2_data_header" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:variable name="Product_Type_Field_Length"     select="10"/>
  <xsl:variable name="L1_Input_Field_Length"         select="12"/>
  <xsl:variable name="L2_Input_Field_Length"         select="$Product_Type_Field_Length"/>
  <xsl:variable name="Reference_System_Field_Length" select="4"/>
  <xsl:variable name="Tide_System_Field_Length"      select="9"/>
  <xsl:variable name="Gravity_Model_Field_Length"    select="64"/>
  <xsl:variable name="Errors_Field_Length"           select="$Product_Type_Field_Length"/>

  <xsl:variable name="Empty_Product_Type_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Product_Type_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_L1_Input_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$L1_Input_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_L2_Input_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$L2_Input_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Reference_System_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Reference_System_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Tide_System_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Tide_System_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Gravity_Model_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Gravity_Model_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Errors_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Errors_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:template match="*/Earth_Explorer_Header/Variable_Header/SPH">

    <xsl:if test="not(starts-with($Mode, 'data_block'))">
      <xsl:apply-templates select="EGG_TRF_2" mode="header"/>
    </xsl:if>

  </xsl:template>

  <xsl:template match="*/Earth_Explorer_Header/Variable_Header/SPH/EGG_TRF_2" mode="header">
  
    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Product_Type"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Product_Type_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Input/L1"/>
      <xsl:with-param name="Empty_Field" select="$Empty_L1_Input_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Input/L2"/>
      <xsl:with-param name="Empty_Field" select="$Empty_L2_Input_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Reference_System"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Reference_System_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Tide_System"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Tide_System_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Gravity_Model"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Gravity_Model_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Errors"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Errors_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

  </xsl:template>

</xsl:stylesheet>
