<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: SST_AUX_2_data_header
Version: 1.1
Date: 14 Jul 2008
-->

<xsl:stylesheet id="SST_AUX_2_data_header" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template match="*/Earth_Explorer_Header/Variable_Header/SPH">

    <xsl:if test="not(starts-with($Mode, 'data_block'))">
      <xsl:apply-templates select="SST_AUX_2" mode="header"/>
    </xsl:if>

  </xsl:template>

  <xsl:template match="*/Earth_Explorer_Header/Variable_Header/SPH/SST_AUX_2" mode="header">

    <xsl:variable name="Product_Type_Keyword"            select="'product_type'"/>
    <xsl:variable name="Model_Name_Keyword"              select="'modelname'"/>
    <xsl:variable name="Earth_Gravity_Constant_Keyword"  select="'earth_gravity_constant'"/>
    <xsl:variable name="Radius_Keyword"                  select="'radius'"/>
    <xsl:variable name="Max_Degree_Keyword"              select="'max_degree'"/>
    <xsl:variable name="Errors_Keyword"                  select="'errors'"/>
    <xsl:variable name="Norm_Keyword"                    select="'norm'"/>
    <xsl:variable name="Tide_System_Keyword"             select="'tide_system'"/>
    <xsl:variable name="Water_Density_Keyword"           select="'water_density'"/>
    <xsl:variable name="End_Of_SST_AUX_2_Header_Keyword" select="'end_of_head'"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Product_Type_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="Product_Type"/>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Model_Name_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="Model_Name"/>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Earth_Gravity_Constant_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="Earth_Gravity_Constant"/>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Radius_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="Spherical_Harmonic_Development/Radius"/>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Max_Degree_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="Spherical_Harmonic_Development/Max_Degree"/>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Errors_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="Errors"/>

    <xsl:text>&#xa;</xsl:text>

    <xsl:if test="Normalization">
    
      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="$Norm_Keyword"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="Normalization"/>

      <xsl:text>&#xa;</xsl:text>

    </xsl:if>

    <xsl:if test="Tide_System">

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="$Tide_System_Keyword"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="Tide_System"/>

      <xsl:text>&#xa;</xsl:text>

    </xsl:if>

    <xsl:if test="Water_Density">

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="$Water_Density_Keyword"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="Water_Density"/>

      <xsl:text>&#xa;</xsl:text>
    
    </xsl:if>

    <xsl:value-of select="$End_Of_SST_AUX_2_Header_Keyword"/>

    <xsl:text>&#xa;</xsl:text>

  </xsl:template>

</xsl:stylesheet>
