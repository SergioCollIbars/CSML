<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: STR_VC__1b_data_records
Version: 1.0
Date: 23 Dec 2010
-->

<xsl:stylesheet id="STR_VC__1b_data_records" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:variable name="GPS_Time_Field_Length"         select="20"/>
  <xsl:variable name="Str_Attitude_Field_Length"     select="15"/>
  <xsl:variable name="Flag_Field_Length"             select="1"/>
  <xsl:variable name="Number3_Field_Length"          select="3"/>

  <xsl:variable name="Empty_GPS_Time_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$GPS_Time_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="Empty_Str_Attitude_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Str_Attitude_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="Empty_Flag_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Flag_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="Empty_Number3_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Number3_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:template match="*/Data_Block/*['STR_VC2_DS' or 'STR_VC3_DS']">
  
    <xsl:for-each select="*['STR_VC2_1i' or 'STR_VC3_1i']">

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Tt_GPS_Str"/>
        <xsl:with-param name="Empty_Field" select="$Empty_GPS_Time_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Str_Attitude/Q1"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Str_Attitude_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Str_Attitude/Q2"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Str_Attitude_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

       <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Str_Attitude/Q3"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Str_Attitude_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

       <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Str_Attitude/Q4"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Str_Attitude_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

     <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Stid"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Cid"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Val_Flag"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Loc_Time"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Bbo_Flag"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Trs_Flag"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Temp_Out_Range_Flag"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Asc_Tc_Flag"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Orb_Cor_Flag"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Seq_Flag"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Est_Conf"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Number3_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Nr_Of_Locks"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Number3_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Nr_Of_Stars"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Number3_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

  </xsl:template>



</xsl:stylesheet>
