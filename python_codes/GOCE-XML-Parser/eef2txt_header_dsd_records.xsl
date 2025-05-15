<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: header_dsd_records
Version: 1.0
Date: 20 Jun 2008
-->

<xsl:stylesheet id="header_dsd_records" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:variable name="Data_Set_Name_Field_Length"      select="28"/>
  <xsl:variable name="Data_Set_Type_Field_Length"      select="1"/>
  <xsl:variable name="Data_Set_File_Name_Field_Length" select="62"/>
  <xsl:variable name="Num_Epochs_Field_Length"         select="11"/>
  <xsl:variable name="Start_GPS_Time_Field_Length"     select="20"/>
  <xsl:variable name="Stop_GPS_Time_Field_Length"      select="$Start_GPS_Time_Field_Length"/>
  <xsl:variable name="Phase_Field_Length"              select="$Data_Set_Type_Field_Length"/>
  <xsl:variable name="Abs_Orbit_Start_Field_Length"    select="6"/>
  <xsl:variable name="Abs_Orbit_Stop_Field_Length"     select="$Abs_Orbit_Start_Field_Length"/>

  <xsl:variable name="Empty_Data_Set_Name_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Data_Set_Name_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Data_Set_Type_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Data_Set_Type_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Data_Set_File_Name_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Data_Set_File_Name_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Num_Epochs_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Num_Epochs_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Start_GPS_Time_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Start_GPS_Time_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Stop_GPS_Time_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Stop_GPS_Time_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Phase_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Phase_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Abs_Orbit_Start_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Abs_Orbit_Start_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Abs_Orbit_Stop_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Abs_Orbit_Stop_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:template match="*/Earth_Explorer_Header/Variable_Header/SPH">
  
    <xsl:variable name="Output_Data_Set_Type" select="'O'"/>

    <xsl:variable name="Unknown_Num_Epochs"      select="'0'"/>
    <xsl:variable name="Unknown_Start_GPS_Time"  select="'0000000000.000000000'"/>
    <xsl:variable name="Unknown_Stop_GPS_Time"   select="'9999999999.999999999'"/>
    <xsl:variable name="Unknown_Phase"           select="'X'"/>
    <xsl:variable name="Unknown_Abs_Orbit_Start" select="$Unknown_Num_Epochs"/>
    <xsl:variable name="Unknown_Abs_Orbit_Stop"  select="$Unknown_Abs_Orbit_Start"/>

    <xsl:variable name="Time" select="Time_Information"/>
    <xsl:variable name="DSDs" select="DSDs/List_of_DSDs"/>

    <xsl:variable name="Start_GPS_Time" select="normalize-space($Time/GPS_Time/Start)"/>
    <xsl:variable name="Stop_GPS_Time" select="normalize-space($Time/GPS_Time/Stop)"/>
    <xsl:variable name="Abs_Orbit_Start" select="normalize-space($Time/Abs_Orbit/Start)"/>
    <xsl:variable name="Abs_Orbit_Stop"  select="normalize-space($Time/Abs_Orbit/Stop)"/>

    <xsl:for-each select="$DSDs/Data_Set_Descriptor">
    
      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Data_Set_Name"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Data_Set_Name_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$HFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Data_Set_Type"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Data_Set_Type_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$HFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="File_Name"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Data_Set_File_Name_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$HFS"/>

      <xsl:choose>
        <xsl:when test="not(normalize-space(Num_Epochs)='')">
          <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="Num_Epochs"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Num_Epochs_Field"/>
            <xsl:with-param name="Justify"     select="$Right_Justified"/>
          </xsl:call-template>
        </xsl:when>
        <xsl:otherwise>
          <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="$Unknown_Num_Epochs"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Num_Epochs_Field"/>
            <xsl:with-param name="Justify"     select="$Right_Justified"/>
          </xsl:call-template>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:value-of select="$HFS"/>

      <xsl:choose>
        <xsl:when test="Data_Set_Type=$Output_Data_Set_Type and not($Start_GPS_Time='')">
          <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="$Start_GPS_Time"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Start_GPS_Time_Field"/>
            <xsl:with-param name="Justify"     select="$Right_Justified"/>
          </xsl:call-template>
        </xsl:when>
	<xsl:otherwise>
          <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="$Unknown_Start_GPS_Time"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Start_GPS_Time_Field"/>
            <xsl:with-param name="Justify"     select="$Right_Justified"/>
          </xsl:call-template>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:value-of select="$HFS"/>
      
      <xsl:choose>
        <xsl:when test="Data_Set_Type=$Output_Data_Set_Type and not($Stop_GPS_Time='')">
	  <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="$Stop_GPS_Time"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Stop_GPS_Time_Field"/>
            <xsl:with-param name="Justify"     select="$Right_Justified"/>
          </xsl:call-template>
        </xsl:when>
	<xsl:otherwise>
	  <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="$Unknown_Stop_GPS_Time"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Stop_GPS_Time_Field"/>
            <xsl:with-param name="Justify"     select="$Right_Justified"/>
          </xsl:call-template>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:value-of select="$HFS"/>

      <xsl:choose>
        <xsl:when test="not(normalize-space(Phase)='')">
          <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="Phase"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Phase_Field"/>
            <xsl:with-param name="Justify"     select="$Left_Justified"/>
          </xsl:call-template>
        </xsl:when>
        <xsl:otherwise>
          <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="$Unknown_Phase"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Phase_Field"/>
            <xsl:with-param name="Justify"     select="$Left_Justified"/>
          </xsl:call-template>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:value-of select="$HFS"/>

      <xsl:choose>
        <xsl:when test="Data_Set_Type=$Output_Data_Set_Type">

          <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="$Abs_Orbit_Start"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Abs_Orbit_Start_Field"/>
            <xsl:with-param name="Justify"     select="$Right_Justified"/>
          </xsl:call-template>

          <xsl:value-of select="$HFS"/>

	  <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="$Abs_Orbit_Stop"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Abs_Orbit_Stop_Field"/>
            <xsl:with-param name="Justify"     select="$Right_Justified"/>
          </xsl:call-template>

        </xsl:when>
        <xsl:otherwise>

	  <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="$Unknown_Abs_Orbit_Start"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Abs_Orbit_Start_Field"/>
            <xsl:with-param name="Justify"     select="$Right_Justified"/>
          </xsl:call-template>

          <xsl:value-of select="$HFS"/>

	  <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="$Unknown_Abs_Orbit_Stop"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Abs_Orbit_Stop_Field"/>
            <xsl:with-param name="Justify"     select="$Right_Justified"/>
          </xsl:call-template>

        </xsl:otherwise>
      </xsl:choose>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

  </xsl:template>

</xsl:stylesheet>
