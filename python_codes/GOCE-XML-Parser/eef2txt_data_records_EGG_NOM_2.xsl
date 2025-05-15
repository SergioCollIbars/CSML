<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: EGG_NOM_2_data_records
Version: 1.0
Date: 20 Jun 2008
-->

<xsl:stylesheet id="EGG_NOM_2_data_records" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:variable name="GPS_Time_Field_Length"         select="20"/>
  <xsl:variable name="Gravity_Gradient_Field_Length" select="15"/>
  <xsl:variable name="Sigma_Field_Length"            select="$Gravity_Gradient_Field_Length"/>
  <xsl:variable name="Flag_Field_Length"             select="1"/>
  <xsl:variable name="Correction_Field_Length"       select="$Gravity_Gradient_Field_Length"/>

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
  <xsl:variable name="Empty_Sigma_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Sigma_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Flag_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Flag_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Correction_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Correction_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:template match="*/Data_Block/EGG_NOM_2/List_of_GG_time_Records">
  
    <xsl:for-each select="GG_time_Record">
    
      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Time_Information/GPS_Time"/>
        <xsl:with-param name="Empty_Field" select="$Empty_GPS_Time_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Gravity_Gradients/XX"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Gravity_Gradient_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Gravity_Gradients/YY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Gravity_Gradient_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Gravity_Gradients/ZZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Gravity_Gradient_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Gravity_Gradients/XY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Gravity_Gradient_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>
      
      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Gravity_Gradients/XZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Gravity_Gradient_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Gravity_Gradients/YZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Gravity_Gradient_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Sigmas/XX"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Sigma_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Sigmas/YY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Sigma_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Sigmas/ZZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Sigma_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>
      
      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Sigmas/XY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Sigma_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Sigmas/XZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Sigma_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Sigmas/YZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Sigma_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Flags/XX"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Flags/YY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Flags/ZZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>
      
      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Flags/XY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Flags/XZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Flags/YZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Direct_Tides/XX"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Direct_Tides/YY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Direct_Tides/ZZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>
      
      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Direct_Tides/XY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Direct_Tides/XZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Direct_Tides/YZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Solid_Earth/XX"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Solid_Earth/YY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Solid_Earth/ZZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>
      
      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Solid_Earth/XY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Solid_Earth/XZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Solid_Earth/YZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Ocean_Tides/XX"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Ocean_Tides/YY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Ocean_Tides/ZZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>
      
      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Ocean_Tides/XY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Ocean_Tides/XZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Ocean_Tides/YZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Pole_Tides/XX"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Pole_Tides/YY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Pole_Tides/ZZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>
      
      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Pole_Tides/XY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Pole_Tides/XZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Tidal/Pole_Tides/YZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Non-Tidal/XX"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Non-Tidal/YY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Non-Tidal/ZZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>
      
      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Non-Tidal/XY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Non-Tidal/XZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Temporal/Non-Tidal/YZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Calibration/XX"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Calibration/YY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Calibration/ZZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>
      
      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Calibration/XY"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Calibration/XZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Corrections/Calibration/YZ"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Quaternions/Q1"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Quaternions/Q2"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Quaternions/Q3"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Quaternions/Q4"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Correction_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

  </xsl:template>

</xsl:stylesheet>
