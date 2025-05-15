<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: EGG_TRF_2_data_records
Version: 1.1
Date: 12 Nov 2010
-->

<xsl:stylesheet id="EGG_TRF_2_data_records" version="1.1" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:variable name="GPS_Time_Field_Length"         select="20"/>
  <xsl:variable name="Phi_Field_Length"              select="15"/>
  <xsl:variable name="Lambda_Field_Length"           select="$Phi_Field_Length"/>
  <xsl:variable name="Radius_Field_Length"           select="13"/>
  <xsl:variable name="Gravity_Gradient_Field_Length" select="$Phi_Field_Length"/>
  <xsl:variable name="Sigma_Field_Length"            select="$Phi_Field_Length"/>
  <xsl:variable name="Flag_Field_Length"             select="1"/>

  <xsl:variable name="Empty_GPS_Time_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$GPS_Time_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Phi_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Phi_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Lambda_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Lambda_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Radius_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Radius_Field_Length"/>
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

  <xsl:template match="*/Data_Block/EGG_TRF_2/List_of_GG_spatial_Records">
  
    <xsl:for-each select="GG_spatial_Record">
    
      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Time_Information/GPS_Time"/>
        <xsl:with-param name="Empty_Field" select="$Empty_GPS_Time_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Position/Phi"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Phi_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Position/Lambda"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Lambda_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Position/Radius_from_Geocenter"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Radius_Field"/>
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

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

  </xsl:template>

</xsl:stylesheet>
