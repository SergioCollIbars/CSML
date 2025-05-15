<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: header_product_record
Version: 1.0
Date: 20 Jun 2008
-->

<xsl:stylesheet id="header_product_record" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:variable name="File_Name_Field_Length"        select="62"/>
  <xsl:variable name="File_Description_Field_Length" select="100"/>
  <xsl:variable name="Mission_Field_Length"          select="4"/>
  <xsl:variable name="File_Class_Field_Length"       select="$Mission_Field_Length"/>
  <xsl:variable name="File_Type_Field_Length"        select="10"/>
  <xsl:variable name="Validity_Start_Field_Length"   select="23"/>
  <xsl:variable name="Validity_Stop_Field_Length"    select="$Validity_Start_Field_Length"/>
  <xsl:variable name="File_Version_Field_Length"     select="$Mission_Field_Length"/>
  <xsl:variable name="System_Field_Length"           select="3"/>
  <xsl:variable name="Creator_Field_Length"          select="12"/>
  <xsl:variable name="Creator_Version_Field_Length"  select="5"/>
  <xsl:variable name="Creation_Date_Field_Length"    select="$Validity_Start_Field_Length"/>
  <xsl:variable name="Ref_Doc_Field_Length"          select="$Validity_Start_Field_Length"/>
  <xsl:variable name="Num_DSD_Field_Length"          select="11"/>

  <xsl:variable name="Empty_File_Name_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$File_Name_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_File_Description_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$File_Description_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Mission_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Mission_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_File_Class_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$File_Class_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_File_Type_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length"    select="$File_Type_Field_Length"/>
      <xsl:with-param name="Empty_Character" select="'_'"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Validity_Start_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Validity_Start_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Validity_Stop_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Validity_Stop_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_File_Version_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length"    select="$File_Version_Field_Length"/>
      <xsl:with-param name="Empty_Character" select="'0'"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_System_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$System_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Creator_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Creator_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Creator_Version_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Creator_Version_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Creation_Date_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Creation_Date_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Ref_Doc_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Ref_Doc_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Num_DSD_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Num_DSD_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:template match="*/Earth_Explorer_Header/Fixed_Header">

    <xsl:variable name="Unknown_Validity_Start" select="'UTC=0000-00-00T00:00:00'"/>
    <xsl:variable name="Unknown_Validity_Stop"  select="'UTC=9999-99-99T99:99:99'"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="File_Name"/>
      <xsl:with-param name="Empty_Field" select="$Empty_File_Name_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="File_Description"/>
      <xsl:with-param name="Empty_Field" select="$Empty_File_Description_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Mission"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Mission_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="File_Class"/>
      <xsl:with-param name="Empty_Field" select="$Empty_File_Class_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="File_Type"/>
      <xsl:with-param name="Empty_Field" select="$Empty_File_Type_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

    <xsl:choose>
      <xsl:when test="not(normalize-space(Validity_Period/Validity_Start)='')">
        <xsl:call-template name="Fixed_Length_Field">
          <xsl:with-param name="Field_Value" select="Validity_Period/Validity_Start"/>
          <xsl:with-param name="Empty_Field" select="$Empty_Validity_Start_Field"/>
          <xsl:with-param name="Justify"     select="$Left_Justified"/>
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <xsl:call-template name="Fixed_Length_Field">
          <xsl:with-param name="Field_Value" select="$Unknown_Validity_Start"/>
          <xsl:with-param name="Empty_Field" select="$Empty_Validity_Start_Field"/>
          <xsl:with-param name="Justify"     select="$Left_Justified"/>
        </xsl:call-template>
      </xsl:otherwise>
    </xsl:choose>

    <xsl:value-of select="$HFS"/>

    <xsl:choose>
      <xsl:when test="not(normalize-space(Validity_Period/Validity_Stop)='')">
        <xsl:call-template name="Fixed_Length_Field">
          <xsl:with-param name="Field_Value" select="Validity_Period/Validity_Stop"/>
          <xsl:with-param name="Empty_Field" select="$Empty_Validity_Stop_Field"/>
          <xsl:with-param name="Justify"     select="$Left_Justified"/>
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <xsl:call-template name="Fixed_Length_Field">
          <xsl:with-param name="Field_Value" select="$Unknown_Validity_Stop"/>
          <xsl:with-param name="Empty_Field" select="$Empty_Validity_Stop_Field"/>
          <xsl:with-param name="Justify"     select="$Left_Justified"/>
        </xsl:call-template>
      </xsl:otherwise>
    </xsl:choose>

    <xsl:value-of select="$HFS"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="File_Version"/>
      <xsl:with-param name="Empty_Field" select="$Empty_File_Version_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

    <xsl:variable name="Tmp_System" select="substring-after(Source/System,'(')"/>
    <xsl:variable name="System"     select="substring-before($Tmp_System,')')"/>
    <xsl:choose>
      <xsl:when test="not($System='')">
        <xsl:call-template name="Fixed_Length_Field">
          <xsl:with-param name="Field_Value" select="$System"/>
          <xsl:with-param name="Empty_Field" select="$Empty_System_Field"/>
          <xsl:with-param name="Justify"     select="$Left_Justified"/>
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <xsl:call-template name="Fixed_Length_Field">
          <xsl:with-param name="Field_Value" select="Source/System"/>
          <xsl:with-param name="Empty_Field" select="$Empty_System_Field"/>
          <xsl:with-param name="Justify"     select="$Left_Justified"/>
        </xsl:call-template>
      </xsl:otherwise>
    </xsl:choose>

    <xsl:value-of select="$HFS"/>

    <xsl:variable name="Tmp_Creator" select="substring-after(Source/Creator,'(')"/>
    <xsl:variable name="Creator"     select="substring-before($Tmp_Creator,')')"/>
    <xsl:choose>
      <xsl:when test="not($Creator='')">
        <xsl:call-template name="Fixed_Length_Field">
          <xsl:with-param name="Field_Value" select="$Creator"/>
          <xsl:with-param name="Empty_Field" select="$Empty_Creator_Field"/>
          <xsl:with-param name="Justify"     select="$Left_Justified"/>
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <xsl:call-template name="Fixed_Length_Field">
          <xsl:with-param name="Field_Value" select="Source/Creator"/>
          <xsl:with-param name="Empty_Field" select="$Empty_Creator_Field"/>
          <xsl:with-param name="Justify"     select="$Left_Justified"/>
        </xsl:call-template>
      </xsl:otherwise>
    </xsl:choose>

    <xsl:value-of select="$HFS"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Source/Creator_Version"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Creator_Version_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Source/Creation_Date"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Creation_Date_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

  </xsl:template>

  <xsl:template match="*/Earth_Explorer_Header/Variable_Header/MPH">

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Ref_Doc"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Ref_Doc_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$HFS"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Num_DSD"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Num_DSD_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

  </xsl:template>

</xsl:stylesheet>
