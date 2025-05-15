<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: EGM_GOC_2_data_records
Version: 1.1
Date: 14 Jul 2008
-->

<xsl:stylesheet id="EGM_GOC_2_data_records" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template match="*/Data_Block/EGM_GOC_2/*[not(local-name()='EGM_GCF_2' or local-name()='EGM_GRP_2')]/List_of_Latitudes">

    <xsl:variable name="Sub_Product" select="local-name(..)"/>

    <xsl:variable name="Latitude_Value_Keyword"              select="'latitude'"/>
    <xsl:variable name="Number_Of_Longitude_Values_Keyword"  select="'number_of_data_values'"/>
    <xsl:variable name="Grid_Format_Default_Max_Line_Length" select="99"/>

    <xsl:variable name="Grid_Possible_Fortran_Format_1" select="'8(F12.8,X)'"/>
    <xsl:variable name="Grid_Possible_Fortran_Format_2" select="'10(F8.3,X)'"/>
    <xsl:variable name="Grid_Possible_Fortran_Format_3" select="'10(F7.4,X)'"/>

    <xsl:variable name="Data_Format">
      <xsl:choose>
        <xsl:when test="$Sub_Product = 'EGM_GEO_2'"><xsl:value-of select="translate($EGM_GEO_2_Data_Format, $LC, $UC)"/></xsl:when>
        <xsl:when test="$Sub_Product = 'EGM_GAN_2'"><xsl:value-of select="translate($EGM_GAN_2_Data_Format, $LC, $UC)"/></xsl:when>
        <xsl:when test="$Sub_Product = 'EGM_GVE_2'"><xsl:value-of select="translate($EGM_GVE_2_Data_Format, $LC, $UC)"/></xsl:when>
        <xsl:when test="$Sub_Product = 'EGM_GVN_2'"><xsl:value-of select="translate($EGM_GVN_2_Data_Format, $LC, $UC)"/></xsl:when>
        <xsl:when test="$Sub_Product = 'EGM_GER_2'"><xsl:value-of select="translate($EGM_GER_2_Data_Format, $LC, $UC)"/></xsl:when>
	<xsl:otherwise>UNKNOWN</xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <xsl:variable name="Total_Number_Of_Fields_Per_Line">
      <xsl:choose>
        <xsl:when test="$Data_Format = $Grid_Possible_Fortran_Format_1">8</xsl:when>
        <xsl:when test="$Data_Format = $Grid_Possible_Fortran_Format_2 or $Data_Format = $Grid_Possible_Fortran_Format_3">10</xsl:when>
        <xsl:otherwise>0</xsl:otherwise>
      </xsl:choose>
    </xsl:variable>
      
    <xsl:variable name="Field_Length">
      <xsl:choose>
        <xsl:when test="$Data_Format = $Grid_Possible_Fortran_Format_1">12</xsl:when>
        <xsl:when test="$Data_Format = $Grid_Possible_Fortran_Format_2">8</xsl:when>
        <xsl:when test="$Data_Format = $Grid_Possible_Fortran_Format_3">7</xsl:when>
        <xsl:otherwise>0</xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <xsl:variable name="Empty_Field">
      <xsl:if test="$Total_Number_Of_Fields_Per_Line &gt; 0 and $Field_Length &gt; 0">
        <xsl:call-template name="Construct_Empty_Field">
          <xsl:with-param name="Field_Length" select="$Field_Length"/>
        </xsl:call-template>
      </xsl:if>
    </xsl:variable>

    <xsl:for-each select="Latitude">

      <xsl:call-template name="Row_Leader_Record">
        <xsl:with-param name="Latitude_Value"                     select="@degree"/>
        <xsl:with-param name="Number_Of_Longitude_Values"         select="List_of_Longitudes/@count"/>
        <xsl:with-param name="Latitude_Value_Keyword"             select="$Latitude_Value_Keyword"/>
        <xsl:with-param name="Number_Of_Longitude_Values_Keyword" select="$Number_Of_Longitude_Values_Keyword"/>
      </xsl:call-template>              

      <xsl:text>&#xa;</xsl:text>

      <xsl:choose>
        <xsl:when test="$Total_Number_Of_Fields_Per_Line &gt; 0 and $Field_Length &gt; 0">
          <xsl:call-template name="Combine_Records_To_Number_Of_Fields">
            <xsl:with-param name="Total_Number_Of_Fields_Per_Line" select="$Total_Number_Of_Fields_Per_Line"/>
            <xsl:with-param name="Empty_Field"                     select="$Empty_Field"/>
            <xsl:with-param name="List"                            select="List_of_Longitudes/Longitude/Value"/>
          </xsl:call-template>
	</xsl:when>
        <xsl:otherwise>
          <xsl:call-template name="Combine_Records_To_Max_Length">
            <xsl:with-param name="Max_Length" select="$Grid_Format_Default_Max_Line_Length"/>
            <xsl:with-param name="List"       select="List_of_Longitudes/Longitude/Value"/>
          </xsl:call-template>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

  </xsl:template>

  <xsl:template name="Row_Leader_Record">
    
    <xsl:param name="Latitude_Value"/>
    <xsl:param name="Number_Of_Longitude_Values"/>
    <xsl:param name="Latitude_Value_Keyword"/>
    <xsl:param name="Number_Of_Longitude_Values_Keyword"/>

    <xsl:variable name="Max_Length_Before_Decimal_Point">
      <xsl:choose>
        <xsl:when test="starts-with($Latitude_Value,'-') or starts-with($Latitude_Value,'+')">4</xsl:when>
        <xsl:otherwise>3</xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Latitude_Value_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:call-template name="Justify_On_Decimal_Point">
      <xsl:with-param name="Value"             select="$Latitude_Value"/>
      <xsl:with-param name="Justify_Character" select="$Value_Prefix"/>
      <xsl:with-param name="Max_Length"        select="$Max_Length_Before_Decimal_Point"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Number_Of_Longitude_Values_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value"        select="$Number_Of_Longitude_Values"/>
      <xsl:with-param name="Empty_Field"        select="$Empty_Number_Of_Values_Field"/>
      <xsl:with-param name="Justify"            select="$Right_Justified"/>
      <xsl:with-param name="Truncate_If_Longer" select="'no'"/>
    </xsl:call-template>

  </xsl:template>
  
  <xsl:template match="*/Data_Block/EGM_GOC_2/EGM_GCF_2/List_of_ICGEM_Records">

    <xsl:variable name="ICGEM_Keyword_Field_Length" select="5"/>
    <xsl:variable name="Degree_Field_Length"        select="4"/>
    <xsl:variable name="Order_Field_Length"         select="$Degree_Field_Length"/>

    <xsl:variable name="Empty_ICGEM_Keyword_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$ICGEM_Keyword_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Degree_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Degree_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Order_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Order_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:for-each select="ICGEM_Record">
      
      <xsl:variable name="ICGEM_Keyword" select="translate(@type, $UC, $LC)"/>

      <xsl:choose>
        <xsl:when test="not($ICGEM_Keyword = 'comments')">
          <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value"        select="$ICGEM_Keyword"/>
            <xsl:with-param name="Empty_Field"        select="$Empty_ICGEM_Keyword_Field"/>
            <xsl:with-param name="Justify"            select="$Left_Justified"/>
            <xsl:with-param name="Truncate_If_Longer" select="'no'"/>
          </xsl:call-template>
        </xsl:when>
	<xsl:otherwise>
          <xsl:value-of select="."/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:apply-templates select="*" mode="process_ICGEM_record">
        <xsl:with-param name="Empty_Degree_Field" select="$Empty_Degree_Field"/>
        <xsl:with-param name="Empty_Order_Field"  select="$Empty_Order_Field"/>
      </xsl:apply-templates>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

  </xsl:template>

  <xsl:template match="*" mode="process_ICGEM_record">

    <xsl:param name="Empty_Degree_Field"/>
    <xsl:param name="Empty_Order_Field"/>

    <xsl:choose>
      <xsl:when test="*">

        <xsl:choose>
          <xsl:when test="C/* and S/*">
            <xsl:value-of select="$DFS"/>
	    <xsl:call-template name="Prepend_Signless_Value">
	      <xsl:with-param name="Value"  select="C/Calibrated"/>
	      <xsl:with-param name="Prefix" select="$Value_Prefix"/>
            </xsl:call-template>
            <xsl:value-of select="$DFS"/>
	    <xsl:call-template name="Prepend_Signless_Value">
	      <xsl:with-param name="Value"  select="S/Calibrated"/>
	      <xsl:with-param name="Prefix" select="$Value_Prefix"/>
            </xsl:call-template>
            <xsl:value-of select="$DFS"/>
	    <xsl:call-template name="Prepend_Signless_Value">
	      <xsl:with-param name="Value"  select="C/Formal"/>
	      <xsl:with-param name="Prefix" select="$Value_Prefix"/>
            </xsl:call-template>
            <xsl:value-of select="$DFS"/>
	    <xsl:call-template name="Prepend_Signless_Value">
	      <xsl:with-param name="Value"  select="S/Formal"/>
	      <xsl:with-param name="Prefix" select="$Value_Prefix"/>
            </xsl:call-template>
          </xsl:when>
          <xsl:otherwise>
            <xsl:apply-templates select="*" mode="process_ICGEM_record">
              <xsl:with-param name="Empty_Degree_Field" select="$Empty_Degree_Field"/>
              <xsl:with-param name="Empty_Order_Field"  select="$Empty_Order_Field"/>
            </xsl:apply-templates>
          </xsl:otherwise>
        </xsl:choose>

      </xsl:when>
      <xsl:otherwise>

        <xsl:value-of select="$DFS"/>

        <xsl:choose>
          <xsl:when test="local-name(.) = 'Degree'">
            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value"        select="."/>
              <xsl:with-param name="Empty_Field"        select="$Empty_Degree_Field"/>
              <xsl:with-param name="Justify"            select="$Right_Justified"/>
              <xsl:with-param name="Truncate_If_Longer" select="'no'"/>
            </xsl:call-template>
	  </xsl:when>
          <xsl:when test="local-name(.) = 'Order'">
            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value"        select="."/>
              <xsl:with-param name="Empty_Field"        select="$Empty_Order_Field"/>
              <xsl:with-param name="Justify"            select="$Right_Justified"/>
              <xsl:with-param name="Truncate_If_Longer" select="'no'"/>
            </xsl:call-template>
	  </xsl:when>
          <xsl:when test="local-name(.) = 'Comments' or local-name(.) = 'Date'">
            <xsl:value-of select="."/>
	  </xsl:when>
          <xsl:otherwise>
	    <xsl:call-template name="Prepend_Signless_Value">
              <xsl:with-param name="Value"  select="."/>
              <xsl:with-param name="Prefix" select="$Value_Prefix"/>
            </xsl:call-template>
          </xsl:otherwise>
        </xsl:choose>

      </xsl:otherwise>
    </xsl:choose>

  </xsl:template>

  <xsl:template match="*/Data_Block/EGM_GOC_2/EGM_GRP_2/List_of_PDF_Files/PDF_File">

    <xsl:variable name="Encoding" select="translate(@encoding, $UC, $LC)"/>

    <xsl:if test="$Encoding = 'base64'">
      <xsl:value-of select="."/>
    </xsl:if>

  </xsl:template>

</xsl:stylesheet>
