<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: EGM_GOC_2_data_header
Version: 1.1
Date: 14 Jul 2008
-->

<xsl:stylesheet id="EGM_GOC_2_data_header" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template match="*/Earth_Explorer_Header/Variable_Header/SPH/EGM_GOC_2">

    <xsl:if test="not(starts-with($Mode, 'data_block'))">
      <xsl:choose>
        <xsl:when test="$Product">
	  <xsl:if test="not($Product = 'EGM_GRP_2')">
            <xsl:apply-templates select="*[local-name()=$Product]" mode="header"/>
	  </xsl:if>
        </xsl:when>
        <xsl:otherwise>
          <xsl:apply-templates mode="header"/>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:if>

  </xsl:template>

  <xsl:template match="*/Earth_Explorer_Header/Variable_Header/SPH/EGM_GOC_2/*[not(local-name()='EGM_GCF_2' or local-name()='EGM_GRP_2')]" mode="header">
 
    <xsl:variable name="Data_Set_Name_Keyword"                    select="'data_set_name'"/>
    <xsl:variable name="Northern_Latitude_Keyword"                select="'northern_latitude'"/>
    <xsl:variable name="Southern_Latitude_Keyword"                select="'southern_latitude'"/>
    <xsl:variable name="Western_Longitude_Keyword"                select="'western_longitude'"/>
    <xsl:variable name="Eastern_Longitude_Keyword"                select="'eastern_longitude'"/>
    <xsl:variable name="Latitude_Cell_Size_Keyword"               select="'latitude_cell_size'"/>
    <xsl:variable name="Longitude_Cell_Size_Keyword"              select="'longitude_cell_size'"/>
    <xsl:variable name="Number_Of_Cells_In_Latitude_Dir_Keyword"  select="'number_of_cells_latitude_dir'"/>
    <xsl:variable name="Number_Of_Cells_In_Longitude_Dir_Keyword" select="'number_of_cells_longitude_dir'"/>
    <xsl:variable name="Mean_Or_Point_Values_Keyword"             select="'mean_(0)_or_point_(1)_values'"/>
    <xsl:variable name="Geocentric_Or_Geodetic_Latitude_Keyword"  select="'geocentric(0)_geodetic(1)_lat'"/>
    <xsl:variable name="Reference_Ellipsoid_Keyword"              select="'reference_ellipsoid'"/>
    <xsl:variable name="Format_Keyword"                           select="'format_of_data'"/>
    <xsl:variable name="Gap_Value_Keyword"                        select="'gap_value'"/>
    <xsl:variable name="Description_Keyword"                      select="'description_of_data'"/>
    <xsl:variable name="Unit_Keyword"                             select="'unit'"/>
    <xsl:variable name="End_Of_Header_Keyword"                    select="'end_of_header'"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Data_Set_Name_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Data_Set_Name" select="Data_Information/Dataset_Name"/>

    <xsl:value-of select="$Data_Set_Name"/>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Northern_Latitude_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Northern_Latitude" select="Coordinate_Information/Latitude/Northern_Border"/>

    <xsl:call-template name="Justify_On_Decimal_Point">
      <xsl:with-param name="Value"             select="$Northern_Latitude"/>
      <xsl:with-param name="Justify_Character" select="$Value_Prefix"/>
      <xsl:with-param name="Max_Length">
        <xsl:choose>
          <xsl:when test="starts-with($Northern_Latitude,'-') or starts-with($Northern_Latitude,'+')">4</xsl:when>
          <xsl:otherwise>3</xsl:otherwise>
        </xsl:choose>
      </xsl:with-param>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Southern_Latitude_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Southern_Latitude" select="Coordinate_Information/Latitude/Southern_Border"/>

    <xsl:call-template name="Justify_On_Decimal_Point">
      <xsl:with-param name="Value"             select="$Southern_Latitude"/>
      <xsl:with-param name="Justify_Character" select="$Value_Prefix"/>
      <xsl:with-param name="Max_Length">
        <xsl:choose>
          <xsl:when test="starts-with($Southern_Latitude,'-') or starts-with($Southern_Latitude,'+')">4</xsl:when>
          <xsl:otherwise>3</xsl:otherwise>
        </xsl:choose>
      </xsl:with-param>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Western_Longitude_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Western_Longitude" select="Coordinate_Information/Longitude/Western_Border"/>

    <xsl:call-template name="Justify_On_Decimal_Point">
      <xsl:with-param name="Value"             select="$Western_Longitude"/>
      <xsl:with-param name="Justify_Character" select="$Value_Prefix"/>
      <xsl:with-param name="Max_Length">
        <xsl:choose>
          <xsl:when test="starts-with($Western_Longitude,'-') or starts-with($Western_Longitude,'+')">4</xsl:when>
          <xsl:otherwise>3</xsl:otherwise>
        </xsl:choose>
      </xsl:with-param>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Eastern_Longitude_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Eastern_Longitude" select="Coordinate_Information/Longitude/Eastern_Border"/>

    <xsl:call-template name="Justify_On_Decimal_Point">
      <xsl:with-param name="Value"             select="$Eastern_Longitude"/>
      <xsl:with-param name="Justify_Character" select="$Value_Prefix"/>
      <xsl:with-param name="Max_Length">
        <xsl:choose>
          <xsl:when test="starts-with($Eastern_Longitude,'-') or starts-with($Eastern_Longitude,'+')">4</xsl:when>
          <xsl:otherwise>3</xsl:otherwise>
        </xsl:choose>
      </xsl:with-param>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Latitude_Cell_Size_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Latitude_Cell_Size" select="Coordinate_Information/Latitude/Cell_Information/Size"/>

    <xsl:call-template name="Prepend_Signless_Value">
      <xsl:with-param name="Value"  select="$Latitude_Cell_Size"/>
      <xsl:with-param name="Prefix" select="$Value_Prefix"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Longitude_Cell_Size_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Longitude_Cell_Size" select="Coordinate_Information/Longitude/Cell_Information/Size"/>

    <xsl:call-template name="Prepend_Signless_Value">
      <xsl:with-param name="Value"  select="$Longitude_Cell_Size"/>
      <xsl:with-param name="Prefix" select="$Value_Prefix"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Number_Of_Cells_In_Latitude_Dir_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Number_Of_Cells_In_Latitude_Dir" select="Coordinate_Information/Latitude/Cell_Information/Number_of_Cells"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value"        select="$Number_Of_Cells_In_Latitude_Dir"/>
      <xsl:with-param name="Empty_Field"        select="$Empty_Number_Of_Values_Field"/>
      <xsl:with-param name="Justify"            select="$Right_Justified"/>
      <xsl:with-param name="Truncate_If_Longer" select="'no'"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Number_Of_Cells_In_Longitude_Dir_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Number_Of_Cells_In_Longitude_Dir" select="Coordinate_Information/Longitude/Cell_Information/Number_of_Cells"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value"        select="$Number_Of_Cells_In_Longitude_Dir"/>
      <xsl:with-param name="Empty_Field"        select="$Empty_Number_Of_Values_Field"/>
      <xsl:with-param name="Justify"            select="$Right_Justified"/>
      <xsl:with-param name="Truncate_If_Longer" select="'no'"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Mean_Or_Point_Values_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Mean_Or_Point_Values"    select="Flags/Mean_or_Point_Values"/>
    <xsl:variable name="LC_Mean_Or_Point_Values" select="translate($Mean_Or_Point_Values, $UC, $LC)"/>
    <xsl:variable name="Mean_Values_Indicator"   select="'mean'"/>
    <xsl:variable name="Point_Values_Indicator"  select="'point'"/>

    <xsl:choose>
      <xsl:when test="$LC_Mean_Or_Point_Values = translate($Mean_Values_Indicator, $UC, $LC)">
	<xsl:value-of select="'0'"/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:if test="$LC_Mean_Or_Point_Values = translate($Point_Values_Indicator, $UC, $LC)">
	  <xsl:value-of select="'1'"/>
	</xsl:if>
      </xsl:otherwise>
    </xsl:choose>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Geocentric_Or_Geodetic_Latitude_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Geocentric_Or_Geodetic_Latitude"    select="Flags/Geocentric_or_Geodetic_Latitudes"/>
    <xsl:variable name="LC_Geocentric_Or_Geodetic_Latitude" select="translate($Geocentric_Or_Geodetic_Latitude, $UC, $LC)"/>
    <xsl:variable name="Geocentric_Latitude_Indicator"      select="'geocentric'"/>
    <xsl:variable name="Geodetic_Latitude_Indicator"        select="'geodetic'"/>

    <xsl:choose>
      <xsl:when test="$LC_Geocentric_Or_Geodetic_Latitude = translate($Geocentric_Latitude_Indicator, $UC, $LC)">
	<xsl:value-of select="'0'"/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:if test="$LC_Geocentric_Or_Geodetic_Latitude = translate($Geodetic_Latitude_Indicator, $UC, $LC)">
	  <xsl:value-of select="'1'"/>
	</xsl:if>
      </xsl:otherwise>
    </xsl:choose>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Reference_Ellipsoid_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Reference_Ellipsoid" select="Reference_Ellipsoid"/>

    <xsl:value-of select="$Reference_Ellipsoid"/>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Format_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:value-of select="$EGM_GAN_2_Data_Format"/>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Gap_Value_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Gap_Value" select="Gap_Value"/>

    <xsl:value-of select="$Gap_Value"/>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Description_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Description" select="Data_Information/Description"/>

    <xsl:value-of select="$Description"/>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Unit_Keyword"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Unit" select="Data_Information/Unit"/>

    <xsl:value-of select="$Unit"/>

    <xsl:text>&#xa;</xsl:text>
    
    <xsl:value-of select="$End_Of_Header_Keyword"/>

    <xsl:text>&#xa;</xsl:text>

  </xsl:template>

  <xsl:template match="*/Earth_Explorer_Header/Variable_Header/SPH/EGM_GOC_2/EGM_GCF_2" mode="header">

    <xsl:variable name="Product_Type_Keyword"            select="'product_type'"/>
    <xsl:variable name="Model_Name_Keyword"              select="'modelname'"/>
    <xsl:variable name="Earth_Gravity_Constant_Keyword"  select="'earth_gravity_constant'"/>
    <xsl:variable name="Radius_Keyword"                  select="'radius'"/>
    <xsl:variable name="Max_Degree_Keyword"              select="'max_degree'"/>
    <xsl:variable name="Errors_Keyword"                  select="'errors'"/>
    <xsl:variable name="Norm_Keyword"                    select="'norm'"/>
    <xsl:variable name="Tide_System_Keyword"             select="'tide_system'"/>
    <xsl:variable name="Water_Density_Keyword"           select="'water_density'"/>
    <xsl:variable name="End_Of_EGM_GCF_2_Header_Keyword" select="'end_of_head'"/>

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

    <xsl:value-of select="$End_Of_EGM_GCF_2_Header_Keyword"/>

    <xsl:text>&#xa;</xsl:text>

  </xsl:template>

</xsl:stylesheet>
