<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: SST_AUX_2_data_records
Version: 1.0
Date: 20 Jun 2008
-->

<xsl:stylesheet id="SST_AUX_2_data_records" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:template match="*/Data_Block/SST_AUX_2/List_of_ICGEM_Records">

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

</xsl:stylesheet>
