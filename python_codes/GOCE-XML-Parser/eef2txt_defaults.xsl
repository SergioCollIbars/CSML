<?xml version="1.0" encoding="US-ASCII"?>

<!--
Defaults Stylesheet
Version: 1.1
Date: 23 Dec 2010
-->

<xsl:stylesheet id="defaults" version="1.1" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:variable name="HFS"             select="'|'"/>
  <xsl:variable name="DFS"             select="' '"/>
  <xsl:variable name="Left_Justified"  select="'Left'"/>
  <xsl:variable name="Right_Justified" select="'Right'"/>
  <xsl:variable name="LC"              select="'abcdefghijklmnopqrstuvwxyz'"/>
  <xsl:variable name="UC"              select="'ABCDEFGHIJKLMNOPQRSTUVWXYZ'"/>

  <xsl:template match="*|/">
    <xsl:apply-templates/>
  </xsl:template>

  <xsl:template match="text()|@*">
    <xsl:apply-templates/>
  </xsl:template>

  <xsl:template name="Fixed_Length_Field">
    
    <xsl:param name="Field_Value"/>
    <xsl:param name="Empty_Field"/>
    <xsl:param name="Justify"/>
    <xsl:param name="Truncate_If_Longer" select="'yes'"/>
    <xsl:param name="Append_If_Shorter"  select="'yes'"/>

    <xsl:variable name="Value_Length" select="string-length($Field_Value)"/>
    <xsl:variable name="Field_Length" select="string-length($Empty_Field)"/>

    <xsl:choose>
      <xsl:when test="$Value_Length &lt; $Field_Length">

        <xsl:choose>
          <xsl:when test="$Justify = $Right_Justified">
            <xsl:value-of select="substring($Empty_Field, 1+$Value_Length)"/>
            <xsl:value-of select="$Field_Value"/>
          </xsl:when>
          <xsl:otherwise>
	    <xsl:value-of select="$Field_Value"/>
	    <xsl:if test="$Append_If_Shorter = 'yes'">
              <xsl:value-of select="substring($Empty_Field, 1+$Value_Length)"/>
	    </xsl:if>
          </xsl:otherwise>
        </xsl:choose>

      </xsl:when>
      <xsl:when test="$Value_Length &gt; $Field_Length">

        <xsl:choose>
          <xsl:when test="$Truncate_If_Longer = 'yes'">
            <xsl:value-of select="substring($Field_Value, 1, $Field_Length)"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="$Field_Value"/>
          </xsl:otherwise>
        </xsl:choose>

      </xsl:when>
      <xsl:otherwise>

        <xsl:value-of select="$Field_Value"/>

      </xsl:otherwise>
    </xsl:choose>

  </xsl:template>

  <xsl:template name="Justify_On_Decimal_Point">
  
    <xsl:param name="Value"/>
    <xsl:param name="Justify_Character" select="' '"/>
    <xsl:param name="Max_Length"        select="0"/>
    <xsl:param name="Length"            select="string-length(substring-before($Value, '.'))"/>
  
    <xsl:choose>
      <xsl:when test="$Length &lt;= $Max_Length + string-length($Justify_Character)">

        <xsl:if test="$Length &lt; $Max_Length">
          <xsl:value-of select="$Justify_Character"/>
	</xsl:if>

        <xsl:call-template name="Justify_On_Decimal_Point">
          <xsl:with-param name="Value"             select="$Value"/>
          <xsl:with-param name="Justify_Character" select="$Justify_Character"/>
          <xsl:with-param name="Max_Length"        select="$Max_Length"/>
          <xsl:with-param name="Length"            select="$Length + string-length($Justify_Character)"/>
        </xsl:call-template>

      </xsl:when>
      <xsl:otherwise>
        
	<xsl:call-template name="Prepend_Signless_Value">
          <xsl:with-param name="Value"  select="$Value"/>
          <xsl:with-param name="Prefix" select="$Value_Prefix"/>
        </xsl:call-template>

      </xsl:otherwise>
    </xsl:choose>
  
  </xsl:template>

  <xsl:template name="Combine_Records_To_Number_Of_Fields">

    <xsl:param name="Total_Number_Of_Fields_Per_Line"/>
    <xsl:param name="Number_Of_Fields" select="1"/>
    <xsl:param name="Empty_Field"/>
    <xsl:param name="List"/>
    <xsl:param name="Line"/>

    <xsl:choose>
      <xsl:when test="$List">

        <xsl:variable name="Current_List" select="$List[$Line] | $List[position() > 1][not($Line)]"/>
	<xsl:variable name="Current_Line">
	  <xsl:choose>
            <xsl:when test="$Line"><xsl:value-of select="$Line"/></xsl:when>
	    <xsl:otherwise>
	      <xsl:call-template name="Fixed_Length_Field">
	        <xsl:with-param name="Field_Value" select="$List[1][not($Line)]"/>
	        <xsl:with-param name="Empty_Field" select="$Empty_Field"/>
	        <xsl:with-param name="Justify"     select="$Right_Justified"/>
              </xsl:call-template>
	    </xsl:otherwise>
	  </xsl:choose>
	</xsl:variable>

        <xsl:choose>
          <xsl:when test="$Number_Of_Fields &lt; $Total_Number_Of_Fields_Per_Line">

	    <xsl:variable name="Current_Value">
	      <xsl:if test="$Current_List[1]">
	        <xsl:call-template name="Fixed_Length_Field">
	          <xsl:with-param name="Field_Value" select="$Current_List[1]"/>
	          <xsl:with-param name="Empty_Field" select="$Empty_Field"/>
	          <xsl:with-param name="Justify"     select="$Right_Justified"/>
                </xsl:call-template>
	      </xsl:if>
            </xsl:variable>
 
	    <xsl:call-template name="Combine_Records_To_Number_Of_Fields">
              <xsl:with-param name="Total_Number_Of_Fields_Per_Line" select="$Total_Number_Of_Fields_Per_Line"/>
              <xsl:with-param name="Number_Of_Fields"                select="$Number_Of_Fields + 1"/>
              <xsl:with-param name="Empty_Field"                     select="$Empty_Field"/>
              <xsl:with-param name="List"                            select="$Current_List[position() > 1]"/>
              <xsl:with-param name="Line">
	        <xsl:choose>
                  <xsl:when test="$Current_Value and not($Current_Value = '')">
		    <xsl:value-of select="concat($Current_Line, $DFS, $Current_Value)"/>
		  </xsl:when>
	          <xsl:otherwise>
		    <xsl:value-of select="$Current_Line"/>
		  </xsl:otherwise>
	        </xsl:choose>
	      </xsl:with-param>
            </xsl:call-template>

          </xsl:when>
          <xsl:otherwise>
            
            <xsl:value-of select="concat($Current_Line, '&#xa;')"/>
            
	    <xsl:call-template name="Combine_Records_To_Number_Of_Fields">
              <xsl:with-param name="Total_Number_Of_Fields_Per_Line" select="$Total_Number_Of_Fields_Per_Line"/>
              <xsl:with-param name="Empty_Field"                     select="$Empty_Field"/>
              <xsl:with-param name="List"                            select="$Current_List"/>
              <xsl:with-param name="Line"                            select="''"/>
            </xsl:call-template>

          </xsl:otherwise>
        </xsl:choose>

      </xsl:when>
      <xsl:otherwise>
        
        <xsl:value-of select="$Line"/>

      </xsl:otherwise>
    </xsl:choose>

  </xsl:template>

  <xsl:template name="Combine_Records_To_Max_Length">

    <xsl:param name="Max_Length"/>
    <xsl:param name="List"/>
    <xsl:param name="Line"/>

    <xsl:choose>
      <xsl:when test="$List">

        <xsl:variable name="Current_List" select="$List[$Line] | $List[position() > 1][not($Line)]"/>
	<xsl:variable name="Current_Line" select="concat($Line, $List[1][not($Line)])"/>
	
        <xsl:choose>
          <xsl:when test="string-length($Current_Line) + string-length($DFS) + string-length($Current_List[1]) &lt;= $Max_Length">
             
            <xsl:variable name="Current_Value">
	      <xsl:call-template name="Prepend_Signless_Value">
                <xsl:with-param name="Value"  select="$Current_List[1]"/>
                <xsl:with-param name="Prefix" select="$Value_Prefix"/>
              </xsl:call-template>
            </xsl:variable>

	    <xsl:call-template name="Combine_Records_To_Max_Length">
              <xsl:with-param name="Max_Length" select="$Max_Length"/>
              <xsl:with-param name="List"       select="$Current_List[position() > 1]"/>
              <xsl:with-param name="Line"       select="concat($Current_Line, $DFS, $Current_Value)"/>
            </xsl:call-template>

          </xsl:when>
          <xsl:otherwise>
            
	    <xsl:call-template name="Prepend_Signless_Value">
              <xsl:with-param name="Value"  select="concat($Current_Line, '&#xa;')"/>
              <xsl:with-param name="Prefix" select="$Value_Prefix"/>
            </xsl:call-template>
            
	    <xsl:call-template name="Combine_Records_To_Max_Length">
              <xsl:with-param name="Max_Length" select="$Max_Length"/>
              <xsl:with-param name="List"       select="$Current_List"/>
              <xsl:with-param name="Line"       select="''"/>
            </xsl:call-template>

          </xsl:otherwise>
        </xsl:choose>

      </xsl:when>
      <xsl:otherwise>
        
	<xsl:call-template name="Prepend_Signless_Value">
          <xsl:with-param name="Value"  select="$Line"/>
          <xsl:with-param name="Prefix" select="$Value_Prefix"/>
        </xsl:call-template>

      </xsl:otherwise>
    </xsl:choose>

  </xsl:template>

  <xsl:template name="Prepend_Signless_Value">

    <xsl:param name="Value"/>
    <xsl:param name="Prefix" select="' '"/>

    <xsl:choose>
      <xsl:when test="starts-with($Value,'-') or starts-with($Value,'+')">
        <xsl:value-of select="$Value"/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="concat($Prefix, $Value)"/>
      </xsl:otherwise>
    </xsl:choose>

  </xsl:template>

  <xsl:template name="Trim_Trailing_Spaces">

    <xsl:param name="String"/>
    <xsl:param name="String_Length" select="string-length($String)"/>

    <xsl:choose>
      <xsl:when test="translate(substring($String,$String_Length,1),' ','')">
        <xsl:value-of select="substring($String,1,$String_Length)"/>
      </xsl:when>
      <xsl:when test="$String_Length &lt; 2"/>
      <xsl:otherwise>
        <xsl:call-template name="Trim_Trailing_Spaces">
          <xsl:with-param name="String"        select="$String"/>
          <xsl:with-param name="String_Length" select="$String_Length - 1"/>
        </xsl:call-template>
      </xsl:otherwise>
    </xsl:choose>

  </xsl:template>

  <xsl:template name="Construct_Empty_Field">
  
    <xsl:param name="Field_Length"    select="0"/>
    <xsl:param name="Empty_Character" select="' '"/>
    <xsl:param name="Empty_Field"     select="$Empty_Character"/>
  
    <xsl:variable name="Current_Length" select="string-length($Empty_Field)"/>

    <xsl:choose>
      <xsl:when test="$Current_Length &lt; $Field_Length">

        <xsl:call-template name="Construct_Empty_Field">
          <xsl:with-param name="Field_Length"    select="$Field_Length"/>
          <xsl:with-param name="Empty_Character" select="$Empty_Character"/>
          <xsl:with-param name="Empty_Field"     select="concat($Empty_Field,$Empty_Character)"/>
        </xsl:call-template>

      </xsl:when>
      <xsl:otherwise>
        
	<xsl:value-of select="$Empty_Field"/>

      </xsl:otherwise>
    </xsl:choose>
  
  </xsl:template>



  <xsl:template name="Recombine_3_Gradients">      
      
     <xsl:param name="XX"/>
     <xsl:param name="YY"/>
     <xsl:param name="ZZ"/>
     <xsl:param name="Nb_Comp"          select="3"/>

     <xsl:variable name="Total_Length"     select="string-length($XX)"/>
     <xsl:variable name="Length"           select="($Total_Length - $Nb_Comp + 1) div $Nb_Comp"/>

     <xsl:variable name="Comp_Length"      select="15"/>

     <xsl:variable name="Empty_Field">
       <xsl:call-template name="Construct_Empty_Field">
         <xsl:with-param name="Field_Length" select="$Comp_Length"/>
       </xsl:call-template>
     </xsl:variable>

     <xsl:choose>
       <xsl:when test="$Nb_Comp &gt;= 1">

         <xsl:variable name="XXdeb" select="substring($XX,1,$Length)"/>
         <xsl:variable name="YYdeb" select="substring($YY,1,$Length)"/>
         <xsl:variable name="ZZdeb" select="substring($ZZ,1,$Length)"/>
          
         <xsl:call-template name="Fixed_Length_Field">
	    <xsl:with-param name="Field_Value" select="$XXdeb"/>
	    <xsl:with-param name="Empty_Field" select="$Empty_Field"/>
	    <xsl:with-param name="Justify"     select="$Left_Justified"/>
         </xsl:call-template>

         <xsl:value-of select="$DFS"/>


         <xsl:call-template name="Fixed_Length_Field">
	    <xsl:with-param name="Field_Value" select="$YYdeb"/>
	    <xsl:with-param name="Empty_Field" select="$Empty_Field"/>
	    <xsl:with-param name="Justify"     select="$Left_Justified"/>
         </xsl:call-template>

         <xsl:value-of select="$DFS"/>

         <xsl:call-template name="Fixed_Length_Field">
	    <xsl:with-param name="Field_Value" select="$ZZdeb"/>
	    <xsl:with-param name="Empty_Field" select="$Empty_Field"/>
	    <xsl:with-param name="Justify"     select="$Left_Justified"/>
         </xsl:call-template>

         <xsl:value-of select="$DFS"/>

        <xsl:call-template name="Recombine_3_Gradients">
           <xsl:with-param name="XX"        select="substring-after($XX,' ')"/>
           <xsl:with-param name="YY"        select="substring-after($YY,' ')"/>
           <xsl:with-param name="ZZ"        select="substring-after($ZZ,' ')"/>
           <xsl:with-param name="Nb_Comp"   select="$Nb_Comp - 1"/>
        </xsl:call-template>
     
      </xsl:when>
     </xsl:choose>     

  </xsl:template>


</xsl:stylesheet>
