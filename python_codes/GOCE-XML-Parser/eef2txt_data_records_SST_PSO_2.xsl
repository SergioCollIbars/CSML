<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: SST_PSO_2_data_records
Version: 1.1
Date: 14 Jul 2008
-->

<xsl:stylesheet id="SST_PSO_2_data_records" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:variable name="Last_File_Indicator" select="'last_file'"/>
  <xsl:variable name="End_Of_File_Keyword" select="'EOF'"/>

  <xsl:template match="*/Data_Block/SST_PSO_2/*[local-name()='SST_PKI_2' or local-name()='SST_PRD_2']/List_of_SP3c_Records">
  
    <xsl:variable name="Unit_Multiplier" select="'b^n'"/>

    <xsl:variable name="Epoch_Header_Record_Symbol"                          select="'*'"/> 
    <xsl:variable name="Pos_And_Clock_Record_Symbol"                         select="'P'"/> 
    <xsl:variable name="Correlation_Record_Symbol"                           select="'E'"/> 
    <xsl:variable name="Pos_And_Clock_Correlation_Record_Symbol"             select="concat($Correlation_Record_Symbol, $Pos_And_Clock_Record_Symbol)"/> 
    <xsl:variable name="Vel_And_Clock_Change_Rate_Record_Symbol"             select="'V'"/> 
    <xsl:variable name="Vel_And_Clock_Change_Rate_Correlation_Record_Symbol" select="concat($Correlation_Record_Symbol, $Vel_And_Clock_Change_Rate_Record_Symbol)"/> 

    <xsl:variable name="Short_Symbol_Field_Length" select="1"/>
    <xsl:variable name="Vehicle_ID_Field_Length"   select="3"/>

    <xsl:variable name="Coordinate_Field_Length" select="14"/>
    <xsl:variable name="Clock_Field_Length"      select="$Coordinate_Field_Length"/>
    
    <xsl:variable name="Coordinates_Std_Deviation_Field_Length"   select="2"/>
    <xsl:variable name="Clock_Std_Deviation_Field_Length"         select="3"/>
    <xsl:variable name="E_Coordinates_Std_Deviation_Field_Length" select="4"/>
    <xsl:variable name="E_Clock_Std_Deviation_Field_Length"       select="7"/>

    <xsl:variable name="Flag_Field_Length"        select="$Short_Symbol_Field_Length"/>
    <xsl:variable name="Correlation_Field_Length" select="8"/>

    <xsl:variable name="Velocity_Field_Length"          select="$Coordinate_Field_Length"/>
    <xsl:variable name="Clock_Change_Rate_Field_Length" select="$Coordinate_Field_Length"/>

    <xsl:variable name="Velocity_Std_Deviation_Field_Length"            select="$Coordinates_Std_Deviation_Field_Length"/>
    <xsl:variable name="Clock_Change_Rate_Std_Deviation_Field_Length"   select="$Clock_Std_Deviation_Field_Length"/>
    <xsl:variable name="E_Velocity_Std_Deviation_Field_Length"          select="$E_Coordinates_Std_Deviation_Field_Length"/>
    <xsl:variable name="E_Clock_Change_Rate_Std_Deviation_Field_Length" select="$E_Clock_Std_Deviation_Field_Length"/>

    <xsl:variable name="Empty_Short_Symbol_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Short_Symbol_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Vehicle_ID_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Vehicle_ID_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="Empty_Coordinate_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Coordinate_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Clock_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Clock_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="Empty_Coordinates_Std_Deviation_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Coordinates_Std_Deviation_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Clock_Std_Deviation_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Clock_Std_Deviation_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_E_Coordinates_Std_Deviation_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$E_Coordinates_Std_Deviation_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_E_Clock_Std_Deviation_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$E_Clock_Std_Deviation_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="Empty_Flag_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Flag_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Correlation_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Correlation_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="Empty_Velocity_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Velocity_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Clock_Change_Rate_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Clock_Change_Rate_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="Empty_Velocity_Std_Deviation_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Velocity_Std_Deviation_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Clock_Change_Rate_Std_Deviation_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Clock_Change_Rate_Std_Deviation_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_E_Velocity_Std_Deviation_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$E_Velocity_Std_Deviation_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_E_Clock_Change_Rate_Std_Deviation_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$E_Clock_Change_Rate_Std_Deviation_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:for-each select="SP3c_Record">

      <xsl:apply-templates select="Time_Information/GPS_Time/Start/Gregorian" mode="process_epoch_header_record">
        <xsl:with-param name="Symbol" select="$Epoch_Header_Record_Symbol"/>
      </xsl:apply-templates>

      <xsl:apply-templates select="List_of_Satellite_IDs" mode="process_satellites">
        <xsl:with-param name="Unit_Multiplier"                                     select="$Unit_Multiplier"/>
        <xsl:with-param name="Pos_And_Clock_Record_Symbol"                         select="$Pos_And_Clock_Record_Symbol"/>
        <xsl:with-param name="Pos_And_Clock_Correlation_Record_Symbol"             select="$Pos_And_Clock_Correlation_Record_Symbol"/>
        <xsl:with-param name="Vel_And_Clock_Change_Rate_Record_Symbol"             select="$Vel_And_Clock_Change_Rate_Record_Symbol"/>
        <xsl:with-param name="Vel_And_Clock_Change_Rate_Correlation_Record_Symbol" select="$Vel_And_Clock_Change_Rate_Correlation_Record_Symbol"/>
        <xsl:with-param name="Empty_Short_Symbol_Field"                            select="$Empty_Short_Symbol_Field"/>
        <xsl:with-param name="Empty_Vehicle_ID_Field"                              select="$Empty_Vehicle_ID_Field"/>
        <xsl:with-param name="Empty_Coordinate_Field"                              select="$Empty_Coordinate_Field"/>
        <xsl:with-param name="Empty_Clock_Field"                                   select="$Empty_Clock_Field"/>
        <xsl:with-param name="Empty_Coordinates_Std_Deviation_Field"               select="$Empty_Coordinates_Std_Deviation_Field"/>
        <xsl:with-param name="Empty_Clock_Std_Deviation_Field"                     select="$Empty_Clock_Std_Deviation_Field"/>
        <xsl:with-param name="Empty_E_Coordinates_Std_Deviation_Field"             select="$Empty_E_Coordinates_Std_Deviation_Field"/>
        <xsl:with-param name="Empty_E_Clock_Std_Deviation_Field"                   select="$Empty_E_Clock_Std_Deviation_Field"/>
        <xsl:with-param name="Empty_Flag_Field"                                    select="$Empty_Flag_Field"/>
        <xsl:with-param name="Empty_Correlation_Field"                             select="$Empty_Correlation_Field"/>
        <xsl:with-param name="Empty_Velocity_Field"                                select="$Empty_Velocity_Field"/>
        <xsl:with-param name="Empty_Clock_Change_Rate_Field"                       select="$Empty_Clock_Change_Rate_Field"/>
        <xsl:with-param name="Empty_Velocity_Std_Deviation_Field"                  select="$Empty_Velocity_Std_Deviation_Field"/>
        <xsl:with-param name="Empty_Clock_Change_Rate_Std_Deviation_Field"         select="$Empty_Clock_Change_Rate_Std_Deviation_Field"/>
        <xsl:with-param name="Empty_E_Velocity_Std_Deviation_Field"                select="$Empty_E_Velocity_Std_Deviation_Field"/>
        <xsl:with-param name="Empty_E_Clock_Change_Rate_Std_Deviation_Field"       select="$Empty_E_Clock_Change_Rate_Std_Deviation_Field"/>
      </xsl:apply-templates>

    </xsl:for-each>

    <xsl:if test="not(starts-with($Mode, $Last_File_Indicator)) and contains($Mode, $Last_File_Indicator)">
      <xsl:value-of select="$End_Of_File_Keyword"/>
      <xsl:text>&#xa;</xsl:text>
    </xsl:if>

  </xsl:template>

  <xsl:template match="*" mode="process_epoch_header_record">

    <xsl:param name="Symbol" select="'  '"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Symbol"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Symbol_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>
  
    <xsl:text> </xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Year"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Year_Start_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>
  
    <xsl:text> </xsl:text>
  
    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Month"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Month_Start_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>
  
    <xsl:text> </xsl:text>
  
    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Day_of_Month"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Day_Of_Month_Start_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>
  
    <xsl:text> </xsl:text>
  
    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Hour"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Hour_Start_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>
  
    <xsl:text> </xsl:text>
  
    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Minute"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Minute_Start_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>
  
    <xsl:text> </xsl:text>
  
    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="Second"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Second_Start_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

  </xsl:template>

  <xsl:template match="*" mode="process_satellites">

    <xsl:param name="Unit_Multiplier"/>
    <xsl:param name="Pos_And_Clock_Record_Symbol"/>
    <xsl:param name="Pos_And_Clock_Correlation_Record_Symbol"/>
    <xsl:param name="Vel_And_Clock_Change_Rate_Record_Symbol"/>
    <xsl:param name="Vel_And_Clock_Change_Rate_Correlation_Record_Symbol"/>
    <xsl:param name="Empty_Short_Symbol_Field"/>
    <xsl:param name="Empty_Vehicle_ID_Field"/>
    <xsl:param name="Empty_Coordinate_Field"/>
    <xsl:param name="Empty_Clock_Field"/>
    <xsl:param name="Empty_Coordinates_Std_Deviation_Field"/>
    <xsl:param name="Empty_Clock_Std_Deviation_Field"/>
    <xsl:param name="Empty_E_Coordinates_Std_Deviation_Field"/>
    <xsl:param name="Empty_E_Clock_Std_Deviation_Field"/>
    <xsl:param name="Empty_Flag_Field"/>
    <xsl:param name="Empty_Correlation_Field"/>
    <xsl:param name="Empty_Velocity_Field"/>
    <xsl:param name="Empty_Clock_Change_Rate_Field"/>
    <xsl:param name="Empty_Velocity_Std_Deviation_Field"/>
    <xsl:param name="Empty_Clock_Change_Rate_Std_Deviation_Field"/>
    <xsl:param name="Empty_E_Velocity_Std_Deviation_Field"/>
    <xsl:param name="Empty_E_Clock_Change_Rate_Std_Deviation_Field"/>

    <xsl:for-each select="*">

      <xsl:variable name="Vehicle_ID" select="local-name(.)"/>

      <xsl:if test="Position">
        
        <xsl:variable name="E_Coordinates_Std_Deviation_Unit" select="'mm'"/>
        <xsl:variable name="E_Clock_Std_Deviation_Unit"       select="'psec'"/>
        <xsl:variable name="Coordinates_Std_Deviation_Unit"   select="concat($Unit_Multiplier, ' ', $E_Coordinates_Std_Deviation_Unit)"/>
        <xsl:variable name="Clock_Std_Deviation_Unit"         select="concat($Unit_Multiplier, ' ', $E_Clock_Std_Deviation_Unit)"/>

	<xsl:if test="Standard_Deviations/Position/@unit = $Coordinates_Std_Deviation_Unit and Standard_Deviations/Clock/@unit = $Clock_Std_Deviation_Unit">

          <xsl:variable name="P_Line">

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="$Pos_And_Clock_Record_Symbol"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Short_Symbol_Field"/>
              <xsl:with-param name="Justify"     select="$Left_Justified"/>
            </xsl:call-template>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="$Vehicle_ID"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Vehicle_ID_Field"/>
              <xsl:with-param name="Justify"     select="$Left_Justified"/>
            </xsl:call-template>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Position/X"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Coordinate_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Position/Y"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Coordinate_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Position/Z"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Coordinate_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Clock"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Clock_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Position[@unit = $Coordinates_Std_Deviation_Unit]/X"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Coordinates_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Position[@unit = $Coordinates_Std_Deviation_Unit]/Y"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Coordinates_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Position[@unit = $Coordinates_Std_Deviation_Unit]/Z"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Coordinates_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Clock[@unit = $Clock_Std_Deviation_Unit]"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Clock_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Flags/Clock/Event"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
              <xsl:with-param name="Justify"     select="$Left_Justified"/>
            </xsl:call-template>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Flags/Clock/Correction_Prediction"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
              <xsl:with-param name="Justify"     select="$Left_Justified"/>
            </xsl:call-template>

            <xsl:text>  </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Flags/Orbit/Maneuver"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
              <xsl:with-param name="Justify"     select="$Left_Justified"/>
            </xsl:call-template>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Flags/Orbit/Prediction"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Flag_Field"/>
              <xsl:with-param name="Justify"     select="$Left_Justified"/>
            </xsl:call-template>

          </xsl:variable>

	  <xsl:if test="$P_Line">

            <xsl:call-template name="Trim_Trailing_Spaces">
              <xsl:with-param name="String" select="$P_Line"/>
            </xsl:call-template>

            <xsl:text>&#xa;</xsl:text>

	  </xsl:if>

        </xsl:if>

	<xsl:if test="Standard_Deviations/Position/@unit = $E_Coordinates_Std_Deviation_Unit and Standard_Deviations/Clock/@unit = $E_Clock_Std_Deviation_Unit">

          <xsl:variable name="EP_Line">

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="$Pos_And_Clock_Correlation_Record_Symbol"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Symbol_Field"/>
              <xsl:with-param name="Justify"     select="$Left_Justified"/>
            </xsl:call-template>

            <xsl:text>  </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Position[@unit = $E_Coordinates_Std_Deviation_Unit]/X"/>
              <xsl:with-param name="Empty_Field" select="$Empty_E_Coordinates_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Position[@unit = $E_Coordinates_Std_Deviation_Unit]/Y"/>
              <xsl:with-param name="Empty_Field" select="$Empty_E_Coordinates_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Position[@unit = $E_Coordinates_Std_Deviation_Unit]/Z"/>
              <xsl:with-param name="Empty_Field" select="$Empty_E_Coordinates_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Clock[@unit = $E_Clock_Std_Deviation_Unit]"/>
              <xsl:with-param name="Empty_Field" select="$Empty_E_Clock_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Correlations/Position/XY"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Correlation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Correlations/Position/XZ"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Correlation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Correlations/Position/XC"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Correlation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Correlations/Position/YZ"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Correlation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Correlations/Position/YC"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Correlation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Correlations/Position/ZC"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Correlation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

          </xsl:variable>

	  <xsl:if test="$EP_Line">

            <xsl:call-template name="Trim_Trailing_Spaces">
              <xsl:with-param name="String" select="$EP_Line"/>
            </xsl:call-template>

            <xsl:text>&#xa;</xsl:text>

	  </xsl:if>

        </xsl:if>

      </xsl:if>

      <xsl:if test="Velocity">

        <xsl:variable name="E_Velocity_Std_Deviation_Unit"          select="'10^-4 mm/s'"/>
        <xsl:variable name="E_Clock_Change_Rate_Std_Deviation_Unit" select="'10^-4 psec/sec'"/>
        <xsl:variable name="Velocity_Std_Deviation_Unit"            select="concat($Unit_Multiplier, ' ', $E_Velocity_Std_Deviation_Unit)"/>
        <xsl:variable name="Clock_Change_Rate_Std_Deviation_Unit"   select="concat($Unit_Multiplier, ' ', $E_Clock_Change_Rate_Std_Deviation_Unit)"/>

	<xsl:if test="Standard_Deviations/Velocity/@unit = $Velocity_Std_Deviation_Unit and Standard_Deviations/Clock_Change_Rate/@unit = $Clock_Change_Rate_Std_Deviation_Unit">

          <xsl:variable name="V_Line">

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="$Vel_And_Clock_Change_Rate_Record_Symbol"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Short_Symbol_Field"/>
              <xsl:with-param name="Justify"     select="$Left_Justified"/>
            </xsl:call-template>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="$Vehicle_ID"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Vehicle_ID_Field"/>
              <xsl:with-param name="Justify"     select="$Left_Justified"/>
            </xsl:call-template>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Velocity/X"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Velocity_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Velocity/Y"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Velocity_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Velocity/Z"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Velocity_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Clock_Change_Rate"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Clock_Change_Rate_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Velocity[@unit = $Velocity_Std_Deviation_Unit]/X"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Velocity_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Velocity[@unit = $Velocity_Std_Deviation_Unit]/Y"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Velocity_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Velocity[@unit = $Velocity_Std_Deviation_Unit]/Z"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Velocity_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Velocity[@unit = $Clock_Change_Rate_Std_Deviation_Unit]"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Clock_Change_Rate_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

          </xsl:variable>

	  <xsl:if test="$V_Line">

            <xsl:call-template name="Trim_Trailing_Spaces">
              <xsl:with-param name="String" select="$V_Line"/>
            </xsl:call-template>

            <xsl:text>&#xa;</xsl:text>

	  </xsl:if>

        </xsl:if>

	<xsl:if test="Standard_Deviations/Velocity/@unit = $E_Velocity_Std_Deviation_Unit and Standard_Deviations/Clock/@unit = $E_Clock_Change_Rate_Std_Deviation_Unit">

          <xsl:variable name="EV_Line">

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="$Vel_And_Clock_Change_Rate_Correlation_Record_Symbol"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Symbol_Field"/>
              <xsl:with-param name="Justify"     select="$Left_Justified"/>
            </xsl:call-template>

            <xsl:text>  </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Velocity[@unit = $E_Velocity_Std_Deviation_Unit]/X"/>
              <xsl:with-param name="Empty_Field" select="$Empty_E_Coordinates_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Velocity[@unit = $E_Velocity_Std_Deviation_Unit]/Y"/>
              <xsl:with-param name="Empty_Field" select="$Empty_E_Coordinates_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Velocity[@unit = $E_Velocity_Std_Deviation_Unit]/Z"/>
              <xsl:with-param name="Empty_Field" select="$Empty_E_Coordinates_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Clock_Change_Rate[@unit = $E_Clock_Change_Rate_Std_Deviation_Unit]"/>
              <xsl:with-param name="Empty_Field" select="$Empty_E_Clock_Std_Deviation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Correlations/Velocity/XY"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Correlation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Correlations/Velocity/XZ"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Correlation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Correlations/Velocity/XC"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Correlation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Correlations/Velocity/YZ"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Correlation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Correlations/Velocity/YC"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Correlation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

            <xsl:text> </xsl:text>

            <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="Standard_Deviations/Correlations/Velocity/ZC"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Correlation_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>

          </xsl:variable>

	  <xsl:if test="$EV_Line">

            <xsl:call-template name="Trim_Trailing_Spaces">
              <xsl:with-param name="String" select="$EV_Line"/>
            </xsl:call-template>

            <xsl:text>&#xa;</xsl:text>

	  </xsl:if>

        </xsl:if>

      </xsl:if>

    </xsl:for-each>

  </xsl:template>

  <xsl:template match="*/Data_Block/SST_PSO_2/SST_PRM_2/List_of_Rotation_Records">
  
    <xsl:variable name="Time_Field_Length" select="9"/>
    <xsl:variable name="Q1_Field_Length"   select="19"/>
    <xsl:variable name="Q2_Field_Length"   select="$Q1_Field_Length"/>
    <xsl:variable name="Q3_Field_Length"   select="$Q1_Field_Length"/>
    <xsl:variable name="Q4_Field_Length"   select="$Q1_Field_Length"/>

    <xsl:variable name="Empty_Time_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Time_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Q1_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Q1_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Q2_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Q2_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Q3_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Q3_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Q4_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Q4_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:for-each select="Rotation_Record">

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Time_Information/Time_Since_Reference_Epoch"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Time_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Quaternions/Q1"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Q1_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Quaternions/Q2"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Q2_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Quaternions/Q3"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Q3_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:value-of select="$DFS"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="Quaternions/Q4"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Q4_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>
  
    <xsl:if test="not(starts-with($Mode, $Last_File_Indicator)) and contains($Mode, $Last_File_Indicator)">
      <xsl:value-of select="concat($Comment_Indicator, ' ', $End_Of_File_Keyword)"/>
      <xsl:text>&#xa;</xsl:text>
    </xsl:if>

  </xsl:template>

  <xsl:template match="*/Data_Block/SST_PSO_2/SST_PCV_2/List_of_Covariance_Records">

    <xsl:variable name="Row_Field_Length"    select="7"/>
    <xsl:variable name="Column_Field_Length" select="$Row_Field_Length"/>
    <xsl:variable name="Value_Field_Length"  select="14"/>

    <xsl:variable name="Empty_Row_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Row_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Column_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Column_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Value_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Value_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:for-each select="Covariance_Record">

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="@row"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Row_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="@column"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Column_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="."/>
        <xsl:with-param name="Empty_Field" select="$Empty_Value_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>
  
    <xsl:if test="not(starts-with($Mode, $Last_File_Indicator)) and contains($Mode, $Last_File_Indicator)">
      <xsl:value-of select="concat($Comment_Indicator, ' ', $End_Of_File_Keyword)"/>
      <xsl:text>&#xa;</xsl:text>
    </xsl:if>

  </xsl:template>

  <xsl:template match="*/Data_Block/SST_PSO_2/SST_PRP_2/List_of_PDF_Files/PDF_File">

    <xsl:variable name="Encoding" select="translate(@encoding, $UC, $LC)"/>

    <xsl:if test="$Encoding = 'base64'">
      <xsl:value-of select="."/>
    </xsl:if>

  </xsl:template>

</xsl:stylesheet>
