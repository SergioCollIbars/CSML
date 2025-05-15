<?xml version="1.0" encoding="US-ASCII"?>

<!--
Stylesheet Name: SST_PSO_2_data_header
Version: 1.1
Date: 14 Jul 2008
-->

<xsl:stylesheet id="SST_PSO_2_data_header" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:variable name="Date_Separator" select="'-'"/>
  <xsl:variable name="Time_Separator" select="':'"/>

  <xsl:variable name="Creation_Program_Keyword" select="'Program that created the file'"/>
  <xsl:variable name="Creation_Date_Keyword"    select="'Date of creation'"/>
  <xsl:variable name="First_Epoch_Keyword"      select="'First epoch'"/>
  <xsl:variable name="End_Of_Header_Keyword"    select="'End of header'"/>

  <xsl:variable name="Creation_Program_Field_Length" select="30"/>
  <xsl:variable name="Creation_Date_Field_Length"    select="19"/>
  <xsl:variable name="First_Epoch_Field_Length"      select="23"/>

  <xsl:variable name="Empty_Creation_Program_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Creation_Program_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_Creation_Date_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$Creation_Date_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="Empty_First_Epoch_Field">
    <xsl:call-template name="Construct_Empty_Field">
      <xsl:with-param name="Field_Length" select="$First_Epoch_Field_Length"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:template match="*/Earth_Explorer_Header/Variable_Header/SPH/SST_PSO_2">

    <xsl:if test="not(starts-with($Mode, 'data_block'))">
      <xsl:choose>
        <xsl:when test="$Product">
	  <xsl:if test="not($Product = 'SST_PRP_2')">
            <xsl:apply-templates select="*[local-name()=$Product]" mode="header"/>
	  </xsl:if>
        </xsl:when>
        <xsl:otherwise>
          <xsl:apply-templates mode="header"/>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:if>

  </xsl:template>

  <xsl:template match="*/Earth_Explorer_Header/Variable_Header/SPH/SST_PSO_2/*[local-name()='SST_PKI_2' or local-name()='SST_PRD_2']" mode="header">

    <xsl:variable name="Two_Character_Field_Length"   select="2"/>
    <xsl:variable name="Three_Character_Field_Length" select="3"/>
    <xsl:variable name="Four_Character_Field_Length"  select="4"/>
    <xsl:variable name="Five_Character_Field_Length"  select="5"/>

    <xsl:variable name="Ten_Column_Float_Field_Length"      select="10"/>
    <xsl:variable name="Twelve_Column_Float_Field_Length"   select="12"/>
    <xsl:variable name="Fourteen_Column_Float_Field_Length" select="14"/>
    <xsl:variable name="Eighteen_Column_Float_Field_Length" select="18"/>

    <xsl:variable name="Four_Column_Integer_Field_Length" select="$Four_Character_Field_Length"/>
    <xsl:variable name="Six_Column_Integer_Field_Length"  select="6"/>
    <xsl:variable name="Nine_Column_Integer_Field_Length" select="9"/>

    <xsl:variable name="Pos_Or_Vel_Flag_Field_Length"    select="1"/>
    <xsl:variable name="Number_Of_Epochs_Field_Length"   select="7"/>
    <xsl:variable name="Data_Used_Field_Length"          select="$Five_Character_Field_Length"/>
    <xsl:variable name="Coordinate_System_Field_Length"  select="$Five_Character_Field_Length"/>
    <xsl:variable name="Orbit_Type_Field_Length"         select="$Three_Character_Field_Length"/>
    <xsl:variable name="Agency_Field_Length"             select="$Four_Character_Field_Length"/>

    <xsl:variable name="GPS_Week_Field_Length"          select="$Four_Character_Field_Length"/>
    <xsl:variable name="Seconds_Of_Week_Field_Length"   select="15"/>
    <xsl:variable name="Epoch_Interval_Field_Length"    select="$Fourteen_Column_Float_Field_Length"/>
    <xsl:variable name="Mod_Jul_Day_Start_Field_Length" select="$Five_Character_Field_Length"/>
    <xsl:variable name="Fractional_Day_Field_Length"    select="$Seconds_Of_Week_Field_Length"/>

    <xsl:variable name="Number_Of_Sats_Field_Length" select="$Two_Character_Field_Length"/>
    <xsl:variable name="Sat_ID_Field_Length"         select="$Three_Character_Field_Length"/>
    <xsl:variable name="Accuracy_Field_Length"       select="$Three_Character_Field_Length"/>

    <xsl:variable name="File_Type_Field_Length"   select="$Two_Character_Field_Length"/>
    <xsl:variable name="Time_System_Field_Length" select="$Three_Character_Field_Length"/>

    <xsl:variable name="Base_for_Pos_or_Vel_Field_Length"  select="$Ten_Column_Float_Field_Length"/>
    <xsl:variable name="Base_for_Clk_or_Rate_Field_Length" select="$Twelve_Column_Float_Field_Length"/>

    <xsl:variable name="Comment_Field_Length" select="57"/>

    <xsl:variable name="Empty_Pos_Or_Vel_Flag_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Pos_Or_Vel_Flag_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Number_Of_Epochs_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Number_Of_Epochs_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Data_Used_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Data_Used_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Coordinate_System_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Coordinate_System_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Orbit_Type_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Orbit_Type_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Agency_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Agency_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="Empty_GPS_Week_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$GPS_Week_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Seconds_Of_Week_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Seconds_Of_Week_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Epoch_Interval_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Epoch_Interval_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Mod_Jul_Day_Start_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Mod_Jul_Day_Start_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Fractional_Day_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Fractional_Day_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="Empty_Number_Of_Sats_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Number_Of_Sats_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Sat_ID_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Sat_ID_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Accuracy_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Accuracy_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="Empty_File_Type_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$File_Type_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Time_System_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Time_System_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Two_Character_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length"    select="$Two_Character_Field_Length"/>
	<xsl:with-param name="Empty_Character" select="'c'"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Three_Character_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length"    select="$Three_Character_Field_Length"/>
	<xsl:with-param name="Empty_Character" select="'c'"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Four_Character_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length"    select="$Four_Character_Field_Length"/>
	<xsl:with-param name="Empty_Character" select="'c'"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Five_Character_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length"    select="$Five_Character_Field_Length"/>
	<xsl:with-param name="Empty_Character" select="'c'"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="Empty_Base_for_Pos_or_Vel_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Base_for_Pos_or_Vel_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Base_for_Clk_or_Rate_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Base_for_Clk_or_Rate_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Ten_Column_Float_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Ten_Column_Float_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Twelve_Column_Float_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Twelve_Column_Float_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Fourteen_Column_Float_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Fourteen_Column_Float_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Eighteen_Column_Float_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Eighteen_Column_Float_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="Empty_Four_Column_Integer_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Four_Column_Integer_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Six_Column_Integer_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Six_Column_Integer_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Nine_Column_Integer_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Nine_Column_Integer_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="Empty_Comment_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Comment_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="Format_Version" select="Original_Source/Format/Version"/>
    <xsl:variable name="Version_Symbol" select="concat($Comment_Indicator, $Format_Version)"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Version_Symbol"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Symbol_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Pos_Or_Vel" select="Pos_or_Vel"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Pos_Or_Vel"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Pos_Or_Vel_Flag_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Year_Start" select="Time_Information/GPS_Time/Start/Gregorian/Year"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Year_Start"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Year_Start_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="Month_Start" select="Time_Information/GPS_Time/Start/Gregorian/Month"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Month_Start"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Month_Start_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="Day_Of_Month_Start" select="Time_Information/GPS_Time/Start/Gregorian/Day_of_Month"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Day_Of_Month_Start"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Day_Of_Month_Start_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="Hour_Start" select="Time_Information/GPS_Time/Start/Gregorian/Hour"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Hour_Start"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Hour_Start_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="Minute_Start" select="Time_Information/GPS_Time/Start/Gregorian/Minute"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Minute_Start"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Minute_Start_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="Second_Start" select="Time_Information/GPS_Time/Start/Gregorian/Second"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Second_Start"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Second_Start_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="Number_Of_Epochs" select="Epoch_Information/Num_Epochs"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Number_Of_Epochs"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Number_Of_Epochs_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="Data_Used" select="Data_Used"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Data_Used"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Data_Used_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="Coordinate_System" select="Coordinate_Sys"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Coordinate_System"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Coordinate_System_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="Orbit_Type" select="Orbit_Type"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Orbit_Type"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Orbit_Type_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="Agency" select="Agency"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Agency"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Agency_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:variable name="Alt_Time_Info_Symbol" select="concat($Comment_Indicator, $Comment_Indicator)"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Alt_Time_Info_Symbol"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Symbol_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="GPS_Week" select="Time_Information/GPS_Time/Start/GPS/Week"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$GPS_Week"/>
      <xsl:with-param name="Empty_Field" select="$Empty_GPS_Week_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="Seconds_Of_Week" select="Time_Information/GPS_Time/Start/GPS/Seconds_of_Week"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Seconds_Of_Week"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Seconds_Of_Week_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="Epoch_Interval" select="Epoch_Information/Interval"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Epoch_Interval"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Epoch_Interval_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="Mod_Jul_Day_Start" select="Time_Information/GPS_Time/Start/Mod_Jul_Day/Day"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Mod_Jul_Day_Start"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Mod_Jul_Day_Start_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text> </xsl:text>

    <xsl:variable name="Fractional_Day" select="Time_Information/GPS_Time/Start/Mod_Jul_Day/Fractional_Day"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Fractional_Day"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Fractional_Day_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:variable name="Sat_ID_Symbol"              select="'+'"/>
    <xsl:variable name="Not_Applicable_Sat_ID"      select="0"/>
    <xsl:variable name="Number_Of_Sat_ID_Lines"     select="5"/>
    <xsl:variable name="Number_Of_Sat_IDs_Per_Line" select="17"/>

    <xsl:variable name="Satellite_Descriptors" select="List_of_Satellite_Descriptors"/>
    <xsl:variable name="Number_Of_Sats"        select="$Satellite_Descriptors/@count"/>

    <xsl:for-each select="(//*)[position() &lt;= $Number_Of_Sat_ID_Lines]">

      <xsl:variable name="Sat_ID_Line_Number" select="position()"/>
      <xsl:variable name="Sat_ID_Count_Start" select="$Number_Of_Sat_IDs_Per_Line * ($Sat_ID_Line_Number - 1)"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="$Sat_ID_Symbol"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Symbol_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:text>  </xsl:text>

      <xsl:choose>
        <xsl:when test="$Sat_ID_Line_Number = 1">
          <xsl:choose>
            <xsl:when test="$Number_Of_Sats">
	      <xsl:call-template name="Fixed_Length_Field">
                <xsl:with-param name="Field_Value" select="$Number_Of_Sats"/>
                <xsl:with-param name="Empty_Field" select="$Empty_Number_Of_Sats_Field"/>
                <xsl:with-param name="Justify"     select="$Right_Justified"/>
              </xsl:call-template>
            </xsl:when>
            <xsl:otherwise>
	      <xsl:call-template name="Fixed_Length_Field">
                <xsl:with-param name="Field_Value" select="'0'"/>
                <xsl:with-param name="Empty_Field" select="$Empty_Number_Of_Sats_Field"/>
                <xsl:with-param name="Justify"     select="$Right_Justified"/>
              </xsl:call-template>
            </xsl:otherwise>
          </xsl:choose>
        </xsl:when>
        <xsl:otherwise>
	  <xsl:value-of select="$Empty_Number_Of_Sats_Field"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:text>   </xsl:text>

      <xsl:for-each select="(//*)[position() &lt;= $Number_Of_Sat_IDs_Per_Line]">

        <xsl:variable name="Sat_ID_Count" select="$Sat_ID_Count_Start + position()"/>
        
	<xsl:choose>
	  <xsl:when test="$Satellite_Descriptors/Satellite_Descriptor[$Sat_ID_Count]/Satellite_ID">
	    <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="$Satellite_Descriptors/Satellite_Descriptor[$Sat_ID_Count]/Satellite_ID"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Sat_ID_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>
	  </xsl:when>
	  <xsl:otherwise>
	    <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="$Not_Applicable_Sat_ID"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Sat_ID_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>
	  </xsl:otherwise>
        </xsl:choose>

      </xsl:for-each>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

    <xsl:variable name="Accuracy_Symbol"         select="concat($Sat_ID_Symbol, $Sat_ID_Symbol)"/>
    <xsl:variable name="Not_Applicable_Accuracy" select="$Not_Applicable_Sat_ID"/>

    <xsl:for-each select="(//*)[position() &lt;= $Number_Of_Sat_ID_Lines]">

      <xsl:variable name="Accuracy_Line_Number" select="position()"/>
      <xsl:variable name="Accuracy_Count_Start" select="$Number_Of_Sat_IDs_Per_Line * ($Accuracy_Line_Number - 1)"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="$Accuracy_Symbol"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Symbol_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:text>       </xsl:text>

      <xsl:for-each select="(//*)[position() &lt;= $Number_Of_Sat_IDs_Per_Line]">

        <xsl:variable name="Accuracy_Count" select="$Accuracy_Count_Start + position()"/>

        <xsl:choose>
	  <xsl:when test="$Satellite_Descriptors/Satellite_Descriptor[$Accuracy_Count]/Accuracy">
	    <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="$Satellite_Descriptors/Satellite_Descriptor[$Accuracy_Count]/Accuracy"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Accuracy_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>
	  </xsl:when>
	  <xsl:otherwise>
	    <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="$Not_Applicable_Accuracy"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Accuracy_Field"/>
              <xsl:with-param name="Justify"     select="$Right_Justified"/>
            </xsl:call-template>
	  </xsl:otherwise>
        </xsl:choose>

      </xsl:for-each>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

    <xsl:variable name="Character_Symbol"          select="'%c'"/>
    <xsl:variable name="Number_Of_Character_Lines" select="2"/>

    <xsl:variable name="File_Type"   select="Original_Source/Format/Type"/>
    <xsl:variable name="Time_System" select="Time_Information/System"/>

    <xsl:for-each select="(//*)[position() &lt;= $Number_Of_Character_Lines]">
      
      <xsl:variable name="Character_Line_Number" select="position()"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="$Character_Symbol"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Symbol_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:text> </xsl:text>

      <xsl:choose>
        <xsl:when test="$Character_Line_Number = 1">
	  <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="$File_Type"/>
            <xsl:with-param name="Empty_Field" select="$Empty_File_Type_Field"/>
            <xsl:with-param name="Justify"     select="$Left_Justified"/>
          </xsl:call-template>
        </xsl:when>
        <xsl:otherwise>
	  <xsl:value-of select="$Empty_Two_Character_Field"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:text> </xsl:text>

      <xsl:value-of select="$Empty_Two_Character_Field"/>

      <xsl:text> </xsl:text>

      <xsl:choose>
        <xsl:when test="$Character_Line_Number = 1">
	  <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="$Time_System"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Time_System_Field"/>
            <xsl:with-param name="Justify"     select="$Left_Justified"/>
          </xsl:call-template>
        </xsl:when>
        <xsl:otherwise>
	  <xsl:value-of select="$Empty_Three_Character_Field"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:text> </xsl:text>

      <xsl:value-of select="$Empty_Three_Character_Field"/>

      <xsl:for-each select="(//*)[position() &lt;= 4]">
        <xsl:text> </xsl:text>
        <xsl:value-of select="$Empty_Four_Character_Field"/>
      </xsl:for-each>

      <xsl:for-each select="(//*)[position() &lt;= 4]">
        <xsl:text> </xsl:text>
        <xsl:value-of select="$Empty_Five_Character_Field"/>
      </xsl:for-each>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

    <xsl:variable name="Float_Symbol"          select="'%f'"/>
    <xsl:variable name="Number_Of_Float_Lines" select="2"/>

    <xsl:variable name="Base_for_Pos_or_Vel"  select="Base_for_Pos_or_Vel"/>
    <xsl:variable name="Base_for_Clk_or_Rate" select="Base_for_Clk_or_Rate"/>

    <xsl:for-each select="(//*)[position() &lt;= $Number_Of_Float_Lines]">

      <xsl:variable name="Float_Line_Number" select="position()"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="$Float_Symbol"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Symbol_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:text> </xsl:text>

      <xsl:choose>
        <xsl:when test="$Float_Line_Number = 1">
	  
	  <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="$Base_for_Pos_or_Vel"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Base_for_Pos_or_Vel_Field"/>
            <xsl:with-param name="Justify"     select="$Right_Justified"/>
          </xsl:call-template>

          <xsl:text> </xsl:text>

	  <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="$Base_for_Clk_or_Rate"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Base_for_Clk_or_Rate_Field"/>
            <xsl:with-param name="Justify"     select="$Right_Justified"/>
          </xsl:call-template>

        </xsl:when>
        <xsl:otherwise>
	  
	  <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="'0.0000000'"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Ten_Column_Float_Field"/>
            <xsl:with-param name="Justify"     select="$Right_Justified"/>
          </xsl:call-template>

          <xsl:text> </xsl:text>

	  <xsl:call-template name="Fixed_Length_Field">
            <xsl:with-param name="Field_Value" select="'0.000000000'"/>
            <xsl:with-param name="Empty_Field" select="$Empty_Twelve_Column_Float_Field"/>
            <xsl:with-param name="Justify"     select="$Right_Justified"/>
          </xsl:call-template>

        </xsl:otherwise>
      </xsl:choose>

      <xsl:text> </xsl:text>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="'0.00000000000'"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Fourteen_Column_Float_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:text> </xsl:text>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="'0.000000000000000'"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Eighteen_Column_Float_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

    <xsl:variable name="Integer_Symbol"                       select="'%i'"/>
    <xsl:variable name="Number_Of_Integer_Lines"              select="2"/>
    <xsl:variable name="Number_Of_Four_Column_Integer_Fields" select="4"/>
    <xsl:variable name="Number_Of_Six_Column_Integer_Fields"  select="$Number_Of_Four_Column_Integer_Fields"/>

    <xsl:for-each select="(//*)[position() &lt;= $Number_Of_Integer_Lines]">

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="$Integer_Symbol"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Symbol_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:for-each select="(//*)[position() &lt;= $Number_Of_Four_Column_Integer_Fields]">

        <xsl:text> </xsl:text>

        <xsl:call-template name="Fixed_Length_Field">
          <xsl:with-param name="Field_Value" select="'0'"/>
          <xsl:with-param name="Empty_Field" select="$Empty_Four_Column_Integer_Field"/>
          <xsl:with-param name="Justify"     select="$Right_Justified"/>
        </xsl:call-template>

      </xsl:for-each>

      <xsl:for-each select="(//*)[position() &lt;= $Number_Of_Six_Column_Integer_Fields]">

        <xsl:text> </xsl:text>

        <xsl:call-template name="Fixed_Length_Field">
          <xsl:with-param name="Field_Value" select="'0'"/>
          <xsl:with-param name="Empty_Field" select="$Empty_Six_Column_Integer_Field"/>
          <xsl:with-param name="Justify"     select="$Right_Justified"/>
        </xsl:call-template>

      </xsl:for-each>

      <xsl:text> </xsl:text>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="'0'"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Nine_Column_Integer_Field"/>
        <xsl:with-param name="Justify"     select="$Right_Justified"/>
      </xsl:call-template>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

    <xsl:variable name="Comment_Symbol"          select="'/*'"/>
    <xsl:variable name="Number_Of_Comment_Lines" select="4"/>

    <xsl:variable name="Comments_Line_One_To_Four"    select="substring-after(Comments, '&#xA;')"/>
    <xsl:variable name="Comments_Line_Two_To_Four"    select="substring-after($Comments_Line_One_To_Four, '&#xA;')"/>
    <xsl:variable name="Comments_Line_Three_And_Four" select="substring-after($Comments_Line_Two_To_Four, '&#xA;')"/>
    <xsl:variable name="Comments_Line_Four"           select="substring-after($Comments_Line_Three_And_Four, '&#xA;')"/>

    <xsl:variable name="Normalized_Comments_Line_One">
      <xsl:choose>
        <xsl:when test="normalize-space($Comments_Line_Two_To_Four) = ''">
          <xsl:value-of select="normalize-space($Comments_Line_One_To_Four)"/>
	</xsl:when>
	<xsl:otherwise>
          <xsl:value-of select="substring-before(normalize-space($Comments_Line_One_To_Four), concat(' ', normalize-space($Comments_Line_Two_To_Four)))"/>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>   
    <xsl:variable name="Normalized_Comments_Line_Two">
      <xsl:choose>
        <xsl:when test="normalize-space($Comments_Line_Three_And_Four) = ''">
          <xsl:value-of select="normalize-space($Comments_Line_Two_To_Four)"/>
	</xsl:when>
	<xsl:otherwise>
          <xsl:value-of select="substring-before(normalize-space($Comments_Line_Two_To_Four), concat(' ', normalize-space($Comments_Line_Three_And_Four)))"/>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>   
    <xsl:variable name="Normalized_Comments_Line_Three">
      <xsl:choose>
        <xsl:when test="normalize-space($Comments_Line_Four) = ''">
          <xsl:value-of select="normalize-space($Comments_Line_Three_And_Four)"/>
	</xsl:when>
	<xsl:otherwise>
          <xsl:value-of select="substring-before(normalize-space($Comments_Line_Three_And_Four), concat(' ', normalize-space($Comments_Line_Four)))"/>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>   
    <xsl:variable name="Normalized_Comments_Line_Four"  select="normalize-space($Comments_Line_Four)"/>

    <xsl:for-each select="(//*)[position() &lt;= $Number_Of_Comment_Lines]">

      <xsl:variable name="Comment_Line_Number" select="position()"/>

      <xsl:call-template name="Fixed_Length_Field">
        <xsl:with-param name="Field_Value" select="$Comment_Symbol"/>
        <xsl:with-param name="Empty_Field" select="$Empty_Symbol_Field"/>
        <xsl:with-param name="Justify"     select="$Left_Justified"/>
      </xsl:call-template>

      <xsl:text> </xsl:text>

      <xsl:choose>
          <xsl:when test="$Comment_Line_Number = 1">
	    <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="$Normalized_Comments_Line_One"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Comment_Field"/>
              <xsl:with-param name="Justify"     select="$Left_Justified"/>
            </xsl:call-template>
	  </xsl:when>
          <xsl:when test="$Comment_Line_Number = 2">
	    <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="$Normalized_Comments_Line_Two"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Comment_Field"/>
              <xsl:with-param name="Justify"     select="$Left_Justified"/>
            </xsl:call-template>
	  </xsl:when>
          <xsl:when test="$Comment_Line_Number = 3">
	    <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="$Normalized_Comments_Line_Three"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Comment_Field"/>
              <xsl:with-param name="Justify"     select="$Left_Justified"/>
            </xsl:call-template>
	  </xsl:when>
          <xsl:when test="$Comment_Line_Number = 4">
	    <xsl:call-template name="Fixed_Length_Field">
              <xsl:with-param name="Field_Value" select="$Normalized_Comments_Line_Four"/>
              <xsl:with-param name="Empty_Field" select="$Empty_Comment_Field"/>
              <xsl:with-param name="Justify"     select="$Left_Justified"/>
            </xsl:call-template>
	  </xsl:when>
      </xsl:choose>

      <xsl:text>&#xa;</xsl:text>

    </xsl:for-each>

  </xsl:template>

  <xsl:template match="*/Earth_Explorer_Header/Variable_Header/SPH/SST_PSO_2/SST_PCV_2" mode="header">

    <xsl:variable name="Var_Cov_Matrix_Keyword"                select="'VARIANCE-COVARIANCE MATRIX'"/>
    <xsl:variable name="Corresponding_Kinematic_Orbit_Keyword" select="'Corresponding kinematic orbit'"/>
    <xsl:variable name="Time_Step_Size_Keyword"                select="'Time step size'"/>
    <xsl:variable name="RMS_Keyword"                           select="'RMS of unit weight'"/>
    <xsl:variable name="Parameters_Keyword"                    select="'Parameters'"/>

    <xsl:variable name="Var_Cov_Matrix_Field_Length"                select="62"/>
    <xsl:variable name="Corresponding_Kinematic_Orbit_Field_Length" select="$Var_Cov_Matrix_Field_Length"/>
    <xsl:variable name="Time_Step_Size_Field_Length"                select="6"/>
    <xsl:variable name="RMS_Field_Length"                           select="10"/>
    <xsl:variable name="Parameters_Field_Length"                    select="$Creation_Program_Field_Length"/>

    <xsl:variable name="Empty_Var_Cov_Matrix_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Var_Cov_Matrix_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Corresponding_Kinematic_Orbit_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Corresponding_Kinematic_Orbit_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Time_Step_Size_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Time_Step_Size_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_RMS_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$RMS_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Parameters_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Parameters_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Var_Cov_Matrix_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Var_Cov_Matrix" select="Var_Cov_Matrix/File_Name"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value"       select="$Var_Cov_Matrix"/>
      <xsl:with-param name="Empty_Field"       select="$Empty_Var_Cov_Matrix_Field"/>
      <xsl:with-param name="Justify"           select="$Left_Justified"/>
      <xsl:with-param name="Append_If_Shorter" select="'no'"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Corresponding_Kinematic_Orbit_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Corresponding_Kinematic_Orbit" select="Corresponding_Kinematic_Orbit/File_Name"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value"       select="$Corresponding_Kinematic_Orbit"/>
      <xsl:with-param name="Empty_Field"       select="$Empty_Corresponding_Kinematic_Orbit_Field"/>
      <xsl:with-param name="Justify"           select="$Left_Justified"/>
      <xsl:with-param name="Append_If_Shorter" select="'no'"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Creation_Program_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Creation_Program" select="Original_Source/Creator"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value"       select="$Creation_Program"/>
      <xsl:with-param name="Empty_Field"       select="$Empty_Creation_Program_Field"/>
      <xsl:with-param name="Justify"           select="$Left_Justified"/>
      <xsl:with-param name="Append_If_Shorter" select="'no'"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Creation_Date_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Creation_Date" select="Original_Source/Creation_Date"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Creation_Date"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Creation_Date_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $First_Epoch_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="CV_First_Epoch_Separator" select="$DFS"/>

    <xsl:variable name="CV_First_Epoch_Year" select="Time_Information/GPS_Time/Start/Gregorian/Year"/>
    <xsl:variable name="CV_First_Epoch_Month">
      <xsl:choose>
        <xsl:when test="string-length(Time_Information/GPS_Time/Start/Gregorian/Month) &lt; 2">
          <xsl:value-of select="concat('0', Time_Information/GPS_Time/Start/Gregorian/Month)"/>
	</xsl:when>
	<xsl:otherwise>
          <xsl:value-of select="Time_Information/GPS_Time/Start/Gregorian/Month"/>
	</xsl:otherwise>
      </xsl:choose>
    </xsl:variable>  
    <xsl:variable name="CV_First_Epoch_Day">
      <xsl:choose>
        <xsl:when test="string-length(Time_Information/GPS_Time/Start/Gregorian/Day_of_Month) &lt; 2">
          <xsl:value-of select="concat('0', Time_Information/GPS_Time/Start/Gregorian/Day_of_Month)"/>
	</xsl:when>
	<xsl:otherwise>
          <xsl:value-of select="Time_Information/GPS_Time/Start/Gregorian/Day_of_Month"/>
	</xsl:otherwise>
      </xsl:choose>    
    </xsl:variable>
    <xsl:variable name="CV_First_Epoch_Date" select="concat($CV_First_Epoch_Year, $Date_Separator, $CV_First_Epoch_Month, $Date_Separator, $CV_First_Epoch_Day)"/>
    <xsl:variable name="CV_First_Epoch_Hour">
      <xsl:choose>
        <xsl:when test="string-length(Time_Information/GPS_Time/Start/Gregorian/Hour) &lt; 2">
          <xsl:value-of select="concat('0', Time_Information/GPS_Time/Start/Gregorian/Hour)"/>
	</xsl:when>
	<xsl:otherwise>
          <xsl:value-of select="Time_Information/GPS_Time/Start/Gregorian/Hour"/>
	</xsl:otherwise>
      </xsl:choose>    
    </xsl:variable>
    <xsl:variable name="CV_First_Epoch_Minute">
      <xsl:choose>
        <xsl:when test="string-length(Time_Information/GPS_Time/Start/Gregorian/Minute) &lt; 2">
          <xsl:value-of select="concat('0', Time_Information/GPS_Time/Start/Gregorian/Minute)"/>
	</xsl:when>
	<xsl:otherwise>
          <xsl:value-of select="Time_Information/GPS_Time/Start/Gregorian/Minute"/>
	</xsl:otherwise>
      </xsl:choose>    
    </xsl:variable>
    <xsl:variable name="CV_First_Epoch_Second">
      <xsl:choose>
        <xsl:when test="string-length(Time_Information/GPS_Time/Start/Gregorian/Second) &lt; 2">
          <xsl:value-of select="concat('0', Time_Information/GPS_Time/Start/Gregorian/Second)"/>
	</xsl:when>
	<xsl:otherwise>
          <xsl:value-of select="Time_Information/GPS_Time/Start/Gregorian/Second"/>
	</xsl:otherwise>
      </xsl:choose>    
    </xsl:variable>
    <xsl:variable name="CV_First_Epoch_Time" select="concat($CV_First_Epoch_Hour, $Time_Separator, $CV_First_Epoch_Minute, $Time_Separator, $CV_First_Epoch_Second)"/>

    <xsl:variable name="CV_Time_System" select="Time_Information/System"/>
    <xsl:variable name="CV_First_Epoch" select="concat($CV_First_Epoch_Date, $CV_First_Epoch_Separator, $CV_First_Epoch_Time, $CV_First_Epoch_Separator, $CV_Time_System)"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$CV_First_Epoch"/>
      <xsl:with-param name="Empty_Field" select="$Empty_First_Epoch_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:variable name="Time_Step_Size_Unit" select="Time_Information/Time_Step_Size/@unit"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Time_Step_Size_Keyword, ' (', $Time_Step_Size_Unit, ')', ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Time_Step_Size" select="Time_Information/Time_Step_Size"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Time_Step_Size"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Time_Step_Size_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $RMS_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="RMS" select="RMS_of_Unit_Weight"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$RMS"/>
      <xsl:with-param name="Empty_Field" select="$Empty_RMS_Field"/>
      <xsl:with-param name="Justify"     select="$Right_Justified"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Parameters_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Parameters" select="Parameters"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value"       select="$Parameters"/>
      <xsl:with-param name="Empty_Field"       select="$Empty_Parameters_Field"/>
      <xsl:with-param name="Justify"           select="$Left_Justified"/>
      <xsl:with-param name="Append_If_Shorter" select="'no'"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>
    
    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $End_Of_Header_Keyword)"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

  </xsl:template>

  <xsl:template match="*/Earth_Explorer_Header/Variable_Header/SPH/SST_PSO_2/SST_PRM_2" mode="header">

    <xsl:variable name="Transformation_Keyword"           select="'TRANSFORMATION'"/>
    <xsl:variable name="Reference_Epoch_Keyword"          select="'Reference epoch'"/>
    <xsl:variable name="Transformation_Direction_Keyword" select="'Transformation direction'"/>
    <xsl:variable name="Pole_File_Keyword"                select="'Pole file'"/>
    <xsl:variable name="Nutation_Model_Keyword"           select="'Nutation model'"/>
    <xsl:variable name="Nutation_Offsets_Keyword"         select="'Nutation offsets'"/>
    <xsl:variable name="Subdaily_Model_Keyword"           select="'Subdaily model'"/>

    <xsl:variable name="Transformation_Field_Length"           select="62"/>
    <xsl:variable name="Reference_Epoch_Field_Length"          select="$First_Epoch_Field_Length"/>
    <xsl:variable name="Transformation_Direction_Field_Length" select="$Creation_Program_Field_Length"/>
    <xsl:variable name="Pole_File_Field_Length"                select="$Creation_Program_Field_Length"/>
    <xsl:variable name="Nutation_Model_Field_Length"           select="$Creation_Program_Field_Length"/>
    <xsl:variable name="Nutation_Offsets_Field_Length"         select="$Creation_Program_Field_Length"/>
    <xsl:variable name="Subdaily_Model_Field_Length"           select="$Creation_Program_Field_Length"/>

    <xsl:variable name="Empty_Transformation_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Transformation_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Reference_Epoch_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Reference_Epoch_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Transformation_Direction_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Transformation_Direction_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Pole_File_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Pole_File_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Nutation_Model_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Nutation_Model_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Nutation_Offsets_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Nutation_Offsets_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="Empty_Subdaily_Model_Field">
      <xsl:call-template name="Construct_Empty_Field">
        <xsl:with-param name="Field_Length" select="$Subdaily_Model_Field_Length"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Transformation_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Transformation" select="Transformation/File_Name"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value"       select="$Transformation"/>
      <xsl:with-param name="Empty_Field"       select="$Empty_Transformation_Field"/>
      <xsl:with-param name="Justify"           select="$Left_Justified"/>
      <xsl:with-param name="Append_If_Shorter" select="'no'"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Creation_Program_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Creation_Program" select="Original_Source/Creator"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value"       select="$Creation_Program"/>
      <xsl:with-param name="Empty_Field"       select="$Empty_Creation_Program_Field"/>
      <xsl:with-param name="Justify"           select="$Left_Justified"/>
      <xsl:with-param name="Append_If_Shorter" select="'no'"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Creation_Date_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Creation_Date" select="Original_Source/Creation_Date"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Creation_Date"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Creation_Date_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Reference_Epoch_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Reference_Epoch" select="Epoch_Information/Reference"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$Reference_Epoch"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Reference_Epoch_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $First_Epoch_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="RM_First_Epoch_Separator" select="$DFS"/>

    <xsl:variable name="RM_First_Epoch_Year" select="Time_Information/GPS_Time/Start/Gregorian/Year"/>
    <xsl:variable name="RM_First_Epoch_Month">
      <xsl:choose>
        <xsl:when test="string-length(Time_Information/GPS_Time/Start/Gregorian/Month) &lt; 2">
          <xsl:value-of select="concat('0', Time_Information/GPS_Time/Start/Gregorian/Month)"/>
	</xsl:when>
	<xsl:otherwise>
          <xsl:value-of select="Time_Information/GPS_Time/Start/Gregorian/Month"/>
	</xsl:otherwise>
      </xsl:choose>
    </xsl:variable>  
    <xsl:variable name="RM_First_Epoch_Day">
      <xsl:choose>
        <xsl:when test="string-length(Time_Information/GPS_Time/Start/Gregorian/Day_of_Month) &lt; 2">
          <xsl:value-of select="concat('0', Time_Information/GPS_Time/Start/Gregorian/Day_of_Month)"/>
	</xsl:when>
	<xsl:otherwise>
          <xsl:value-of select="Time_Information/GPS_Time/Start/Gregorian/Day_of_Month"/>
	</xsl:otherwise>
      </xsl:choose>    
    </xsl:variable>
    <xsl:variable name="RM_First_Epoch_Date" select="concat($RM_First_Epoch_Year, $Date_Separator, $RM_First_Epoch_Month, $Date_Separator, $RM_First_Epoch_Day)"/>
    <xsl:variable name="RM_First_Epoch_Hour">
      <xsl:choose>
        <xsl:when test="string-length(Time_Information/GPS_Time/Start/Gregorian/Hour) &lt; 2">
          <xsl:value-of select="concat('0', Time_Information/GPS_Time/Start/Gregorian/Hour)"/>
	</xsl:when>
	<xsl:otherwise>
          <xsl:value-of select="Time_Information/GPS_Time/Start/Gregorian/Hour"/>
	</xsl:otherwise>
      </xsl:choose>    
    </xsl:variable>
    <xsl:variable name="RM_First_Epoch_Minute">
      <xsl:choose>
        <xsl:when test="string-length(Time_Information/GPS_Time/Start/Gregorian/Minute) &lt; 2">
          <xsl:value-of select="concat('0', Time_Information/GPS_Time/Start/Gregorian/Minute)"/>
	</xsl:when>
	<xsl:otherwise>
          <xsl:value-of select="Time_Information/GPS_Time/Start/Gregorian/Minute"/>
	</xsl:otherwise>
      </xsl:choose>    
    </xsl:variable>
    <xsl:variable name="RM_First_Epoch_Second">
      <xsl:choose>
        <xsl:when test="string-length(Time_Information/GPS_Time/Start/Gregorian/Second) &lt; 2">
          <xsl:value-of select="concat('0', Time_Information/GPS_Time/Start/Gregorian/Second)"/>
	</xsl:when>
	<xsl:otherwise>
          <xsl:value-of select="Time_Information/GPS_Time/Start/Gregorian/Second"/>
	</xsl:otherwise>
      </xsl:choose>    
    </xsl:variable>
    <xsl:variable name="RM_First_Epoch_Time" select="concat($RM_First_Epoch_Hour, $Time_Separator, $RM_First_Epoch_Minute, $Time_Separator, $RM_First_Epoch_Second)"/>

    <xsl:variable name="RM_Time_System" select="Time_Information/System"/>
    <xsl:variable name="RM_First_Epoch" select="concat($RM_First_Epoch_Date, $RM_First_Epoch_Separator, $RM_First_Epoch_Time, $RM_First_Epoch_Separator, $RM_Time_System)"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="$RM_First_Epoch"/>
      <xsl:with-param name="Empty_Field" select="$Empty_First_Epoch_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Transformation_Direction_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Transformation_Direction" select="Transformation/Direction"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value"       select="$Transformation_Direction"/>
      <xsl:with-param name="Empty_Field"       select="$Empty_Transformation_Direction_Field"/>
      <xsl:with-param name="Justify"           select="$Left_Justified"/>
      <xsl:with-param name="Append_If_Shorter" select="'no'"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Pole_File_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Pole_File" select="Pole_File"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value"       select="$Pole_File"/>
      <xsl:with-param name="Empty_Field"       select="$Empty_Pole_File_Field"/>
      <xsl:with-param name="Justify"           select="$Left_Justified"/>
      <xsl:with-param name="Append_If_Shorter" select="'no'"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Nutation_Model_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Nutation_Model" select="Nutation/Model"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value"       select="$Nutation_Model"/>
      <xsl:with-param name="Empty_Field"       select="$Empty_Nutation_Model_Field"/>
      <xsl:with-param name="Justify"           select="$Left_Justified"/>
      <xsl:with-param name="Append_If_Shorter" select="'no'"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>
    
    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Nutation_Offsets_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Nutation_Offsets" select="Nutation/Offsets"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value"       select="$Nutation_Offsets"/>
      <xsl:with-param name="Empty_Field"       select="$Empty_Nutation_Offsets_Field"/>
      <xsl:with-param name="Justify"           select="$Left_Justified"/>
      <xsl:with-param name="Append_If_Shorter" select="'no'"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $Subdaily_Model_Keyword, ':')"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:variable name="Subdaily_Model" select="Subdaily_Model"/>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value"       select="$Subdaily_Model"/>
      <xsl:with-param name="Empty_Field"       select="$Empty_Subdaily_Model_Field"/>
      <xsl:with-param name="Justify"           select="$Left_Justified"/>
      <xsl:with-param name="Append_If_Shorter" select="'no'"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

    <xsl:call-template name="Fixed_Length_Field">
      <xsl:with-param name="Field_Value" select="concat($Comment_Indicator, ' ', $End_Of_Header_Keyword)"/>
      <xsl:with-param name="Empty_Field" select="$Empty_Header_Keyword_Field"/>
      <xsl:with-param name="Justify"     select="$Left_Justified"/>
    </xsl:call-template>

    <xsl:text>&#xa;</xsl:text>

  </xsl:template>

</xsl:stylesheet>
