%--------------------------------------------------------------------------
% Script file for reading consecutive STR_VC3_1b products.
%--------------------------------------------------------------------------
% 
% CODED: GOCE PDGS Team 
% CONTACTS: eohelp@esa.int 
%           
%    COPYRIGHT (c) 2011 ESA/ESRIN 
%    This is free software; you can redistribute it and/or modify it 
%    under the terms of the GNU General Public License, version 2, as 
%    published by the Free Software Foundation. 
% 
%    The software is distributed in the hope that it will be useful, but 
%    WITHOUT ANY WARRANTY; without even the implied warranty of 
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
%    General Public License for more details. 
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Date: 2011/03/14 13:15:04 $
% $Revision: 1.2 $
%
% $Log: Read_STR_VC3_1b_dir.m,v $
% Revision 1.2  2011/03/14 13:15:04  cmfuser
% header fixed
%
% Revision 1.1.1.1  2011/03/14 12:56:38  cmfuser
% SpuriousTracks repository  creation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: The routine reads all STR_VC3_1b products in current directory.
%
% INPUTS: One or more GOCE *STR_VC3_1b*.EEF files
%
% OUTPUTS: - Time GPS
%          - Star Tracker Attitude quaternions (Q1, Q2, Q3, Q4) 
%          - Identifier of the Start Sensor currently used, valid range 1 to 3 (STR_ID)
%          - Camera identifier, from which camera the attitude is measured
%            (CID)
%          - Flag indicating if the current attitude is valid 0=invalid,
%            1=valid (Val_Flag)
%          - Flag set to 1 if the ACS timestamps are given as free running
%            time since boot. It is set to 0 if the ASC has locked on to the
%            time received in the timing packets (Loc_Time)
%          - Flag indicating if the ASC has detected any Big Bright objects
%            0=no detection, 1=detection (BBO)
%          - Flag set to 1 if the ASC has received a tome reference signal.
%            Otherwise set to 0 (TRS)
%          - Flag set to 1 if one or more of the HK temperature
%            measurements show a temperature out of range, otherwise set to
%            0 (Temp_Range)
%          - Flag set to 1 after the ASC have received an AscTimeTC packet,
%            otherwise set to 0 (Asc_Flag)
%          - Flag set to 1 if the attitude determination SW uses any kind
%            of orbit correction, otherwise set to 0 (Orb_Flag)
%          - Flag set to 0 if initial and finetuning attitude determination
%            s are performed. Set to 1 if only a fine-tuning attitude det. is
%            performed (Seq_Flag)
%          - Estimate confidence field. Value from 0 to 255 (Est_Conf)
%          - Number of Locks. It is used only if Sequence field is 0. It
%            indicates how many stars are used to acquire the initial
%            attitude (Loc_Num)
%          - Actual number of stars found in the image (STR_Stars)
%
% REQUIREMENTS: Read_STR_L1b_ds.m Matlab routine in the same directory 
%               of Read_STR_VC3_1b_dir.m routine.
%
% NOTES ON USAGE: launch in a MATLAB shell the Read_STR_VC3_1b_dir.m
%                 routine.

d=dir('*STR_VC3*_1b*.EEF');


if ispc==1
%Windows OS    
    p='.\';
else
%Linux/Mac OS
    p='./';
end

for i=1:length(d)
    disp(cat(2,'Processing file ',num2str(i), ' of ',num2str(length(d))));
    [GPS_Time,Q1,Q2,Q3,Q4,Val_Flag,STR_ID,STR_Stars,CID,Loc_Time,BBO,TRS,Temp_Range,Asc_Flag,Orb_Flag,Seq_Flag,Est_Conf,Loc_Num]= Read_STR_L1b_ds(p,d(i).name);
    
    if i==1
        GPS_Time_final=GPS_Time;
        Q1_final=Q1;
        Q2_final=Q2;
        Q3_final=Q3;
        Q4_final=Q4;
        Val_Flag_final=Val_Flag;
        STR_ID_final=STR_ID;
        STR_Stars_final=STR_Stars;
        CID_final=CID;
        Loc_Time_final=Loc_Time;
        BBO_final=BBO;
        TRS_final=TRS;
        Temp_Range_final=Temp_Range;
        Asc_Flag_final=Asc_Flag;
        Orb_Flag_final=Orb_Flag;
        Seq_Flag_final=Seq_Flag;
        Est_Conf_final=Est_Conf;
        Loc_Num_final=Loc_Num;
    else
        GPS_Time_final=[GPS_Time_final; GPS_Time];
        Q1_final=[Q1_final; Q1];
        Q2_final=[Q2_final; Q2];
        Q3_final=[Q3_final; Q3];
        Q4_final=[Q4_final; Q4];
        Val_Flag_final=[Val_Flag_final; Val_Flag];
        STR_ID_final=[STR_ID_final; STR_ID];
        STR_Stars_final=[STR_Stars_final; STR_Stars];
        CID_final=[CID_final; CID];
        Loc_Time_final=[Loc_Time_final; Loc_Time];
        BBO_final=[BBO_final; BBO];
        TRS_final=[TRS_final; TRS];
        Temp_Range_final=[Temp_Range_final; Temp_Range];
        Asc_Flag_final=[Asc_Flag_final; Asc_Flag];
        Orb_Flag_final=[Orb_Flag_final; Orb_Flag];
        Seq_Flag_final=[Seq_Flag_final; Seq_Flag];
        Est_Conf_final=[Est_Conf_final; Est_Conf];
        Loc_Num_final=[Loc_Num_final; Loc_Num];
    end
    clear GPS_Time Q1 Q2 Q3 Q4 Val_Flag STR_ID STR_Stars CID Loc_Time BBO TRS Temp_Range Asc_Flag Orb_Flag Seq_Flag Est_Conf Loc_Num;


end
