%--------------------------------------------------------------------------
% Script file for reading consecutive EGG_NOM_1b products.
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
% $Log: Read_EGG_NOM_1b_dir.m,v $
% Revision 1.2  2011/03/14 13:15:04  cmfuser
% header fixed
%
% Revision 1.1.1.1  2011/03/14 12:56:38  cmfuser
% SpuriousTracks repository  creation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: The routine reads all EGG_NOM_1b products in current directory.
%
% INPUTS: One or more GOCE *EGG_NOM_1b*.EEF files
%
% OUTPUTS: - Time GPS
%          - Common (CCM) and Differential (CDM) accelerations for each accelerometer pair 
%            and each degree of freedom.
%          - Gravity Gradient Tensor components xx, yy, zz, xy, xz, yz (EGG)
%          - Gradiometer Inertial Attitude Quaternions (IAQ)
%          - Quality Flags on proof mass Control Voltages (Q_FLAGS)
%
% REQUIREMENTS: Read_EGG_L1b_ds.m Matlab routine in the same directory 
%               of Read_EGG_NOM_1b_dir.m routine.
%
% NOTES ON USAGE: launch in a MATLAB shell the Read_EGG_NOM_1b_dir.m
%                 routine.


d=dir('*EGG_NOM_1b*.EEF');

if ispc==1
%Windows OS    
    p='.\';
else
%Linux/Mac OS
    p='./';
end

for i=1:length(d)
    disp(cat(2,'Processing file ',num2str(i), ' of ',num2str(length(d))));
    
    [TT_GPS, EGG, CCM, CDM, IAQ, Q_FLAGS]=Read_EGG_L1b_ds(p,d(i).name);
    
    eval(cat(2,'TT_GPS_',num2str(i),'=TT_GPS;'));
    eval(cat(2,'CCM_',num2str(i),'=CCM;'));
    eval(cat(2,'CDM_',num2str(i),'=CDM;'));
    eval(cat(2,'EGG_',num2str(i),'=EGG;'));
    eval(cat(2,'IAQ_',num2str(i),'=IAQ;'));
    eval(cat(2,'Q_FLAGS_',num2str(i),'=Q_FLAGS;'));
    
    if i==1
        EGG_final=EGG;
        TT_GPS_final=TT_GPS;
        CCM_final=CCM;
        CDM_final=CDM;
        IAQ_final=IAQ;
        Q_FLAGS_final=Q_FLAGS;
        
        eval(cat(2,'clear ','EGG_',num2str(i)));
        eval(cat(2,'clear ','TT_GPS_',num2str(i)));
        eval(cat(2,'clear ','CCM_',num2str(i)));
        eval(cat(2,'clear ','CDM_',num2str(i)));
        eval(cat(2,'clear ','IAQ_',num2str(i)));
        eval(cat(2,'clear ','Q_FLAGS_',num2str(i)));
        eval(cat(2,'clear ','TT_GPS_GGT_',num2str(i)));
    else
        eval(cat(2,'EGG_final=','[EGG_final;EGG_',num2str(i),'];'))
        eval(cat(2,'TT_GPS_final=','[TT_GPS_final;TT_GPS_',num2str(i),'];'))
        eval(cat(2,'CCM_final=','[CCM_final;CCM_',num2str(i),'];'))
        eval(cat(2,'CDM_final=','[CDM_final;CDM_',num2str(i),'];'))
        eval(cat(2,'IAQ_final=','[IAQ_final;IAQ_',num2str(i),'];'))
        eval(cat(2,'Q_FLAGS_final=','[Q_FLAGS_final;Q_FLAGS_',num2str(i),'];'))
        
        eval(cat(2,'clear ','EGG_',num2str(i)));
        eval(cat(2,'clear ','TT_GPS_',num2str(i)));
        eval(cat(2,'clear ','CCM_',num2str(i)));
        eval(cat(2,'clear ','CDM_',num2str(i)));
        eval(cat(2,'clear ','IAQ_',num2str(i)));
        eval(cat(2,'clear ','Q_FLAGS_',num2str(i)));
        
    end

    clear EGG TT_GPS CCM CDM IAQ Q_FLAGS;


end
