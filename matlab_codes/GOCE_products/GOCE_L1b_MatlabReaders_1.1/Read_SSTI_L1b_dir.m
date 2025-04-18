%--------------------------------------------------------------------------
% Script file for reading consecutive SST_NOM_1b products.
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
% $Log: Read_SSTI_L1b_dir.m,v $
% Revision 1.2  2011/03/14 13:15:04  cmfuser
% header fixed
%
% Revision 1.1.1.1  2011/03/14 12:56:38  cmfuser
% SpuriousTracks repository  creation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: The routine reads all all SST_NOM_1b products in current directory.
%
% INPUTS: One or more GOCE *SST_NOM_1b*.EEF files
%
% OUTPUTS: - Time GPS
%          - GOCE Position x,y,z (POS_PVT) 
%          - GOCE Velocity x,y,z (VEL_PVT)
%          - Position Covariance Matrix (COV_POS)
%          - Velocity Covariance Matrix (COV_VEL)
%
% REQUIREMENTS: Read_SSTI_L1b_ds.m Matlab routine in the same directory 
%               of Read_SSTI_NOM_1b_dir.m routine.
%
% NOTES ON USAGE: launch in a MATLAB shell the Read_SST_NOM_1b_dir.m
%                 routine.
clear;
clc;
close all;

% constants
Re = 6378.1370E3; % [m]

addpath("data/")
clear;
clc;

d=dir('*SST_NOM_1b*.EEF');

if ispc==1
%Windows OS    
    p='.\';
else
%Linux/Mac OS
    p='./';
end


for i=1:length(d)
    
    disp(cat(2,'Processing file ',num2str(i), ' of ',num2str(length(d))));
    
    [TT_GPS_PVT, POS_PVT, VEL_PVT, COV_POS, COV_VEL]=Read_SSTI_L1b_ds(p,d(i).name);
    
    eval(cat(2,'TT_GPS_PVT_',num2str(i),'=TT_GPS_PVT;'));
    eval(cat(2,'POS_PVT_',num2str(i),'=POS_PVT;'));
    eval(cat(2,'VEL_PVT_',num2str(i),'=VEL_PVT;'));
    eval(cat(2,'COV_POS_',num2str(i),'=COV_POS;'));
    eval(cat(2,'COV_VEL_',num2str(i),'=COV_VEL;'));
    
    if i==1
        TT_GPS_PVT_final=TT_GPS_PVT;
        POS_PVT_FINAL=POS_PVT;
        VEL_PVT_FINAL=VEL_PVT;
        COV_POS_FINAL=COV_POS;
        COV_VEL_FINAL=COV_VEL;
        
        eval(cat(2,'clear ','TT_GPS_PVT_',num2str(i)));
        eval(cat(2,'clear ','POS_PVT_',num2str(i)));
        eval(cat(2,'clear ','VEL_PVT_',num2str(i)));
        eval(cat(2,'clear ','COV_POS_',num2str(i)));
        eval(cat(2,'clear ','COV_VEL_',num2str(i)));
        
    else
        eval(cat(2,'TT_GPS_PVT_final=','[TT_GPS_PVT_final;TT_GPS_PVT_',num2str(i),'];'))
        eval(cat(2,'POS_PVT_FINAL=','[POS_PVT_FINAL; POS_PVT_',num2str(i),'];'))
        eval(cat(2,'VEL_PVT_FINAL=','[VEL_PVT_FINAL; VEL_PVT_',num2str(i),'];'))
        eval(cat(2,'COV_POS_FINAL=','[COV_POS_FINAL; COV_POS_',num2str(i),'];'))
        eval(cat(2,'COV_VEL_FINAL=','[COV_VEL_FINAL; COV_VEL_',num2str(i),'];'))
        
        eval(cat(2,'clear ','TT_GPS_PVT_',num2str(i)));
        eval(cat(2,'clear ','POS_PVT_',num2str(i)));
        eval(cat(2,'clear ','VEL_PVT_',num2str(i)));
        eval(cat(2,'clear ','COV_POS_',num2str(i)));
        eval(cat(2,'clear ','COV_VEL_',num2str(i)));
       
    end
    clear TT_GPS_PVT POS_PVT VEL_PVT COV_POS COV_VEL;
      
end

% remove outliers from data level 1
Nt = length(TT_GPS_PVT_final);
for j = 1:Nt
    radius = vecnorm(POS_PVT_FINAL(j ,:));
    H = (radius - Re)./1E3; % [Km]
    if (H > 500)
        POS_PVT_FINAL(j, :) = NaN * ones(3, 1);
        VEL_PVT_FINAL(j, :) = NaN * ones(3, 1);
    end
end

% Save data
ch = "16_Nov_2012";
matName = ["pos", "vel", "posCov", "velCov", "time"];
var     = ["POS_PVT_FINAL", "VEL_PVT_FINAL", "COV_POS_FINAL", "COV_VEL_FINAL", "TT_GPS_PVT_final"]; 
for j = 1:5
    name = matName(j) + "_" + ch + ".mat";
    save("data/"+name, var(j));
end

