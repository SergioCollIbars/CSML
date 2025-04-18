%--------------------------------------------------------------------------
% Function to read in GOCE SST_NOM_1b products
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
% $Log: Read_SSTI_L1b_ds.m,v $
% Revision 1.2  2011/03/14 13:15:04  cmfuser
% header fixed
%
% Revision 1.1.1.1  2011/03/14 12:56:38  cmfuser
% SpuriousTracks repository  creation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: The routine reads the SST_NOM_1b datasets of the products
%              in the current directory.
%
% INPUTS: One or more GOCE *SST_NOM_1b*.EEF files
%
% OUTPUTS: - Time GPS
%          - GOCE Position x,y,z (POS_PVT) 
%          - GOCE Velocity x,y,z (VEL_PVT)
%          - Position Covariance Matrix (COV_POS)
%          - Velocity Covariance Matrix (COV_VEL)
%
%
% NOTES ON USAGE: The routine is called by Read_SST_NOM_1b_dir.m routine.  
%                [TT_GPS_PVT, POS_PVT, VEL_PVT, COV_POS,COV_VEL]=Read_SSTI_L1b_ds(p,d(i).name);




function [TT_GPS_PVT_out, POS_PVT_out, VEL_PVT_out, COV_POS_out, COV_VEL_out]= Read_SSTI_L1b_ds(p,f)


fid=fopen(cat(2,p,f),'r');

k=1;
l=1;
n=1;
o=1;
p=1;
z=1;

while 1
    tline=fgetl(fid);

    if tline==-1
        break
    end

    if (~isempty(strfind(tline,'<SST_NAV_1i>'))) %get time tags
        tline=fgetl(fid);
        time_string=tline; %get time
        Time_nav(o,:)=str2num(time_string(1,26:end-9));
        tline=fgetl(fid);
        tline=fgetl(fid);
        tline=fgetl(fid);
        pos_string=tline; %get X NAV pos
        Pos_nav1=str2num(pos_string(1,14:end-4));
        tline=fgetl(fid);
        pos_string=tline;%get Y NAV pos
        Pos_nav2=str2num(pos_string(1,14:end-4));
        tline=fgetl(fid);
        pos_string=tline; %get Z NAV pos
        Pos_nav3=str2num(pos_string(1,14:end-4));
        Pos_nav_final(p,:)=[Pos_nav1 Pos_nav2 Pos_nav3];
        p=p+1;
        o=o+1;
    end

    if (~isempty(strfind(tline,'<SST_PVT_1i>'))) %get time tags
        tline=fgetl(fid);
        time_string=tline;
        Time(n,:)=str2num(time_string(1,26:end-9));
        n=n+1;
    end
    if (~isempty(strfind(tline,'<Position unit="m">')))
        pos_string=tline;
        pos(k,:)=str2num(pos_string(1,30:end-11));
        k=k+1;
    end

    if (~isempty(strfind(tline,'<SST_Vel unit="m/s">')))
        vel_string=tline;
        vel(l,:)=str2num(vel_string(1,29:end-10));
        l=l+1;
    end
    
    if (~isempty(strfind(tline,'<SST_Cov_Pt unit="m^2">')))
       pos_string=tline;
       tline=fgetl(fid);
       COV_POS(z,:,1)=str2num(tline(17:end-7));
       tline=fgetl(fid);
       COV_POS(z,:,2)=str2num(tline(17:end-7));
       tline=fgetl(fid);
       COV_POS(z,:,3)=str2num(tline(17:end-7));
       tline=fgetl(fid);
       COV_POS(z,:,4)=str2num(tline(17:end-7));
       tline=fgetl(fid);
       tline=fgetl(fid);
       tline=fgetl(fid);
       COV_VEL(z,:,1)=str2num(tline(17:end-7));
       tline=fgetl(fid);
       COV_VEL(z,:,2)=str2num(tline(17:end-7));
       tline=fgetl(fid);
       COV_VEL(z,:,3)=str2num(tline(17:end-7));
       z=z+1; 
    end
    

end

TT_GPS_PVT_out=Time;
POS_PVT_out=pos;
VEL_PVT_out=vel;
COV_POS_out=COV_POS;
COV_VEL_out=COV_VEL;
fclose(fid);
