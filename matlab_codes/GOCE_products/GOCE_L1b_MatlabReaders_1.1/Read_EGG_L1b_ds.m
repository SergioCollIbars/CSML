%--------------------------------------------------------------------------
% Function to read in GOCE EGG_NOM_1b products
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
% $Log: Read_EGG_L1b_ds.m,v $
% Revision 1.2  2011/03/14 13:15:04  cmfuser
% header fixed
%
% Revision 1.1.1.1  2011/03/14 12:56:38  cmfuser
% SpuriousTracks repository  creation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: The routine reads the EGG_NOM_1b datasets of the products
%              in the current directory.
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
% NOTES ON USAGE: The routine is called by Read_EGG_NOM_1b_dir.m routine.  
%                [TT_GPS, EGG, CCM, CDM, IAQ,Q_FLAGS]=Read_EGG_L1b_ds(p,d(i).name);



function [TT_GPS_out, EGG_out, CCM_out, CDM_out, IAQ_out, Q_FLAGS_out]= Read_EGG_L1b_ds(p,f)


fid=fopen(cat(2,p,f),'r');

k=1;
l=1;
n=1;
g=1;
h=1;
READ=0;

while 1
    tline=fgetl(fid);

    if tline==-1
        break
    end

    if (~isempty(strfind(tline,'<EGG_CCD_DS>')))
        READ=1;
    end
    
    if (~isempty(strfind(tline,'</EGG_CCD_DS>')))
        READ=0;
    end
    
    
    if (~isempty(strfind(tline,'<Tt_GPS unit="s">')) && READ)
        pos_string=tline;
        TT_GPS(n,:)=str2num(pos_string(22:end-9));
        tline=fgetl(fid);
        tline=fgetl(fid);
        CCM(n,1,:)=str2num(tline(21:end-4));
        tline=fgetl(fid);
        CCM(n,2,:)=str2num(tline(21:end-4));
        tline=fgetl(fid);
        CCM(n,3,:)=str2num(tline(21:end-4));
        tline=fgetl(fid);
        tline=fgetl(fid);
        tline=fgetl(fid);
        CDM(g,1,:)=str2num(tline(21:end-4));
        tline=fgetl(fid);
        CDM(g,2,:)=str2num(tline(21:end-4));
        tline=fgetl(fid);
        CDM(g,3,:)=str2num(tline(21:end-4));
        
        n=n+1;
        g=g+1;
    end
    
    if (~isempty(strfind(tline,'<EGG_GGT_DS>')))
        READ_GGT=1;
    end
    
    if (~isempty(strfind(tline,'</EGG_GGT_DS>')))
        READ_GGT=0;
    end
       
    
    if (~isempty(strfind(tline,'<U_G unit="1/s^2">')))
        pos_string=tline;
        EGG(k,:)=str2num(pos_string(23:end-6));
        tline=fgetl(fid);
        tline=fgetl(fid);
        tline=fgetl(fid);
        Q_FLAGS(h,:)=str2num(tline(11:end-7));
        k=k+1;
        h=h+1;
    end
    
    if (~isempty(strfind(tline,'<Q_Grad>')))
        pos_string=tline;
        IAQ(l,:)=str2num(pos_string(13:end-9));
        l=l+1;
    end
    
end

TT_GPS_out=TT_GPS;
EGG_out=EGG;
CCM_out=CCM;
CDM_out=CDM;
IAQ_out=IAQ;
Q_FLAGS_out=Q_FLAGS;

fclose(fid);

