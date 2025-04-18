%--------------------------------------------------------------------------
% Function to read in GOCE STR_VC2/3_1b products
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
% $Log: Read_STR_L1b_ds.m,v $
% Revision 1.2  2011/03/14 13:15:04  cmfuser
% header fixed
%
% Revision 1.1.1.1  2011/03/14 12:56:38  cmfuser
% SpuriousTracks repository  creation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: The routine reads the STR_VC2/3_1b datasets of the products
%              in the current directory.
%
% INPUTS: One or more GOCE *STR_VC2/3_1b*.EEF files
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
% NOTES ON USAGE: The routine is called by Read_STR_VC2/3_1b_dir.m routine.  
%                 [GPS_Time,Q1,Q2,Q3,Q4,Val_Flag,STR_ID,STR_Stars,CID,Loc_Time,BBO,TRS,Temp_Range,Asc_Flag,Orb_Flag,Seq_Flag,Est_Conf,Loc_Num]= Read_STR_L1b_ds(p,d(i).name);

function [GPS_Time,Q1,Q2,Q3,Q4,Val_Flag,STR_ID,STR_Stars,CID,Loc_Time,BBO,TRS,Temp_Range,Asc_Flag,Orb_Flag,Seq_Flag,Est_Conf,Loc_Num]= Read_STR_L1b_ds(p,f)
fid=fopen(cat(2,p,f),'r');

k=1;
l=1;
m=1;
n=1;
j=1;
g=1;
p=1;
q=1;
r=1;
s=1;
t=1;
u=1;
v=1;
z=1;
x=1;
w=1;
h=1;
e=1;
READ=0;
while 1
    tline=fgetl(fid);

    if tline==-1
        break
    end

    if (~isempty(strfind(tline,'<Q1>')))
        pos_string=tline;
        Q1(k,:)=str2num(pos_string(9:end-5));
        k=k+1;
    end

    if (~isempty(strfind(tline,'<Q2>')))
        vel_string=tline;
        Q2(l,:)=str2num(vel_string(9:end-5));
        l=l+1;
    end

    if (~isempty(strfind(tline,'<Q3>')))
        sigma_string=tline;
        Q3(m,:)=str2num(sigma_string(9:end-5));
        m=m+1;
    end

    if (~isempty(strfind(tline,'<Q4>')))
        sigma_string=tline;
        Q4(n,:)=str2num(sigma_string(9:end-5));
        n=n+1;
    end
    
    if (~isempty(strfind(tline,'<Val_Flag>')))
        sigma_string=tline;
        Val_Flag(j,:)=str2num(sigma_string(14:end-11));
        j=j+1;
    end
    
    if (~isempty(strfind(tline,'<Stid>')))
        sigma_string=tline;
        STR_ID(g,:)=str2num(sigma_string(10:end-7));
        g=g+1;
    end
    
    if (~isempty(strfind(tline,'<Nr_Of_Stars>')))
        sigma_string=tline;
        STR_Stars(p,:)=str2num(sigma_string(17:end-14));
        p=p+1;
    end
    
    if (~isempty(strfind(tline,'<Cid>')))
       sigma_string=tline;
       CID(r,:)=str2num(sigma_string(9:end-6));
       r=r+1;
    end
    
    if (~isempty(strfind(tline,'<Loc_Time>')))
        sigma_string=tline;
        Loc_Time(s,:)=str2num(sigma_string(14:end-11));
        s=s+1;
    end

    if (~isempty(strfind(tline,'<Bbo_Flag>')))
        sigma_string=tline;
        BBO(t,:)=str2num(sigma_string(14:end-11));
        t=t+1;
    end
    
    if (~isempty(strfind(tline,'<Trs_Flag>')))
        sigma_string=tline;
        TRS(u,:)=str2num(sigma_string(14:end-11));
        u=u+1;
    end
    
    if (~isempty(strfind(tline,'<Temp_Out_Range_Flag>')))
        sigma_string=tline;
        Temp_Range(v,:)=str2num(sigma_string(25:end-22));
        v=v+1;
    end
    
    if (~isempty(strfind(tline,'<Asc_Tc_Flag>')))
        sigma_string=tline;
        Asc_Flag(z,:)=str2num(sigma_string(17:end-14));
        z=z+1;
    end
    
    if (~isempty(strfind(tline,'<Orb_Cor_Flag>')))
        sigma_string=tline;
        Orb_Flag(x,:)=str2num(sigma_string(18:end-15));
        x=x+1;
    end 
    
    if (~isempty(strfind(tline,'<Seq_Flag>')))
        sigma_string=tline;
        Seq_Flag(w,:)=str2num(sigma_string(14:end-11));
        w=w+1;
    end
    
    if (~isempty(strfind(tline,'<Est_Conf>')))
        sigma_string=tline;
        Est_Conf(h,:)=str2num(sigma_string(14:end-11));
        h=h+1;
    end
    
    if (~isempty(strfind(tline,'<Nr_Of_Locks>')))
        sigma_string=tline;
        Loc_Num(e,:)=str2num(sigma_string(17:end-14));
        e=e+1;
    end
    
    if (~isempty(strfind(tline,'<STR_VC2_1i>')) || ~isempty(strfind(tline,'<STR_VC3_1i>')))
        READ=1;
    end
    if (~isempty(strfind(tline,'<Tt_GPS_Str unit="s">')) && READ)
        sigma_string=tline;
        GPS_Time(q,:)=str2num(sigma_string(25:end-13));
        q=q+1;
    end
end
fclose(fid);
