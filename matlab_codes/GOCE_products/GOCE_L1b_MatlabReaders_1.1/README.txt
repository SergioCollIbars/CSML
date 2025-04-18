GOCE Level-1b MATLAB READERS
COPYRIGHT (c) 2011 ESA/ESRIN 


This is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 2, as published by the Free Software Foundation. 

The software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

The readers have been developed by GOCE PDGS Team in ESA-ESRIN.

--------------------------------------------------------------------------------
 $Date: 2011/03/14 13:07:09 $
 $Revision: 1.3 $

 $Log: README.txt,v $
 Revision 1.3  2011/03/14 13:07:09  cmfuser
 init fixed

 Revision 1.2  2011/03/14 13:06:15  cmfuser
 dos2unix applied

 Revision 1.1.1.1  2011/03/14 12:56:38  cmfuser

--------------------------------------------------------------------------------

PACKAGE CONTENT: 
Read_EGG_NOM_1b_dir.m
Read_SSTI_L1b_dir.m
Read_STR_VC2_1b_dir.m
Read_STR_VC3_1b_dir.m
Read_EGG_L1b_ds.m
Read_SSTI_L1b_ds.m
Read_STR_L1b_ds.m
UserManual.pdf
README.txt - This file. 

DESCRIPTION:
The enclosed routines read one or more standard GOCE .EEF files. 

USAGE:
The EEF Files to be read must be located in the same directory as the routines. 

Read_EGG_NOM_1b_dir.m: 
Reads one or more EGG_NOM_1b products located in the local directory. 
Requires :  Read_EGG_L1b_ds.m

Read_SSTI_L1b_dir.m: 
Reads one or more SST_NOM_1b products located in the local directory. 
Requires :  Read_SSTI_L1b_ds.m 

Read_STR_VC2_1b_dir.m: 
Reads one or more STR_VC2_1b products located in the local directory. 
Requires Read_STR_L1b_ds.m 

Read_STR_VC3_1b_dir.m: 
Reads STR_VC3_1b products located in the local directory. 
Requires Read_STR_L1b_ds.m 

Notes

The *_ds.m routines cannot be used as stand alone readers. They are called by the related *_dir.m routines.


