#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PARNAM_LNGTH    10
#define VERSION_LNGTH   30

#define TRUE   1

#define FILENAME_LNGTH  35

/* These are constants for REGRES/R/DUZ header files */
#define PROGNAM_LNGTH   10
#define RUNDES_LNGTH    280
#define DATIM_LNGTH     10
#define STANAM_LNGTH    10

#define SIZEINT		sizeof(long long)
#define SIZEDOUBLE	sizeof(double)
#define SIZECHAR	sizeof(char)

struct resid_data {
			char                    orig_name[PARNAM_LNGTH+1];
                        char                    new_name[PARNAM_LNGTH+1];
                        struct resid_data       *nextitem;
                  };
typedef struct resid_data       RESID;
typedef RESID       *RESIDLINK;


void main(int argc, char *argv[])
{
  /*------------------------------------------------------------------------------|

    Name:       rename_params_in_reg

    Purpose:    Rename parameters in a given reg file

    Input:      A reg file and a file containing old and new names
		Default filetype is REGRES

    Output:     A new reg file with changed parameter names

    Return:

  |-------------------------------------------------------------------------------|

    History:
    Date		Author			Comment
    --------    -------------   ----------------------------------|
    05/04/2017	Furqan Ahmed	Created by changing split_reg


  |------------------------------------------------------------------------------*/

  int 		changeCount, s, m,e, i, j, k, npar, sescnt, sestot, pntr, extra,
		cnt, npar_arc, rowsize, sizem, rows;

  long long 	narc, nsta, idsta, idsat; 		/* File header variables */
  double 	re, flat, stacor[3];
  char          inpFlag;
  char          *filein, *fileout, *parnamefile, *difffilename, newpar[PARNAM_LNGTH+1];
  char          tmppar[PARNAM_LNGTH+1];
  char 		prognam[PROGNAM_LNGTH], 
		rundes[RUNDES_LNGTH], 
		date[DATIM_LNGTH], 
		time[DATIM_LNGTH], 
		stanam[STANAM_LNGTH];
  
  long long 	ns, nf, np, nga, iarc, igrps;		/* Arc header variables */
  long long 	*igrpcnt, *igrptot, ityp, ista;
  double	stasat[6],
  		*aprval_loc, *parsig_loc, *parsca_loc,
  		*aprval_loc_all, *parsig_loc_all, *parsca_loc_all;
  double 	utc, *H, *tmp, *y, sigma, arctim,
		diff, int_part, frac_part; 

  char 		*parnam_loc, *parnam_loc_all, plus[2], line[RUNDES_LNGTH];
  FILE          *parfile, *datafilein, *datafileout, *difffile;
  //RESIDLINK	difftail, difflist = NULL;
  RESIDLINK     tempstruct, partail, partail_prev, parlist = NULL;

  int nTotalObs = 0;    //Number of total observations read from the input file
  int nWrittenObs = 0;	//Number of total observations written to the output file
  char singleSpace = " ";
  /*-----------------------------------------------------------------------------*/


 /* Help output if input argument is "-h", "--help" or empty */
 
  if(argc==1 || strcmp(argv[1],"-h")==0 || strcmp(argv[1],"--help")==0){
	printf("\nProgram: rename_params_in_reg\n");
	printf("==================\n");
	printf("Created: 2017-05-04 by Furqan Ahmed\n");
	printf("\nPURPOSE: Rename parameters in a given REG file and create a new REG files containing new parameter names.\n");
	printf("\nUSAGE: rename_params_in_reg <file_in> <file_out> <param_name_file>\n\n");
	printf("\twhere\n\n");
	printf("\t<file_in> is the full path to the input REG file.\n");
	printf("\t<file_out> is the full path to the output REG file.\n");
	printf("\t<param_name_file> is the full path of a text file containing the list of parameter names to be changed. The file should consist of two columns that are each exactly 10 characters in width, separated by a space.\n\tThe first column should be the original name of the parameter to be changed, and the second column should" 
" contain the new parameter name. For example:\n"
"\tSST1  AC1X SST1  XXXX\n"
"\tSST1  AC1Y SST1  YYYY\n"
"\tSST1  AC1Z SST1  ZZZZ\n"
"\twould change the parameter \'SST1  AC1X\' to \'SST1  XXXX\', and so on.\n");

	printf("\nOUTPUT: a REG file with path and name specified by <file_out>, containing parameters with new names. All the parameters which are not specified in <param_name_file> are written with their original names in <file_out>.\n\n");

	printf("EXAMPLE:\n\n");
	printf("\tTo run: rename_params_in_reg grcA_16-08-08_RL00_GPSRL_RL05UD.354092.reg newReg.reg par2change.txt\n\n");
	printf("\tOutput file obtained: newReg.reg (containing parameters with new names).\n\n");


	exit(1);

  }


  filein = argv[1];


  /*Read input file into datafilein*/
  if ( (datafilein = fopen64(filein, "rb")) == NULL ){
    fprintf(stderr, "Error opening file %s!\n", filein);
    exit(1);
  }
  else{
    fprintf(stdout, "Processing file %s\n", filein); 
  }

 /*Create output file name based on input argument file_out*/
 fileout = argv[2];
 parnamefile = argv[3]; 

 /* Read file containing the new parameter name map */
  parfile = fopen( parnamefile, "r" );


 while(1){

    /* Get next line and trim leading whitespace */
    fgets( line, RUNDES_LNGTH, parfile );
    if ( feof(parfile) ) break;

    if (parlist == NULL){
      parlist = (struct resid_data *) malloc(sizeof(RESID));
      partail = parlist;
    }
    else if ( partail -> nextitem == NULL ){
      partail -> nextitem = (struct resid_data *) malloc(sizeof(RESID));
      partail = partail -> nextitem;
    }

    sscanf( line, "%10c%1c%10c", &(partail->orig_name), &tmppar, &(partail->new_name) );

    /* Trim trailing white space */
    /* NOTE: have to do it this way b/c '\0' is not recognized as white space */

    for ( s=PARNAM_LNGTH-1 ; s>=0 ; s-- ){
    /* printf("\n   orig s = %i",s); */
      if ( partail->orig_name[s] != '\0' && !isspace( partail->orig_name[s] ) ){
        partail->orig_name[s+1] = '\0';
        break;
      }
    }

    for ( s=PARNAM_LNGTH-1 ; s>=0 ; s-- ){
    /* printf("\n   new  s = %i",s); */
      if ( partail->new_name[s] != '\0' && !isspace( partail->new_name[s] ) ){
        partail->new_name[s+1] = '\0';
        break;
      }
    }


  }

 /*Open fileout binary file in write mode*/
 if ( (datafileout = fopen64(fileout, "wb")) == NULL ){
    fprintf(stderr, "Error opening file %s!\n", fileout);
    exit(1);
  }

  /* Read File Header */
  fread( prognam, SIZECHAR, PROGNAM_LNGTH, datafilein );
  fwrite( prognam, SIZECHAR, PROGNAM_LNGTH, datafileout );

/*printf("\n\nProgram Name:\t%s\n",prognam); */
  printf("\n\nProgram Name:\t%-.10s\n",prognam);

  fread( rundes, SIZECHAR, RUNDES_LNGTH, datafilein );
  fwrite( rundes, SIZECHAR, RUNDES_LNGTH, datafileout );
  printf("Run Descr:\t%-.70s\n",&rundes[0]); 
  printf("\t\t%-.70s\n",&rundes[70]);
  printf("\t\t%-.70s\n",&rundes[140]);
  printf("\t\t%-.70s\n",&rundes[210]);
  fread( date, SIZECHAR, DATIM_LNGTH, datafilein );
  fwrite( date, SIZECHAR, DATIM_LNGTH, datafileout );
  fread( time, SIZECHAR, DATIM_LNGTH, datafilein );
  fwrite( time, SIZECHAR, DATIM_LNGTH, datafileout );
/*printf("\nDate/Time:\t%s,\t%s\n", date, time ); */
  printf("\nDate/Time:\t%-.10s,\t%-.10s\n", date, time );

  fread( &narc, SIZEINT, 1, datafilein );
  fwrite( &narc, SIZEINT, 1, datafileout );
  fread( &nsta, SIZEINT, 1, datafilein );
  fwrite( &nsta, SIZEINT, 1, datafileout );

  printf("narc: %i\n", (int) narc );
  printf("nsta: %i\n", (int) nsta );

  if (nsta != 0){
    for (i=0;i<nsta;i++) fread( stanam, SIZECHAR, STANAM_LNGTH, datafilein );
    for (i=0;i<nsta;i++) fwrite( stanam, SIZECHAR, STANAM_LNGTH, datafileout );
    for (i=0;i<nsta;i++) fread( &idsta, SIZEINT, 1, datafilein );
    for (i=0;i<nsta;i++) fwrite( &idsta, SIZEINT, 1, datafileout );
    fread( &re, SIZEDOUBLE, 1, datafilein );
    fwrite( &re, SIZEDOUBLE, 1, datafileout );
    fread( &flat, SIZEDOUBLE, 1, datafilein );
    fwrite( &flat, SIZEDOUBLE, 1, datafileout );
    printf("Re: %le\tFlat: %le\n", re,flat );
    for (i=0;i<nsta;i++) fread( &stacor, SIZEDOUBLE, 3, datafilein );
    for (i=0;i<nsta;i++) fwrite( &stacor, SIZEDOUBLE, 3, datafileout );
  }

  /* Read Arc Header */
  fread( &ns, SIZEINT, 1, datafilein );
  fwrite( &ns, SIZEINT, 1, datafileout );
  printf("\nNumber of satellites (ns): %i\n", (int) ns );
  fread( &nf, SIZEINT, 1, datafilein );
  fwrite( &nf, SIZEINT, 1, datafileout );
  printf("Number of force model parameters (ns) excl state: %i\n", (int) nf );
  fread( &np, SIZEINT, 1, datafilein );
  fwrite( &np, SIZEINT, 1, datafileout );
  printf("Number of parameters (np) excl state: %i\n", (int) np );
  fread( &nga, SIZEINT, 1, datafilein );
  fwrite( &nga, SIZEINT, 1, datafileout );
  printf("Number of measurement model params (nga): %i\n", (int) nga );
  fread( &iarc, SIZEINT, 1, datafilein );
  fwrite( &iarc, SIZEINT, 1, datafileout );
  fread( &utc, SIZEDOUBLE, 1, datafilein );
  fwrite( &utc, SIZEDOUBLE, 1, datafileout );
  printf("Arc number: %i\tArc epoch: %le\n", (int) iarc, utc );
  

  fread( &igrps, SIZEINT, 1, datafilein );
  fwrite( &igrps, SIZEINT, 1, datafileout );
  printf("Number of session groups: %i\n", (int) igrps );

  sescnt = 0;
  sestot = 0;

  if ( igrps > 0 ){

    igrpcnt = (long long*) calloc( (int) igrps , SIZEINT ); 
    igrptot = (long long*) calloc( (int) igrps , SIZEINT ); 

    for ( i=0 ; i<igrps ; i++ ){
      fread( &igrpcnt[i], SIZEINT, 1, datafilein );
      fwrite( &igrpcnt[i], SIZEINT, 1, datafileout );
      sescnt = sescnt + (int) igrpcnt[i];
    }
    for ( i=0 ; i<igrps ; i++ ){
      fread( &igrptot[i], SIZEINT, 1, datafilein );
      fwrite( &igrptot[i], SIZEINT, 1, datafileout );
      sestot = sestot + (int) igrptot[i];
    }
    
  }

  printf("Total number of session parameters: %i\n", sestot );

  for ( i=0 ; i<ns ; i++ ){
    fread( &idsat, SIZEINT, 1, datafilein );
    fwrite( &idsat, SIZEINT, 1, datafileout );
  }
  for ( i=0 ; i<ns ; i++ ){
    fread( &stasat, SIZEDOUBLE, 6, datafilein );
    fwrite( &stasat, SIZEDOUBLE, 6, datafileout );
  }

  /* npar is the number of params physically stored in the file */
  npar = 6*ns+np;
  /* npar_arc is the total number of params after the sessions parameters are unpacked */
  npar_arc = npar - sescnt + sestot;

  printf("Total number of parameters incl state (npar_arc): %i\n", npar_arc );

  aprval_loc = (double*)calloc(npar,SIZEDOUBLE );
  parnam_loc = (char*)calloc(npar,PARNAM_LNGTH );
  parsig_loc = (double*)calloc(npar,SIZEDOUBLE );
  parsca_loc = (double*)calloc(npar,SIZEDOUBLE );

  fread( aprval_loc, SIZEDOUBLE, npar, datafilein );
  fwrite( aprval_loc, SIZEDOUBLE, npar, datafileout );
  fread( parnam_loc, SIZECHAR, PARNAM_LNGTH*(npar), datafilein );

  cnt=0;

  tempstruct = parlist;
  char *oldN, *newN, tempN[PARNAM_LNGTH+1], tempN2[PARNAM_LNGTH+1];
  while(tempstruct != NULL){ //WhileX starts
  oldN = tempstruct->orig_name;
  newN = tempstruct->new_name;

  //Now search parnam_loc for each oldN and replace it with corresponding 10character newN
  for ( i=0; i<npar ; i++ ){ //i loop begins
    for (j=0 ; j<PARNAM_LNGTH ; j++ ) tempN[j] = parnam_loc[i * PARNAM_LNGTH + j];
  
  for ( s=PARNAM_LNGTH-1 ; s>=0 ; s-- ){
      if ( oldN[s] != '\0' && !isspace( oldN[s] ) ){
        oldN[s+1] = '\0';
        break;
      }
    }
  for ( s=PARNAM_LNGTH-1 ; s>=0 ; s-- ){
      if ( newN[s] != '\0' && !isspace( newN[s] ) ){
        newN[s+1] = '\0';
        break;
      }
    }
  for ( s=PARNAM_LNGTH-1 ; s>=0 ; s-- ){
      if ( tempN[s] != '\0' && !isspace( tempN[s] ) ){
        tempN[s+1] = '\0';
        break;
      }
    } 

   pad(newN, 10); //Append spaces if length < 10 characters

   if (strcmp(tempN,oldN)==0)
   { //replace element by element in this loop
      
       changeCount++;
       //printf("Match found %d\n",changeCount);
       for (j=0; j<PARNAM_LNGTH; j++) parnam_loc[i * PARNAM_LNGTH + j] = newN[j];

   }

   for (j=0 ; j<PARNAM_LNGTH ; j++ ) tempN2[j] = parnam_loc[i * PARNAM_LNGTH + j];

 

   } //i loop ends
  

  tempstruct=tempstruct->nextitem;
  } //WhileX ends


  fwrite( parnam_loc, SIZECHAR, PARNAM_LNGTH*(npar), datafileout );
  fread( parsig_loc, SIZEDOUBLE, npar, datafilein );
  fwrite( parsig_loc, SIZEDOUBLE, npar, datafileout );
  fread( parsca_loc, SIZEDOUBLE, npar, datafilein );
  fwrite( parsca_loc, SIZEDOUBLE, npar, datafileout );

  /* extra is the number of elements in addition to the partials in each record of the data file */
  extra = 5;
  rowsize=npar+extra+igrps;
  cnt = 0;


  tmp = ( double * ) calloc( rowsize, SIZEDOUBLE );


  /* The following is needed to round off diff so that the (long long) conversion
     further down works properly.  When C does a (long long) conversion of a float
     it simply drops the digits past the decimal.  So the value 29.999 turns into
     29, not 30.  This can cause the time comparison below to be off by one second. */

  diff = (utc - 2451544.5)*86400;
  frac_part = modf( diff, &int_part );

  if ( fabs( frac_part ) < 0.5 )
    diff = floor( diff );
  else
    diff = ceil( diff );

  while ( TRUE ){

    fread( tmp, SIZEDOUBLE, rowsize, datafilein );

    /* Test to see if we are at end of file */
    if ( feof(datafilein) ) break;


    else { 

                nTotalObs+=1;

                fwrite( tmp, SIZEDOUBLE, rowsize, datafileout );
		nWrittenObs+=1;
                                   


	} 

  }


  free(aprval_loc);
  free(parsig_loc);
  free(parsca_loc);
  free(parnam_loc);

  printf("\n");
  printf("Number of observations in the input file: %d\n",nTotalObs);
  printf("Number of observations written to the output file: %d\n\n",nWrittenObs);

  printf("A total of %d parameter names changed in this run.\n",changeCount);
 
}

void pad(char *s, int length)
{
  int l;
  l = strlen(s); /* its length */
  while(l<length) {
    s[l] = ' '; /* insert a space */
    l++;
  }
  s[l]= '\0'; /* strings need to be terminated in a null */
}

