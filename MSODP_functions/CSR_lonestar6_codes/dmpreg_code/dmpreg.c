#include <stdio.h>
#include <stdlib.h> 
#include <string.h>

#define PARNAM_LNGTH    10
#define VERSION_LNGTH   30
#define FILENAME_LNGTH  35

/* These are constants for REGRES/R/DUZ header files */
#define PROGNAM_LNGTH   10
#define RUNDES_LNGTH    280
#define DATIM_LNGTH     10
#define STANAM_LNGTH    10

#define SIZEINT   sizeof(long long)
#define SIZEDOUBLE  sizeof(double)
#define SIZECHAR  sizeof(char)

/* VERSION MODIFIED TO RUN IN MACOS. Sergio Coll-Ibars 03/05/2026 */

int main(int argc, char *argv[])
{
  int     me, i, j, k, npar, sescnt, sestot, pntr, extra,
    cnt, npar_arc, rowsize, sizem, rows, mode, Hrows;

  long long   narc, nsta, idsta, idsat;     
  double  re, flat, stacor[3];
  char    *filename, newpar[PARNAM_LNGTH+1];
  char    prognam[PROGNAM_LNGTH], 
    rundes[RUNDES_LNGTH], 
    date[DATIM_LNGTH], 
    time[DATIM_LNGTH], 
    stanam[STANAM_LNGTH];
  
  long long   ns, nf, np, nga, iarc, igrps;   
  long long   *igrpcnt, *igrptot, ityp, ista;
  double  stasat[6],
      *aprval_loc, *parsig_loc, *parsca_loc,
      *aprval_loc_all, *parsig_loc_all, *parsca_loc_all;
  double  utc, *H, *tmp, *y, sigma, arctim; 

  char    *parnam_loc, *parnam_loc_all;
  FILE          *datafile;

  if ( argc == 1 || argc > 4 || strcmp( argv[1] , "--help") == 0 ){
    printf("\nUsage:  dmpreg <REGRES file> <mode (optional)> <# of H rows (optional)>\n");
    printf("\nPurpose: Dump parameters and partials from regres file\n");
    printf("\nOutput: To screen or redirect to file with >\n");
    exit(1);
  }

  filename = argv[1];
  
  mode = 1;
  if (argc > 2) mode = atoi(argv[2]);

  Hrows = 10;
  if (mode == 2) Hrows = -1;
  else if (argc == 4 && mode == 3) Hrows = atoi(argv[3]);

  /* CHANGE 1: Changed fopen64 to fopen for macOS compatibility */
  if ( (datafile = fopen(filename, "rb")) == NULL ){
    fprintf(stderr, "Error opening file %s!\n", filename);
    exit(1);
  }
  else{
    fprintf(stdout, "Processing file %s\n", filename); 
  }
   
  /* Read File Header */
  fread( prognam, SIZECHAR, PROGNAM_LNGTH, datafile );
  printf("\n\nProgram Name:\t%-.10s\n",prognam);
  fread( rundes, SIZECHAR, RUNDES_LNGTH, datafile );
  printf("Run Descr:\t%-.70s\n",&rundes[0]); 
  fread( date, SIZECHAR, DATIM_LNGTH, datafile );
  fread( time, SIZECHAR, DATIM_LNGTH, datafile );
  printf("\nDate/Time:\t%-.10s,\t%-.10s\n", date, time );

  fread( &narc, SIZEINT, 1, datafile );
  fread( &nsta, SIZEINT, 1, datafile );

  if (nsta != 0){
    for (i=0;i<nsta;i++) fread( stanam, SIZECHAR, STANAM_LNGTH, datafile );
    for (i=0;i<nsta;i++) fread( &idsta, SIZEINT, 1, datafile );
    fread( &re, SIZEDOUBLE, 1, datafile );
    fread( &flat, SIZEDOUBLE, 1, datafile );
    for (i=0;i<nsta;i++) fread( &stacor, SIZEDOUBLE, 3, datafile );
  }

  /* Read Arc Header */
  fread( &ns, SIZEINT, 1, datafile );
  fread( &nf, SIZEINT, 1, datafile );
  fread( &np, SIZEINT, 1, datafile );
  fread( &nga, SIZEINT, 1, datafile );
  fread( &iarc, SIZEINT, 1, datafile );
  fread( &utc, SIZEDOUBLE, 1, datafile );
  fread( &igrps, SIZEINT, 1, datafile );

  if ( igrps > 0 ){
    igrpcnt = (long long*) calloc( (int) igrps , SIZEINT ); 
    igrptot = (long long*) calloc( (int) igrps , SIZEINT ); 

    for ( i=0 ; i<igrps ; i++ ) {
      fread( &igrpcnt[i], SIZEINT, 1, datafile );
      sescnt = sescnt + (int) igrpcnt[i];
    }
    for ( i=0 ; i<igrps ; i++ ) {
      fread( &igrptot[i], SIZEINT, 1, datafile );
      sestot = sestot + (int) igrptot[i];
    }
  }

  for ( i=0 ; i<ns ; i++ ) fread( &idsat, SIZEINT, 1, datafile );
  for ( i=0 ; i<ns ; i++ ) fread( &stasat, SIZEDOUBLE, 6, datafile );

  npar = 6*ns+np;
  npar_arc = npar - sescnt + sestot;

  aprval_loc = (double *) calloc( npar , SIZEDOUBLE );
  aprval_loc_all = (double *) calloc( npar_arc , SIZEDOUBLE );
  parnam_loc = (char *) calloc( npar , PARNAM_LNGTH );
  parnam_loc_all = (char *) calloc( npar_arc , PARNAM_LNGTH );
  parsig_loc = (double *) calloc( npar , SIZEDOUBLE );
  parsig_loc_all = (double *) calloc( npar_arc , SIZEDOUBLE );
  parsca_loc = (double *) calloc( npar , SIZEDOUBLE );
  parsca_loc_all = (double *) calloc( npar_arc , SIZEDOUBLE );

  fread( aprval_loc, SIZEDOUBLE, npar, datafile );
  fread( parnam_loc, SIZECHAR, PARNAM_LNGTH*(npar), datafile );
  fread( parsig_loc, SIZEDOUBLE, npar, datafile );
  fread( parsca_loc, SIZEDOUBLE, npar, datafile );

  sescnt = 0;
  sestot = 0;
  cnt = 0;

  if (igrps > 0 ){
    for ( i=0 ; i<igrps ; i++ ){
      for ( j=1 ; j<=igrptot[i] ; j++ ){
        sprintf(newpar, "%-6.6s%4i", &parnam_loc[sescnt * PARNAM_LNGTH], j);
        for ( k=0 ; k<PARNAM_LNGTH ; k++ )  parnam_loc_all[cnt * PARNAM_LNGTH + k] = newpar[k];
        aprval_loc_all[cnt] = aprval_loc[sescnt];
        parsig_loc_all[cnt] = parsig_loc[sescnt];
        parsca_loc_all[cnt] = parsca_loc[sescnt];
        cnt++;
      }
      sescnt = sescnt + (int) igrpcnt[i];
      sestot = sestot + (int) igrptot[i];
    }  
  }

  for ( i=sescnt ; i<npar ; i++ ){
    for (j=0 ; j<PARNAM_LNGTH ; j++ ) parnam_loc_all[cnt * PARNAM_LNGTH + j] = parnam_loc[i * PARNAM_LNGTH + j];
    aprval_loc_all[cnt] = aprval_loc[i];
    parsig_loc_all[cnt] = parsig_loc[i];
    parsca_loc_all[cnt] = parsca_loc[i];
    cnt++;
  }

  for ( i=0 ; i<npar_arc ; i++ ){
    if ( mode == 1){
      if ( i < 10 || npar_arc - i < 10 ){
        printf("%-.10s          %22.15le %10.5le %10.5le\n", &parnam_loc_all[i*PARNAM_LNGTH], 
                 aprval_loc_all[i], parsig_loc_all[i], parsca_loc_all[i] );
      }
    }
    else if ( mode == 2 || mode == 3){
      printf("%-.10s          %22.15le %10.5le %10.5le\n", &parnam_loc_all[i*PARNAM_LNGTH], 
                 aprval_loc_all[i], parsig_loc_all[i], parsca_loc_all[i] );
    }
  }

  extra = 5;
  rowsize=npar+extra+igrps;
  npar_arc = npar - sescnt + sestot;

  tmp = ( double * ) calloc( rowsize, SIZEDOUBLE );
  H = ( double * ) calloc( npar_arc, SIZEDOUBLE );
  y = ( double * ) calloc( 1 , SIZEDOUBLE );

  if ( H == NULL || y == NULL ){
   printf("Error in dmpreg.c: not enough memory to allocate H and y!\n");
   return 1;
  }

  sizem = Hrows; 

  for ( rows=0 ; sizem<0 || rows<sizem ; rows++ ){
    for ( i=0 ; i<npar_arc ; i++ ) H[i] = 0.0;
    for ( i=0 ; i<rowsize ; i++ ) tmp[i] = 0.0;

    fread( tmp, SIZEDOUBLE, rowsize, datafile );

    arctim = tmp[0];
    ityp = (long long) tmp[1];
    ista = (long long) tmp[2];
    *y = tmp[3];
    sigma = tmp[4];

    if ( feof(datafile) ) break;

    sestot = 0;
    sescnt = 0;
    if ( igrps > 0 ){
      for ( i=0 ; i<igrps ; i++ ){
        pntr = (int) tmp[extra+i] + sestot - 1;
        for ( j=0 ; j<igrpcnt[i] ; j++ ) H[ pntr+j ] = tmp[extra + igrps + sescnt + j];
        sestot = sestot + (int) igrptot[i];
        sescnt = sescnt + (int) igrpcnt[i];
      }
    }

    for ( i = sescnt ; i < npar ; i++ ) H[ sestot + (i-sescnt) ] = tmp[extra + igrps + i ];
    
    if (mode == 1 || mode == 2) printf("% 9.1lf  % 6.5le % 6.5le % 6.5le % 6.5le % 6.5le % 6.5le % 6.5le % 6.5le\n",arctim,H[0],H[1],H[2],H[npar_arc-3],H[npar_arc-2],H[npar_arc-1],*y,sigma);
    else if (mode == 3){ 
      printf("% 9.1lf  ",arctim);
      for ( j=0 ; j<npar_arc ; j++ ) printf("% 6.5le ",H[j]);
      printf("% 6.5le % 6.5le\n",*y,sigma);
    }
  }

  free(aprval_loc);
  free(parsig_loc);
  free(parsca_loc);
  free(parnam_loc);

  /* CHANGE 2: Added 'return 0;' so main correctly returns an int */
  return 0;
}