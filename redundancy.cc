//
// redundancy.cc
//
// C++ Main Program for Redundancy Calculations 
//
//  $Id: redundancy.cc,v 1.1.1.1 1997/08/23 17:50:33 jak Exp $
//
//  Author: John Kassebaum
//
// $Log: redundancy.cc,v $
// Revision 1.1.1.1  1997/08/23 17:50:33  jak
// BoxCounter disappeared from CVS ... this is a replacement. -jak
//
// Revision 1.4  1995/11/08 05:28:57  jak
// Added back lost RCS variables. -jak
//
// Revision 1.1  1994/05/05  04:55:44  jak
// New files - added a much faster boxcounter, and the redundancy calculator. -jak
//
//
//

static char rcsid_redundancy_cc[] = "$Id: redundancy.cc,v 1.1.1.1 1997/08/23 17:50:33 jak Exp $";

#include "FastBoxCounter.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <octave/Matrix.h>

#define USAGE    "Usage: %s [-h] [-v] [-s] < ascii_data_file \n"
#define DESCR_1    " -h  prints this help message and exits.\n"
#define DESCR_2    " -v  verbosely prints out data\n"
#define DESCR_3    " -s  prints statistics\n"
#define DESCR_4_1  " -eo #  sets starting embedding dimension for search (default = 1)\n"
#define DESCR_4_2  " -ef #  sets ending embedding dimension for search (default = 8)\n"
#define DESCR_4_3  " -ed #  sets delta for embedding dimension search (default = 1)\n"
#define DESCR_5    " -do #  sets start delay in embedding search (default = 1)\n"
#define DESCR_6    " -df #  sets final delay in embedding search (default = 16)\n"
#define DESCR_7    " -dd #  sets delay delta in embedding search (default = 1)\n"
#define DESCR_8    " -q #  sets highest quantization level (default = 8)\n"
#define DESCR_9    " ascii_data_file  -  A file of floating point numbers\n"
#define DESCR_10   "     in ascii format, with one number per line. \n example:\n"
#define DESCR_11   " -47.52\n -39.59\n -27.72\n -15.84\n"

#define BLOCKSIZE    1024
#define Null(A)             ((A *) 0)
#define New(A)              ((A *) malloc( sizeof(A) ) )
#define NewBlock(A,N)       ((A *) malloc( sizeof(A) * (N)) )
#define BiggerBlock(A,B,N)  ((A *) realloc( (void *)(B), sizeof(A) * (N)))

main(int argc,char **argv)
{
    int c,i,j,k,m;
    int verbose, stats, flag;
    long int data_count;
    long int space_alloc;
    float min_value, max_value;
    double *tempf_p, epsilon;
    int step_count;
    float step;
    int num_of_levels;
    int embedding_start, embedding_finish, embedding_delta;
    int delay_start, delay_finish, delay_delta;
    float *data;
    double **redundancy;
	int red_maxdim, red_maxqlevel;
    FastBoxCounter *myBoxCounterp;

    verbose = 0;
    stats = 0;
    num_of_levels = 8;
    embedding_start = 1;    embedding_finish = 8;   embedding_delta = 1;
    delay_start = 1;        delay_finish = 16;      delay_delta = 1;
    for (c=1; c< argc; c++) {
        if ( !strcmp( argv[ c ],"-h") || !strcmp( argv[ c ],"-help")){
            fprintf(stderr,USAGE,argv[0]);
            fprintf(stderr,DESCR_1);
            fprintf(stderr,DESCR_2);
            fprintf(stderr,DESCR_3);
            fprintf(stderr,DESCR_4_1);
            fprintf(stderr,DESCR_4_2);
            fprintf(stderr,DESCR_4_3);
            fprintf(stderr,DESCR_5);
            fprintf(stderr,DESCR_6);
            fprintf(stderr,DESCR_7);
            fprintf(stderr,DESCR_8);
            fprintf(stderr,DESCR_9);
            fprintf(stderr,DESCR_10);
            fprintf(stderr,DESCR_11);
            exit(0);
        } else if (!strcmp( argv[ c ],"-v")){
            verbose = 1;
        } else if (!strcmp( argv[ c ],"-s")){
            stats = 1;
        } else if (!strcmp( argv[ c ],"-eo")){
            embedding_start = atoi( argv[++c] );
        } else if (!strcmp( argv[ c ],"-ef")){
            embedding_finish = atoi( argv[++c] );
        } else if (!strcmp( argv[ c ],"-ed")){
            embedding_delta = atoi( argv[++c] );
        } else if (!strcmp( argv[ c ],"-q")){
            num_of_levels = atoi( argv[++c] );;
        } else if (!strcmp( argv[ c ],"-do")){
            delay_start = atoi( argv[++c] );;
        } else if (!strcmp( argv[ c ],"-df")){
            delay_finish = atoi( argv[++c] );;
        } else if (!strcmp( argv[ c ],"-dd")){
            delay_delta = atoi( argv[++c] );;
        }
    }
    
    data_count = 0;
    space_alloc = 0;
    data = (float *)0;
    
    if( data = NewBlock( float, BLOCKSIZE ) ){
        space_alloc = BLOCKSIZE;
    } else {
        perror(argv[0]);
        abort();
    };
        
  // Get space for redundancy temporary storage
	red_maxdim = 1 + (embedding_finish - embedding_start + 1 ) / embedding_delta;
	red_maxqlevel = num_of_levels;
    if( redundancy = NewBlock( double *, red_maxdim ) ){
        for(i=0; i<red_maxdim; i++){
		    if( ( redundancy[i] = NewBlock( double, red_maxqlevel ) ) == 0){
			    fprintf(stderr,"%s: can't get space for redundancy[%d]! \n",argv[0],i);
				fprintf(stderr,"red_maxdim = %d, red_maxqlevel = %d\n", red_maxdim, red_maxqlevel);
				perror(argv[0]);
				abort();
			}
		}
    } else {
		fprintf(stderr,"%s: can't get space for redundancy! \n",argv[0]);
        perror(argv[0]);
        abort();
    };
	
    flag = 0;
    while ( scanf( "%f", &data[ data_count ]) != EOF ){
        if (!flag) {
            min_value = max_value = data[ data_count ];
            flag = 1;
        } else {
            if( data[ data_count ] < min_value )
                min_value = data[ data_count ];
            if( data[ data_count ] > max_value )
                max_value = data[ data_count ];
        }
        if (verbose) fprintf(stdout, "%f\n", data[ data_count ]);
        data_count++;
        if ( data_count >= space_alloc ){
            if( data = BiggerBlock( float, data, (space_alloc + BLOCKSIZE) ) ){
                space_alloc += BLOCKSIZE;
            } else {
                perror(argv[0]);
                abort();
            };
        };
    }
    
    if (stats) {
	    fprintf(stdout, "Redundancy calculations:\n");
	    fprintf(stdout, "Embedding delay search: from %d to %d by %d\n",delay_start, delay_finish, delay_delta);
	    fprintf(stdout, "Quantization Levels up to: %d bits\n", num_of_levels);
	    fprintf(stdout, "Embedding Dimensions searched: %d to %d by %d \n", embedding_start, embedding_finish, embedding_delta);
        fprintf(stdout, "Space Allocated for %d floats.\n", space_alloc);
        fprintf(stdout, "Space Used for %d floats.\n", data_count);
        fprintf(stdout, "Minimum value = %f.\n", min_value);
        fprintf(stdout, "Maximum value = %f.\n", max_value);
        fflush(stdout);
    }
        
    myBoxCounterp = new FastBoxCounter( num_of_levels );
	for(i=delay_start; i<= delay_finish; i += delay_delta){
	  //
	  //Calculate Redundancies for each delay time
	  //
		if(embedding_start != 1)
		    m = 1;
		else
		    m = embedding_start;
	    for(j=0; m<= embedding_finish; m += embedding_delta, j++){
			myBoxCounterp->use_time_series( data, data_count, m, i, min_value, max_value);
			tempf_p = myBoxCounterp->dimension(1);
			for( k=0; k< num_of_levels; k++ ){
                //epsilon = ( max_value - min_value)/pow(2.0,(double)k);
			    if( j != 0 ){
				    redundancy[j][k] =  m * redundancy[0][k] - tempf_p[k];
				} else {
				    redundancy[j][k] = tempf_p[k] ;
				}
			}
			if(embedding_start != 1)
			    if( !j ) m = embedding_start-embedding_delta;
		}
	  //
	  //Print Redundancy Data
	  //
		fprintf(stdout,"\ndelay = %d: ", i);
		for(m=embedding_start; m<= embedding_finish; m += embedding_delta){
			fprintf(stdout," m = %3d",m);
			if( m < embedding_finish ) fprintf(stdout,",\t",m);
		}
		fprintf(stdout,"\n");
		for( k=0; k< num_of_levels; k++ ){
		    fprintf(stdout,"%d:\t",k+1);
			if(embedding_start == 1)
				j=0;
			else 
				j=1;
			for(m=embedding_start; m<= embedding_finish; m += embedding_delta, j++){
				fprintf(stdout,"%f",redundancy[j][k]);
				if( m < embedding_finish ) fprintf(stdout,",\t",m);
			}
			fprintf(stdout,"\n");
		}
		fprintf(stdout,"\n");
		fflush(stdout);
	}
    
	for(i=0; i<red_maxdim; i++)
	    free( redundancy[i] );
    free( redundancy );
    free( data );
    delete myBoxCounterp;
}
