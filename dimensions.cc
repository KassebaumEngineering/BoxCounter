//
// dimensions.cc
//
// C++ Main Program for Generalized Dimensions Calculations 
//
//  $Id: dimensions.cc,v 1.1 1997/08/23 17:50:32 jak Exp $
//
//  Author: John Kassebaum
//
// $Log: dimensions.cc,v $
// Revision 1.1  1997/08/23 17:50:32  jak
// Initial revision
//
// Revision 1.9  1995/11/08 05:23:16  jak
// Fixed loss of RCS variables! -jak
//
// Revision 1.5  1994/05/05  04:54:40  jak
// Made fixes and additions.  -jak
//
// Revision 1.4  1994/05/03  17:30:41  jak
// Changed the stderrs to stdouts and added an fflush().   -jak
//
// Revision 1.3  1994/05/03  17:11:31  jak
// Removed all tabs from the formatting of the file. -jak
//
// Revision 1.2  1994/05/03  17:07:41  jak
// Changed the way the commenting around the log lines looks. -jak
//
// Revision 1.1  1994/05/03  16:56:20  jak
// First checkin of working boxcounter software  -jak
//
//

static char rcsid_dimensions_cc[] = "$Id: dimensions.cc,v 1.1 1997/08/23 17:50:32 jak Exp $";

#include "FastBoxCounter.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define USAGE    "Usage: %s [-h] [-v] [-s] < ascii_data_file \n"
#define DESCR_1    " -h  prints this help message and exits.\n"
#define DESCR_2    " -v  verbosely prints out data\n"
#define DESCR_3    " -s  prints statistics\n"
#define DESCR_4    " -e #  sets embedding dimension (default = 1)\n"
#define DESCR_5    " -do #  sets start delay in embedding search (default = 1)\n"
#define DESCR_6    " -df #  sets final delay in embedding search (default = 16)\n"
#define DESCR_7    " -dd #  sets delay delta in embedding search (default = 1)\n"
#define DESCR_8    " -q #  sets highest quantization level (default = 8)\n"
#define DESCR_8a   " -min #  sets the data minimum value arbitrarily\n"
#define DESCR_8b   " -max #  sets the data maximum value arbitrarily\n"
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
    int c,i,j,k;
    int verbose, stats, flag;
    long int data_count;
    long int space_alloc;
    float min_value, max_value;
    float SETmin_value, SETmax_value;
    double tempf, epsilon;
    int step_count;
    float step;
    int num_of_levels;
    int embedding;
    int delay_start, delay_finish, delay_delta;
    float *data;
    FastBoxCounter *myBoxCounterp;

    SETmin_value = SETmax_value = 0.0;
    verbose = 0;
    stats = 0;
    num_of_levels = 8;
    embedding = 1;
    delay_start = 1;
    delay_finish = 16;
    delay_delta = 1;
    for (c=1; c< argc; c++) {
        if ( !strcmp( argv[ c ],"-h") || !strcmp( argv[ c ],"-help")){
            fprintf(stderr,USAGE,argv[0]);
            fprintf(stderr,DESCR_1);
            fprintf(stderr,DESCR_2);
            fprintf(stderr,DESCR_3);
            fprintf(stderr,DESCR_4);
            fprintf(stderr,DESCR_5);
            fprintf(stderr,DESCR_6);
            fprintf(stderr,DESCR_7);
            fprintf(stderr,DESCR_8);
            fprintf(stderr,DESCR_8a);
            fprintf(stderr,DESCR_8b);
            fprintf(stderr,DESCR_9);
            fprintf(stderr,DESCR_10);
            fprintf(stderr,DESCR_11);
            exit(0);
        } else if (!strcmp( argv[ c ],"-v")){
            verbose = 1;
        } else if (!strcmp( argv[ c ],"-s")){
            stats = 1;
        } else if (!strcmp( argv[ c ],"-e")){
            embedding = atoi( argv[++c] );
        } else if (!strcmp( argv[ c ],"-q")){
            num_of_levels = atoi( argv[++c] );;
        } else if (!strcmp( argv[ c ],"-min")){
            SETmin_value = atof( argv[++c] );;
        } else if (!strcmp( argv[ c ],"-max")){
            SETmax_value = atof( argv[++c] );;
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
	    fprintf(stdout, "Generalized Dimensions calculations:\n");
	    fprintf(stdout, "Embedding delay search: from %d to %d by %d\n",delay_start, delay_finish, delay_delta);
	    fprintf(stdout, "Quantization Levels up to: %d bits\n", num_of_levels);
	    fprintf(stdout, "Embedding Dimension searched: %d \n", embedding);
        fprintf(stdout, "Space Allocated for %d floats.\n", space_alloc);
        fprintf(stdout, "Space Used for %d floats.\n", data_count);
        fprintf(stdout, "Minimum value = %f, using %f.\n", min_value, (SETmin_value == 0.0)?min_value:SETmin_value );
        fprintf(stdout, "Maximum value = %f, using %f.\n", max_value, (SETmax_value == 0.0)?max_value:SETmax_value );
        fprintf(stdout, "Scale Time Series by %lf\n", ((0x01 << num_of_levels) - 1)/(max_value-min_value) );
        fflush(stdout);
    }

    if( SETmin_value != 0.0 ) 
        min_value = SETmin_value;

    if( SETmax_value != 0.0 ) 
        max_value = SETmax_value;

        
    myBoxCounterp = new FastBoxCounter( num_of_levels );
    for(i=delay_start; i<= delay_finish; i += delay_delta){
        myBoxCounterp->use_time_series( data, data_count, embedding, i, min_value, max_value);
        fprintf(stdout,"\ndelay = %d\n\tepsilon  \tln(epsilon)\t Q = 0  \t Q == 1  \t Q == 2",i);
        for( k=0; k< num_of_levels; k++ ){
//            epsilon = (max_value - min_value)/pow(2.0,(double)(k));
            epsilon = ((0x01 << num_of_levels) - 1)/pow(2.0,(double)(k));
            fprintf(stdout,"\n");
//            fprintf(stdout,"%d:\t",k+1);
            fprintf(stdout,"%d:\t%lf\t%lf\t",k+1,epsilon,log(epsilon));
            for( j=0; j<3; j++ ) {
                if( j==1 )
                    tempf = myBoxCounterp->dimension(j,k);
                else
                    tempf = log( myBoxCounterp->dimension(j,k) );
//                tempf = tempf / (double) (k+1);
                tempf = tempf / log( epsilon );
//                tempf = tempf;
                if( j != 1)
                    tempf /=  (1.0 - (double)j);
                fprintf(stdout,"%lf\t", tempf);
            }
        }
        fprintf(stdout,"\n");
        fflush(stdout);
    }
    
    free( data );
    delete myBoxCounterp;
}
