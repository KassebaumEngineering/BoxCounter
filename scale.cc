//
// dimensions.cc
//
// C++ Main Program for Generalized Dimensions Calculations 
//
//  $Id: scale.cc,v 1.1 1997/08/23 17:50:32 jak Exp $
//
//  Author: John Kassebaum
//
// $Log: scale.cc,v $
// Revision 1.1  1997/08/23 17:50:32  jak
// Initial revision
//
//
//

static char rcsid_dimensions_cc[] = "$Id: scale.cc,v 1.1 1997/08/23 17:50:32 jak Exp $";

#include "FastBoxCounter.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define USAGE    "Usage: %s [-h] [-v] [-s] < ascii_data_file \n"
#define DESCR_1    " -h  prints this help message and exits.\n"
#define DESCR_2    " -v  verbosely prints out data\n"
#define DESCR_3    " -s  prints statistics\n"
#define DESCR_4    " -e #  sets embedding dimension (default = 1)\n"
#define DESCR_5    " -d #  sets delay for embedding (default = 1)\n"
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
    double temp0f,  temp1f, temp2f, epsilon0, epsilon1, epsilon2, tempval;
    int step_count;
    float step;
    int num_of_levels;
    int embedding;
    int delay;
    float *data;
    FastBoxCounter *corrIntegral_1_p, *corrIntegral_2_p;

    SETmin_value = SETmax_value = 0.0;
    verbose = 0;
    stats = 0;
    num_of_levels = 8;
    embedding = 1;
    delay = 1;
    for (c=1; c< argc; c++) {
        if ( !strcmp( argv[ c ],"-h") || !strcmp( argv[ c ],"-help")){
            fprintf(stderr,USAGE,argv[0]);
            fprintf(stderr,DESCR_1);
            fprintf(stderr,DESCR_2);
            fprintf(stderr,DESCR_3);
            fprintf(stderr,DESCR_4);
            fprintf(stderr,DESCR_5);
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
        } else if (!strcmp( argv[ c ],"-d")){
            delay = atoi( argv[++c] );;
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
	    fprintf(stdout, "Scale calculations:\n");
	    fprintf(stdout, "Embedding delay: %d\n",delay);
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

        
    corrIntegral_1_p = new FastBoxCounter( num_of_levels );
    corrIntegral_2_p = new FastBoxCounter( num_of_levels );

    corrIntegral_1_p->use_time_series( data, data_count, embedding, delay, min_value, max_value);
    corrIntegral_1_p->pre_calc();

    tempval = (( (max_value - min_value) * sqrt(2.0) ) - (max_value - min_value)) / 2.0;
    corrIntegral_2_p->use_time_series( data, data_count, embedding, delay, min_value-tempval, max_value+tempval);
    corrIntegral_2_p->pre_calc();

	epsilon1 = ((0x01 << num_of_levels) - 1)/pow(2.0,(double)(0));
	epsilon2 = ((0x01 << num_of_levels) - 1)/pow(2.0,(double)(0)) / sqrt(2.0);

	temp1f = corrIntegral_2_p->dimension(2,0);   // correlation dimension
	temp2f = corrIntegral_1_p->dimension(2,0); 
//fprintf(stderr, "\ntemp1f = %lf, temp2f = %lf, \nepsilon1 = %lf, epsilon2 = %lf\n", temp1f, temp2f, epsilon1, epsilon2);
	fprintf( stdout, "%lf	%lf\n", log( (epsilon1+epsilon2)/2.0 ), 
		(log(temp1f)-log(temp2f))/(log(epsilon1)-log(epsilon2)) );

	for( k=1; k< num_of_levels-1; k++ ){

		epsilon0 = epsilon1;
		epsilon1 = epsilon2;
		epsilon2 = ((0x01 << num_of_levels) - 1)/pow(2.0,(double)(k));
        temp0f = temp1f;   
        temp1f = temp2f;   
        temp2f = corrIntegral_2_p->dimension(2,k);   
//fprintf(stderr, "\ntemp0f = %lf, temp1f = %lf, temp2f = %lf, \nepsilon0 = %lf, epsilon1 = %lf, epsilon2 = %lf\n", temp0f, temp1f, temp2f, epsilon0, epsilon1, epsilon2);
        fprintf( stdout, "%lf	%lf\n", log( epsilon1 ), 
            ((log(temp0f)-log(temp1f))/(log(epsilon0)-log(epsilon1)) + 
             (log(temp0f)-log(temp2f))/(log(epsilon0)-log(epsilon2))*2.0 +
             (log(temp1f)-log(temp2f))/(log(epsilon1)-log(epsilon2)) ) / 4.0 );

		epsilon0 = epsilon1;
		epsilon1 = epsilon2;
		epsilon2 = ((0x01 << num_of_levels) - 1)/pow(2.0,(double)(k)) / sqrt(2.0);
        temp0f = temp1f;   
        temp1f = temp2f;   
        temp2f = corrIntegral_1_p->dimension(2,k); 
//fprintf(stderr, "\ntemp0f = %lf, temp1f = %lf, temp2f = %lf, \nepsilon0 = %lf, epsilon1 = %lf, epsilon2 = %lf\n", temp0f, temp1f, temp2f, epsilon0, epsilon1, epsilon2);
        fprintf( stdout, "%lf	%lf\n", log( epsilon1 ), 
            ( (log(temp0f)-log(temp1f))/(log(epsilon0)-log(epsilon1)) + 
              (log(temp0f)-log(temp2f))/(log(epsilon0)-log(epsilon2))*2.0 +
              (log(temp1f)-log(temp2f))/(log(epsilon1)-log(epsilon2)) ) / 4.0 );

	}

	epsilon0 = epsilon1;
	epsilon1 = epsilon2;
    temp0f = temp1f;   
    temp1f = temp2f;   
//fprintf(stderr, "\ntemp0f = %lf, temp1f = %lf, \nepsilon0 = %lf, epsilon1 = %lf\n", temp0f, temp1f, epsilon0, epsilon1);
	fprintf( stdout, "%lf	%lf\n", log( (epsilon0+epsilon1)/2.0 ), 
			(log(temp0f)-log(temp1f))/(log(epsilon0)-log(epsilon1)) );

	fprintf(stdout,"\n");
	fflush(stdout);
    
    free( data );
    delete corrIntegral_1_p;
    delete corrIntegral_2_p;
}
