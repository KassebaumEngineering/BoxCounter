//
// FastBoxCounter.cc
//
// C++ Implementation for the FastBoxCounter Class 
//
//  $Id: FastBoxCounter.cc,v 1.1 1997/08/23 17:50:31 jak Exp $
//
//  Author: John Kassebaum
//
// $Log: FastBoxCounter.cc,v $
// Revision 1.1  1997/08/23 17:50:31  jak
// Initial revision
//
// Revision 1.4  1995/11/08 05:28:53  jak
// Added back lost RCS variables. -jak
//
// Revision 1.1  1994/05/05  04:55:41  jak
// New files - added a much faster boxcounter, and the redundancy calculator. -jak
//
//
//

static char rcsid_FastBoxCounter_cc[] = "$Id: FastBoxCounter.cc,v 1.1 1997/08/23 17:50:31 jak Exp $";

#include "FastBoxCounter.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define DEBUG

#define BLOCKSIZE    1024
#define Null(A)             ((A *)0)
#define New(A)              ((A *) malloc( sizeof(A) ) )
#define NewBlock(A,N)       ((A *) malloc( sizeof(A) * (N)) )
#define BiggerBlock(A,B,N)  ((A *) realloc( (void *)(B), sizeof(A) * (N)))
#define Free(A)             free( (void *)A )

#define MIN(A,B)            ( (A) < (B) ? (A) : (B) )
#define MAX(A,B)            ( (A) > (B) ? (A) : (B) )

#define bzero(A,B)    memset( (A), 0, (B) )
extern "C" {
    extern void memset(char *, char, int);
}

int compare_fast(const void *a_in , const void *b_in ){
    return ((*(unsigned long long *)a_in) - (*(unsigned long long *)b_in));
};

FastBoxCounter:: FastBoxCounter( unsigned int quants ): 
    num_of_quant_levels(quants), pre_calc_is_done(0),
    bit_interleave_is_done(0), q_of_dimension( -1.0 )
{
    int j,k;

    if( !(bitmask = NewBlock(unsigned long long, num_of_quant_levels ))){
        perror("FastBoxCounter:: FastBoxCounter() -> bitmask");
        abort();
    };
    if( !(boxmask = NewBlock(unsigned long long, num_of_quant_levels ))){
        perror("FastBoxCounter:: FastBoxCounter() -> boxmask");
        abort();
    };
    if( !(gen_dimension = NewBlock( double, num_of_quant_levels )) ){
        perror("FastBoxCounter:: FastBoxCounter() -> gen_dimension");
        abort();
    };
    if( !(precalc_gen_dimension[0] = NewBlock( double, num_of_quant_levels )) ){
        perror("FastBoxCounter:: FastBoxCounter() -> precalc_gen_dimension[0]");
        abort();
    };
    if( !(precalc_gen_dimension[1] = NewBlock( double, num_of_quant_levels )) ){
        perror("FastBoxCounter:: FastBoxCounter() -> precalc_gen_dimension[1]");
        abort();
    };
    if( !(precalc_gen_dimension[2] = NewBlock( double, num_of_quant_levels )) ){
        perror("FastBoxCounter:: FastBoxCounter() -> precalc_gen_dimension[2]");
        abort();
    };

};

FastBoxCounter:: ~FastBoxCounter()
{
    if (bitmask) Free( bitmask );
    if (boxmask) Free( boxmask );
    if (bit_interleaved_series) Free( bit_interleaved_series );
    if (gen_dimension) Free( gen_dimension );
    if (precalc_gen_dimension[0]) Free( precalc_gen_dimension[0] );
    if (precalc_gen_dimension[1]) Free( precalc_gen_dimension[1] );
    if (precalc_gen_dimension[2]) Free( precalc_gen_dimension[2] );
};

void FastBoxCounter::use_time_series( const float *aseries, unsigned int N, unsigned int dim, unsigned int delay )
{
    int i;
    float min_value, max_value;

  // calculate min and max
    for( i=0; i< N; i++){
        if (!i) {
            min_value = max_value = aseries[ i ];
        } else {
            if( aseries[ i ] < min_value )
                min_value = aseries[ i ];
            if( aseries[ i ] > max_value )
                max_value = aseries[ i ];
        }
    }

  // call other use_time_series() method
    use_time_series(aseries,N,dim,delay,min_value, max_value);
};

void FastBoxCounter::use_time_series( const float *aseries, unsigned int N, unsigned int dim, unsigned int adelay, float min, float max)
{
    int j,k;
    pre_calc_is_done = 0;
    bit_interleave_is_done = 0;

    original_time_series = aseries;
    embedding_dimension = dim;
    delay = adelay;
    num_of_data_points = N - (embedding_dimension-1)*delay;
    min_val = min;
    max_val = max;
    working_precision = num_of_quant_levels * embedding_dimension;
//    working_precision = 12 * embedding_dimension;
	
    if( num_of_quant_levels * embedding_dimension > sizeof(unsigned long long)*8 ){
        fprintf(stderr,"FastBoxCounter::use_time_series() -> not enough bits for %d dimension\n", embedding_dimension);
        fprintf(stderr,"time series at %d quantization levels.\n", num_of_quant_levels);
        fprintf(stderr,"I suggest you use fewer quantization levels and try again.\n", num_of_quant_levels);
        fprintf(stderr,"max capability -> num_of_quant_levels(%d) * embedding_dimension(%d) (= %d) < %d bits\n", num_of_quant_levels, num_of_quant_levels,  num_of_quant_levels*num_of_quant_levels,  sizeof(unsigned long long)*8);
        abort();
    }

    if( !(bit_interleaved_series = NewBlock(unsigned long long, num_of_data_points ))){
        perror("FastBoxCounter:: use_time_series() -> bit_interleaved_series");
        abort();
    };

    for(k=0; k < num_of_quant_levels; k++){
		bitmask[num_of_quant_levels-k-1] = ( ((unsigned long long) 0x01 ) <<  k );
//		bitmask[k] = (unsigned long long) pow(2.0,(double)k);
#ifdef DEBUG
        fprintf(stderr, "bitmask[ %d ] = %lx\n",num_of_quant_levels-k-1,  (long)bitmask[num_of_quant_levels-k-1]);
#endif
    }
    for(k=0; k < num_of_quant_levels; k++){
        boxmask[k] = ((unsigned long long) 0x00 );
        for( j=embedding_dimension*(num_of_quant_levels-k-1); j < embedding_dimension*num_of_quant_levels; j++){
            boxmask[k] = boxmask[k] + ( ((unsigned long long) 0x01 ) << j );
        }
#ifdef DEBUG
        fprintf(stderr, "%d: %x - %x \n", k+1, (long) bitmask[k], (long)boxmask[k]);
#endif
    }
};


void FastBoxCounter::bit_interleave()
{
    int  i,k,N;
    long j, scaler;
    unsigned int *quantized_input;

    if( !(quantized_input = NewBlock( unsigned int, embedding_dimension )) ){
        perror("FastBoxCounter:: bit_interleave()");
        abort();
    };

    scaler = (0x01 << num_of_quant_levels) - 1 ;
//    scaler = 1 ;

    N = num_of_data_points ; 
    for(j=0; j< N; j++){
        bit_interleaved_series[j] = 0;
        for(k=0 ; k< num_of_quant_levels; k++){
//fprintf(stderr,"working_precision = %d, shift left by %d\n", working_precision, (working_precision-(k+1)));
            for(i=0; i< embedding_dimension; i++){
                quantized_input[i] = (unsigned long long) rint( 
                     (double) scaler * 
                    ((double) original_time_series[j + i*delay] - min_val) / 
                    (max_val - min_val) 
                );

                bit_interleaved_series[j] = (bit_interleaved_series[j] << 1) + 
                    ((bitmask[k] & (unsigned long long)quantized_input[i]) >> (num_of_quant_levels-k-1));

#ifdef DEBUG
        fprintf(stderr,"quantized_input # %d = %lx\n", i, (long) quantized_input[i]);
        fprintf(stderr,"bit_interleaved_series[%d] << 1 = %lx",j,(long) bit_interleaved_series[j]); 
        fprintf(stderr,"( %lx & %lx ) >> %d = %lx \n", (long) bitmask[k], (long) quantized_input[i], (num_of_quant_levels-k-1), (long) ((bitmask[k] & (unsigned long long)quantized_input[i]) >> (num_of_quant_levels-k-1))); 
#endif                

#ifdef DEBUG
//        fprintf(stderr,"( %x & %x ) = %x \n", bitmask[k], quantized_input[i], (bitmask[k] & quantized_input[i]));
//        fprintf(stderr,"interleaved #%d = %x", j, bit_interleaved_series[j]); 
#endif
            }
        }
    }

    qsort((void *)bit_interleaved_series, N, sizeof(unsigned long long), compare_fast);

#ifdef DEBUG
fprintf(stderr,"qsort complete!\n");
#endif
    bit_interleave_is_done = 1;
    Free( quantized_input );
};

void FastBoxCounter::pre_calc()
{
    int j,k;
    unsigned long long  currentBox, lastBox;
    unsigned int   ni, N;
    double         temp, temp1;
    double         sum[3];
    double         range;

    range = max_val - min_val; 

    if (!bit_interleave_is_done) 
        bit_interleave();

    N = num_of_data_points;
	
    for( k=0 ; k< num_of_quant_levels; k++){
        ni = 1;
        sum[0] = 0.0;
        sum[1] = 0.0;
        sum[2] = 0.0;
        lastBox = ( bit_interleaved_series[0] & boxmask[k] );
        for( j=0; j < N ; j++ ){
            currentBox = ( bit_interleaved_series[j] & boxmask[k]);
            if(lastBox == currentBox) {
                ni++; 
            } else {
                lastBox = currentBox;
                temp = ( (double) ni / (double) (N) ); // probability of state being in this volume element
                if ( temp > 0.0 )  
//                    temp1 = - temp *  log( temp ) / log( 2.0 ); 
                    temp1 = temp *  log( temp ); 
                else 
                    temp1 = 0;
                sum[0] += 1.0;
                sum[1] += temp1;
                sum[2] += temp*temp;
                ni = 1;
            }
        }
        temp = ( (double) ni / (double) (N) ); // probability of state being in this volume element
        if ( temp > 0.0 )  
//            temp1 = - temp * log( temp ) / log( 2.0 );
            temp1 = temp * log( temp );
        else 
            temp1 = 0;
        sum[0] += 1.0;
        sum[1] += temp1;
        sum[2] += temp*temp;
	// Return Generalized (Renyi) Information                           // Generalized Dimensions:
//        precalc_gen_dimension[0][k] =  ( log( sum[0] ) / log( 2.0 ) );  //  ( log( sum[0] ) / log( 2.0 ) ) / (k+1);
        precalc_gen_dimension[0][k] =  sum[0];
        precalc_gen_dimension[1][k] =  -sum[1];                          //  sum[1] / (k+1);
//        precalc_gen_dimension[2][k] =  -( log( sum[2] ) / log( 2.0 ) ); // -( log( sum[2] ) / log( 2.0 ) ) / (k+1);
        precalc_gen_dimension[2][k] =  sum[2];
#ifdef DEBUG
{
    double temp0, temp1, temp2;
    temp0 = - log( sum[0] );
//    temp0 = sum[0] ;
    temp0 /= (double) -(k+1);
    temp1 = sum[1];
//    temp1 = sum[1];
    temp1 /= (double) -(k+1);
    temp2 = log( sum[2] );
//    temp2 = sum[2] ;
    temp2 /= (double) -(k+1);
    fprintf(stderr,"%d: sum[0] = %lf, sum[1] = %lf, sum[2] = %lf\n",k, sum[0], -sum[1], sum[2]);
}
#endif
    }

    pre_calc_is_done = 1;
};

double *FastBoxCounter::dimension( double a_q )
{
    int j,k,N;
    unsigned long long  currentBox, lastBox;
    unsigned int   ni;
    double         temp, temp1;
    double         sum;

    double        *rtn;

    N = num_of_data_points; 
    if (!pre_calc_is_done)
        pre_calc();

    if (a_q == 0.0){
        rtn =  precalc_gen_dimension[0];

    } else if (a_q == 1.0){
        rtn =  precalc_gen_dimension[1];

    } else if (a_q == 2.0){
        rtn =  precalc_gen_dimension[2];

    } else {

        if (q_of_dimension != a_q){
            q_of_dimension = a_q;
            for( k=0 ; k< num_of_quant_levels; k++){
                ni = 1;
                sum = 0.0;
                lastBox = ( bit_interleaved_series[0] & boxmask[k] );
                for( j=0; j < num_of_data_points; j++ ){
                    currentBox = ( bit_interleaved_series[j] & boxmask[k]);
                    if(lastBox == currentBox) {
                        ni++;
                    } else {
                        lastBox = currentBox;
                        temp = ( (double) ni / (double) (N) );
                        sum += pow( temp, q_of_dimension);
                        ni = 1;
                    }
                }
                temp = ( (double) ni / (double) (N) );
                sum += pow( temp, q_of_dimension);
                gen_dimension[k] = ( log( sum ) / log( 2.0 ) ) / ( 1.0 - q_of_dimension );
				
            }
        }

        rtn =  gen_dimension;
    }

    return rtn;
};

double *FastBoxCounter::dimension( int a_q )
{
    return dimension( (double) a_q );
};

double FastBoxCounter::dimension( int a_q, int a_quantlevel)
{
    double *ptr;
    ptr = dimension( (double) a_q );

    return ptr[a_quantlevel];
};

