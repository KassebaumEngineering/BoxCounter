
//
// BoxCounter.cc
//
// C++ Implementation for the BoxCounter Class 
//
//  $Id: BoxCounter.cc,v 1.1.1.1 1997/08/23 17:50:33 jak Exp $
//
//  Author: John Kassebaum
//
// $Log: BoxCounter.cc,v $
// Revision 1.1.1.1  1997/08/23 17:50:33  jak
// BoxCounter disappeared from CVS ... this is a replacement. -jak
//
// Revision 1.5  1995/11/08 05:28:50  jak
// Added back lost RCS variables. -jak
//
// Revision 1.3  1994/05/03  17:11:30  jak
// Removed all tabs from the formatting of the file. -jak
//
// Revision 1.2  1994/05/03  17:07:39  jak
// Changed the way the commenting around the log lines looks. -jak
//
// Revision 1.1  1994/05/03  16:56:14  jak
// First checkin of working boxcounter software  -jak
//
//

static char rcsid_BoxCounter_cc[] = "$Id: BoxCounter.cc,v 1.1.1.1 1997/08/23 17:50:33 jak Exp $";

#include "BoxCounter.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define DEBUG

#define BLOCKSIZE    1024
#define Null(A)             ((A *)0)
#define New(A)                ((A *) malloc( sizeof(A) ) )
#define NewBlock(A,N)        ((A *) malloc( sizeof(A) * (N)) )
#define BiggerBlock(A,B,N)  ((A *) realloc( (void *)(B), sizeof(A) * (N)))
#define Free(A)             free( (void *)A )

#define MIN(A,B)            ( (A) < (B) ? (A) : (B) )
#define MAX(A,B)            ( (A) > (B) ? (A) : (B) )

#define bzero(A,B)    memset( (A), 0, (B) )
extern "C" {
    extern void memset(char *, char, int);
}

ArbitraryPrecisionDatum::ArbitraryPrecisionDatum( unsigned int val, unsigned int precision )
{
    float size;
    size = (float)(precision) / (float) (sizeof(unsigned int) * 8);
    if( !(datum = NewBlock( unsigned int,  (length = (int) ceil(size)) )) ){
        perror("ArbitraryPrecisionDatum:: ArbitraryPrecisionDatum() -> datum");
        abort();
    };
    memset( (char *)datum, 0, length * sizeof(unsigned int) );
    datum[0] = val;
};

ArbitraryPrecisionDatum::ArbitraryPrecisionDatum( unsigned long long val, unsigned int precision )
{
    float size;
    size = (float)(precision) / (float) (sizeof(unsigned int) * 8);
    if( !(datum = NewBlock( unsigned int,  (length = (int) ceil(size)) )) ){
        perror("ArbitraryPrecisionDatum:: ArbitraryPrecisionDatum() -> datum");
        abort();
    };
    memset( (char *)datum, 0, length * sizeof(unsigned int) );
    datum[0] = val & 0xffffffff;
    datum[1] = val >> 32;
};

ArbitraryPrecisionDatum::ArbitraryPrecisionDatum( ArbitraryPrecisionDatum& adatum )
{
    length = adatum.length;
    if( !(datum = NewBlock( unsigned int,  length )) ){
        perror("ArbitraryPrecisionDatum:: ArbitraryPrecisionDatum() -> datum");
        abort();
    };
    memcpy( (char *)datum, (char *)adatum.datum, sizeof(unsigned int)*length );
};

ArbitraryPrecisionDatum::~ArbitraryPrecisionDatum()
{
    if(datum) Free( datum );
};

ArbitraryPrecisionDatum ArbitraryPrecisionDatum::operator<<( int shift )
{
    ArbitraryPrecisionDatum rtn( *this );
    int i, bigshift;
    unsigned int carry, nextcarry;

//fprintf(stderr,"(<<: shift = %d ) ", shift);
    bigshift = 0;
    if( abs(shift) > sizeof(unsigned int)*8 ){
        bigshift = shift / sizeof(unsigned int)*8;
        shift = shift % sizeof(unsigned int)*8;
//fprintf(stderr,"(<<: bigshift = %d ) ", bigshift);
//fprintf(stderr,"(<<: shift = %d ) ", shift);
    }

    if (shift > 0){ // shift left
        if( bigshift ){
            if( !(rtn.datum = BiggerBlock( unsigned int, rtn.datum,(rtn.length+bigshift)) ) ){
                perror("ArbitraryPrecisionDatum:: operator <<() -> datum overflow");
                abort();
            };
            rtn.length += bigshift;
            for( i= rtn.length-1; i>= bigshift; i-- ){
                rtn.datum[i] = rtn.datum[i-bigshift];
                rtn.datum[i-bigshift] = 0;
            }
        }
        carry = 0;
        for(i=0; i< rtn.length; i++){
            nextcarry = rtn.datum[i] >> (31 - shift);
            rtn.datum[i] = (rtn.datum[i] << shift) + carry;
            carry = nextcarry;
        }
        if (carry > 0) {
//fprintf(stderr,"(<<: carry = %d (> 0)) ", carry);
            if( !(rtn.datum = BiggerBlock( unsigned int, rtn.datum,(rtn.length+1)) ) ){
                perror("ArbitraryPrecisionDatum:: operator <<() -> datum overflow");
                abort();
            };
            rtn.datum[rtn.length] = carry;
            rtn.length += 1;
        }
    } else if( shift < 0 ) { // shift right (shift is negative)
        if( bigshift ){  // bigshift should be negative
            for( i= rtn.length-1; i>= bigshift; i-- ){
                rtn.datum[i] = rtn.datum[i-bigshift];
                rtn.datum[i-bigshift] = 0;
            }
        }
        carry = 0;
        for(i=rtn.length-1; i>=0 ; i--){
            nextcarry = rtn.datum[i] << (31 + shift);
            rtn.datum[i] = (rtn.datum[i] >> -shift) + carry;
            carry = nextcarry;
        }
    }

    return rtn;
}

ArbitraryPrecisionDatum ArbitraryPrecisionDatum::operator &( ArbitraryPrecisionDatum& adatum)
{
    int i;

    if ( length <= adatum.length){
        ArbitraryPrecisionDatum rtn( *this );
        for(i=0; i< length; i++){
            rtn.datum[i] = rtn.datum[i] & adatum.datum[i];
        }
        return rtn;
    } else {
        ArbitraryPrecisionDatum rtn( adatum );
        for(i=0; i< adatum.length; i++){
            rtn.datum[i] = rtn.datum[i] & adatum.datum[i];
        }
        return rtn;
    }
};

ArbitraryPrecisionDatum ArbitraryPrecisionDatum::operator+( ArbitraryPrecisionDatum& adatum )
{
    ArbitraryPrecisionDatum rtn( *this );
    int i;
    unsigned long long result;
    unsigned int carry, nextcarry;

    if( rtn.length > adatum.length ){
        carry = 0;
        for(i=0; i< adatum.length; i++){
            result = rtn.datum[i] + adatum.datum[i] + carry;
            carry = (unsigned int) (result >> 32);
            rtn.datum[i] = result & 0xffffffff;
        }
        for( i=adatum.length; i< rtn.length; i++ ){
            result = rtn.datum[i] + carry;
            carry = (unsigned int) (result >> 32);
            rtn.datum[i] = result & 0xffffffff;
        }
    } else if ( rtn.length <= adatum.length ) {
        int oldlength;
        if( !(rtn.datum = BiggerBlock( unsigned int, rtn.datum, adatum.length ) ) ){
            perror("ArbitraryPrecisionDatum:: operator+() -> datum overflow");
            abort();
        };
        oldlength = rtn.length;
        rtn.length = adatum.length;
        memset( (char *) &rtn.datum[oldlength], 0, (adatum.length-oldlength) * sizeof(unsigned int) );
        carry = 0;
        for(i=0; i< oldlength; i++){
            result = rtn.datum[i] + adatum.datum[i] + carry;
            carry = (unsigned int) (result >> 32);
            rtn.datum[i] = result & 0xffffffff;
        }
        for( i=oldlength; i< adatum.length; i++ ){
            result = adatum.datum[i] + carry;
            carry = (unsigned int) (result >> 32);
            rtn.datum[i] = result & 0xffffffff;
        }
    } 

    if (carry > 0) {
        if( !(rtn.datum = BiggerBlock( unsigned int, rtn.datum, (rtn.length+1) ) ) ){
            perror("ArbitraryPrecisionDatum:: operator+() -> datum overflow");
            abort();
        };
        rtn.datum[ rtn.length ] = carry;
        rtn.length += 1;
    }

    return rtn;
}

ArbitraryPrecisionDatum ArbitraryPrecisionDatum::operator +( unsigned int anum )
{
    ArbitraryPrecisionDatum rtn( *this );
    int i;
    unsigned long long result;
    unsigned int carry, nextcarry;

    carry = 0;
    for(i=0; i< rtn.length; i++){
        if (i==0)
            result = rtn.datum[i] + anum + carry;
        else
            result = rtn.datum[i] + carry;
        carry = (unsigned int) (result >> 32);
        rtn.datum[i] = result & 0xffffffff;
    }

    if (carry > 0) {
        if( !(rtn.datum = BiggerBlock( unsigned int, rtn.datum,(rtn.length+1)) ) ){
            perror("ArbitraryPrecisionDatum:: operator+() -> datum overflow");
            abort();
        };
        rtn.datum[rtn.length] = carry;
        rtn.length += 1;
    }
                
    return rtn;
};

ArbitraryPrecisionDatum ArbitraryPrecisionDatum::operator +( unsigned long long anum )
{
    ArbitraryPrecisionDatum rtn( *this );
    int i;
    unsigned long long result;
    unsigned int carry, nextcarry;

    carry = 0;
    for(i=0; i< rtn.length; i++){
        if (i==0)
            result = rtn.datum[i] + anum + carry;
        else
            result = rtn.datum[i] + carry;
        carry = (unsigned int) (result >> 32);
        rtn.datum[i] = result & 0xffffffff;
    }

    if (carry > 0) {
        if( !(rtn.datum = BiggerBlock( unsigned int, rtn.datum,(rtn.length+1)) ) ){
            perror("ArbitraryPrecisionDatum:: operator+() -> datum overflow");
            abort();
        };
        rtn.datum[rtn.length] = carry;
        rtn.length += 1;
    }
                
    return rtn;
};

ArbitraryPrecisionDatum ArbitraryPrecisionDatum:: operator =( unsigned int anum )
{
    memset( (char *) datum, 0, sizeof( unsigned int ) * length );

    return (*this + anum);
};

ArbitraryPrecisionDatum ArbitraryPrecisionDatum:: operator =( unsigned long long anum )
{
    memset( (char *) datum, 0, sizeof( unsigned int ) * length );

    return (*this + anum);
};

ArbitraryPrecisionDatum ArbitraryPrecisionDatum::operator =( ArbitraryPrecisionDatum &adatum )
{
    if (length != adatum.length){
        length = adatum.length;
        if( !(datum = BiggerBlock(unsigned int,datum, length) ) ){
            perror("ArbitraryPrecisionDatum:: operator =() -> datum");
            abort();
        };
    }
    memcpy( (char *) datum, (char *) adatum.datum , sizeof(unsigned int)*length );

    return *this;
};

void ArbitraryPrecisionDatum::print( FILE *fp )
{
    int i;
    for( i= (signed int)length - 1; i>=0; i--){
        fprintf( fp, "%8.8x", datum[i]);
        if( i ) fprintf( fp,", ");
    }
};

char *ArbitraryPrecisionDatum::print( char *astring )
{
    char elem[20];
    int i,j;

    for(i=length-1,j=0; i>=0; i--){
        sprintf( elem, "%8.8x", datum[i]);
        if( i ) {
            sprintf( &astring[j],"%s, ", elem);
        } else {
            sprintf( &astring[j], "%s ", elem);
        }
        j += (strlen(elem) + 2);
   }

    return astring;
};

int ArbitraryPrecisionDatum::operator ==( ArbitraryPrecisionDatum &adatum )
{
    return (compare( (void *)this, (void *) &adatum ) == 0);
};

int compare(const void *a_in , const void *b_in ){
    int i,rtn;
    ArbitraryPrecisionDatum *a, *b;

    a = (ArbitraryPrecisionDatum *)a_in;
    b = (ArbitraryPrecisionDatum *)b_in;

    if( a->length > b->length ){
        for( i = a->length-1; i >= b->length; i--)
            if( rtn = a->datum[i]  )
                return rtn;
    } else if ( a->length < b->length ){
        for( i = b->length-1; i >= a->length; i--)
            if( rtn = b->datum[i]  )
                return rtn;
    }

    for( i =MIN(a->length, b->length)-1; i >= 0; i-- ){
        if( rtn = a->datum[i] - b->datum[i] ){ 
            break;
        }
    }

    return rtn;
};

BoxCounter:: BoxCounter( unsigned int quants ): 
    num_of_quant_levels(quants), pre_calc_is_done(0),
    bit_interleave_is_done(0), q_of_dimension( -1.0 )
{
    int j,k;

    if( !(bitmask = new ArbitraryPrecisionDatum[ num_of_quant_levels ]) ){
        perror("BoxCounter:: BoxCounter() -> bitmask");
        abort();
    };
    if( !(boxmask = new ArbitraryPrecisionDatum[ num_of_quant_levels ]) ){
        perror("BoxCounter:: BoxCounter() -> boxmask");
        abort();
    };
    if( !(gen_dimension = NewBlock( double, num_of_quant_levels )) ){
        perror("BoxCounter:: BoxCounter() -> gen_dimension");
        abort();
    };
    if( !(precalc_gen_dimension[0] = NewBlock( double, num_of_quant_levels )) ){
        perror("BoxCounter:: BoxCounter() -> precalc_gen_dimension[0]");
        abort();
    };
    if( !(precalc_gen_dimension[1] = NewBlock( double, num_of_quant_levels )) ){
        perror("BoxCounter:: BoxCounter() -> precalc_gen_dimension[1]");
        abort();
    };
    if( !(precalc_gen_dimension[2] = NewBlock( double, num_of_quant_levels )) ){
        perror("BoxCounter:: BoxCounter() -> precalc_gen_dimension[2]");
        abort();
    };

};

BoxCounter:: ~BoxCounter()
{
    if (bitmask) delete [] bitmask;
    if (boxmask) delete [] boxmask;
    if (bit_interleaved_series) delete [] bit_interleaved_series ;
    if (gen_dimension) Free( gen_dimension );
    if (precalc_gen_dimension[0]) Free( precalc_gen_dimension[0] );
    if (precalc_gen_dimension[1]) Free( precalc_gen_dimension[1] );
    if (precalc_gen_dimension[2]) Free( precalc_gen_dimension[2] );
};

void BoxCounter::use_time_series( const float *aseries, unsigned int N, unsigned int dim, unsigned int delay )
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

void BoxCounter::use_time_series( const float *aseries, unsigned int N, unsigned int dim, unsigned int adelay, float min, float max)
{
    int j,k;
    pre_calc_is_done = 0;
    bit_interleave_is_done = 0;

    original_time_series = aseries;
    embedding_dimension = dim;
    delay = adelay;
    num_of_data_points = N - embedding_dimension*delay;
    min_val = min;
    max_val = max;
    working_precision = num_of_quant_levels * embedding_dimension;
#ifdef LIMIT
    if( num_of_quant_levels * embedding_dimension > 64 ){
        fprintf(stderr,"BoxCounter::use_time_series() -> not enough bits for %d dimension\n", embedding_dimension);
        fprintf(stderr,"time series at %d quantization levels.\n", num_of_quant_levels);
        fprintf(stderr,"I suggest you use fewer quantization levels and try again.\n", num_of_quant_levels);
        fprintf(stderr,"max capabiltiy -> num_of_quant_levels * embedding_dimension < 64 bits\n", num_of_quant_levels);
        abort();
    }
#endif

#ifndef BOGUS
    if( !(bit_interleaved_series = new ArbitraryPrecisionDatum[ num_of_data_points ])){
#else
    if( !(bit_interleaved_series = new ArbitraryPrecisionDatum( 0, working_precision )[ num_of_data_points ] )){
#endif
        perror("BoxCounter:: use_time_series() -> bit_interleaved_series");
        abort();
    };

    for(k=0; k < num_of_quant_levels; k++){
        bitmask[k] = ArbitraryPrecisionDatum( 0x01 ) <<  k;
        boxmask[k] = ArbitraryPrecisionDatum( 0x00 );
        for(j=0; j < embedding_dimension*(k+1); j++){
            boxmask[k] = boxmask[k] + ( ArbitraryPrecisionDatum( 0x01 ) << j );
        }
#ifdef DEBUG
//        fprintf(stderr, "%d: %x - %x \n", k, bitmask[k], boxmask[k]);
        fprintf(stderr, "%d: ", k); bitmask[k].print(stderr); fprintf(stderr, " - "); boxmask[k].print(stderr); 
        fprintf(stderr, "\n");
#endif
    }
};


void BoxCounter::bit_interleave()
{
    int  i,k;
    long j, scaler;
    unsigned int *quantized_input;

    if( !(quantized_input = NewBlock( unsigned int, embedding_dimension )) ){
        perror("BoxCounter:: bit_interleave()");
        abort();
    };

    scaler = (0x01 << num_of_quant_levels) - 1 ;

    for(j=0; j< num_of_data_points - (embedding_dimension-1)*delay; j++){
        bit_interleaved_series[j] = 0;
        for(k=0 ; k< num_of_quant_levels; k++){
//fprintf(stderr,"working_precision = %d, shift left by %d\n", working_precision, (working_precision-1-k));
            for(i=0; i< embedding_dimension; i++){
                quantized_input[i] = (unsigned long long) rint( (double) scaler * ((double)original_time_series[j + i*delay] - min_val) / (max_val - min_val) );

#ifdef DEBUG
//        fprintf(stderr,"quantized_input #%d = %x\n", i, quantized_input[i]);
//        fprintf(stderr,"bit_interleaved_series[%d] >> 1 = ",j); 
//        (bit_interleaved_series[j] >> 1).print( stderr ); 
//        fprintf(stderr,"\n");
//        fprintf(stderr,"( "); bitmask[k].print(stderr); fprintf(stderr," & %x ) << %d = ",quantized_input[i],(working_precision-1-k) ); 
//        ( ArbitraryPrecisionDatum((bitmask[k] & quantized_input[i])) << (working_precision-1-k)).print( stderr ); 
//        fprintf(stderr,"\n");
#endif                

                bit_interleaved_series[j] = (bit_interleaved_series[j] >> 1) + ( ArbitraryPrecisionDatum((bitmask[k] & quantized_input[i])) << (working_precision-1-k));

#ifdef DEBUG
//        fprintf(stderr,"("); bitmask[k].print(stderr); fprintf(stderr," & %x) = ", quantized_input[i]); 
//        (bitmask[k] & quantized_input[i]).print( stderr ); fprintf(stderr,"\n");
//        fprintf(stderr,"interleaved #%d = ", j); 
//        bit_interleaved_series[j].print( stderr ); 
//        fprintf(stderr,"\n");
#endif
            }
        }
    }

    qsort((void *)bit_interleaved_series, (num_of_data_points - (embedding_dimension-1)*delay), sizeof(ArbitraryPrecisionDatum), compare);

#ifdef DEBUG
fprintf(stderr,"qsort complete!\n");
#endif
    bit_interleave_is_done = 1;
    Free( quantized_input );
};

void BoxCounter::pre_calc()
{
    int j,k;
    ArbitraryPrecisionDatum  currentBox, lastBox;
    unsigned int   ni;
    double         temp, temp1;
    double         sum[3];
    double         range;

    range = max_val - min_val; 

    if (!bit_interleave_is_done) 
        bit_interleave();

    for( k=0 ; k< num_of_quant_levels; k++){
        ni = 1;
        sum[0] = 0.0;
        sum[1] = 0.0;
        sum[2] = 0.0;
        lastBox = ( bit_interleaved_series[0] & boxmask[k]);
//fprintf(stderr,"\nbit_interleaved_series[0] = "); bit_interleaved_series[0].print( stderr ); fprintf(stderr,"\n");
//fprintf(stderr,"boxmask[%d] = ",k); boxmask[k].print( stderr ); fprintf(stderr,"\n");
//fprintf(stderr,"lastBox = "); lastBox.print( stderr ); fprintf(stderr,"\n");
        for( j=0; j < (num_of_data_points - (embedding_dimension-1)*delay); j++ ){
            currentBox = ( bit_interleaved_series[j] & boxmask[k]);
//fprintf(stderr,"lastBox = "); lastBox.print( stderr ); fprintf(stderr,"\n");
//fprintf(stderr,"currentBox = "); currentBox.print( stderr ); fprintf(stderr,"\n");
            if(lastBox == currentBox) {
                ni++; 
//fprintf(stderr,"+");
            } else {
//fprintf(stderr,".");
                lastBox = currentBox;
                temp = ( (double) ni / (double) (num_of_data_points - (embedding_dimension-1)*delay - 1) );
                if (temp != 0.0)  
                    temp1 =  temp  * log( temp ); 
                else 
                    temp1 = 0;
                sum[0] += 1;
                sum[1] += temp1;
                sum[2] += temp*temp;
                ni = 1;
            }
        }
        temp = ( (double) ni / (double) (num_of_data_points - (embedding_dimension-1)*delay - 1) );
        if (temp != 0.0)  
            temp1 =  temp * log( temp );
        else 
            temp1 = 0;
        sum[0] += 1;
        sum[1] += temp1;
        sum[2] += temp*temp;
        precalc_gen_dimension[0][k] = sum[0];
        precalc_gen_dimension[1][k] = sum[1];
        precalc_gen_dimension[2][k] = sum[2];
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
    fprintf(stderr,"%d: sum[0] = %lf, sum[1] = %lf, sum[2] = %lf\n",k, temp0, temp1, temp2);
}
#endif
    }

    pre_calc_is_done = 1;
};

double *BoxCounter::dimension( double a_q )
{
    int j,k;
    ArbitraryPrecisionDatum  currentBox, lastBox;
    unsigned int   ni;
    double         temp, temp1;
    double         sum;

    double        *rtn;

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
                        temp = ( (double) ni / (double) (num_of_data_points - (embedding_dimension-1)*delay - 1) );
                        sum += pow( temp, q_of_dimension);
                        ni = 1;
                    }
                }
                temp = ( (double) ni / (double) (num_of_data_points - (embedding_dimension-1)*delay - 1) );
                sum += pow( temp, q_of_dimension);
                gen_dimension[k] = sum;
            }
        }

        rtn =  gen_dimension;
    }

    return rtn;
};

double *BoxCounter::dimension( int a_q )
{
    return dimension( (double) a_q );
};

double BoxCounter::dimension( int a_q, int a_quantlevel)
{
    double *ptr;
    ptr = dimension( (double) a_q );

    return ptr[a_quantlevel];
};

