//
// BoxCounter.h
//
// C++ Interface for the BoxCounter Class 
//
//  $Id: BoxCounter.h,v 1.1 1997/08/23 17:50:31 jak Exp $
//
//  Author: John Kassebaum
//
// $Log: BoxCounter.h,v $
// Revision 1.1  1997/08/23 17:50:31  jak
// Initial revision
//
// Revision 1.3  1995/11/08 05:28:51  jak
// Added back lost RCS variables. -jak
//
// Revision 1.1  1994/05/03  16:56:16  jak
// First checkin of working boxcounter software  -jak
// 
//

static char rcsid_BoxCounter_h[] = "$Id: BoxCounter.h,v 1.1 1997/08/23 17:50:31 jak Exp $";

#ifndef _BoxCounter_h
#define _BoxCounter_h

#include <stdio.h>

//
// InterleavedDatum Class - tool for BoxCounter
//
class ArbitraryPrecisionDatum {
private:
public:
    unsigned int  length;
    unsigned int *datum;
public:
    // ArbitraryPrecisionDatum( Precision_in_bits );
    ArbitraryPrecisionDatum( unsigned int = 0, unsigned int = 32 );
    ArbitraryPrecisionDatum( unsigned long long, unsigned int = 64 );
    ArbitraryPrecisionDatum( ArbitraryPrecisionDatum & );
    ~ArbitraryPrecisionDatum( void );

    ArbitraryPrecisionDatum operator <<( int );
    inline ArbitraryPrecisionDatum operator >>( int ashift ){
        return (*this << -ashift );
    };
    ArbitraryPrecisionDatum operator &( ArbitraryPrecisionDatum& );
    inline ArbitraryPrecisionDatum operator &( unsigned int anint ){
        return (*this & ArbitraryPrecisionDatum( anint ) );
    };
    inline ArbitraryPrecisionDatum operator &( unsigned long long along ){
        return (*this & ArbitraryPrecisionDatum( along ) );
    };
    ArbitraryPrecisionDatum operator +( ArbitraryPrecisionDatum& );
    ArbitraryPrecisionDatum operator +( unsigned int );
    ArbitraryPrecisionDatum operator +( unsigned long long );
    ArbitraryPrecisionDatum operator =( unsigned int );
    ArbitraryPrecisionDatum operator =( unsigned long long );
    ArbitraryPrecisionDatum operator =( ArbitraryPrecisionDatum& );
    int operator ==( ArbitraryPrecisionDatum& );

    void print( FILE * = stdout );
    char *print( char * );

    inline int precision(){
        return (length * sizeof( unsigned int ) * 8);
    };

    friend int compare(const void *, const void *);
};

//
// BoxCounter Class for calculating generalized dimensions of
// an E-dimensional time delay embedding or a scalar time series
//

class BoxCounter {
private:

// inputs
    unsigned int  embedding_dimension;
    unsigned int  num_of_data_points;
    unsigned int  delay;
    unsigned int  num_of_quant_levels;
    const float  *original_time_series; // array length = num_of_data_points

// temporaries
    ArbitraryPrecisionDatum *bitmask; // array length = num_of_quant_levels
    ArbitraryPrecisionDatum *boxmask; // array length = num_of_quant_levels
    ArbitraryPrecisionDatum *bit_interleaved_series; 
                            // array length = num_of_data_points
    int working_precision;

    int           pre_calc_is_done;
    int           bit_interleave_is_done;
    double         min_val, max_val;

// outputs
    double       *gen_dimension; // array length = num_of_quant_levels
    double        q_of_dimension;
    double       *precalc_gen_dimension[3]; // array length = num_of_quant_levels
                                            // for 0, 1, and 2 dimensions
    
public:
  //
  // Constructor for a general BoxCounter:
  // BoxCounter( quantization );
  //
    BoxCounter( unsigned int = 16 );
   ~BoxCounter( void );

  // Queries
    inline unsigned int  get_embedding_dimension( void ) const { return embedding_dimension; };
    inline unsigned int  get_num_of_data_points( void ) const { return num_of_data_points; };
    inline unsigned int  get_delay( void ) const { return delay; };
    inline unsigned int  get_num_of_quant_levels( void ) const { return num_of_quant_levels; };
    inline int           is_bit_interleave_done( void ) const { return bit_interleave_is_done; };
    inline int           is_pre_calc_done( void ) const { return pre_calc_is_done; };

  //
  // use_time_series( pointer to a time series, number of data points, embedding dimension, delay );
  //
    void use_time_series( const float *, unsigned int , unsigned int, unsigned int);

  //
  // same as above except with pre-calculated min and max values
  //
    void use_time_series( const float *, unsigned int , unsigned int, unsigned int, float , float );

  //
  // construct sorted bit-interleaved series
  //
    void bit_interleave( void );

  //
  // construct Generalized dimensions
  //
    void pre_calc( void );

  //
  // calculate q-dimension
  //
    double *dimension( double );
    double *dimension( int );
    double dimension( int, int );

};


#endif
