//
// FastBoxCounter.h
//
// C++ Interface for the FastBoxCounter Class 
//
//  $Id: FastBoxCounter.h,v 1.1 1997/08/23 17:50:32 jak Exp $
//
//  Author: John Kassebaum
//
/* $Log: FastBoxCounter.h,v $
/* Revision 1.1  1997/08/23 17:50:32  jak
/* Initial revision
/*
 * Revision 1.5  1995/11/08 05:32:31  jak
 * Hmmm .. not a mistake after all ... CVS thinks this is straight C - oh well. -jak
 *
 *  Revision 1.4  1995/11/08  05:30:33  jak
 *  oops - a small mistake in the header comment - now fixed. -jak
 *
 *  Revision 1.3  1995/11/08  05:28:54  jak
 *  Added back lost RCS variables. -jak
 *  
 *  Revision 1.1  1994/05/05  04:55:42  jak
 *  New files - added a much faster boxcounter, and the redundancy calculator. -jak
 */
//

static char rcsid_FastBoxCounter_h[] = "$Id: FastBoxCounter.h,v 1.1 1997/08/23 17:50:32 jak Exp $";

#ifndef _FastBoxCounter_h
#define _FastBoxCounter_h

#include <stdio.h>

extern int compare(const void *, const void *);

//
// BoxCounter Class for calculating generalized dimensions of
// an E-dimensional time delay embedding or a scalar time series
//

class FastBoxCounter {
private:

// inputs
    unsigned int  embedding_dimension;
    unsigned int  num_of_data_points;
    unsigned int  delay;
    unsigned int  num_of_quant_levels;
    const float  *original_time_series; // array length = num_of_data_points

// temporaries
    unsigned long long *bitmask; // array length = num_of_quant_levels
    unsigned long long *boxmask; // array length = num_of_quant_levels
    unsigned long long *bit_interleaved_series; 
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
  // FastBoxCounter( quantization );
  //
    FastBoxCounter( unsigned int = 8 );
   ~FastBoxCounter( void );

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
