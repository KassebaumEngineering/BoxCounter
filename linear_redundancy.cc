//
// linear_redundancy.cc
//
// C++ Main Program for Linear Redundancy Calculations 
//
//  $Id: linear_redundancy.cc,v 1.1 1997/08/23 17:50:31 jak Exp $
//
//  Author: John Kassebaum
//
// $Log: linear_redundancy.cc,v $
// Revision 1.1  1997/08/23 17:50:31  jak
// Initial revision
//
// Revision 1.4  1995/11/08 05:28:56  jak
// Added back lost RCS variables. -jak
//
// Revision 1.1  1994/05/09  07:39:38  jak
// Added a new tool for linear_redundancy calculations! -jak
//
//
//

static char rcsid_linear_redundancy_cc[] = "$Id: linear_redundancy.cc,v 1.1 1997/08/23 17:50:31 jak Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define USAGE    "Usage: %s [-h] [-v] [-s] < ascii_data_file \n"
#define DESCR_1    " -h  prints this help message and exits.\n"
#define DESCR_2    " -v  verbosely prints out data\n"
#define DESCR_3    " -s  prints statistics\n"
#define DESCR_4_1  " -eo #  sets starting embedding dimension for search (default = 2)\n"
#define DESCR_4_2  " -ef #  sets ending embedding dimension for search (default = 5)\n"
#define DESCR_4_3  " -ed #  sets delta for embedding dimension search (default = 1)\n"
#define DESCR_5    " -do #  sets start delay in embedding search (default = 1)\n"
#define DESCR_6    " -df #  sets final delay in embedding search (default = 16)\n"
#define DESCR_7    " -dd #  sets delay delta in embedding search (default = 1)\n"
#define DESCR_8    " ascii_data_file  -  A file of floating point numbers\n"
#define DESCR_9    "     in ascii format, with one number per line. \n example:\n"
#define DESCR_10   " -47.52\n -39.59\n -27.72\n -15.84\n"

#define BLOCKSIZE    1024
#define Null(A)             ((A *) 0)
#define New(A)              ((A *) malloc( sizeof(A) ) )
#define NewBlock(A,N)       ((A *) malloc( sizeof(A) * (N)) )
#define BiggerBlock(A,B,N)  ((A *) realloc( (void *)(B), sizeof(A) * (N)))

/* modified from numerical recipes */
static void tred2(double *a, int n, double *d, double *e)
{
  int l,k,j,i;
  double hh,h,g,f;  

  for(i=n-1; i>0; i--) {
    l=i-1;
    h=0.0;
    if (l > 0) {
      double scale;

      scale=0.0;

      for (k=0;k<i;k++)	scale += fabs( *(a + (i * n) + k) );

      if ( scale == 0.0) 
	e[i]=*(a + (i * n) + l);
      else {

	for (k=0;k<i;k++) {
	  double val;
	  val = *(a + (i * n) + k) / scale;
	  *(a + (i * n) + k) = val;
	  h += val*val;
	} /* for k */ 

	f=*(a + (i * n) + l);
	g = (f>0.0) ? -sqrt(h) : sqrt(h); 
	e[i]=scale*g;
	h -= f*g;
	*(a + (i * n) + l)=f-g;
	f=0.0;

	for (j=0;j<i;j++) {
	  /* Next statement can be omitted if eigenvectors not wanted */
//	  *(a + (j * n) + i)=*(a + (i * n) + j)/h;
	  g=0.0;
	  for (k=0;k<=j;k++) g += *(a + (j * n) + k) * *(a + (i * n) + k);
	  for (k=j+1;k<i;k++) g += *(a + (k * n) + j) * *(a + (i * n) + k);
	  e[j]=g/h;
	  f += e[j] * *(a + (i * n) + j);
	} /* for j */

	hh=f/(h+h);
	for (j=0;j<i;j++) {
	  f= *(a + (i * n) + j);
	  e[j]=g=e[j] - hh * f;
	  for (k=0;k<=j;k++) *(a + (j * n) + k) -= (f * e[k] + g * *(a + (i * n) + k) );
	} /* for j */
      }
    } else
      e[i]=*(a + (i * n) + l);
    
    d[i]=h;

  } /* for i */
  
  /* Next statement can be omitted if eigenvectors not wanted */
//  d[0]=0.0;
  e[0]=0.0;
  
  /* Contents of this loop can be omitted if eigenvectors not
     wanted except for statement d[i]=a[i][i]; */
  for (i=0;i<n;i++) {
//
//    if (d[i]) {
//      for (j=0;j<i;j++) {
//	g=0.0;
//	for (k=0;k<i;k++)  g += *(a + (i * n) + k) * *(a + (k * n) + j);
//	for (k=0;k<i;k++)  *(a + (k * n) + j) -= g * *(a + (k * n) + i);
//      } /* for j */
//    }
//
    d[i]=*(a + (i * n) + i);
//    *(a + (i * n) + i)=1.0;
//    for (j=0;j<i;j++) *(a + (j * n) + i) = *(a + (i * n) + j) = 0.0;
  } /* for i */
}

/* modified from numerical recipes */
static void tqli(double *d, double *e, int n, double *z)
{
  int m,l,iter,i; //k;
  double s,r,p,g,f,dd,c,b;

#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
  
  for (i=1;i<n;i++) e[i-1]=e[i];
  e[n-1]=0.0;
  
  l=0;

  for (l=0;l<n;l++) {
    iter=0;
    do {
      for (m=l;m<n-1;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	if (fabs(e[m])+dd == dd) break;
      }
      if (m != l) {
	if (iter++ > 30) {
	  fprintf(stderr, "\nEigenvector calculation does not seem to be converging...\n");
	}
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=sqrt((g*g)+1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  if (fabs(f) >= fabs(g)) {
	    c=g/f;
	    r=sqrt((c*c)+1.0);
	    e[i+1]=f*r;
	    c *= (s=1.0/r);
	  } else {
	    s=f/g;
	    r=sqrt((s*s)+1.0);
	    e[i+1]=g*r;
	    s *= (c=1.0/r);
	  }
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  p=s*r;
	  d[i+1]=g+p;
	  g=c*r-b;
	  /* Next loop can be omitted if eigenvectors not wanted */
//	  for (k=0;k<n;k++) {
//	    f = *(z + (k * n) + i + 1);
//	    *(z + (k * n) + i + 1) = s * *(z + (k * n) + i) + c * f;
//	    *(z + (k * n) + i) = c * *(z + (k * n) + i) - s * f;
//	  }
	}
	d[l]=d[l]-p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m != l);
  }
}


/* modified from numerical recipes */
static void eigsrt(double *d, double *v, int n)
{
  int k,j,i;
  double  p;
  
  for (i=0;i<n-1;i++) {
    p = d[k=i];
    for (j=i+1;j<n;j++)
      if (d[j] >= p) p = d[k=j];
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j=0;j<n;j++) {
	p = *(v + (j * n) + i);
	*(v + (j * n) + i) = *(v + (j * n) + k);
	*(v + (j * n) + k) = p;
      }
    }
  }
}

void calculate_eigen_vectors_symmetric(double *eigen_vector_matrix, double *eigen_value_vector, int n)
{
  double *offdiag;

  if( (offdiag = NewBlock( double, n )) == (double *)0 ){
    perror("calculate_eigen_vectors_symmetric()");
    abort();
  };

  /* Call two functions to tridiagonalize the eigen_vector_matrix matrix and then
     return the eigenvalues (in eigen_value_vector) and the eigenvectors (in eigen_vector_matrix).
     The k-th column of eigen_vector_matrix will be the normalized eigenvector corresponding
     to eigen_value_vector[k] */
  tred2( eigen_vector_matrix, n, eigen_value_vector, offdiag );
  tqli( eigen_value_vector, offdiag, n, eigen_vector_matrix );
  
  /* Sort the eigenvalues into descending order and rearrange the columns of 
     eigen_vector_matrix accordingly. */
// don't need to waste time with this - I only need the eigenvaues - order is unimportant
//  eigsrt( eigen_value_vector, eigen_vector_matrix, n );

  free( (char *)offdiag );
  
}

void matrixcopy(double **dest, double **source, int rows, int cols)
{
    register int r,c;
    register double *dptr, *sptr;

//
// Assumes an array of pointers to columns
//
    dptr = dest[0];
    for(r=0; r<rows; r++){
        dest[r] = dptr;
        sptr = source[r];
        for(c=0; c<cols; c++)
            *dptr++ = *sptr++;
    }
}

main(int argc,char **argv)
{
    int c,i,j,k,m;
    int verbose, stats, flag;
    long int data_count;
    long int space_alloc;
    double min_value, max_value, mean_value;
    double *tempf_p;
    int step_count;
    double step;
    int embedding_start, embedding_finish, embedding_delta;
    int delay_start, delay_finish, delay_delta;
    double *data;

    double **covariance_mat;
    double **eigenvector_mat;
    double *eigenvectors; // traditional C - space reserve for double[][]
    double *eigenvalue_vec;

    verbose = 0;
    stats = 0;
    embedding_start = 2;    embedding_finish = 6;   embedding_delta = 1;
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
            exit(0);
        } else if (!strcmp( argv[ c ],"-v")){
            verbose = 1;
        } else if (!strcmp( argv[ c ],"-s")){
            stats = 1;
        } else if (!strcmp( argv[ c ],"-eo")){
            embedding_start = atoi( argv[++c] );
        } else if (!strcmp( argv[ c ],"-ef")){
            embedding_finish = atoi( argv[++c] ) + 1;
        } else if (!strcmp( argv[ c ],"-ed")){
            embedding_delta = atoi( argv[++c] );
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
    data = (double *)0;
    
    if( data = NewBlock( double, BLOCKSIZE ) ){
        space_alloc = BLOCKSIZE;
    } else {
        perror(argv[0]);
        abort();
    };
        
    flag = 0;
    mean_value = 0.0;
    while ( scanf( "%lf", &data[ data_count ]) != EOF ){
        mean_value += data[ data_count ];
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
            if( data = BiggerBlock( double, data, (space_alloc + BLOCKSIZE) ) ){
                space_alloc += BLOCKSIZE;
            } else {
                perror(argv[0]);
                abort();
            };
        };
    }
    mean_value = mean_value / data_count;

    if (stats) {
	    fprintf(stdout, "Linear Redundancy calculations:\n");
	    fprintf(stdout, "Embedding delay search: from %d to %d by %d\n",delay_start, delay_finish, delay_delta);
	    fprintf(stdout, "Embedding Dimensions searched: %d to %d by %d \n", embedding_start, embedding_finish-1, embedding_delta);
        fprintf(stdout, "Space Allocated for %d doubles.\n", space_alloc);
        fprintf(stdout, "Space Used for %d doubles.\n", data_count);
        fprintf(stdout, "Minimum value = %f.\n", min_value);
        fprintf(stdout, "Maximum value = %f.\n", max_value);
        fprintf(stdout, "Mean value = %f.\n", mean_value);
        fflush(stdout);
    }

  //
  // get space and initialize double **covariance_mat; 
  //
	if( (covariance_mat = NewBlock( double *, embedding_finish )) != 0 ){
        for(i=0; i<embedding_finish; i++){
            if( (covariance_mat[i] = NewBlock( double, embedding_finish )) == 0 ){
                char msg[80];
                sprintf( msg,"Covariance Matrix[%d]",i);
				perror(msg);
				abort();
            } else {
                for(j=0; j<embedding_finish; j++) 
                    covariance_mat[i][j] = 0.0;
            }
        }
	} else {
        perror("Covariance Matrix");
        abort();
    }

  //
  // get space and initialize double **eigenvector_mat;
  //
	if( (eigenvector_mat = NewBlock( double *, embedding_finish )) == 0 ){
        perror("Eigen Vector Matrix");
        abort();
    }
	if( (eigenvectors = (double *)new double[embedding_finish][embedding_finish]) != 0 ){
        for(i=0; i<embedding_finish; i++){
            eigenvector_mat[i] = &(eigenvectors[ i * embedding_finish ]);
			for(j=0; j<embedding_finish; j++) 
			    eigenvector_mat[i][j] = 0.0;
         }
	} else {
        perror("Eigen Vector Space");
        abort();
    }

  //
  // get space and initialize double *eigenvalue_vec;
  //
	if( (eigenvalue_vec = NewBlock( double, embedding_finish )) == 0 ){
		char msg[80];
		sprintf( msg,"Eigen Value Vector");
		perror(msg);
		abort();
	} else {
		for(j=0; j<embedding_finish; j++) 
			eigenvalue_vec[j] = 0.0;
	}
  
  //
  // Perform Calculations for Linear Redundancy
  //
	fprintf(stdout, "delay\t");
	for(k=embedding_start; k< embedding_finish; k += embedding_delta)
		fprintf(stdout, "m = %d    \t", k);
	fprintf(stdout, "\n");
    fflush(stdout);

    for(i=delay_start; i<delay_finish; i += delay_delta){
        for(k=0; k<embedding_finish; k++)
			for(j=0; j<embedding_finish; j++) 
				covariance_mat[k][j] = 0.0;
        for(j=0; j<(data_count - (embedding_finish-1)*i ); j++){
            int r,c;
            for(r=0, k=j; k< j + (embedding_finish-1)*i; k += i, r++ ){
                for(c=r, m=k; m < j + (embedding_finish-1)*i; m += i, c++){
                    register double tempval;
                    tempval = (data[k] - mean_value) * (data[m] - mean_value);
                    covariance_mat[r][c] += tempval;
                    covariance_mat[c][r] += tempval;
                }
            }
        }
        for(k=0; k<embedding_finish; k++)
			for(j=0; j<embedding_finish; j++) 
				covariance_mat[k][j] /= (data_count - (embedding_finish-1)*i - 1);
      //
      // Now I can calculate for all embeddings of interest for this delay
      //
        fprintf(stdout, "%d\t", i);
        for(k=embedding_start; k< embedding_finish; k += embedding_delta){
            double LinearRedundancyValue;
          //
          // Calculate Eigenvalues of the Covariance Matrix
          //
            matrixcopy(eigenvector_mat, covariance_mat, k, k);
            calculate_eigen_vectors_symmetric(eigenvectors, eigenvalue_vec, k);
          //
          // Ln = -0.5 * sum [ (m from 0 to k ) of Log( Eigenvalue m )]
          //
            for( LinearRedundancyValue = 0.0, m=0; m < k; m++){
                LinearRedundancyValue += (log(covariance_mat[m][m]) - log( eigenvalue_vec[m] ));
            }
            fprintf(stdout, "%f\t", 0.5 * LinearRedundancyValue);
        }
        fprintf(stdout, "\n");
        fflush(stdout);
    }
	
  //
  // Free Dynamic Data Space
  //
	for(j=0; j<embedding_finish; j++) 
       free( covariance_mat[j] );
    free( covariance_mat );
    free( eigenvectors );
    free( eigenvector_mat );
    free( eigenvalue_vec );
    free( data );
}
