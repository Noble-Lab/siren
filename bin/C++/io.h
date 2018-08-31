/**
 *  * \file regression.h
 *   * AUTHOR: Alex Hu
 *    * CREATE DATE: February 12, 2016
 *     * PROJECT: DIA Regression
 *      * \brief Functions to read in and process data
 *       * matrices
 *        *****************************************************************************/

#ifndef IO_H
#define IO_H

#include "dbl_matrix_rowmajor.h"
#include "getRSS.cpp"
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "dbl_matrix.h"

using namespace std;

vector<dbl_matrix*> get_y_rowmajorx( const char * xfile, const char * yfile, bool SQRT = true, bool rescale = true ){
        
	printf("Reading X: %s\n", xfile);
        clock_t begin = clock();	
	dbl_matrix_rowmajor *x = new dbl_matrix_rowmajor( xfile );
        printf("%ix%i\n", x->nrow, x->ncol );
        printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000);


	printf("Reading Y: %s\n", yfile);
	dbl_matrix *y = new dbl_matrix( yfile );
        printf("%ix%i\n", y->nrow, y->ncol );
        printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000);

	// Change the dimensions so that y-nrow is the same as x-> nrow	
	if( y-> nrow == x-> nrow ){
          printf("X and Y dimensions match\n");
	}
	else if( y-> nrow == x-> ncol ){
          printf("X and Y dimensions don't match, transposing X\n");
	  x->t(); 
	}
	else if( y->ncol == x->nrow ){
          printf("X and Y dimensions don't match, transposing Y\n");
	  y->t();
         printf("New Y: %ix%i\n", y->nrow, y->ncol );
	}
	else if( y->ncol == x->ncol ){
          printf("X and Y dimensions match but in the wrong way, transposing X and Y\n");
          x->t();
	  y->t();
	}
	else{
	  printf("X and Y's dimensions don't match and nothing can be done about it :(\n");
	}

	if( SQRT ){
		x->colnorm();
		printf("Square rooting x and y\n");
		x->sqrt();
		y->sqrt();
	}
	else{ printf("\n\n\nNot square rooting.\n"); }

	x->max_norm();
	y->max_norm();

	if ( rescale ){
		printf("Normalizing columns of x.\n");
		x-> col_magnorm();
	}
	else{ printf("\n\n\nNot normalizing.\n"); }



	printf("Creating vector\n");
	vector<dbl_matrix*> * xy = new vector<dbl_matrix*>();
	xy->push_back(x); xy->push_back(y);
	printf("Created:!\n");	

        printf("X: %ix%i\n", x->nrow, x->ncol );
        printf("Y: %ix%i\n\n", y->nrow, y->ncol );

	return(*xy);
}

vector<dbl_matrix*> get_xy_rowmajor( const char * xfile, const char * yfile, bool SQRT = true, bool rescale = true ){
        
	printf("Reading X: %s\n", xfile);
        clock_t begin = clock();	
	dbl_matrix_rowmajor *x = new dbl_matrix_rowmajor( xfile );
        printf("%ix%i\n", x->nrow, x->ncol );
        printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000);
	

	printf("Reading Y: %s\n", yfile);
	dbl_matrix_rowmajor *y = new dbl_matrix_rowmajor( yfile );
        printf("%ix%i\n", y->nrow, y->ncol );
        printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000);

	// Change the dimensions so that y-nrow is the same as x-> nrow	
	if( y-> nrow == x-> nrow ){
          printf("X and Y dimensions match\n");
	}
	else if( y-> nrow == x-> ncol ){
          printf("X and Y dimensions don't match, transposing X\n");
	  x->t(); 
	}
	else if( y->ncol == x->nrow ){
          printf("X and Y dimensions don't match, transposing Y\n");
	  y->t();
         printf("New Y: %ix%i\n", y->nrow, y->ncol );
	}
	else if( y->ncol == x->ncol ){
          printf("X and Y dimensions match but in the wrong way, transposing X and Y\n");
          x->t();
	  y->t();
	}
	else{
	  printf("X and Y's dimensions don't match and nothing can be done about it :(\n");
	}

	if( SQRT ){
		x->colnorm();
		printf("Square rooting x and y\n");
		x->sqrt();
		y->sqrt();
	}
	else{ printf("\n\n\nNot square rooting.\n"); }

	x->max_norm();
	y->max_norm();

	if ( rescale ){
		printf("Normalizing columns of x.\n");
		x-> col_magnorm();
	}
	else{ printf("\n\n\nNot normalizing.\n"); }



	printf("Creating vector\n");
	vector<dbl_matrix*> * xy = new vector<dbl_matrix*>();
	xy->push_back(x); xy->push_back(y);
	printf("Created:!\n");	

        printf("X: %ix%i\n", x->nrow, x->ncol );
        printf("Y: %ix%i\n\n", y->nrow, y->ncol );

	return(*xy);
}

vector<dbl_matrix*> get_xy( const char * xfile, const char * yfile, bool SQRT = true, bool rescale = true ){
        
	printf("\nReading X: %s\n", xfile);
        clock_t begin = clock();	
	dbl_matrix *x = new dbl_matrix( xfile );
        printf("%ix%i\n", x->nrow, x->ncol );
        printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000);

	printf("\nReading Y: %s\n", yfile);
	dbl_matrix *y = new dbl_matrix( yfile );
        printf("%ix%i\n", y->nrow, y->ncol );
        printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000);

	// Change the dimensions so that y-nrow is the same as x-> nrow	
	if( y-> nrow == x-> nrow ){
          printf("X and Y dimensions match\n");
	}
	else if( y-> nrow == x-> ncol ){
          printf("X and Y dimensions don't match, transposing X\n");
	  x->t(); 
	}
	else if( y->ncol == x->nrow ){
          printf("X and Y dimensions don't match, transposing Y\n");
	  y->t();
         printf("New Y: %ix%i\n", y->nrow, y->ncol );
	}
	else if( y->ncol == x->ncol ){
          printf("X and Y dimensions match but in the wrong way, transposing X and Y\n");
          x->t();
	  y->t();
	}
	else{
	  printf("X and Y's dimensions don't match and nothing can be done about it :(\n");
	}

	if( SQRT ){
		colnorm( x->m, x->nrow, x->ncol );
		printf("Square rooting x and y\n");
		x->sqrt();
		y->sqrt();
	}
	else{ printf("\n\n\nNot square rooting.\n"); }

	if( true ){
		x->max_norm();
		y->max_norm();
	}
	else{
		printf("THIS IS ABNORMAL!!! NOT MAX NORMING!\n");
	}
	
	if ( rescale ){
		printf("Normalizing columns of x.\n");
		col_magnorm( x->m, x->nrow, x->ncol );
	}
	else{ printf("\n\n\nNot normalizing.\n"); }



	printf("Creating vector\n");
	vector<dbl_matrix*> * xy = new vector<dbl_matrix*>();
	xy->push_back(x); xy->push_back(y);
	printf("Created:!\n");	

        printf("X: %ix%i\n", x->nrow, x->ncol );
        printf("Y: %ix%i\n\n", y->nrow, y->ncol );

	return(*xy);
}

vector<dbl_matrix*> get_xy_withms1( const char * xfile1,const char * xfile2, const char * yfile1, const char * yfile2, double ms1weight = 1.0, bool SQRT = true, bool rescale = true ){

	printf("\n\nSQRT is %i\n\n\n",SQRT);

	// Read them in        
	vector<dbl_matrix*> xy1 = get_xy( xfile1, yfile1, SQRT, rescale );
	vector<dbl_matrix*> xy2 = get_xy( xfile2, yfile2, SQRT, rescale );

	// Normalize the columns of y1 to the columns of y2
	// with an optional scaling factor of ms1weight
	// xy[1] is y, xy[0] is x
	// A negative ms1weight means we don't do any rescaling
	if( ms1weight >= 0 ){
		double a = sqrt( ms1weight );
		magnorm( xy2[1]->m, xy2[1]->nrow, xy2[1]->ncol, xy1[1]->m, xy1[1]->nrow, xy1[1]->ncol, a );
	}

	// Now bind the columns of 
	xy1[0]->rbind( xy2[0] );
	delete xy2[0];
	xy1[1]->rbind( xy2[1] );
	delete xy2[1];
	

	return(xy1);

}

#endif 
