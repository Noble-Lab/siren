/**
 * \file regression.h
 * AUTHOR: Alex Hu
 * CREATE DATE: March 5, 2015
 * PROJECT: DIA Regression
 * \brief Returns peptide abundances in DIA spectra that optimize
 * the least-squares regression objective function
 *****************************************************************************/
#ifndef REGRESSION_H
#define REGRESSION_H
#include "io.h"
#include "scorePeptides.h"
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

class regression {
public:
	int M; // Number of m/z bins
	int T; // Number of scans
	int N; // Number of candidate peptides
	int MT; // M*T
	int NT;
	int maxit;
	double stepsize;
	double avg_resid;
	int reset_epsilon_period;
	double lambda; // L1 regularization parameter
	double percent_sparsity;
	
	int M1;
	int M2;

	// Input
	dbl_matrix Y; // Observed spectra, MxT
	dbl_matrix X; // Theoretical spectra, MxN
	dbl_matrix B; // Learned abundances, NxT

	dbl_matrix_rowmajor * Yrm;
	dbl_matrix_rowmajor * Xrm;

	bool sparsexy; // True if both X and Y are given as sparse rowmajor matrices // False if not, and X is consequently made sparse.

	char const * outdir;
	char const * algorithm;
	virtual void setB( dbl_matrix * b );

	// If we have a strict prior on where elution profiles
	// can and cannot appear
	bool * E; // NxT for elution time priors. 
		
	// Output
	
	// Initialization
	regression(){};
	regression(dbl_matrix * y, dbl_matrix * x, char const * outdir, double psparsity);
	regression(dbl_matrix * y, dbl_matrix * x, char const * outdir, double psparsity, bool * e);
	void reset();
	void initialize_X2(); // Use this when initializing the theoretical spectra when learning it
	virtual void estimate_lambda();
	bool debias;  
	bool compute_sse;

	// Computing the objective
	virtual double ls();
	double ls(int mstart, int mend);	
	virtual double total_objective();


	// Optimization
	virtual void learn();	
	virtual void write_B(const char * fn);	
	virtual void set_stepsize( double ss );
	virtual ~regression();
	
};

regression::regression( dbl_matrix * y, dbl_matrix * x, char const * od, double psparsity ){

	Y = *y;
	X = *x;

	//Yrm = dynamic_cast<dbl_matrix_rowmajor * >(y);
	Xrm = dynamic_cast<dbl_matrix_rowmajor * >(x);
	sparsexy = false;
	//if( Yrm != NULL && Xrm != NULL ){
	if( Xrm != NULL ){
		printf("X are rowmajor sparse!\n");
		sparsexy = true;
	}	
	else{
		printf("X is not sparse\n");
	}

	outdir = od;	
	percent_sparsity = psparsity;


        M = Y.nrow;
        T = Y.ncol;
        N = X.ncol;
	NT = N*T;
	MT = M*T;
	B = *(new dbl_matrix( T, N, Y.scale/X.scale ));	
	E = NULL;
	
	compute_sse = false;
	debias = false;
//	estimate_lambda();	
	maxit = 100;
	reset_epsilon_period = -1;
}

// e is a boolean matrix that denotes
// in which scans a peptide's inferred abundance
// is allowed to be non-zero.
regression::regression( dbl_matrix * y, dbl_matrix * x, char const * od, double psparsity, bool * e = NULL){

	Y = *y;
	X = *x;

	//Yrm = dynamic_cast<dbl_matrix_rowmajor * >(y);
	Xrm = dynamic_cast<dbl_matrix_rowmajor * >(x);
	sparsexy = false;
	if(  Xrm != NULL ){
	//if( Yrm != NULL && Xrm != NULL ){
		printf("X is rowmajor sparse!\n");
		sparsexy = true;
	}	
	else{	
		printf("X is not sparse\n");
	}

	outdir = od;	
	percent_sparsity = psparsity;

        M = Y.nrow;
        T = Y.ncol;
        N = X.ncol;
	NT = N*T;
	MT = M*T;
	B = *(new dbl_matrix( T, N, Y.scale/X.scale ));	
	E = e;

	compute_sse = false;	
	maxit = 100;	
	debias = false;
	//estimate_lambda();	
	reset_epsilon_period = -1;
}

// If percent_sparsity is positive, treat it as percent sparsity
// If it is negative, treat it as negative lambda
void regression::estimate_lambda(){
        clock_t begin = clock();
	printf("Estimating Lambda\n");
        if( percent_sparsity == 0 ){
                lambda = 0;
        }
	else if( percent_sparsity < 0 ){
		printf("Negative sparsity given\n");
		printf("Interpreting as negative lambda\n");
		printf( "X scale: %f, Y scale: %f\n", X.scale, Y.scale );		
		lambda = -percent_sparsity/X.scale/Y.scale;
		printf("given lambda: %f, scaled lambda: %f\n",-percent_sparsity, lambda);
	}
        else{
		printf( "X: (%ix%i), Y:(%ix%i)\n", X.nrow, X.ncol, Y.nrow, Y.ncol );
                double * half_derivs0 = Tmultiply( X.m, M, N, Y.m, M, T, 1 );
		printf("Computed half_derivs0\n");
		printf("N=%i, T=%i\n", N, T );
		lambda = 2*percentile( half_derivs0, N, T, percent_sparsity );
                delete [] half_derivs0;
        }
        printf("lambda for %f%: %f\n", percent_sparsity, lambda );
        printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000);
}

double regression::ls(){
	double sum = 0;
	if( compute_sse ){
	//	printf("ls()\n");
		double * model = multiply( X.m, M, N, B.m, N, T, 1.0 );
		double r = 0;
	//	clock_t begin = clock();
		for( int i=0; i<MT; i++ ){
			r = model[i]-Y[i];
			sum += r*r;
		}
		delete [] model;
	}
	else{
		sum = -1;
	}
//	printf("least-squares computation: %f (s)\n", double(clock()-begin)/CLOCKS_PER_SEC); 
	return(sum);
};

double regression::ls( int mstart, int mend ){
//	printf("ls()\n");
	int newM = mend-mstart;
	double * model = multiplyT( &(X.m)[mstart*N], newM, N, B.m, T, N, 1.0 );
	double sum = 0;
	double r = 0;
//	clock_t begin = clock();
	int newend = mend*T;
	for( int i=mstart*T; i<MT; i++ ){
		r = model[i]-Y[i];
		sum += r*r;
	}
	delete [] model;
//	printf("least-squares computation: %f (s)\n", double(clock()-begin)/CLOCKS_PER_SEC); 
	return(sum);
};

double regression::total_objective(){
	double objective = ls() + norm1(B.m,T,N)*lambda;
	return( objective );
};

void regression::learn(){
	printf("Regression learn() does nothing.\n");
}

void regression::set_stepsize( double ss ){
	printf("Regression set_stepsize() does nothing.\n");
}

void regression::write_B( const char* fn = "BOO!"){
	printf("Regression write_B() does nothing\n");
}


void regression::setB( dbl_matrix * b ){
	printf("Regression setB() does nothing\n");
}

regression::~regression(){
	printf("Deleting X\n");
	delete &X;
	printf("Deleting B\n");
	delete &B;
	printf("Deleting Y\n");
	delete &Y;
//	printf("Deleting Model\n");
//	delete [] model;
}

#endif
