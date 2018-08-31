/**
 * \file bitonicBGD.h
 * AUTHOR: Alex Hu
 * CREATE DATE: March 5, 2015
 * PROJECT: DIA Regression
 * \brief Returns peptide abundances in DIA spectra that optimize
 * the bitonic regression objective function
 *****************************************************************************/

/*
#include "alternateBW.h"
#include "scorePeptides.h"
#include <stdio.h>
#include <stdlib.h>
*/
#include "lassoSGD.h"
//#include "lassoADAM.h"
int rssA;
int rssB;

class lassoBGD: public regression {

public:	
	// Initialization
        lassoBGD(): regression() {};
        lassoBGD(dbl_matrix * y, dbl_matrix * x, char const * outdir, double percent_sparsity, bool * e);
	lassoBGD(lassoSGD * sgd);

	void reset();
	int xnnz;
	int maxnz;

	dbl_matrix * derivs0;
	dbl_matrix * XvX2;
//	double * XvX2csr;
//	int * XvX2nzix;
//	int XvX2N;
	dbl_matrix * dls;
	void setB( dbl_matrix * b );

	// Optimization
	void computeL( double l);
	dbl_matrix * Bprev;
	dbl_matrix * By;
	void pl();
	void update_dls();
	double step();	
	virtual void learn();
	virtual void set_stepsize( double stepsize );
	void iterate_track(); // This will iterate and time how long each
	// iteration takes as well as the ls() objective
	double L; // Lipschitz constant	
	double epsilon; // Lipschitz constant	
	void update_epsilon();

	double stept;
	double steptprev;
	int itnum;

	virtual void write_B( const char * fn );	
	~lassoBGD();
	
};
/*
void lassoBGD::update_epsilon(){
	//epsilon = 1.0/( (epochnum) + 2 );
	epsilon = 1.0/( (epochnum*epochnum) + 2 );
	//epsilon = 1.0/( sqrt(epochnum) + 2 );
	printf(" epsilon: %f\n", epsilon);
	epsilons.push_back(epsilon);
};*/

void lassoBGD::computeL( double l = 100.0 ){
	//This should be 2*maximum eigenvalue( XT dot X );	
	//But this is just a constance that is hopefully
	//adequate for general use.
	L = l;
}

void lassoBGD::set_stepsize( double ss=100.0 ){
	L = ss;
}

// Assume that b is the correct dimension
void lassoBGD::setB( dbl_matrix * b ){

	if( b->nrow != B.ncol || b->ncol != B.nrow ){
		printf("B dimensions do not match!\n");
		return;
	}	

	int size = (b->nrow)*(b->ncol);
	for( int i=0; i<size; i++ ){
		B.m[i] = b->m[i];
		Bprev->m[i] = b->m[i];
		By->m[i] = b->m[i];
	}

}

lassoBGD::lassoBGD(dbl_matrix * y, dbl_matrix * x, char const * outdir, double percent_sparsity, bool * e=NULL): regression(y, x, outdir, percent_sparsity, e){
	maxnz=0;
	printf("1. This is X: %ix%i\n", X.nrow, X.ncol );
	
	clock_t begin = clock();

	if( X.oldMIndices == NULL ){
		X.sparsify( &Y );
		M = X.nrow;
	}
	printf("2. This is X: %ix%i\n", X.nrow, X.ncol );
	MT = M*T;
	NT = N*T;	
	itnum = 0;
	
	printf("Initializing...\n");
	rssA = getCurrentRSS();
	X.print();
	Y.print();
	printf("3. This is X: %ix%i\n", X.nrow, X.ncol );
	printf("Computing Derivs 0...\n");
	derivs0 = X.Tmultiply( &Y, -2.0 );
	derivs0->print();
	printf("Tmultiplied!\n");
	E = e;
	if( E != NULL ){
		write_csr_matrix( E, N, T, "etest.txt");
		// Set the original derivatives to be 0
		// where we think the peptides should not be
		for( int i=0; i<NT; i++ ){
			if( E[i] == false ){
				derivs0->m[i] = 0;
			}
		}
	}

	printf("Derivs0: %imb\n", (getCurrentRSS()-rssA)/1000000 );
	//print_matrix( derivs0, N, T, "derivs0");
	rssA = getCurrentRSS();
	printf("This is X: %ix%i\n", X.nrow, X.ncol );
	XvX2 = X.Tmultiply( &X, 2.0 );

	printf("XvX2: %imb\n", (getCurrentRSS()-rssA)/1000000 );
	//print_matrix( XvX2, N, N, "XvX2");
	computeL();	

	rssA = getCurrentRSS();
	printf("B: %imb\n", (getCurrentRSS()-rssA)/1000000 );
	By = new dbl_matrix(N,T); // T rows, N columns
	Bprev = new dbl_matrix(N,T); // T rows, N columns
	dls = new dbl_matrix(N,T);
	
	printf("M: %i, T: %i, N: %i, MT: %i\n", M, T, N, MT );
	debias = false;

	printf("Done constructing!\n");
};

void lassoBGD::write_B( const char * fn = "bB.txt" ){
	char filename[200];
        sprintf(filename,"%s/%s",outdir,fn);
	printf("Printing B (%ix%i)\n", N, T );
	int nnz = write_csr_matrix( B.m, N, T, filename );
}

lassoBGD::lassoBGD(lassoSGD * sgd){

	Y = sgd->Y;
	X = sgd->X;
	outdir = sgd->outdir;
	percent_sparsity = sgd-> percent_sparsity;
	
	M = sgd->M;
	T = sgd->T;
	N = sgd->N;
	NT = sgd->NT;
	MT = sgd->MT;
	stepsize = sgd->stepsize;	

	B = *(new dbl_matrix( N,T ));
	transpose( sgd->B.m, T, N, B.m );
	E = sgd->E;

	maxit = sgd->maxit;
	compute_sse = sgd->compute_sse;
	debias = false;
	lambda = sgd->lambda;

	maxnz=sgd->maxnz;
	printf("1. This is X: %ix%i\n", X.nrow, X.ncol );
	printf("Prev Algorithm: %s\n", sgd->algorithm );
	algorithm = sgd->algorithm;
	printf("Later Algorithm: %s\n", algorithm );
	
	clock_t begin = clock();

	if( X.oldMIndices == NULL ){
		X.sparsify( &Y );
	}
	printf("2. This is X: %ix%i\n", X.nrow, X.ncol );
	itnum = 0;
	
	printf("Initializing...\n");
	rssA = getCurrentRSS();
	X.print();
	Y.print();
	printf("3. This is X: %ix%i\n", X.nrow, X.ncol );
	printf("Computing Derivs 0...\n");
	derivs0 = X.Tmultiply( &Y, -2.0 );
	derivs0->print();
	printf("Tmultiplied!\n");
	printf("Derivs0: %imb\n", (getCurrentRSS()-rssA)/1000000 );
	
	// Enforce the elution prior 
	if( E != NULL ){
		E = transpose(sgd->E, T, N );
		write_csr_matrix( E, N, T, "etest.txt");
		// Set the original derivatives to be 0
		// where we think the peptides should not be
		for( int i=0; i<NT; i++ ){
			if( E[i] == false ){
				derivs0->m[i] = 0;
			}
		}
	}
	
	//print_matrix( derivs0, N, T, "derivs0");
	rssA = getCurrentRSS();
	printf("This is X: %ix%i\n", X.nrow, X.ncol );
	XvX2 = X.Tmultiply( &X, 2.0 );

	printf("XvX2: %imb\n", (getCurrentRSS()-rssA)/1000000 );
	//print_matrix( XvX2, N, N, "XvX2");
	computeL();	

	rssA = getCurrentRSS();
	printf("B: %imb\n", (getCurrentRSS()-rssA)/1000000 );
	// Is it ok if I start with a non-zero B??
	// Yes, I think so!
	By = B.copy_dbl_matrix(); // T rows, N columns
	Bprev = B.copy_dbl_matrix(); // T rows, N columns
	dls = new dbl_matrix(N,T);
	
	printf("M: %i, T: %i, N: %i, MT: %i\n", M, T, N, MT );
	debias = false;

	printf("Done constructing!\n");
};

// Computes delta f(By) and puts it into dls
void lassoBGD::update_dls(){
	// dls is NxT, same as By
	multiply( XvX2->m, N, N, By->m, N, T, 1.0, dls->m);
	//multiply_csr( XvX2csr, XvX2nzix, N, XvX2N, By, N, T, dls );
	add( dls->m, derivs0->m, N, T, dls->m ); 
};

// Sets elements in B equal to pL(By)
void lassoBGD::pl(){

	//clock_t begin = clock();
//	rssA = getCurrentRSS();
	update_dls();		
	for( int i=0; i<NT; i++ ){ (*dls)[i] += lambda; }
	linearcombo( By->m, dls->m, N, T, B.m, 1.0, -1.0/L );	
//	printf("time lapsed: %f (s)\n", double(clock()-begin)/CLOCKS_PER_SEC); 
}

double lassoBGD::step(){
	double val;
	clock_t begin = clock();
	printf("ALGORITHM: %s\n", algorithm);
	if( strcmp( algorithm, "fista") == 0 || strcmp( algorithm, "both") == 0 ){
		double * temp = B.m;
		B.m = Bprev->m;
		Bprev->m = temp;
		pl();
		double a = (steptprev-1)/stept;
		itnum += 1;
		steptprev = stept;
		stept = (1+sqrt(1.0+4.0*steptprev*steptprev ))/2.0;
		// Val is the l1 norm of the absolute value of the difference 
		// between the two iterations
		val = linearcombo( B.m, Bprev->m, N, T, By->m, 1.0+a, -1.0*a );	
        	printf("\n fista in %f seconds. nnz in B: %i\n",double(clock()-begin)/CLOCKS_PER_SEC,B.nnz());
	}
	else{
		printf("Batch\n");	
		itnum += 1;
		//Update derivatives
		multiply( XvX2->m, N, N, B.m, N, T, 1.0, dls->m);
		add( dls->m, derivs0->m, N, T, dls->m ); 
		multiply( XvX2->m, N, N, B.m, N, T, 1.0, dls->m);
		add( dls->m, derivs0->m, N, T, dls->m ); 
		// Val is the l1 norm of the absolute value of the difference 
		// between the two iterations
		double scalar = -1.0/stepsize/(itnum+2);
		printf("stepsize = %f, SCALAR = %f\n", stepsize, scalar );
		val = linearcombo( B.m, dls->m, N, T, B.m, 1.0, scalar );	
        	printf("\n batch in %f seconds. nnz in B: %i\n",double(clock()-begin)/CLOCKS_PER_SEC,B.nnz());
	}
	return(val);
}

void lassoBGD::learn(){
	printf("\nIterating!\n");
	steptprev = 1.0;
	stept = 1.0;
	itnum = 0;
	printf("maxit=%i\n", maxit);
	double delta;
	char filename[500];
	printf("SSE at %i: %f\n", itnum, ls() );
	Tmultiply( X.m, M, N, Y.m, M, T, -2.0, dls->m );	
	printf("computed dls\n");
	delta = step();
	while(  itnum < maxit ){
		printf("%i:%f\n",itnum,delta);
		delta = step();
		printf("SSE at %i: %f\n", itnum, ls() );
		printf("norm1: %f\n", norm1(B.m,B.nrow,B.ncol) );
		printf("lasso objective: %f\n", total_objective() );
		if(  itnum % 200 == 0 ){
		//if( true ){
			sprintf(filename, "%s/%i_bgd_B.txt", outdir, itnum );
			write_csr_matrix( B.m, N, T, filename );
		}
	}
	int peakrss = getPeakRSS();
	printf("Num batch Its: %i\nPeak RSS: %i mb\n", itnum, peakrss/1000000 );
//	sprintf(filename, "%s/%i_bgd_B.txt", outdir, itnum );
//	write_csr_matrix( B.m, N, T, filename );
		
}

lassoBGD::~lassoBGD(){
	printf("bgd destructing\n");
	delete derivs0;
	printf("2\n");
	delete XvX2;
	printf("3\n");
	delete dls;
	printf("4\n");
	delete Bprev;
	printf("5\n");
	delete By;
	printf("6\n");
}

/*
void bitonicBGD::iterate_track(){
	
	clock_t begin;
	//print_matrix( y, m, t);

//	reset();
	steptprev = 1.0;
	stept = 1.0;
	itnum = 0;
	double lspen;
	char filename[999];
	double delta = step();
	while( delta > .010 && itnum < maxit ){
		printf("%i:%f\n",itnum,delta);
		
		begin = clock();
		delta = step();
		printf("time lapsed: %f (s)\n", double(clock()-begin)/CLOCKS_PER_SEC); 	
		lspen = ls();
		printf("ls: %f \n", lspen); 	
		if( itnum % 15 == 0 ){
			sprintf(filename, "%s/%i_bgd_B.txt", outdir, itnum );
			write_csr_matrix( B, N, T, filename );
		}
	}
	int peakrss = getPeakRSS();
	printf("Num Its: %i\nPeak RSS: %i mb\n", itnum, peakrss/1000000 );
	sprintf(filename, "%s/%i_bgd_B.txt", outdir, itnum );
	write_csr_matrix( B, N, T, filename );
		
}


*/
