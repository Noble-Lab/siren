/**
 * \file regression.h
 * AUTHOR: Alex Hu
 * CREATE DATE: March 5, 2015
 * PROJECT: DIA Regression
 * \brief Returns peptide abundances in DIA spectra that optimize
 * the least-squares regression objective function
 *****************************************************************************/
#ifndef lassoSGD_H
#define lassoSGD_H

#include "regression.h"
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

class lassoSGD: public regression {
public:
	void reset();
	int xnnz;
	void sparsifyXandY();	
	int * nzx;// M x maxnz matrix that lists the columns of X where row m is nonzero. Each list ends at -1.	
	double* Xcsr; // Compressed row matrix whose non-zero locations are in nzx.
	int maxnz;
	int * oldMIndices;
	void initialize_X2(); // Use this when initializing the theoretical spectra when learning it

	lassoSGD();
	lassoSGD(dbl_matrix * y, dbl_matrix * x, char const * outdir, double percent_sparsity, bool * e );


	void estimate_lambda_sample( int nsamples, double percent_sparsity );	
	virtual void estimate_lambda();	 
	double epsilon;
	vector<double> epsilons;
	double * Bt;
	virtual void setB( dbl_matrix * b ); // Here, this b should be NxT


	// SGD optimization
	int itnum;
	int epochnum;	
	double gradient;
	int * peaksperpeptide;
	double get_resid( int m, int t );
	virtual	void update(int m, int t ); 
	int* randsequence; // This is for sampling Y without replacement
	virtual double epoch();
	void update_epsilon();
	virtual void learn();
	virtual void set_stepsize( double ss );
	void learn_track();
	bool eprior;

	int nzix;
	int nzend;	
	int n;
	int Bix;
	double diff;

	double ls();
	double oldsse;
	double newsse;
	double ssediff;
	int pollperiod;
	
	// For batch
	void get_resids( );
	double get_resids_fast( );
	double * resids;
	double * gradients;

	// Variables for the function update( int m, int t, double epsilon );
	double residU;
	double deltaU;
	double max_deltaU; // Contains the maximum delta encountered for any Bn,t in the update

	int nzixU;
	int nzendU;
	int BixU;

	// Output stuff
	virtual void write_B( const char * fn );
	~lassoSGD();
};

// If percent_sparsity is positive, treat it as percent sparsity
// If it is negative, treat it as negative lambda
void lassoSGD::estimate_lambda(){
        clock_t begin = clock();
	printf("Estimating Lambda\n");
        if( percent_sparsity == 0 ){
                lambda = 0;
        }
	else if( percent_sparsity < 0 ){
		printf("Negative sparsity given\n");
		printf("Interpretting it as the negative value of lambda\n");
		printf( "X scale: %f, Y scale: %f\n", X.scale, Y.scale );		
		lambda = -percent_sparsity/X.scale/Y.scale;
		printf("given lambda: %f, scaled lambda: %f\n",-percent_sparsity, lambda);
	}
        else{
		printf( "X: (%ix%i), Y:(%ix%i)\n", X.nrow, X.ncol, Y.nrow, Y.ncol );
                //double * half_derivs0 = Tmultiply( X.m, M, N, Y.m, M, T, 1 );
		printf("Sparse multiplying\n");
		double * half_derivs0 = sparseTmultiply( Xcsr, nzx, Y.m , M, T, N, maxnz, 1.0 );
		printf("Computed half_derivs0\n");
		printf("N=%i, T=%i\n", N, T );
		lambda = 2*percentile( half_derivs0, N, T, percent_sparsity );
                delete [] half_derivs0;
        }
        printf("lambda for %f%: %f\n", percent_sparsity, lambda );
        printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000);
}

lassoSGD::lassoSGD(){
	
}

lassoSGD::lassoSGD(dbl_matrix * y, dbl_matrix * x, char const * outdir, double percent_sparsity, bool * e = NULL ): regression(y, x, outdir, percent_sparsity, e){

	printf("constructing!\n");
	maxnz = 0;
	printf("X: (%ix%i)\n",x->nrow,x->ncol);
	printf("Y: (%ix%i)\n",y->nrow,y->ncol);
	
	clock_t begin = clock();
	peaksperpeptide = zeros_int(N,1);
	
	if( ! sparsexy ){
		printf("Sparsifying M=%i to...",M); fflush(stdout);
		sparsifyXandY();
	}
	else{
		peaksperpeptide = Xrm->nz_per_column;
		Xcsr = Xrm->X;
		nzx = Xrm->nzx;
		maxnz = Xrm->maxnz;	
	}

	MT = M*T;	
	NT = N*T;	
	Bt = new double[NT];
	eprior = false;
	if( E != NULL ){	
		bool * Et = transpose(E, N, T );
		printf("!!!! (N,T) = (%ix%i)\n", N, T );
		delete [] E;
		E = Et;
		eprior = true; 
	}	

	estimate_lambda();
	pollperiod = 1;
	oldsse = 0;
	for( int i=0; i<MT; i++ ){
		if( Y[i] != 0 ){
			oldsse += Y[i]*Y[i];
		}
	}	
	printf("Initial SSE: %f\n", oldsse );
	
	printf("M: %i, T: %i, N: %i,  MT: %i\n", M, T, N, MT );
	epsilon = 1;
	epsilons.push_back(epsilon);

	itnum = 1;
	epochnum = 0;
	debias = false;

	//srand ( unsigned ( std::time(NULL) ) );
	srand ( unsigned (0 ) );
	randsequence = new int[MT];
	for( int i=0; i<MT; i++ ){ randsequence[i] = i; } 
};

void lassoSGD::setB( dbl_matrix * b ){

	if( b->nrow != B.ncol || b->ncol != B.nrow ){
		printf("B dimensions do not match!\n");
		printf("Given b: %ix%i\n", b->nrow, b->ncol);
		printf("Target B: %ix%i\n", B.ncol, B.nrow);
		return;
	}	
	printf("Setting B!\n");
	int size = b->nrow*b->ncol;
	for( int i=0; i<size; i++ ){
		Bt[i] = b->m[i];
	}
	transpose( b->m, b->nrow, b->ncol, B.m );
}

void lassoSGD::write_B( const char * fn = "sB.txt" ){
	char filename[200];
        sprintf(filename,"%s/%s",outdir,fn);
	transpose(B.m,T,N,Bt);
	printf("Printing Bt (%ix%i), scale: %f\n", N, T, B.scale );
	int nnz = write_csr_matrix( Bt, N, T, filename, B.scale );
	printf("Printed!\n");
}

// Looks for m/z bins that are unpopulated by X and removes them from
// both X and Y
void lassoSGD::sparsifyXandY(){

	clock_t begin = clock();
	// Identify the m's to keep
	//bool occupied[M] = {false};
	bool * occupied = new bool[M];
	for( int m=0; m<M; m++ ){ occupied[m] = false; };
	int newM = 0;
	int offset;
	int nz;
	maxnz = 0;
	xnnz = 0;

	for( int m=0; m<M; m++ ){
		nz=0;
		for( int n=0; n<N; n++ ){
			if( X[m*N+n] > 0 ){
				if( ! occupied[m] ){
					newM++;
				}
				occupied[m] = true;
				peaksperpeptide[n] += 1;
				nz+=1;		
			}
		}
		xnnz += nz;
		if( nz > maxnz ){ maxnz = nz; }
	}
	printf(" M=%i\n", newM);
	
	// Create the new X and Y
	dbl_matrix *newX = new dbl_matrix(newM,N,X.scale);
	dbl_matrix *newY = new dbl_matrix(newM,T,Y.scale);
	nzx = zeros_int(newM,maxnz); for(int i=0; i<newM*maxnz; i++ ){ nzx[i] = -1; }
	oldMIndices = zeros_int(newM,1);
	Xcsr = zeros(newM,maxnz); 
	// Copy into them only the m's to keep
	int newm = 0;
	int xnoffset;
	int nzix;

	for( int m=0; m<M; m++ ){
		nzix = newm*maxnz;
		if( !occupied[m] ){ continue; }
		for( int n=0; n<N; n++ ){
			if( X[m*N+n] > 0.0 ){
				newX->m[newm*N+n] = X[m*N+n];
				nzx[nzix] = n;
				Xcsr[nzix] = X[m*N+n];
				nzix++;
			}
		}
		// Copy Ys
		for( int t=0; t<T; t++ ){
			if( Y[m*T+t] == 0 ){ continue; }
			newY->m[newm*T+t] = Y[m*T+t];
		}
		oldMIndices[newm] = m;
		newm++;
	}
	if( DELETE ){
		X.~dbl_matrix();
		Y.~dbl_matrix();
	}
	X = *newX;
	Y = *newY;
	M = newM;
	MT = M*T;	
	printf("maxnz=%i\n",maxnz);
	printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000); 
};

void lassoSGD::initialize_X2(){
	return;	
}

void lassoSGD::reset(){
	for( int i=0; i<N*T; i++ ){
		//B[i] = 0;
	}
	itnum=1;
	epochnum=0;
};

// Optimization
double lassoSGD::get_resid( int m, int t ){
	double resid = 0;	
	/*
	if( sparsexy ){
		// Here, the t actually refers to a column number within the sparse matrix, not
		// the literal scan number. Then each column refers to a specific scan number.
		nzix = m*(Yrm->maxnz);
		nzend = nzix + Yrm->maxnz;
		while( nzix < nzend ){
			if( Yrm->nzx[nzix] == t ){
				resid = -1*Yrm->X[nzix];
				break;
			}
			nzix++;
		}

	}
	else{
		resid = -1*Y[m*T+t];	
	}*/

	resid = -1*Y[m*T+t];	
	int nzix = m*maxnz;
	int n=nzx[nzix];
	int nzend = nzix+maxnz;	
	int Bix;

	while( n != -1 ){
		Bix = t*N+n;
		if( B[Bix] != 0 ){
			resid += B[Bix]*Xcsr[nzix];
		}
		nzix++;
		if( nzix == nzend ){ break; }
		n=nzx[nzix];
	}
	return( resid );
};



void lassoSGD::update(int m, int t ){ // Returns the value of greatest update
	residU = get_resid( m, t );
	avg_resid += residU*residU;
	max_deltaU = 0;
	
	nzix = m*maxnz;
	nzend = nzix+maxnz;	
	n=nzx[nzix];

	//nzx is a M*maxnz matrix that, which each row m contains
	//the peptides that have a theoretical peak at m/z bin m
	while( n != -1 ){
		Bix = t*N+n;
		if( (eprior && E[Bix]) || !eprior ){	
//		if( true ){	
			gradient = 2*Xcsr[nzix]*residU;
			if( debias ){ 
				if( B[Bix] > 0.0 ){
					deltaU = min( gradient*epsilon, B[Bix] ); 
				}
			}
			else{
				deltaU = min( (gradient + lambda/peaksperpeptide[n])*epsilon, B[Bix] );
			}
	
			B[Bix] = B[Bix] - deltaU;
			if( deltaU < 0 ){ deltaU = -1*deltaU; }
			if( deltaU > max_deltaU ){
				max_deltaU = deltaU;
			}
		}
		nzix++;
		if( nzix == nzend ){ break; }
		n=nzx[nzix];
	}
	itnum++;
	//update_epsilon();
};


// Optimization
// Performs XB-Y and puts the result in resids

// Optimization
// Performs XB-Y and puts the result in resids
double lassoSGD::get_resids_fast(){
	// Initialize the residuals with the values in Y times -1
	if( resids == NULL ){
		resids = new double[MT];
	}
	for( int i=0; i<MT; i++ ){
		resids[i] = -Y[i];
	}
	printf("Getting resids fast!\n");
	int nzix;
	int n;
	int nzend;
	int Bix;
	int rix;
	int tN;
	for( int m=0; m<M; m++ ){
		rix = m*T;	
		for( int t=0; t<T; t++ ){	
			nzix = m*maxnz;
			n=nzx[nzix];
			nzend = nzix+maxnz;	
			while( n != -1 ){
				//printf("-m:%i, t:%i, n:%i\n",m,t,n);
				Bix = t*N+n;
				if( B[Bix] != 0 ){
					resids[rix] += B[Bix]*Xcsr[nzix];
				}
				nzix++;
				if( nzix == nzend ){ break; }
				n=nzx[nzix];
			}
			rix++;
		}
	}
	printf("-done\n");
	//return(0);//?return( norm2(resids,M,T) );
	return( norm2(resids,M,T) );
};
// Optimization
void lassoSGD::get_resids(){

	// Initialize the residuals with the values in Y times -1
	if( resids == NULL ){
		resids = new double[MT];
	}
	multiplyT( X.m, M, N, B.m, T, N, 1.0, resids );
	for( int i=0; i<MT; i++ ){
		resids[i] = resids[i]-Y[i];
	}
};

double lassoSGD::epoch(){
        double max_delta=0;
	avg_resid=0;
        double delta;
        int m;
        int t;
        clock_t begin = clock();
        begin = clock();
	random_shuffle(&randsequence[0], &randsequence[MT]);
	update_epsilon();
        for( int i=0; i<MT; i++ ){
		m = randsequence[i]/T;
                t = randsequence[i]-m*T;
		update( m, t );
                if( max_deltaU > max_delta ){
                        max_delta = max_deltaU;
                }
        }
	avg_resid /= MT;
        printf("\n in %f seconds. Max delta: %f, nnz in B: %i\n",double(clock()-begin)/CLOCKS_PER_SEC,max_delta,B.nnz());
        printf("total_resid: %f\n",avg_resid*MT);
        epochnum++;
        return(max_delta);
};


void lassoSGD::update_epsilon(){
	epsilon = 1.0/( (stepsize*epochnum*2) + 2.0 );
	//epsilon = 1.0/( (stepsize*max(epochnum,15)) + 2.0 );
//	epsilon = 1.0/( stepsize*(epochnum*epochnum) + 2 );
	//epsilon = 1.0/( stepsize*sqrt(epochnum) + 2 );
	printf("stepsize = %f, epochnum: %i, epsilon: %f\n", stepsize, epochnum, epsilon);
	epsilons.push_back(epsilon);
};

void lassoSGD::learn(){
	printf("Outdir: %s\n",outdir);
	/*
	if( sparsexy ){
		printf("Magnitude of Y: %f\n", Yrm->norm2());
	}
	else{
		printf("Magnitude of Y: %f\n", Y.norm2());
	}*/
	printf("Magnitude of Y: %f\n", Y.norm2());
	itnum = 1;
	epsilon = 1;
	double max_delta;
	clock_t begin;
	int nnz;
	double ratio;
	double lasso_objective;
	
	for( int i=0; i<maxit; i++ ){
		max_delta = epoch();
	        printf("Peak RSS: %i mb\n", getPeakRSS()/1000000 );
		if( i < 2 ){
			newsse = ls();		
			lasso_objective = total_objective();
			ssediff = oldsse - newsse;
			printf("SSE at %i epochs: %f\n", i+1, newsse );
			printf("norm1 at %i epochs: %f\n", i+1, norm1(B.m, B.nrow, B.ncol) );
			printf("objective at %i epochs: %f\n", i+1, lasso_objective );
		}
		else if( i % (pollperiod) == 0 ){
			oldsse = newsse;
			newsse = ls();
			lasso_objective = total_objective();
			printf("SSE at %i epochs: %f\n", i+1, newsse );
			printf("norm1 at %i epochs: %f\n", i+1, norm1(B.m, B.nrow, B.ncol) );
			printf("objective at %i epochs: %f\n", i+1, lasso_objective );
			ratio = (oldsse-newsse)/ssediff;
//			printf(" ratio: %f\n", ratio );
		}
		/*
		if(  itnum % 200 == 0 ){
		//if( true ){
			sprintf(filename, "%s/%i_bgd_B.txt", outdir, itnum );
			write_csr_matrix( B.m, N, T, filename );
		}*/
	}
	write_B();
	printf("Num epochs: %i\n", epochnum );	

};	

double lassoSGD::ls(){
	double sum = 0;
	if( compute_sse ){
	//if( false ){
		double * model = multiplyT( X.m, M, N, B.m, T, N, 1.0 );
		double r = 0;
		for( int i=0; i<MT; i++ ){
			r = model[i]-Y[i];
			sum += r*r;
		}
		delete [] model;
	}
	else{
		sum = -1;
	}
	return(sum);
};

void lassoSGD::set_stepsize( double ss=1.0 ){
	stepsize = ss;
}

lassoSGD::~lassoSGD(){
	printf("--Destructing a lassoSGD class--\n");
	delete [] Bt;
	delete [] randsequence;
	delete [] nzx;
	delete [] Xcsr;
	delete [] oldMIndices;
	delete [] peaksperpeptide;	
	if( resids != NULL ){
		delete [] resids;
	}
	if( gradients != NULL ){
		delete [] gradients;
	}
}


#endif
