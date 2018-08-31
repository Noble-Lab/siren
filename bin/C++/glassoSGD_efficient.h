/**
 * \file regression.h
 * AUTHOR: Alex Hu
 * CREATE DATE: March 5, 2015
 * PROJECT: DIA Regression
 * \brief Returns peptide abundances in DIA spectra that optimize
 * the least-squares regression objective function
 *****************************************************************************/
#ifndef glassoSGD_H
#define glassoSGD_H

#include "produce_fragments_efficient.h"
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <vector>


class glassoSGD: public regression {
public:
	// Initialization
	glassoSGD(): regression() {};
	glassoSGD(dbl_matrix * y, dbl_matrix * x, char const * outdir, double percent_sparsity, dbl_matrix * x1, double o1, double bw1, double o2, double bw2, bool * e, char const * precursor_annfile );
	//glassoSGD(dbl_matrix * y, dbl_matrix * x, char const * outdir, double percent_sparsity, double o2, double bw2, bool * e, char const * precursor_annfile );
	void reset();
	int xnnz;
	void sparsifyXandY();	
	void estimate_lambda();
	int * nzx;// M x maxnz matrix that lists the columns of X where row m is nonzero. Each list ends at -1.	
	double* Xcsr; // Compressed row matrix whose non-zero locations are in nzx.
	int maxnz;
	int * oldMIndices;
	void initialize_X2(); // Use this when initializing the theoretical spectra when learning it

	void estimate_lambda_sample( int nsamples, double percent_sparsity );	
	double epsilon;
	vector<double> epsilons;
	double * Bt;
	virtual void setB( dbl_matrix * b ); // Here, this b should be NxT

	// SGD optimization
	int itnum;
	int epochnum;	
	int batchnum;
	double lsgradient; // with respect to just the least-squares loss
	int * peaksperpeptide;
	double get_resid( int m, int t );
	void get_resids( );
	double get_resids_fast();
	virtual	void update(int m, int t ); 
	int* randsequence; // This is for sampling Y without replacement
	virtual double epoch();
	double batch();
	int batchits;
	void update_epsilon();
	void update_epsilon_batch();
	virtual void learn();
	virtual void set_stepsize( double ss );
	void learn_track();
	bool eprior;
	double total_resid;
	int nzp;

	double ls();
	double total_objective();
	double oldsse;
	double newsse;
	double ssediff;
	int pollperiod;
	double * grouppenalty_per_scan();
	double * resids;

	// Variables for the function update( int m, int t, double epsilon );
	double gp; // Group penalty for the variable in question
	double dgp; // Gradient of the group penalty for the variable in question
	int gpnorm; // Decides the relationship between dgp and abundance. 0 is for none, 1 is for dgp*abundance, 
	// 2 is for dgp*abundance*abundance, etc.
	double residU;
	double deltaU;
	double max_deltaU; // Contains the maximum delta encountered for any Bn,t in the update
	double * peptide_magnitudes;

	// Group Lasso Variables
	// that stores the penalty associated with peptide n and group pn.
	fragmentgroup * fg;
	fragmentgroup::groupcsrmatrix * gm;	
	double * Bcsr; // Switched from glassoSGD.h
	// row in Bcsr corresponds to a scan/time point, and each column is a peptide
	// whose identity depends on bzx	
	int * nzb;	

	int nzixU;
	int nzendU;
	int BixU;

	// Output stuff
	virtual void write_B( const char * fn );
	~glassoSGD();
};

glassoSGD::glassoSGD(dbl_matrix * y, dbl_matrix * x, char const * outdir, double percent_sparsity, dbl_matrix * x1, double o1, double bw1, double o2, double bw2, bool * e = NULL, char const * precursor_annfile = NULL ): regression(y, x, outdir, percent_sparsity, e){

	printf("constructing!\n");
	maxnz = 0;
	batchits = 0;	
	batchnum = 0;
	resids = NULL;
	printf("Sparsifying M=%i to...",M); fflush(stdout);
	clock_t begin = clock();
	peaksperpeptide = zeros_int(N,1);
	sparsifyXandY();
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
	
	stepsize = 1.0;
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
	peptide_magnitudes = column_magnitude( X.m, X.nrow, X.ncol  );

	// Build the group lasso stuff
	if( precursor_annfile == NULL ){
		fg = new fragmentgroup( x1, T, o1, bw1, o2, bw2, outdir );
	}
	else{
		fg = new fragmentgroup( precursor_annfile, N, T, o1, bw1, o2, bw2, outdir );
	}
	estimate_lambda();

	Bcsr = Xcsr; // Switched from glassoSGD.h
	// row in Bcsr corresponds to a scan/time point, and each column is a peptide
	// whose identity depends on bzx	
	//print_matrix( peptide_magnitudes, x->ncol, 1, "Magnitudes" );
	nzb = nzx;	
};

void glassoSGD::setB( dbl_matrix * b ){

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

void glassoSGD::write_B( const char * fn = "sB.txt" ){
	char filename[200];
        sprintf(filename,"%s/%s",outdir,fn);
	transpose(B.m,T,N,Bt);
	printf("Printing Bt (%ix%i)\n", N, T );
	int nnz = write_csr_matrix( Bt, N, T, filename );
}

// If percent_sparsity is positive, treat it as percent sparsity
// If it is negative, treat it as negative lambda
void glassoSGD::estimate_lambda(){
        clock_t begin = clock();
	printf("Estimating Lambda\n");
        if( percent_sparsity == 0 ){
                lambda = 0;
        }
	else if( percent_sparsity < 0 ){
		printf("Negative lambda given\n");
		lambda = -percent_sparsity;
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

// Looks for m/z bins that are unpopulated by X and removes them from
// both X and Y
void glassoSGD::sparsifyXandY(){

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
	dbl_matrix *newX = new dbl_matrix(newM,N);
	dbl_matrix *newY = new dbl_matrix(newM,T);
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
};

void glassoSGD::initialize_X2(){
	return;	
}

void glassoSGD::reset(){
	for( int i=0; i<N*T; i++ ){
		//B[i] = 0;
	}
	itnum=1;
	epochnum=0;
};

// Optimization
double glassoSGD::get_resid( int m, int t ){
	double resid = -1*Y[m*T+t];	
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

// Optimization
void glassoSGD::get_resids(){

	// Initialize the residuals with the values in Y times -1
	if( resids == NULL ){
		resids = new double[MT];
	}
	/*/
	for( int i=0; i<MT; i++ ){
		resids[i] = -Y[i];
	}

	int nzix;
	int n;
	int nzend;
	int Bix;
	int rix;

	for( int m=0; m<M; m++ ){
		nzix = m*maxnz;
		n=nzx[nzix];
		nzend = nzix+maxnz;	

		rix = m*T;	
		for( int t=0; t<T; t++ ){	
			nzix = m*maxnz;
			while( n != -1 ){
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
	*/
	multiplyT( X.m, M, N, B.m, T, N, 1.0, resids );
	for( int i=0; i<MT; i++ ){
		resids[i] = resids[i]-Y[i];
	}
};

// Optimization
// Performs XB-Y and puts the result in resids
double glassoSGD::get_resids_fast(){
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
	double norm2resids = norm2(resids,M,T);
	return( norm2resids );
};

// Here, Y is TxM
// t and m switched places from lassoSGD.h
double glassoSGD::batch(){ // Returns the value of greatest update

        clock_t begin = clock();
        begin = clock();
	update_epsilon_batch();
	double total_resid = get_resids_fast();
	printf("Get fast resid magnitude: %f\n", total_resid );
	int max_delta = 0;
	int nzix;
	int nzend;	
	int n;
	int Xix; // Replaces Bix
	double diff;
	
	dbl_matrix& X2 = B; // Switched from LassoSGD.h, it is MxN

	int p; // This is the index of the fragmentation group in Fn
	// F is a list of pointers to matrices
	// F is a list of pointers to matrices
	double * Fnmcsr;// 
	double * Fnmcsr_gradients;// 
	int Fnix; //
	int * nzFnm; // This will be the nonzero fragmentation groups associated with a particular
	// peptide n and mz bin m
	max_deltaU = 0;	
	
	int residix = 0;
	for( int t=0; t<M; t++ ){ // t indexes the rows
		for( int m=0; m<T; m++ ){ // m indexes the columns
			nzix = t*maxnz;
			nzend = nzix+maxnz;	
			n=nzb[nzix];	
			//nzx is a M*maxnz matrix that, which each row m contains
			//the peptides that have a theoretical peak at m/z bin m
			//For each peptide that is non-zero at a time point
			// Here, lambda will only be for the group lasso. There is no 
			// L1 regularization for X2 as a whole.
			// If precursor n is present at time t
			while( n != -1 ){
				gm = fg->precursorgroups[n];
				// If precursor n possesses groups that can contribute to m
				if( m < gm->maxm ){
					//printf("m: %i, maxm: %i, n: %i, t: %i\n",m,gm->maxm,n,t);
					Xix = m*N+n; // Find the index in X2 corresponding to the mass bin m in peptide n
					lsgradient = 2*Bcsr[nzix]*resids[residix];
			//		if( lsgradient != 0 ){ printf("\tlsgradient: %f\n",glsradient); }
					Fnmcsr = &(gm->csr[m*gm->maxp]); //FIXME
					Fnmcsr_gradients = &(gm->csr_gradients[m*gm->maxp]); //FIXME 
					nzFnm = &(gm->nz[m*gm->maxp]); // FIXMEmes Grimmelman
	
					Fnix = 0;
					p = nzFnm[Fnix];
					while( p != -1 && Fnix < gm->maxp ){
						// Add to the total gradient
						Fnmcsr_gradients[Fnix] += lsgradient;			
						Fnix++;
						p = nzFnm[Fnix];
					}
				}
				nzix++;
				if( nzix == nzend ){ break; }
				n=nzx[nzix];
			}
			residix++;	
		}
	}
	// Compute the updates and then perform them
	// Then reinitialize the lsgradients to 0
	for( int n=0; n<N; n++ ){
		gm = fg->precursorgroups[n];
		for( int m=0; m<gm->maxm; m++ ){
			Xix = m*N+n; // Find the index in X2 corresponding to the mass bin m in peptide n
			Fnmcsr = &(gm->csr[m*gm->maxp]); //FIXME
			Fnmcsr_gradients = &(gm->csr_gradients[m*gm->maxp]); //FIXME
			nzp = gm->nzperm[m];
			nzFnm = &(gm->nz[m*gm->maxp]); // FIXME
	
			Fnix = 0;
			p = nzFnm[Fnix];
			while( p != -1 && Fnix < gm->maxp ){
				// Add to the total gradient somehow
				gp = gm->penalties[p];	
				if( gp == 0 ){
					dgp = lambda;
				}
				else{
					dgp = lambda*Fnmcsr[Fnix]/gp;	
				}	
				// scale dgp by the abundance
				for( int gpn = 0; gpn<gpnorm; gpn++ ){
					dgp *= peptide_magnitudes[n];
				}
				deltaU = min( (Fnmcsr_gradients[Fnix]+dgp)*epsilon, Fnmcsr[Fnix] );
				deltaU /= nzp;
				Fnmcsr_gradients[Fnix] = 0;			
                                gm->penalties[p] = sqrt( max(gp*gp + deltaU*(deltaU-2*Fnmcsr[Fnix]),0.0) ); // Update penalty
                                Fnmcsr[Fnix] = Fnmcsr[Fnix] - deltaU; // Update Fnmcsr
                                X2[Xix] = X2[Xix] - deltaU; // Update X2; 
                                if( fabs(deltaU) > fabs(max_deltaU) ){
                                        max_deltaU = deltaU;
                                } 
				Fnix++;
				p = nzFnm[Fnix];
			}
		}
	}
	itnum++;
	batchnum++;
        printf("in %f seconds. Max delta: %f, total_resid: %f , nnz in B: %i\n", double(clock()-begin)/CLOCKS_PER_SEC, max_deltaU, total_resid, B.nnz());
	return( max_deltaU );
};

// Here, Y is TxM
// t and m switched places from lassoSGD.h
void glassoSGD::update(int t, int m ){ // Returns the value of greatest update
	residU = get_resid( t, m );
	total_resid += residU*residU;
	max_deltaU = 0;
	int nzix = t*maxnz;
	int nzend = nzix+maxnz;	
	int n=nzb[nzix];
	int Xix; // Replaces Bix
	double diff;
	
	dbl_matrix& X2 = B; // Switched from LassoSGD.h, it is MxN
	
	int p; // This is the index of the fragmentation group in Fn
	// F is a list of pointers to matrices
	// F is a list of pointers to matrices
	double * Fnmcsr;// This should be replaced with moop.
	int Fnix; //
	int * nzFnm; // This will be the nonzero fragmentation groups associated with a particular
	// peptide n and mz bin m
	//nzx is a M*maxnz matrix that, which each row m contains
	//the peptides that have a theoretical peak at m/z bin m
	//For each peptide that is non-zero at a time point
	// Here, lambda will only be for the group lasso. There is no 
	// L1 regularization for X2 as a whole.
	// If precursor n is present at time t
	while( n != -1 ){
		//printf(" n: %i\n", n);
		//printf(" fg->precursorgroups.size(): %i\n", fg->precursorgroups.size() );
		gm = fg->precursorgroups[n];
		nzp = gm->nzperm[m]; // FIXME
		// If precursor n possesses groups that can contribute to m
		if( m < gm->maxm ){
			//printf("m: %i, maxm: %i, n: %i, t: %i\n",m,gm->maxm,n,t);
			Xix = m*N+n; // Find the index in X2 corresponding to the mass bin m in peptide n
			lsgradient = 2*Bcsr[nzix]*residU;
			Fnmcsr = &(gm->csr[m*gm->maxp]); //FIXME
			nzFnm = &(gm->nz[m*gm->maxp]); // FIXME
	
			Fnix = 0;
			p = nzFnm[Fnix];
			while( p != -1 && Fnix < gm->maxp ){
				// Compute dgp, the gradient w.r.t. the group penalty
				gp = gm->penalties[p];	
				if( gp == 0 ){
					dgp = lambda;
				}
				else{
					dgp = lambda*Fnmcsr[Fnix]/gp;	
				}	
				// scale dgp by the abundance
				for( int gpn = 0; gpn<gpnorm; gpn++ ){
					dgp *= peptide_magnitudes[n];
				}

				// Compute deltaU
				deltaU = min( (lsgradient + dgp/peaksperpeptide[n])*epsilon, Fnmcsr[Fnix] );	
				deltaU /= nzp;
				if( deltaU > 10 ){
					//printf("resid: %f, Fnmcsr[Fnix]:%f, t:%i, m:%i, p:%i, fnix: %i, n:%i\n", residU, Fnmcsr[Fnix], t, m,p,Fnix, n );
				}
				// Update the penalty, the intensities in Fn, and X2	
				gm->penalties[p] = sqrt( max(gp*gp + deltaU*(deltaU-2*Fnmcsr[Fnix]),0.0) ); // Update penalty
				Fnmcsr[Fnix] = Fnmcsr[Fnix] - deltaU; // Update Fnmcsr
				X2[Xix] = X2[Xix] - deltaU; // Update X2; 

				if( fabs(deltaU) > fabs(max_deltaU) ){
					max_deltaU = deltaU;
				}

				Fnix++;
				p = nzFnm[Fnix];
			}
		}
		nzix++;
		if( nzix == nzend ){ break; }
		n=nzx[nzix];
	}
	itnum++;
	//update_epsilon();
};

double glassoSGD::epoch(){
        double max_delta=0;
        double delta;
        int m;
        int t;
        clock_t begin = clock();
        begin = clock();
	random_shuffle(&randsequence[0], &randsequence[MT]);
	update_epsilon();
	total_resid = 0;
        for( int i=0; i<MT; i++ ){
		m = randsequence[i]/T;
                t = randsequence[i]-m*T;
		update( m, t );
                if( max_deltaU > max_delta ){
                        max_delta = max_deltaU;
                }
        }
        printf("\n in %f seconds. Max delta: %f, total_resid: %f , nnz in B: %i\n",double(clock()-begin)/CLOCKS_PER_SEC,max_delta,total_resid, B.nnz());
        epochnum++;
        return(max_delta);
};

void glassoSGD::update_epsilon(){
	int stepvar = (epochnum+batchnum);
	if( reset_epsilon_period > 0 ){
		stepvar = stepvar % reset_epsilon_period;
	}
	epsilon = 1.0/( stepsize*stepvar + 3.0 );
	//epsilon = 1.0/( stepsize*(epochnum*epochnum) + 2 );
	//epsilon = 1.0/( stepsize*sqrt(epochnum+batchnum) + 2 );
	printf("stepsize = %f, epochnum: %i, epsilon: %f\n", stepsize, epochnum, epsilon);
	epsilons.push_back(epsilon);
};

void glassoSGD::update_epsilon_batch(){
	int stepvar = (epochnum+batchnum)/2;
//	if( reset_epsilon_period > 0 ){
//		stepvar = stepvar % reset_epsilon_period;
//	}
	epsilon = 1.0/( stepsize*stepvar + 3.0 );
	printf("stepsize = %f, epochnum: %i, epsilon: %f\n", stepsize, epochnum, epsilon);
	epsilons.push_back(epsilon);
};

void glassoSGD::learn(){
	printf("Outdir: %s\n",outdir);
	printf("Magnitude of Y: %f\n", Y.norm2());
	
	itnum = 1;
	epsilon = 1;
	double max_delta;
	clock_t begin;
	int nnz;
	double ratio;
	double glasso_objective;
	char filename[500];	
	//maxit = 0;
	printf("\n\nMAXIT:%i\n",maxit);
	for( int i=0; i<maxit; i++ ){
		max_delta = epoch();
	        printf("Peak RSS: %i mb\n", getPeakRSS()/1000000 );
		if( i < 2 ){
			newsse = ls();		
			glasso_objective = total_objective();
			ssediff = oldsse - newsse;
			printf("SSE at %i epochs: %f\n", i+1, newsse );
			printf("norm1 at %i epochs: %f\n", i+1, norm1(B.m, B.nrow, B.ncol) );
			printf("objective at %i epochs: %f\n", i+1, glasso_objective );
		}
		else if( i % (pollperiod) == 0 ){
			oldsse = newsse;
			newsse = ls();
			glasso_objective = total_objective();
			printf("SSE at %i epochs: %f\n", i+1, newsse );
			printf("norm1 at %i epochs: %f\n", i+1, norm1(B.m, B.nrow, B.ncol) );
			printf("objective at %i epochs: %f\n", i+1, glasso_objective );
			ratio = (oldsse-newsse)/ssediff;
//			printf(" ratio: %f\n", ratio );
		}
		//if(  itnum % 200 == 0 ){
	}
	if( true ){
		sprintf(filename, "%s%i_sgd_B.txt", outdir, maxit );
		printf("The filename is %s\n", filename );
		write_csr_matrix( B.m, T, N, filename );
	}
	printf("Num epochs: %i\n", epochnum );	
	
	for( int i=0; i<batchits; i++ ){
		max_delta = batch();
	        printf("Peak RSS: %i mb\n", getPeakRSS()/1000000 );
		if( i < 2 ){
			newsse = ls();		
			glasso_objective = total_objective();
			ssediff = oldsse - newsse;
			printf("SSE at %i batches: %f\n", i+1, newsse );
			printf("norm1 at %i batches: %f\n", i+1, norm1(B.m, B.nrow, B.ncol) );
			printf("objective at %i batches: %f\n", i+1, glasso_objective );
		}
		else if( i % (pollperiod) == 0 ){
			oldsse = newsse;
			newsse = ls();
			glasso_objective = total_objective();
			printf("SSE at %i batches: %f\n", i+1, newsse );
			printf("norm1 at %i batches: %f\n", i+1, norm1(B.m, B.nrow, B.ncol) );
			printf("objective at %i batches: %f\n", i+1, glasso_objective );
			ratio = (oldsse-newsse)/ssediff;
//			printf(" ratio: %f\n", ratio );
		}
//		if(  i >= (batchits-6)  ){
//			sprintf(filename, "%s%i_bgd_B.txt", outdir, i );
//			printf("The filename is %s\n", filename );
//			write_csr_matrix( B.m, T, N, filename );
//		}
		
	}
	
	if( true ){
		sprintf(filename, "%s%i_bgd_B.txt", outdir, batchits );
		printf("The filename is %s\n", filename );
		write_csr_matrix( B.m, T, N, filename );
	}



};	

double glassoSGD::ls(){
	double sum = 0;
	if( compute_sse ){
//	if( true ){
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

void glassoSGD::set_stepsize( double ss=1.0 ){
	stepsize = ss;
}

double glassoSGD::total_objective(){
	double objective = ls() + fg->group_penalty()*lambda;//norm1(B.m,T,N)*lambda;
	return( objective );
};

glassoSGD::~glassoSGD(){
	delete [] Bt;
	delete [] randsequence;
	delete [] nzx;
	delete [] Xcsr;
	delete [] oldMIndices;
	delete [] peaksperpeptide;	
	delete [] peptide_magnitudes;
	if( resids != NULL ){
		delete [] resids;
	}	
}

#endif
