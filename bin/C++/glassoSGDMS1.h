/**
 * \file regression.h
 * AUTHOR: Alex Hu
 * CREATE DATE: March 5, 2015
 * PROJECT: DIA Regression
 * \brief Returns peptide abundances in DIA spectra that optimize
 * the least-squares regression objective function
 *****************************************************************************/
#include "lassoSGD.h"
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

//class glassoSGD: public regression {
class glassoSGDMS1: public lassoSGD{
	public:
		// Variables that store group information
		int G; // number of groups
		int GT;
		double * gmagnitudes; // GxT Stores the magnitude of each group ||Bg||2 over each time point
		int * gbyn; // Stores which group each n belongs to 

		glassoSGDMS1( dbl_matrix * y, dbl_matrix * x, char const * outdir, double percent_sparsity, char const * groupfilename,  bool * e );

		// Initialization of group information
		void initializeGroups( char const * groupfilename );

		// Optimization
		double dgp;
		void update(int m, int t ); // Returns the value of greatest update
		void learn_batch();
		double batch(); // Returns the value of greatest update

		int nzix;
		int nzend;	
		int n;
		int Bix;
		double diff;
		int gix; // the index of the group in consideration
		double gmag;
		double sse;

		~glassoSGDMS1();

};


glassoSGDMS1::glassoSGDMS1(dbl_matrix * y, dbl_matrix * x, char const * outdir, double percent_sparsity, char const * groupfilename, bool * e = NULL ): lassoSGD( y, x, outdir, percent_sparsity, e ) {
	printf("Constructed lassoSGD, constructing glassoSGDMS1\n");
	initializeGroups( groupfilename );	
	printf("Done!\n");
}

// This is formatted such tter jug math riddle solutionhat each row corresponds to one of the
// features in X, and the third column corresponds to which group
// that feature belongs to, and the groups are consecutive integers
// starting at 0.
void glassoSGDMS1::initializeGroups( char const * groupfilename ){ 

	
	G = 0;
	gbyn = new int[N];

	ifstream file( groupfilename );
	string line;
	getline( file, line);
	istringstream iss(line);
	
	printf("N=%i\n",N);
	
	int g;
	int n=0;
	char meh[50];
	//G = 0;
	while(file.good()){
		getline( file, line );
		if( line.length() < 4 ){
			continue;
		}	
		istringstream s(line);
		s >> meh; 
		s >> meh; 
		s >> meh;
		g = atoi(meh);
		
		gbyn[n]=g;
		if( (g+1)>G ){
			G = g+1;	
		}
		n++;
	}
	file.close();

	printf("G=%i\n",G);	
	gmagnitudes = zeros( G, T );
	GT = G*T;
}



// This update assumes the gmagnitude and B were already in agreement
// And this will preserve the agreement.
void glassoSGDMS1::update(int m, int t ){ // Returns the value of greatest update
	//printf("Updating (%i,%i)\n", m, t );
	residU = get_resid( m, t );
	avg_resid += residU*residU;
	max_deltaU = 0;
	
	nzix = m*maxnz;
	nzend = nzix+maxnz;	
	n=nzx[nzix];

	//nzx is a M*maxnz matrix that, which each row m contains
	//the peptides that have a theoretical peak at m/z bin m
	while( n != -1 ){
		//printf("-n:%i\n",n);
		Bix = t*N+n;
		gix = gbyn[n]*T+t; // gbyn[n] is the row number, t is the column number
		gmag = gmagnitudes[gix];
		if( (eprior && E[Bix]) || !eprior ){	
//		if( true ){	
			gradient = 2*Xcsr[nzix]*residU;
			if( debias ){ 
				if( B[Bix] > 0.0 ){
					deltaU = min( gradient*epsilon, B[Bix] ); 
				}
			}
			else{
				if( gmag == 0 ){
					dgp = lambda;
				}
				else{
					dgp = lambda*B[Bix]/gmag;
				}
				deltaU = min( (gradient+dgp)*epsilon, B[Bix] ); //peaksperpeptide[n])*epsilon, B[Bix] );
			}
			//printf("N: %i, %f %f\n", Bix, deltaU, gradient );
			// Update the group magnitude using a different formulation
			gmagnitudes[gix] = max(0.0,sqrt( deltaU*(deltaU-2*B[Bix]) + gmag*gmag ));
			// Update B			
			B[Bix] = B[Bix] - deltaU;
			// Update the group magnitude
			//gmagnitude[gix] = sqrt( deltaU*(-B[Bix]) + gmag*gmag );

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

// Here, Y is MxT
// t and m switched places from lassoSGD.h
double glassoSGDMS1::batch(){ // Returns the value of greatest update
	update_epsilon();

        clock_t begin = clock();
        begin = clock();
	if( gradients == NULL ){
		gradients = zeros(N,T);
	}
	sse = get_resids_fast(); // Resids is MxT
	int max_delta = 0;
	int nzix;
	int nzend;	
	int n;
	double diff;
		
	max_deltaU = 0;	
	
	int residix = 0;
	// Compute the gradient with respect to the
	// least-squares penalty
	printf("batching\n");
	for( int m=0; m<M; m++ ){ // m indexes the rows of resids, and the rows of x
//		printf("m: %i\n", m);
		for( int t=0; t<T; t++ ){ // t indexes the columns resids, rows of B
//			printf("-m: %i", m);
			nzix = m*maxnz;
			nzend = nzix+maxnz;	
			n=nzx[nzix];	
			//nzx is a M*maxnz matrix that, which each row m contains
			//the peptides that have a theoretical peak at m/z bin m
			//For each peptide that is non-zero at a time point
			// Here, lambda will only be for the group lasso. There is no 
			// L1 regularization for X2 as a whole.
			// If precursor n is present at time t
			while( n != -1 ){
				Bix = t*N+n;
				gradients[Bix] += 2*Xcsr[nzix]*resids[residix];
				nzix++;
				if( nzix == nzend ){ break; }
				n=nzx[nzix];
			}
			residix++;
		}
	}
	// Compute the updates and then perform them
	Bix = 0;
	int nnz = 0;
	for( int t=0; t<T; t++ ){
		for( n=0; n<N; n++ ){
			gix = gbyn[n]*T+t; // gbyn[n] is the row number, t is the column number
			gmag = gmagnitudes[gix];
//			if( (eprior && E[Bix]) || !eprior ){	
			if( true ){
				gradient = gradients[Bix];
				if( debias ){ 
					if( B[Bix] > 0.0 ){
						deltaU = min( gradient*epsilon, B[Bix] ); 
					}
				}
				else{
					if( gmag == 0 ){
						dgp = lambda;
					}
					else{
						dgp = lambda*B[Bix]/gmag;
					}
					deltaU = min( (gradient+dgp)*epsilon, B[Bix] ); //peaksperpeptide[n])*epsilon, B[Bix] );
				}
				//printf("N: %i, %f %f\n", Bix, deltaU, gradient );
				// Update B	
				B[Bix] = B[Bix] - deltaU;
				if( B[Bix] != 0 ){
					nnz++;				
				}	
				if( deltaU < 0 ){ deltaU = -1*deltaU; }
				if( deltaU > max_deltaU ){
					max_deltaU = deltaU;
				}
			}
			Bix++;
		}
	}
	printf("Actual nnz for B: %i, size: %i\n", nnz, B.nrow*B.ncol);

	// Then reinitialize the lsgradients to 0
	// and update the group magnitudes
	// Zero out the group magnitudes
	for( Bix=0; Bix<GT; Bix++ ){
		gmagnitudes[Bix] = 0;
	}
	Bix = 0;	
	for( int t=0; t<T; t++ ){
		for( n=0; n<N; n++ ){
			gradients[Bix] = 0;
			gix = gbyn[n]*T+t;
			gmagnitudes[gix] += B[Bix]*B[Bix];	
			Bix++;		
		}
	}
	for( Bix=0; Bix<GT; Bix++ ){
		gmagnitudes[Bix] = pow(gmagnitudes[Bix],0.5);
	}
        printf("\n in %f seconds. Max delta: %f, nnz in B: %i\n",double(clock()-begin)/CLOCKS_PER_SEC,max_delta,B.nnz());
	itnum++;
	return(max_deltaU);
}

void glassoSGDMS1::learn_batch(){
	printf("Outdir: %s\n",outdir);
	printf("Magnitude of Y: %f\n", Y.norm2());
	
	itnum = 1;
	epsilon = 1;
	double max_delta;
	clock_t begin;
	int nnz;
	double ratio;
	double lasso_objective;
	
	for( int i=0; i<maxit; i++ ){
		max_delta = batch();
	        printf("Peak RSS: %i mb\n", getPeakRSS()/1000000 );
		if( i < 2 ){
			lasso_objective = total_objective();
			ssediff = oldsse - newsse;
			printf("SSE at %i batches: %f\n", i+1, sse );
			printf("norm1 at %i batches: %f\n", i+1, norm1(B.m, B.nrow, B.ncol) );
			printf("objective at %i batches: %f\n", i+1, lasso_objective );
		}
		else if( i % (pollperiod) == 0 ){
			oldsse = newsse;
			lasso_objective = total_objective();
			printf("SSE at %i batches: %f\n", i+1, sse );
			printf("norm1 at %i batches: %f\n", i+1, norm1(B.m, B.nrow, B.ncol) );
			printf("objective at %i batches: %f\n", i+1, lasso_objective );
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
	char filename[200];
	sprintf(filename, "%s/%i_bgd_B.txt", outdir, itnum );
	write_B( filename ){
	printf("Num epochs: %i\n", epochnum );	

};	
 
glassoSGDMS1::~glassoSGDMS1(){
//	delete [] gmagnitudes;









//	delete [] gbyn;
//	delete [] gbyn;
}



