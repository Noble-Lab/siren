/**
 * \file regression.h
 * AUTHOR: Alex Hu
 * CREATE DATE: March 5, 2015
 * PROJECT: DIA Regression
 * \brief Returns peptide abundances in DIA spectra that optimize
 * the least-squares regression objective function
 *****************************************************************************/
#include "lassoSGD.h"
//#include "lassoADAM.h"
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

class bitonicSGD: public lassoSGD{
public:
	// Initialization/destruction
	bitonicSGD(): lassoSGD() {};
	bitonicSGD(dbl_matrix * y, dbl_matrix * x, char const * outdir, double percent_sparsity);
	~bitonicSGD();

	// Variables for bitonic optimization
	// Up and down
	dbl_matrix * Bt;
	double * W; // Stores elution profile weights
	double * P; // Stores elution profile penalties for optimization of B
	double * up;
	double * down;
	void compute_bitonic_weights();

	// Computing the objective
	double bitonic_objective();
	double bitonic_penalties();
	//double total_objective();

	// Optimization of B
	void update(int m, int t);
	double get_bipenalty( int Bix );
	void learn(int npercentiles);

	// Optimization of W
	double compute_Ws(double percent_annotated_peaks ); // This one finds just one non-zero W per peptide.
	void compute_unweighted_penalties();

	// Variables for the function update
	double bipenaltyU;

};

bitonicSGD::bitonicSGD(dbl_matrix * y, dbl_matrix * x, char const * outdir, double percent_sparsity): lassoSGD(y, x, outdir, percent_sparsity){

        Bt = new dbl_matrix( B.ncol, B.nrow );
	W = zeros( N, T ); // Stores elution profile weights
	P = zeros( N, T ); // Stores elution profile penalties for optimization of B
	up = zeros( N, T );
	down = zeros (N, T );

};

// Here B is NxT!
// Output is also N,T!
// Filles P with the penalties assiated with 
// B
void bitonicSGD::compute_unweighted_penalties(){

	// Zero out penalties
	for( int i=0; i<NT; i++ ){ P[i] = 0; }
	int Bix;
	int Bend;
	int Begin;
	double diff;
	for( int n=0; n<N; n++ ){
		Begin = n*T;
		Bix = Begin;
		Bend = Bix+T-1;
		while( Bix < Bend ){
			diff = (*Bt)[Bix]-(*Bt)[Bix+1];
			if( diff > 0 ){
				for(int pix=(Bix+1); pix<=Bend; pix++){
					P[pix] += diff;
				}
			}
			else if( diff < 0 ){
				for(int pix=Begin; pix<=Bix; pix++){
					P[pix] -= diff;
				}
			}
			Bix++;	
		}	
	}
}

double bitonicSGD::get_bipenalty(int Bix){
	int t = Bix/N;
	int n = Bix-t*N;
	int penix = n*(T-1);
	double penalty = 0;
	if( t > 0 ){
		if(B[Bix] > B[Bix-N]){ // B[n,t] > B[n,t-1]
			penalty += down[penix+t-1];
		}
		else if(B[Bix] < B[Bix-N]){// B[n,t] < B[n,t-1]
			penalty -= up[penix+t-1];
		}
	}
	if( t < T-1 ){
		if(B[Bix] > B[Bix+N]){ // B[n,t] > B[n,t+1]
			penalty += up[penix+t];
		}
		else if(B[Bix] < B[Bix+N]){// B[n,t] < B[n,t+1]
			penalty -= down[penix+t];
		}
	}
	return(penalty);
};	

void bitonicSGD::update(int m, int t ){ // Returns the value of greatest update
	residU = get_resid( m, t );
	max_deltaU = 0;
	
	int nzix = m*maxnz;
	int nzend = nzix+maxnz;	
	int n=nzx[nzix];
	int Bix;
	double diff;

	//nzx is a M*maxnz matrix that, which each row m contains
	//the peptides that have a theoretical peak at m/z bin m
	while( n != -1 ){
		Bix = t*N+n;
		bipenaltyU = get_bipenalty( Bix );
		if( debias && B[Bix] > 0.0 ){ 
			deltaU = min( (2*Xcsr[nzix]*residU + bipenaltyU )*epsilon, B[Bix] ); 
		}
		else{
			deltaU = min( (2*Xcsr[nzix]*residU + lambda/peaksperpeptide[n] + bipenaltyU )*epsilon, B[Bix] );
		}

		//// This is to make sure it doesn't fluctuate up and down too much.
		// If the update crosses a bitonic boundary, it'll force the update to just stay
		// at the boundary.
		if( t > 0 ){
			diff = B[Bix]-B[Bix-N];
			// If weight is greater than previous weight and the delta makes it less than the previous weight
			// Or, if weight is less than previous weight and the delta makes it greater than the previous weiht
			if( (diff > 0 && deltaU > diff) || (  diff < 0 && deltaU < diff ) ){
				deltaU = diff;	
			}
		}
		if( t < (T-1) ){
			diff = B[Bix]-B[Bix+N];
			// If weight is greater than the succeeding weight and the delta makes it less than the previous weight
			// Or, if weight is less than the succeeding weight and the delta makes it greater than the previous weiht
			if( (diff > 0 && deltaU > diff) || (  diff < 0 && deltaU < diff ) ){
				deltaU = diff;	
			}
		}

		B[Bix] = B[Bix] - deltaU;
		if( deltaU < 0 ){ deltaU = -1*deltaU; }
		if( deltaU > max_deltaU ){
			max_deltaU = deltaU;
		}

		nzix++;
		if( nzix == nzend ){ break; }
		n=nzx[nzix];
	}
	itnum++;
	//update_epsilon();
};

void bitonicSGD::compute_bitonic_weights(){
        if( ! up ){
                up = zeros(N,T-1);
        }
        if( ! down ){
                down = zeros(N,T-1);
        }
        zero( up, N, T-1 );
        zero( down, N, T-1 );
        int offset;
        int Woffset;

        for(int n=0; n<N; n++){
                offset = n*(T-1);
                Woffset = n*T;
                down[offset] = W[Woffset];
                for(int t=1; t<(T-1); t++ ){
                        down[offset+t] = down[offset+t-1] + W[Woffset+t];
                }
                up[offset+T-2] = W[Woffset+T-1];
                for(int t=(T-3); t > -1; t-- ){
                        up[offset+t] = up[offset+t+1] + W[Woffset+t+1];
                }
        }
};

// This is the historical
// version of compute Ws that necessarily finds the optimal, single
// non-zero W and finds its optimal value.
double bitonicSGD::compute_Ws( double percent_annotated_peaks ){

	// Zero out old Ws
        for( int i=0; i<NT; i++ ){ if( W[i] != 0 ) { W[i] = 0; } }
        compute_unweighted_penalties();
	int * mintperpep = mincolperrow( P, N, T );
	double * minPperpep = zeros(N,1);
	for( int n=0; n<N; n++ ){ minPperpep[n] = P[n*T+mintperpep[n]]; }
	double omega = percentile( minPperpep, N, 1, percent_annotated_peaks )/2.0;
        double wnew;
        for( int n=0; n<N; n++ ){
                wnew = omega-P[n*T+mintperpep[n]]/2.0;
                if( wnew > 0 ){ W[n*T+mintperpep[n]] = wnew; }
        }
	printf("Percent with peaks: %f\t, omega: %f\n", percent_annotated_peaks, omega );

        delete [] mintperpep;
	delete [] minPperpep;
	return( omega );
}

void bitonicSGD::learn( int npercentiles ){
	printf("Outdir: %s\n",outdir);
	itnum = 1;
	epochnum = 0;
	epsilon = 1;
	double max_delta;

        double * exW = zeros(N,T);
        double ** scores;
        clock_t begin;
        char filename[200];

	if( npercentiles <= 2 ){
		printf("Minimum # of percentiles is 2!\n");
		npercentiles = 2;
	}
	double * percentiles = new double[npercentiles];
	for( int i=0; i<npercentiles; i++ ){
		percentiles[i] = i*100/(npercentiles-1);
	}

	//double percentiles[4] = { 0, 100.0/3.0, 200.0/3.0, 100.0 };
	//Do an epoch without bitonic regression
	for( int i=0; i<npercentiles; i++ ){
		printf("Epoch #%i\n",i);
		// Compute new Ws and do the learning
		if( percentiles[i] > 0 ){
			compute_Ws( percentiles[i] );
		}
		printf("Percentile with peaks: %f\n", percentiles[i]) ;
		max_delta = epoch();
	
		// Write sresults
		transpose(B.m,T,N, Bt->m);
		sprintf(filename,"%s/%i_B.txt",outdir,i);
		double true_sparsity = write_csr_matrix( Bt->m, N, T, filename );
		sprintf(filename,"%s/%i_W.txt",outdir,i);
		write_sparse_matrix( W, N, T, filename );	
		sprintf(filename,"%s/%i_P.txt",outdir,i);
		write_sparse_matrix( P, N, T, filename );	

		// Scoring
		printf("Scoring..\n");
		// FIXME!
		compute_Ws( 100 );
		sprintf(filename,"%s/%i_%s",outdir,i,"scores.txt");

		begin = clock();
		scores = corrscore_W( Y.m, X.m, Bt->m, W, M, T, N );
		printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000); 
		write_scores( scores, 5, N, filename );
		deldoublestarstar( scores, 5);
	        printf("Peak RSS: %i mb\n", getPeakRSS()/1000000 );
	}
	
	delete exW;

};	

bitonicSGD::~bitonicSGD(){

	delete Bt;

//STUFF??
}

