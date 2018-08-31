#include "io.h"

double * rowpens;

int compare2( const void *a, const void *b ){

	int x = *(int*)a;
	int y = *(int*)b;

	if( rowpens[x] > rowpens[y] ){
		return 1;
	}
	if( rowpens[x] < rowpens[y] ){
		return -1;
	} 
	return(0);
}

// Returns a matrix of omegas of size Nxnumnonzero that
// contains the maximum omega for which there are numnonzero-1
// nonzero omegas
double * compute_min_omegas( double * penalties, int N, int T, double sigma, int numnonzero  ){
	int * indices = zeros_int(1,T);
	for( int i=0; i<T; i++ ){ indices[i] = i; } 
	double * omegas = zeros(N,numnonzero);
	double psum;
	double lambda;
	for( int n=0; n<N; n++ ){
		rowpens = &penalties[n*T];
		qsort( &indices[0], T, sizeof(int), compare2 );
		psum = rowpens[ indices[0] ];
		for( int nnz=2; nnz<=numnonzero; nnz++ ){
			psum += rowpens[indices[nnz-1]];
			omegas[n*numnonzero+nnz-1] = (rowpens[indices[nnz-1]]*nnz-psum)/2.0/sigma; 		
		}	
	}
	delete [] indices;
	return omegas;	
}

// Relies on rowpens, which is the bitonic(B) penalties
// W is NxT
double optimize_row( double*W, int N,int T,double*penalties, int * indices, int n, double omega, double sigma ){

	rowpens = &penalties[n*T];
	for( int i=0; i<T; i++ ){ indices[i] = i; W[n*T+i] = 0.0; }
	qsort( &indices[0], T, sizeof(int), compare2 );
	double psum = rowpens[indices[0]];
	W[n*T+indices[0]] = sigma;

	double c = sigma/T;
	double unipen = c*c; //Penalty for a Wn,t=0
	double current_pen;
	double minpen = (unipen*(T-1) + (sigma-c)*(sigma-c))*omega + rowpens[indices[0]]; //Penalty if there is one nonzero value 1
	double lambda;
	double w;
	int ix;
	
	//printf("nnz=1 Penalty: %f\n", minpen);
	for( int nnz=2; nnz<(T+1); nnz++ ){
		//printf("Computing for nnz %i\n", nnz );
		current_pen = unipen*(T-nnz)*omega;
		psum += rowpens[indices[nnz-1]];
		//Compute lambda
		lambda = (2*omega*(nnz*c-sigma) - psum)/nnz;
		// Compute new weights and total penalty
		for( int i=0; i<nnz; i++ ){
			w = c-(rowpens[indices[i]]+lambda)/2.0/omega;
			current_pen += (c-w)*(c-w)*omega + rowpens[indices[i]]*w;	
			if( w < 0.0 || w > sigma ){
		//		printf("Done!\n");
				return(minpen);
			}
		}
		//print_matrix( &W[n*T], 1, T, "W\n");
		//printf("Penalty: %f\n", current_pen);
		if( minpen > current_pen ){
			for( int i=0; i<nnz; i++ ){
				ix = indices[i];
				W[n*T+ix] = c-(rowpens[ix]+lambda)/2.0/omega;
			}
			minpen = current_pen;
		}
		//print_matrix( &W[n*T], 1, T, "W\n");
	}
	return( minpen );
}


void compute_Ws( double * W, int N, int T, double * penalties, double omega, double sigma ){

	int * indices = zeros_int(1,T);
	for( int t=0; t<T; t++ ){ indices[t] = t; }
	for(int n=0; n<N; n++ ){
		optimize_row( W, N,T,penalties, indices, n, omega, sigma );			
	}
	delete [] indices;
}

void compute_Ws_v1( double * W, int N, int T, double * penalties, double omega ){

	//Reset Ws to be full of 0s
	int size = N*T;
	for( int i=0; i<size; i++ ){ W[i] = 0; }
	// Finds the scan t for each peptide that has the smallest penalty
	// And computes its corresponding Wn,t
	int mint;
	double w;
	for( int n=0; n<N; n++ ){
		mint = argmin(&penalties[n*T],1,T); 
		w = omega-penalties[n*T+mint]/2.0;
		if( w > 0 ){ W[n*T+mint] = w; }
	}
}


// Here B is NxT!
// Output is also N,T!
double * compute_unweighted_penalties( double * B, int N, int T ){

	double * penalties = zeros(N,T);
	int Bix;
	int Bend;
	int Begin;
	double diff;
	for( int n=0; n<N; n++ ){
		Begin = n*T;
		Bix = Begin;
		Bend = Bix+T-1;
		while( Bix < Bend ){
			diff = B[Bix]-B[Bix+1];
			if( diff > 0 ){
				for(int pix=(Bix+1); pix<=Bend; pix++){
					penalties[pix] += diff;
				}
			}
			else if( diff < 0 ){
				for(int pix=Begin; pix<=Bix; pix++){
					penalties[pix] -= diff;
				}
			}
			Bix++;	
		}	
	}
	return(penalties);
}


// Here B is NxT!
// Output is also N,T!
void compute_unweighted_penalties( double * B, int N, int T, double * penalties ){

	int NT = N*T;
	for( int i=0; i<NT; i++ ){ penalties[i] = 0; }
	int Bix;
	int Bend;
	int Begin;
	double diff;
	for( int n=0; n<N; n++ ){
		Begin = n*T;
		Bix = Begin;
		Bend = Bix+T-1;
		while( Bix < Bend ){
			diff = B[Bix]-B[Bix+1];
			if( diff > 0 ){
				for(int pix=(Bix+1); pix<=Bend; pix++){
					penalties[pix] += diff;
				}
			}
			else if( diff < 0 ){
				for(int pix=Begin; pix<=Bix; pix++){
					penalties[pix] -= diff;
				}
			}
			Bix++;	
		}	
	}
}


/* 
//This one is for a truncated normal distribution of weights
double * initialize_Ws( int N, int T ){

	double * Ws = zeros( N, T );
	int size = N*T;
	for( int i=0; i<size; i++ ){
		while( Ws[i] > 0.0 ){
			
		}		
	}	
	return( Ws );
}*/

//This one is for randomly cutting 1 into T pieces uniformly
double * initialize_Ws( int N, int T, double sigma ){
	int start;
	double * Ws = zeros(N,T);
	double * samples = zeros(1,T-1);
	srand(time(NULL));
	for( int n=0; n<N; n++ ){
		for( int i=0; i<(T-1); i++ ){
			samples[i] = (float)rand()/RAND_MAX;	
		}
		sort(samples,samples+(T-1));
		start = n*T;
		Ws[start] = samples[0]*sigma;
		for( int i=0; i<(T-2); i++ ){
			Ws[start+i+1] = (samples[i+1]-samples[i])*sigma;
		}
		Ws[start+T-1] = (1-samples[T-2])*sigma;
	}
	delete [] samples;
	return(Ws);
}

double compute_Ws_original( double percent_annotated_peaks, double * B, double * W, int N, int T, double * P ){

        int size = N*T;
        for( int i=0; i<size; i++ ){ W[i] = 0; }

        compute_unweighted_penalties( B, N, T, P );
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

	int NT = N*T;
        delete [] mintperpep;
	return( omega );
}



