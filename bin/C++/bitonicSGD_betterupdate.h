/**
 * \file bitonicSGD.h
 * AUTHOR: Alex Hu
 * CREATE DATE: March 5, 2015
 * PROJECT: DIA Regression
 * \brief Returns peptide abundances in DIA spectra that optimize
 * the bitonic regression objective function
 *****************************************************************************/
#include "alternateBW.h"
#include "scorePeptides.h"
#include "getRSS.cpp"
#include <ctime>
#include <stdio.h>
#include <stdlib.h>

class bitonicSGD {

public:

	int M; // Number of m/z bins
	int T; // Number of scans
	int N; // Number of candidate peptides
	int MT; // M*T
	int NT;
	double lambda; // L1 regularization parameter
	double percent_sparsity;
	int * peaksperpeptide;

	// Input
	double* Y; // Observed spectra, MxT
	double* resid;
	double* X; // Theoretical spectra, MxN
	double* W; // Elution peak weights, NxT
	double *P; // bitonic penalties
	char const * outdir;
	
	// Output
	double* B; // Learned abundances, NxT
	
	// Initialization
	bitonicSGD(double y[], double x[], double w[], int m, int t, int n, char const * outdir, double percent_sparsity);
	void reset();
	int xnnz;
	void sparsifyXandY();	
	int * nzx;// M x maxnz matrix that lists the columns of X where row m is nonzero. Each list ends at -1.	
	double* Xcsr; // Compressed row matrix whose non-zero locations are in nzx.
	int maxnz;

	void estimate_lambda_sample( int nsamples, double percent_sparsity );	
	void estimate_lambda( double percent_sparsity );
	double* up; // Bitonic penalty for (Bn,t - Bn,t+1)+
	double* down; // Bitonic penalty for (Bn,t - Bn,t+1)-
	double epsilon;
	void compute_bitonic_weights();
	bool debias;  

	// Computing the objective
	double ls();
	double bitonic_objective();
	double bitonic_penalties();
	double total_objective();
	int itnum;
	int maxepochs;
	int epochnum;	

	// Optimization
	double get_bipenalty(int Bix);
	void update(int m, int t, double epsilon); 
	int* randsequence; // This is for sampling Y with replacement
	double epoch();
	void update_epsilon();
	void learn();
	void learn_bitonic();

	// Variables for the function update( int m, int t, double epsilon );
	double residU;
	double bipenaltyU;
	double deltaU;
	double max_deltaU; // Contains the maximum delta encountered for any Bn,t in the update

	int nzixU;
	int nzendU;
	int nU;
	int BixU;

	// Output stuff
	
	~bitonicSGD();
	
};

void bitonicSGD::reset(){
	for( int i=0; i<N*(T-1); i++){
		up[i] = 0;
		down[i]=0;
	}	
	compute_bitonic_weights();
	for( int i=0; i<N*T; i++ ){
		//B[i] = 0;
	}
	itnum=1;
	epochnum=0;
};

bitonicSGD get_bitonicSGD(char const * yfile, char const * xfile, char const * wfile, int maxT, char const * outdir, double percent_sparsity){
	
	int m;
	int t;
	int n;
	
	printf("Reading Y: %s\nmaxT=%i\n", yfile,maxT); 
	clock_t begin = clock();
	double * y = read_rowmajor_dblmatrix( yfile, &m, &t, maxT );
	printf("%ix%i\n", m, t );
	//print_matrix( y, m, t);
	printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000); 
	
	printf("Reading X: %s\n", xfile); 
	begin = clock();
	double * x = read_rowmajor_dblmatrix( xfile, &m, &n, -1 );
	printf("%ix%i\n", m, n );
	printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000); 
	//print_matrix( x, m, n);
	
	printf("Reading W: %s\n", wfile); 
	begin = clock();
	//double * w = read_dblmatrix( wfile, &n, &t, maxT );
	double * w = zeros( n, t );
	printf("%ix%i\n", n, t );
	printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000); 
	//print_matrix( w, n, t);
	
	return( bitonicSGD( y, x, w, m, t, n, outdir, percent_sparsity) );
};

// Initialization
bitonicSGD::bitonicSGD( double* y, double* x, double* w, int m, int t, int n, char const * od, double psparsity ){

	Y = y; 
	X = x;
	W = w;
	outdir = od;

	M = m;
	T = t;
	N = n;
	percent_sparsity = psparsity;
	maxnz=0;
	maxepochs=5;
	
	if( SQRT ){
		sqrt_mat(Y,M,T);
		sqrt_mat(X,M,N);
	}
	
	printf("Sparsifying M=%i to...",M); fflush(stdout);
	clock_t begin = clock();
	peaksperpeptide = zeros_int(N,1);
	sparsifyXandY();
	MT = M*T;	
	NT = N*T;	
	
	up = zeros(N,T-1);
	down = zeros(N,T-1);


	B = zeros(T,N); // T rows, N columns
	P = zeros(N,T);
	printf("M: %i, T: %i, N: %i,  MT: %i\n", M, T, N, MT );
	epsilon = 1;
	itnum = 1;
	epochnum = 0;
	debias = false;

	srand ( unsigned (0 )  );
	compute_bitonic_weights();
	randsequence = new int[MT];
	for( int i=0; i<MT; i++ ){ randsequence[i] = i; } 
	lambda = 0.0;
	printf("Estimating lambda at sparsity %f...\n", percent_sparsity);
	/*
	for( int i=0; i<1; i++ ){
		estimate_lambda_sample(0.1*MT, percent_sparsity );	
		printf("%i: %f\n",i,lambda);
		if( lambda > 0.0 ){
			printf("Lambda is actually non-zero.\n");
		}
		else{ printf("Lambda is actually zero :(\n"); } 
	}*/
	//estimate_lambda( percent_sparsity );
	estimate_lambda( percent_sparsity );	
};

// Looks for m/z bins that are unpopulated by X and removes them from
// both X and Y
void bitonicSGD::sparsifyXandY(){
	

	clock_t begin = clock();

	// Identify the m's to keep
	bool occupied[M] = {false};
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
	double *newX = zeros(newM,N);
	double *newY = zeros(newM,T);
	nzx = zeros_int(newM,maxnz); for(int i=0; i<newM*maxnz; i++ ){ nzx[i] = -1; }
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
				newX[newm*N+n] = X[m*N+n];
				nzx[nzix] = n;
				Xcsr[nzix] = X[m*N+n];
				nzix++;
			}
		}
		// Copy Ys
		for( int t=0; t<T; t++ ){
			if( Y[m*T+t] == 0 ){ continue; }
			newY[newm*T+t] = Y[m*T+t];
		}
		newm++;
	}
	if( DELETE ){
		delete [] X;
		delete [] Y;
	}
	X = newX;
	Y = newY;
	M = newM;
	MT = M*T;
	resid = zeros(M,T); for( int i=0; i<MT; i++ ){ if( Y[i] > 0 ){ resid[i] = -Y[i]; } }	
	printf("maxnz=%i\n",maxnz);
	printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000); 
};


void bitonicSGD::estimate_lambda( double percent_sparsity ){
	if( xnnz == 0 ){
		printf("Estimating lambda error: not sparsified yet.\n");
		lambda = 0;
		return;
	}
	clock_t begin = clock();
	
	double * half_derivs0 = Tmultiply( X, M, N, Y, M, T, 1 );
	lambda = 2*percentile( half_derivs0, N, T, percent_sparsity );
	printf("lambda for %f%: %f\n", percent_sparsity, lambda );
	printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000); 
	delete [] half_derivs0;
}

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

double bitonicSGD::ls(){
//	printf("ls()\n");
	double * model = multiplyT( X, M, N, B, T, N, 1.0 );
	double sum = 0;
	double r = 0;
//	clock_t begin = clock();
	for( int i=0; i<MT; i++ ){
		r = model[i]-Y[i];
		sum += r*r;
	}
	delete [] model;
//	printf("least-squares computation: %f (s)\n", double(clock()-begin)/CLOCKS_PER_SEC); 
	return(sum);
};

double bitonicSGD::bitonic_objective(){
//	printf("bitonic_objective()\n");
//	This B has dimensions TxN
	double penalty=0;
	double diff;
	int end = (T-1)*N;
	for( int n=0; n<N; n++ ){ // For each peptide n
		for(int i=n; i<end; i++ ){ // For each time for peptide n
			diff = B[i]-B[i+N];// B[i+N] is the subsequent time point
			if( diff > 0 ){
				penalty = penalty + diff*up[(i/N)+(i%N)]; // THIS IS PROBABLY WRONG
			}
			else if( diff < 0 ){
				penalty = penalty - diff*down[(i/N)+(i%N)]; // THIS IS PROBABLY WRONG
			}
		}
	}
	return(penalty);
};
/*
// This is a NxT matrix of penalties associated with each
// peptide/possible elution peak pair
double * bitonicSGD::bitonic_penalties(){
//	printf("bitonic_objective()\n");
	double * penalties = zeros( N, T);
	double diff;
	for(int i=0; i<(T*N-N); i++ ){
		diff = B[i]-B[i+N];
		if( diff > 0 ){
			penalty = penalty + diff*up[(i/N)+(i%N)];
		}
		else if( diff < 0 ){
			penalty = penalty - diff*down[(i/N)+(i%N)];
		}
	}
	return(penalties);
};*/

double bitonicSGD::total_objective(){
	double objective = ls() + bitonic_objective() + norm1(B,T,N)*lambda;
	return( objective );
};
	
void bitonicSGD::update(int m, int t, double epsilon){ // Returns the value of greatest update
	int rix = m*T+t;
	residU = resid[rix];
	max_deltaU = 0;
	
	nzixU = m*maxnz;
	int nzend = nzixU+maxnz;	
	int n=nzx[nzixU];
	int Bix;
	double diff;
	
	//nzx is a M*maxnz matrix that, which each row m contains
	//the peptides that have a theoretical peak at m/z bin m
	while( n != -1 ){
		Bix = t*N+n;
//		printf("C\n");
		bipenaltyU = get_bipenalty(Bix);
		
		if( debias && B[Bix] > 0.0 ){ 
			deltaU = min( (2*Xcsr[nzixU]*residU + bipenaltyU)*epsilon, B[Bix] ); 
		}
		else{
			deltaU = min( (2*Xcsr[nzixU]*residU + lambda/peaksperpeptide[n] + bipenaltyU)*epsilon, B[Bix] );
		}
		
		//// This is to make sure it doesn't fluctuate up and down too much.
		// If the update crosses a bitonic boundary, it'll force the update to just stay
		// at the boundary.
		if( t > 0 ){
			diff = B[Bix]-B[Bix-T];
			// If weight is greater than previous weight and the delta makes it less than the previous weight
			// Or, if weight is less than previous weight and the delta makes it greater than the previous weiht
			if( (diff > 0 && deltaU > diff) || (  diff < 0 && deltaU < diff ) ){
				deltaU = diff;	
			}
		}
		if( t < (T-1) ){
			diff = B[Bix]-B[Bix+T];
			// If weight is greater than the succeeding weight and the delta makes it less than the previous weight
			// Or, if weight is less than the succeeding weight and the delta makes it greater than the previous weiht
			if( (diff > 0 && deltaU > diff) || (  diff < 0 && deltaU < diff ) ){
				deltaU = diff;	
			}
		}
		////////

		B[Bix] = B[Bix] - deltaU;
		resid[rix] -= deltaU*Xcsr[nzixU];	
		if( deltaU < 0 ){ deltaU = -1*deltaU; }
		if( deltaU > max_deltaU ){
			max_deltaU = deltaU;
		}
		nzixU++;
		if( nzixU == nzend ){ break; }
		n=nzx[nzixU];
	}
	itnum++;
	//update_epsilon();
};

double bitonicSGD::epoch(){
        //printf("epoch %i, e=%f...", epsilon, epochnum ); fflush(stdout);
        double max_delta=0;
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
                update( m, t, epsilon );
                if( max_deltaU > max_delta ){
                        max_delta = max_deltaU;
                }
        }
        printf(" in %f seconds. Max delta: %f, nnz in B: %i\n",double(clock()-begin)/CLOCKS_PER_SEC,max_delta,nnz(B,T,N));
        epochnum++;
        return(max_delta);
};

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

void bitonicSGD::update_epsilon(){
	double maxbitonicpen = max(up,N,T-1); 	
	epsilon = 1.0/( (epochnum) + 2 + maxbitonicpen*2 );
};

void bitonicSGD::learn(){
	printf("Outdir: %s\n",outdir);
	itnum = 1;
	epochnum = 0;
	epsilon = 1;
	double max_delta;
	//printf("%f\n",total_objective());
	double * Bt;
	double * penalties;
	double * exW = zeros(N,T);
	clock_t begin;
	for( int i=0; i<maxepochs; i++ ){
		max_delta = epoch();
	//	printf("%f\n",total_objective());
		char filename[200];
		sprintf(filename,"%s/%i_B.txt",outdir,i);
		Bt = transpose(B,T,N);
		//hardthresh( Bt, N, T, max_delta/2.0 ); 
		double true_sparsity = write_csr_matrix( Bt, N, T, filename );
		//print_matrix( Bt, N, T, "B");	
		if( i > 0 ){
			printf("Scoring..\n");
			penalties = compute_unweighted_penalties(Bt, N, T);
			compute_Ws( exW, N, T, penalties, 0, 1 );
			sprintf(filename,"%s/%i_%s",outdir,i,"scores.txt");
			
			begin = clock();
			double ** scores = corrscore_W( Y, X, Bt, exW, M, T, N );
			printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000); 
			write_scores( scores, 5, N, filename );
			deldoublestarstar( scores, 5);
			delete [] penalties;		
			//double ** scores = auc_score( B, N, T );
			//write_scores( scores, 3, N, filename );
		}
	        printf("Peak RSS: %i mb\n", getPeakRSS()/1000000 );
		delete [] Bt;

	}
	printf("Num epochs: %i\n", epochnum );
	//print_matrix(transpose(B,T,N),N,T);
	
};	

void bitonicSGD::learn_bitonic(){
	printf("Outdir: %s\n",outdir);
	itnum = 1;
	epochnum = 0;
	epsilon = 1;
	double max_delta;
	//printf("%f\n",total_objective());
	double * Bt = zeros(N,T);
	double * penalties;
	double * exW = zeros(N,T);
	clock_t begin;
	char filename[200];

	//double percentiles[4] = { 0, 100.0/3.0, 200.0/3.0, 100.0 };
	double percentiles[5] = { 0, 100/4.0, 200.0/4.0, 300.0/4.0, 100.0 };
	int npercentiles = 5;
	//Do an epoch without bitonic regression
	for( int i=0; i<npercentiles; i++ ){
		printf("Epoch #%i\n",i);
		// Compute new Ws and do the learning
		if( percentiles[i] > 0 ){
			compute_Ws_original( percentiles[i], Bt, W, N, T, P );
			compute_bitonic_weights();
		}
		max_delta = epoch();
		
		// Write sresults
		transpose(B,T,N, Bt);
		sprintf(filename,"%s/%i_B.txt",outdir,i);
		double true_sparsity = write_csr_matrix( Bt, N, T, filename );
		sprintf(filename,"%s/%i_W.txt",outdir,i);
		write_sparse_matrix( W, N, T, filename );	
		sprintf(filename,"%s/%i_P.txt",outdir,i);
		write_matrix( P, N, T, filename );	
		sprintf(filename,"%s/%i_up.txt",outdir,i);
		write_matrix( up, N, T-1, filename );	
		sprintf(filename,"%s/%i_down.txt",outdir,i);
		write_matrix( down, N, T-1, filename );	

		// Scoring
		printf("Scoring..\n");
		penalties = compute_unweighted_penalties(Bt, N, T);
		compute_Ws( exW, N, T, penalties, 0, 1 );
		sprintf(filename,"%s/%i_%s",outdir,i,"scores.txt");

		begin = clock();
		double ** scores = corrscore_W( Y, X, Bt, exW, M, T, N );
		printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000); 
		write_scores( scores, 5, N, filename );
		deldoublestarstar( scores, 5);
		delete [] penalties;		
	        printf("Peak RSS: %i mb\n", getPeakRSS()/1000000 );
	}
	
};	

bitonicSGD::~bitonicSGD(){
	if( !up ){ delete [] up; }
	if( !down ){delete [] down; }
	if( !B ){delete [] B; }
	if( !Y ){delete [] Y; }
	if( !X ){delete [] X; }
	if( !W ){delete [] W; }
	if( !randsequence ){delete [] randsequence; }
}

