using namespace std;

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "alternateBW.h"
#include "deconvolution.h"
#include "bitonicSGD.h"
#include "lassoBGD.h" 
#include "lassoSGD.h" 
#include "io.h" 

extern int opterr;


int main( int argc, char *argv[] ){

	char * x1file = NULL;
	char * x2file = NULL;
	char * y1file = NULL;
	char * y2file = NULL;
	char * Bfile = NULL;
	char * outdir = NULL;
	char * efile = NULL;
	char * algorithm = NULL;
	bool squareroot = true;

	double psparsity = 99.5; // Default value
	double ms1weight = 1.0; // Default value
	int maxits = 100; // Default value
	bool bitonic = false; // Default value
	bool deconvolution = false; // Deconvolution skips the scoring
	bool rescaleX = true; // 
	bool compute_sse = false;
	bool rescore = true;
	double stepsize = -1.0;
	bool score_only = false;

	// Parse the command line
	int c=0;
	opterr = 0;
	while( 1 ){	
		// Here are all the options
		static struct option long_options[] = {
			{"x1", optional_argument, NULL, 'p' }, // MS1 theoretical matrix file
			{"x2", optional_argument, NULL, 'f' }, // MS2 theoretical matrix file
			{"y1", optional_argument, NULL, 'q' }, // MS1 observed matrix file
			{"y2", optional_argument, NULL, 'g' }, // MS2 observed matrix file
			{"B", optional_argument, NULL, 'B' }, // Initial solution for peptide abundance
			{"stepsize", optional_argument, NULL, 'z' }, // stepsize for gradient descent
			{"efile", optional_argument, NULL, 'e' }, // Boolean elution prior file
			{"ms1weight", optional_argument, NULL, 'w' }, // float weight for ms1
			{"outdir", required_argument, NULL, 'o' }, // Outdir
			{"algorithm", optional_argument, NULL, 'a' }, // BGD or SGD
			{"psparsity", required_argument, NULL, 's' }, // percent sparsity of solution
			// For lassoBGD and lassoSGD
			{"maxits", optional_argument, NULL, 'i' }, // maximum number of BGD iterations or SGD epochs	
			// For bitonic, sets a flag on whether to do bitonic
			{"bitonic", no_argument, NULL, 'b' }, // MS2 observed matrix file
			{"unsquared", no_argument, NULL, 'u' }, // whether or not to square root input data
			{"rescore_bfile", no_argument, NULL, 'r' }, // NOT SURE WHAT THIS SHOULD DO RN
			{"score_only", no_argument, NULL, 'R' }, // NOT SURE WHAT THIS SHOULD DO RN
			{"noscale", no_argument, NULL, 'n' }, // whether or not to square root input data
			{"compute_sse", no_argument, NULL, 'c' }, // compute the sse
			{"deconvolution", no_argument, NULL, 'd' }, // MS2 observed matrix file
			{NULL, 0, NULL, 0},
		};
		int option_index = 0;
		c = getopt_long( argc, argv, "p:f:q:g:B:e:w:o:a:s:i:z:burncdR", long_options, &option_index );
		if( c == -1 ){
			break;
		}
		printf("option_index: %i\n", option_index );
		switch(c){
			case 'B':
				Bfile = optarg;
				printf("Bfile: %s\n", Bfile );
				break;
			case 'p':
				x1file = optarg;
				printf("x1file: %s\n", x1file );
				break;
			case 'w':
				ms1weight = atof(optarg);
				printf("ms1weight: %f\n", ms1weight );
				break;
			case 'u':
				squareroot = false;
				break;
			case 'r':
				rescore = true;
				break;
			case 'R':
				score_only = true;
				break;
			case 'n':
				rescaleX = false;
				break;
			case 'f':
				x2file = optarg;
				printf("x2file: %s\n", x2file );
				break;
			case 'q':
				y1file = optarg;
				printf("y1file: %s\n", y1file );
				break;
			case 'g':
				y2file = optarg;
				printf("y2file: %s\n", y2file );
				break;
			case 'z':
				stepsize = atof(optarg);
				printf("stepsize: %f\n", stepsize );
				break;
			case 'e':
				efile = optarg;
				printf("efile: %s\n", efile );
				break;
			case 'o':	
				outdir = optarg;
				printf("outdir: %s\n", outdir );
				break;
			case 'a':
				algorithm = optarg;
				printf("algorithm: %s\n", algorithm );
				break;
			case 's':
				psparsity = atof(optarg);
				printf("psparsity: %f\n", psparsity );
				break;
			case 'i':
				maxits = atoi(optarg);
				printf("maxits: %i\n", maxits );
				break;
			case 'b':
				bitonic = true;
				printf("bitonic: %i\n", bitonic );
				break;
			case 'c':
				compute_sse = true;
				printf("compute_sse: %i\n", compute_sse );
				break;
			case 'd':
				deconvolution = true;
				printf("deconvolution: %i\n", deconvolution );
				break;
			case 0:
				printf("x1file: %s\n", x1file );
				printf("Arg 0\n");
				break;
			default:
				abort();
				break;
		}
	}

	printf("Loading the data\n");
	if( !deconvolution ){
		// Load the data!
		vector<dbl_matrix*> xy;
		//printf("\n\n\nsquareroot is %i\n\n\n",squareroot);
		if( x1file != NULL && y1file != NULL ){
			// If both MS1 and MS2 data is included
			if( x2file != NULL && y2file != NULL ){			
				xy = get_xy_withms1( x1file, x2file, y1file, y2file, ms1weight, squareroot, rescaleX );
			}
			// If just MS1 data is included
			else{
				xy = get_xy( x1file, y1file, squareroot, rescaleX );
			}
		}
		// If just MS2 data is included
		else if( x2file != NULL && y2file != NULL  ){		
				xy = get_xy( x2file, y2file, squareroot, rescaleX );	
		}
		else{
			printf("Insufficient data\n");
			return(0); // FIXME Throw some kind of error!;
		}

		bool * E = NULL;
		int N;
		int T;
		int yT = xy[1]->ncol;

		if( efile != NULL ){
			E = read_rowmajor_bool( efile, &N, &T, yT ); 
			printf("E: (%i,%i)\n", N, T );
		//	E = falses( xy[0]->nrow, xy[1]->ncol );
		}

		// Choose the algorithm
		regression * r;

		//print_matrix( xy[0]->m, xy[0]->nrow, xy[0]->ncol, "X (b) matrix!\n");

		if( strcmp( algorithm, "sgd" ) == 0 || strcmp( algorithm, "both" ) == 0 ){	
			printf("Making SGD\n");	
			if( bitonic ){
				r = new bitonicSGD( xy[1], xy[0], outdir, psparsity );
			}
			else{
				r = new lassoSGD( xy[1], xy[0], outdir, psparsity, E );
			}
		}
		else if( strcmp( algorithm, "bgd" ) == 0 || strcmp( algorithm, "fista" ) == 0 ){
			printf("Making BGD\n");
			r = new lassoBGD( xy[1], xy[0], outdir, psparsity, E );
		}
		else{
			printf("Algorithm does not match 'sgd', 'bgd', 'fista', nor 'both'\n");
			return(0);
		}

		dbl_matrix * B;
		if( Bfile != NULL ){
			B = new dbl_matrix( Bfile ); 
			r->setB( B );
			delete B;
		}

		r->algorithm = algorithm;
		r->maxit = maxits;
		r->stepsize = stepsize;
		r->compute_sse = compute_sse;
		printf("\n\n\nIS SOMETHING HAPPENING?\n\n\n\n");

		// If B was given, don't try to learn it, just score it
		if( ! score_only ){
			r->learn();	
			r->write_B();
		}

		if( strcmp( algorithm, "both" ) == 0 ){	
			printf("\nMaking BGD\n");
			lassoSGD * rs = static_cast< lassoSGD * >(r);
			lassoBGD * rb = new lassoBGD( rs );

			if( ! score_only ){
				rb->learn();
				rb->write_B();
			}
			else{
				r->write_B();
			}
			r = rb;
		}

		// Scoring
		if( rescore ){
			printf("Scoring..\n");

			N = r->N;
			T = r->T;

			if( strcmp( algorithm, "sgd" ) == 0 ){		
				double * penalties;
				lassoSGD * rs = static_cast< lassoSGD * >(r);
				penalties  = compute_unweighted_penalties(rs->Bt, N, T);
				double * exW  = zeros(N,T);
		
				char filename[500];
				compute_Ws( exW, N, T, penalties, 0, 1 );
				sprintf(filename,"%s/final_scores.txt",outdir);
		
				//double ** scores = corrscore_W( r->Y.m, r->X.m, rs->Bt, exW, r->M, T, N );
				//double ** scores = corrscore_W_backgroundsubtracted( r->Y.m, r->X.m, rs->Bt, exW, r->M, T, N );
				//write_scores( scores, 6, N, filename );
				//deldoublestarstar( scores, 6);
			//double ** scores = corrscore_W_variable( r->Y.m, r->X.m, r->B.m, W, int M, int T, int N, int win, char* wintype, int K, int ms1size, int ms2size ){
				//write_scores( scores, 5, N, filename );
				//deldoublestarstar( scores, 5);
				delete [] penalties;		
			}
			else{
				double * penalties;
				penalties  = compute_unweighted_penalties(r->B.m, N, T);
				double * exW  = zeros(N,T);
	
				char filename[500];
				compute_Ws( exW, N, T, penalties, 0, 1 );
				sprintf(filename,"%s/final_scores.txt",outdir);
	
				double ** scores = corrscore_W_backgroundsubtracted( r->Y.m, r->X.m, r->B.m, exW, r->M, T, N );
				write_scores( scores, 6, N, filename );
				deldoublestarstar( scores, 6);
				//double ** scores = corrscore_W_variable( r->Y.m, r->X.m, r->B.m, W, int M, int T, int N, int win, char* wintype, int K, int ms1size, int ms2size ){
				//write_scores( scores, 5, N, filename );
				//deldoublestarstar( scores, 5);
				delete [] penalties;		
			}
		}
	}

	// Deconvolution
	else{


		// Requires x1, y1, and y2
		// Load the data!
		vector<dbl_matrix*> xy;
		//printf("\n\n\nsquareroot is %i\n\n\n",squareroot);
		if( x1file == NULL || y1file == NULL || y2file == NULL){
			printf("Insufficient data for deconvolution\n");
			return(0);
		}
		xy = get_xy( x1file, y1file, squareroot, rescaleX );

		bool * E = NULL;
		int N;
		int T;
		int yT = xy[1]->ncol;

		if( efile != NULL ){
			E = read_rowmajor_bool( efile, &N, &T, yT ); 
			printf("E: (%i,%i)\n", N, T );
		}

		// Choose the algorithm
		regression * r;
		// The MS1 elution profile step should have no L1 regularization
		if( strcmp( algorithm, "sgd" ) == 0 || strcmp( algorithm, "both" ) == 0 ){	
			printf("Making SGD\n");	
			if( bitonic ){
				r = new bitonicSGD( xy[1], xy[0], outdir, psparsity );
			}
			else{
				r = new lassoSGD( xy[1], xy[0], outdir, psparsity, E );
			}
		}
		else if( strcmp( algorithm, "bgd" ) == 0 || strcmp( algorithm, "fista" ) == 0 ){
			printf("Making BGD\n");
			r = new lassoBGD( xy[1], xy[0], outdir, psparsity, E );
		}
		else{
			printf("Algorithm does not match 'sgd', 'bgd', nor 'both'\n");
			return(0);
		}
		
		dbl_matrix * B;

		r->algorithm = algorithm;
		r->maxit = maxits;
		r->stepsize = stepsize;
		r->compute_sse = compute_sse;
		if( Bfile != NULL ){
			printf("Reading and setting B: %s\n", Bfile );
			B = new dbl_matrix( Bfile ); 
			r->setB( B );
			delete B;
		}
		// Only learn a B if a B wasn't given
		else{
			r->learn();	
		}

		if( strcmp( algorithm, "both" ) == 0 ){	
			printf("Making BGD\n");
			lassoSGD * rs = static_cast< lassoSGD * >(r);
			lassoBGD * rb = new lassoBGD( rs );
			r->stepsize = stepsize;
			// Only learn a B if a B wasn't given
			if( ! score_only ){
				rb->learn();
			}
			r = rb;
		}
		r->write_B("ms1_elutionB.txt");

		printf("Segmenting elution profiles\n");
		char filename[500];
		sprintf(filename,"%s/segment_annotations.txt",outdir);
		dbl_matrix *dB;
		if( strcmp( algorithm, "sgd" ) == 0 ){	
			lassoSGD * rs = static_cast< lassoSGD * >(r);
			// This bullshit is to get the correct form of B from lasso SGD,
			// whose B element is the transpose of what you want it to be.
			int T = rs->B.nrow;
			int N = rs->B.ncol;
			rs->B.m = rs->Bt;
			rs->B.nrow = N;
			rs->B.ncol = T;
		 	dB = segment_elution_profiles( &rs->B, filename );
		}
		else{
		 	dB = segment_elution_profiles( &r->B, filename );
		}
		

		sprintf(filename,"%s/ms1_elutionB_segmented.txt",outdir);
		dB->write_csr(filename);	
		col_magnorm( dB->m, dB->nrow, dB->ncol );

		//return(0);

		//delete r;
		printf("Reading y2 matrix\n");
		dbl_matrix *y2 = new dbl_matrix( y2file );
		y2->t();

		printf("Initializing new regression object\n");
		// Choose the algorithm
		if( strcmp( algorithm, "sgd" ) == 0 || strcmp( algorithm, "both" ) == 0 ){	
			printf("Making SGD\n");	
			if( bitonic ){
				r = new bitonicSGD( y2, dB, outdir, 90 );
			}
			else{
				r = new lassoSGD( y2, dB, outdir, 90, NULL );
			}
		}
		else if( strcmp( algorithm, "bgd" ) == 0 ){
			printf("Making BGD\n");
			r = new lassoBGD( y2, dB, outdir, 90, NULL );
		}
		else{
			printf("Algorithm does not match 'sgd', 'bgd', nor 'both'\n");
			return(0);
		}

		r->algorithm = algorithm;
		r->stepsize = stepsize;
		r->maxit = maxits;
		r->compute_sse = compute_sse;
		r->learn();	

		if( strcmp( algorithm, "both" ) == 0 ){	
			r->write_B("deconvolved_sX2.txt");
			printf("Making BGD\n");
			lassoSGD * rs = static_cast< lassoSGD * >(r);
			lassoBGD * rb = new lassoBGD( rs );
			rb->learn();
			r = rb;
		}
		r->write_B("deconvolved_X2.txt");

	}
	printf("Peak RSS: %i mb\n", getPeakRSS()/1000000 );
	printf("Done!\n");
}













