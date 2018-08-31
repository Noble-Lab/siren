/**
 * \file bitonicSGD.cpp
 * AUTHOR: Alex Hu
 * CREATE DATE: March 5, 2015
 * PROJECT: DIA Regression
 * \brief Returns peptide abundances in DIA spectra that optimize
 * the bitonic regression objective function
 *****************************************************************************/
//#include "bitonicSGD.h"
#include "bitonicSGD_betterupdate_old.h"
//#include "bitonicSGD_betterupdate.h"
using namespace std;

//The arguments are: Yfilename Xfilename Wfilename omega max number of scans to use

int main( int argc, char* argv[] ){

	int N = 2;
	int T = 8;
	int M = 3;
	double X[6] = {2, 1, 2, 1, .2, .1};
	double X2[6] = { 1, 0, 0, 1, 0, 0 };
	print_matrix(X,M,N,"Raw");
	colnorm(X,M,N);
	print_matrix( X,M,N, "Normalized");
	magnorm(X2,M,N,X,M,N);
	print_matrix( X,M,N, "Mag Normalized");
	return(0);
	
	double B[16] = {1, 2, 3, 2, 1, 0, 0, 0, 0, 0, 0, 1, 2, 3, 2, 1};
	double W[16] = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0};
	double * Y = multiply( X, M, N, B, N, T );
	double  noise[24] = {0.02773116, -0.09167278, -0.15704929,  0.20243044, -0.05496895, -0.11395801, -0.28393765, -0.13150016,  0.14477702, 0.12675951 -0.07915461, -0.12324123,  0.26817194, -0.03470999, -0.08055519,  0.09053798,  0.14047286, -0.15911021, -0.03714298,  0.02480380, -0.03469844,  0.17957805,  0.21522673, -0.10971076 };
	for( int i=0; i<(M*T); i++ ){
		Y[i] += noise[i];
	}

	double * resid = noise;	

	print_matrix( X, M, N, "X");
	print_matrix( B, N, T, "B");
	print_matrix( Y, M, T, "Y");
	print_matrix( noise, M, T, "noise");

	int start = 0;
	int end = 4;
	
	double ** residscores = corrscore_W_backgroundsubtracted( Y, X, B, W, M, T, N );

	write_scores( residscores, 6, 2, "residscores.txt");

	return (0);
};

/*
//The arguments are: Yfilename Xfilename Wfilename omega max number of scans to use
int main( int argc, char* argv[] ){
	int omega = atoi( argv[4] );
	int maxT = atoi( argv[5] );
	char * outdir = argv[6];
	double psparsity = atof( argv[7] );

	char paramfile[200];
	sprintf(paramfile,"%s/params.txt",outdir);	
	FILE * file;
        file = fopen(paramfile,"w");
        for( int i=1; i<argc; i++ ){ fprintf(file,"%s\n", argv[i] ); }
	fclose(file);

	bitonicSGD sgd = get_bitonicSGD(argv[1], argv[2], argv[3], omega, maxT, outdir, psparsity);
	//print_matrix( sgd.peaksperpeptide, sgd.N, 1, "Peaks per peptide" );
	char filename[200];
        double true_sparsity;
	double l1;
        double ls;
	double * Bt;

	
        for( double sparsity=95; sparsity>= 0; sparsity -= 5 ){	
		printf("Iterating for sparsity %f\n", sparsity);
        	sgd.estimate_lambda( sparsity );      
        	sgd.learn();  
		printf("Computing ls\n",psparsity);
        	ls = sgd.ls();        
        	Bt = transpose( sgd.B, sgd.T, sgd.N );
		l1 = sgd.lambda*norm1(Bt,sgd.N,sgd.T);
                sprintf(filename,"/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2015-04-29/%f_sgd_B.txt",sparsity);
		true_sparsity = write_csr_matrix( Bt, sgd.N, sgd.T, filename );      
      	  	printf("ls_sparsity_truesparsity\t%f\t%f\t%f\t%f\n", ls, sparsity, true_sparsity, l1 );
		sgd.reset();
		delete [] Bt;
	}
	return (0);
};*/


