/**
 * \file bitonicSGD_prevB.cpp
 * AUTHOR: Alex Hu
 * CREATE DATE: March 5, 2015
 * PROJECT: DIA Regression
 * \brief Returns peptide abundances in DIA spectra that optimize
 * the bitonic regression objective function
 *****************************************************************************/
#include "bitonicSGD_betterupdate_old.h"
using namespace std;

//The arguments are: Yfilename Xfilename Wfilename omega max number of scans to use
int main( int argc, char* argv[] ){
	
        for( int i=1; i<argc; i++ ){ printf("%i: %s\n", i, argv[i] );}

	bitonicSGD sgd;
	int maxT;
	char * outdir;
	double psparsity;
	char paramfile[200];
	
	int epochnum = atoi( argv[5] );
	maxT = atoi( argv[6] );
	outdir = argv[7];
	psparsity = atof( argv[8] );
	
	sprintf(paramfile,"%s/params.txt",outdir);	
	FILE * file;
        file = fopen(paramfile,"w");
        for( int i=1; i<argc; i++ ){ fprintf(file,"%s\n", argv[i] ); }
	fclose(file);
	sgd = get_bitonicSGD_prevB(argv[1], argv[2], argv[3], argv[4], epochnum, maxT, outdir, psparsity);
	
	char filename[200];
        double true_sparsity=0;
	double l1;
        double ls;
	double * Bt;

        sgd.learn();  
	printf("Computing ls\n",psparsity);
        ls = sgd.ls();        
        Bt = transpose( sgd.B, sgd.T, sgd.N );
	l1 = sgd.lambda*norm1(Bt,sgd.N,sgd.T);
        sprintf(filename,"%s/_B.txt",outdir);
	for( int i=0; i<sgd.NT; i++ ){ if( Bt != 0 ){ true_sparsity++; }}
      	true_sparsity = true_sparsity/sgd.NT;
	printf("SSE\tlambda\tl1penalty\tpercent non-zero\n");
      	printf("%f\t%f\t%f\t%f\n", ls, sgd.lambda, l1, true_sparsity );
	
	delete [] Bt;	
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


