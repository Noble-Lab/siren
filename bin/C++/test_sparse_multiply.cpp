#include "lassoBGD.h"

int main( int argc, char* argv[] ){
        
	/*
	const char* y2file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-07-29/SimulatedData/y2_N400T1000_0.txt"; 
        const char* bfile = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-07-29/SimulatedData/profs_N400T1000_0.txt"; 
        //const char* x1file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/x1_test.txt"; 
	const char* x1file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-07-29/SimulatedData/targets1_N400T1000_0.txt";
        const char* outdir = "/net/noble/vol3/user/alexhu/proj/DIARegression/bin/C++/2016-03-03/test_pf";
	*/

	//const char* y2file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/100groups_2peps_y2.txt"; 
        //const char* bfile = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/100groups_2peps_profs.txt"
        //const char* x1file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/100groups_2peps_x1.txt";

	/*
	const char* y2file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/2groups_2peps_y2.txt"; 
        const char* bfile = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/2groups_2peps_profs.txt"; 
	const char* x1file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/2groups_2peps_x1.txt";
        const char* outdir = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/TestGroupLasso/";
	*/
	
	for( int c = 0; c<argc; c++ ){
		printf("Arg %i: %s\n", c, argv[c] );
	}

	/*
	const char* yfile = argv[1];
        const char* bfile = argv[2];
	const char* xfile = argv[3];
        const char* outdir = argv[6];
	double sparsity = atof( argv[7] );
	bool * E;*/

	const char* yfile = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-09-26/SimulatedData/y1_N100T500_0.txt";
        const char* bfile = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-09-26/SimulatedData/profs_N100T500_0.txt";
	const char* xfile = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-09-26/SimulatedData/targets1_N100T500_0.txt";
        const char* outdir = "./";
	double sparsity = 0.0;
	bool * E = NULL;

	printf("Reading Y\n");
        dbl_matrix *y = new dbl_matrix( yfile );
	printf("Reading Profiles B\n");
        dbl_matrix *b = new dbl_matrix( bfile );
	printf("Reading X \n");
        dbl_matrix *x = new dbl_matrix( xfile );

	printf("Making glassoSGD.\n");

	lassoSGD * sgd = new lassoSGD( y, x, outdir, sparsity, E );
	printf("YB is:\n");
	sgd->Y.print();


	int M = sgd->M;

	//double * dmult = zeros( sgd->M, sgd->T );
	//double * smult = zeros( sgd->M, sgd->T );
	double * dmult = zeros( sgd->N, sgd->T );
	double * smult = zeros( sgd->N, sgd->T );
	
	int T = sgd->T;
	int N = sgd->N;
	printf("X is:\n");
	sgd->X.print();
	printf("M: %i, N: %i\n", sgd->M, sgd->N );

	printf("dense multiplying\n");
	Tmultiply( sgd->X.m, sgd->M, sgd->N, sgd->Y.m, M, T, 1.0, dmult );
	//multiply( sgd->X.m, sgd->M, sgd->N, b->m, N, T, 1.0, dmult );
	printf("sparse multiplying\n");
	sparseTmultiply( sgd->Xcsr, sgd->nzx, sgd->Y.m, M, T, N, sgd->maxnz, 1.0, smult ); 
	//sparseMultiply( sgd->Xcsr, sgd->nzx, b->m, M, T, sgd->maxnz, smult ); 

	write_csr_matrix( dmult, N, T, "dense_multiplied.txt");
	write_csr_matrix( smult, N, T, "sparse_multiplied.txt");

	//int MT = M*T;
	int NT = N*T;
	printf("NT = %i\n", NT );
	//printf("MT = %f\n", MT );
	bool same = true;
	for( int i=0; i<NT; i++ ){
		if( smult[i] != dmult[i] ){
			printf("unequal %i: %f %f\n", i, smult[i], dmult[i]);
			same = false;
			break;
		}
	}
	if(same){
		printf("it worked!\n");
	}

	printf("Done!\n");
	return(0);
}





