#include "produce_fragments.h"

int main( int argc, char *argv[] ){
	printf("Reading data\n");

	const char* y2file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-07-29/SimulatedData/y2_N400T1000_0.txt"; 
	const char* bfile = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-07-29/SimulatedData/profs_N400T1000_0.txt"; 
	const char* x1file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-07-29/SimulatedData/targets1_N400T1000_0.txt"; 
	//const char* x1file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/x1_test.txt"; 
	const char* outdir = "/net/noble/vol3/user/alexhu/proj/DIARegression/bin/C++/2016-03-03/test_pf";

	dbl_matrix *y2 = new dbl_matrix( y2file );
	dbl_matrix *b = new dbl_matrix( bfile );
	dbl_matrix *x1 = new dbl_matrix( x1file );

	int N = x1->ncol;
	int M = y2->nrow;


	double pm = 749;
	double bw1 = 0.1;
	//double pm = 741.237;
	//double bw1 = 0.10005079;
	double o2 = 0.4;
	double bw2 = 0.10005079*10;


	printf("y2: %ix%i\n", y2->nrow, y2->ncol);
	fragmentgroup * fg = new fragmentgroup( x1, M, pm, bw1, o2, bw2 );
        //write_csr_matrix( bool *M, int nrow, int ncol, "bmasses.txt" );
        write_csr_matrix( fg->bmasses, M, M, "fg.txt" );
	int n = N/2+1;
	fragmentgroup::groupcsrmatrix * gm = new fragmentgroup::groupcsrmatrix( fg, n );
	write_csr_matrix( gm->nz, gm->maxm, gm->maxp, "csrints.txt" );
}
