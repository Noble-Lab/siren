/**
 * \file alternateBW.cpp
 * AUTHOR: Alex Hu
 * CREATE DATE: March 5, 2015
 * PROJECT: DIA Regression
 * \brief Returns peptide abundances in DIA spectra that optimize
 * the bitonic regression objective function
 *****************************************************************************/
//#include "bitonicSGD.h"
#include "bitonicSGD_betterupdate_old.h"
//include "bitonicSGD_betterupdate.h"
//#include "scorePeptides.h"
using namespace std;

//The arguments are: Yfilename Xfilename Bfilename Width Scorefilename Sparsity
int main( int argc, char* argv[] ){

	int m;
	int t;
	int n;

	char * s = "moo";
	if( strcmp( s, "moo") == 0 ){
		printf("%s == %s\n", s, "moo" );
	}
	
	printf("Reading B: %s\n", argv[3]); 
	clock_t begin = clock();
	double * b = read_rowmajor_dblmatrix( argv[3], &n, &t, -1 );
	printf("%ix%i\n", n, t );
	printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000); 
	//print_matrix( w, n, t);

        for( int i=1; i<argc; i++ ){ printf("%i: %s\n", i, argv[i] );}
	
	printf("Reading Y: %s\nmaxT=%i\n", argv[1],t); 
	begin = clock();
	double * y = read_rowmajor_dblmatrix( argv[1], &m, &t, t );
	printf("%ix%i\n", m, t );
	//print_matrix( y, m, t);
	printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000); 
	
	printf("Reading X: %s\n", argv[2]); 
	begin = clock();
	double * x = read_rowmajor_dblmatrix( argv[2], &m, &n, -1 );
	printf("%ix%i\n", m, n );
	printf("time lapsed: %f (s), peak RSS: %f(mb)\n", double(clock()-begin)/CLOCKS_PER_SEC, getPeakRSS()/1000000); 
	//print_matrix( x, m, n);

	printf("Computing W\n");
        double * exW = zeros(n,t);
	double * penalties = compute_unweighted_penalties(b, n, t);
	compute_Ws( exW, n, t, penalties, 0, 1 );

	printf("Computing scores\n");
	int fixedwin = 3;
	double ** scores = corrscore_W_minwin( y, x, b, exW, m, t, n, fixedwin );	
	write_scores( scores, 5, n, argv[4] );
	return (0);
};
