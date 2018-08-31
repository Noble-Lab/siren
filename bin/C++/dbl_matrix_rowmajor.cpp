#include "dbl_matrix_rowmajor.h"

int main( int argc, char* argv[] ){
        
	/*
	const char* y2file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-07-29/SimulatedData/y2_N400T1000_0.txt"; 
        const char* bfile = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-07-29/SimulatedData/profs_N400T1000_0.txt"; 
        //const char* x1file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/x1_test.txt"; 
	const char* x1file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-07-29/SimulatedData/targets1_N400T1000_0.txt";
        const char* outdir = "/net/noble/vol3/user/alexhu/proj/DIARegression/bin/C++/2016-03-03/test_pf";
	*/

	//const char* y2file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/100groups_2peps_y2.txt"; 
        //const char* x1file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/100groups_2peps_x1.txt";

	/*
	const char* y2file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/2groups_2peps_y2.txt"; 
        const char* bfile = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/2groups_2peps_profs.txt"; 
	const char* x1file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/2groups_2peps_x1.txt";
        const char* outdir = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/TestGroupLasso/";
	*/

	printf("Building y\n");	
        //dbl_matrix_rowmajor y = dbl_matrix_rowmajor("/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-08-17/SimulatedData/100groups_2peps_profs.txt");
        //dbl_matrix_rowmajor y = dbl_matrix_rowmajor("/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-07-29/SimulatedData/profs_N400T1000_0.txt");
        dbl_matrix* y = new dbl_matrix_rowmajor("/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-07-29/SimulatedData/profs_N400T1000_0.txt");
        //dbl_matrix *y = new dbl_matrix("/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-07-29/SimulatedData/profs_N400T1000_0.txt");

	dbl_matrix_rowmajor * moo = dynamic_cast<dbl_matrix_rowmajor * >(y);
	if( moo == NULL ){
		printf("moo is null\n");
	}
	else{
		printf("moo is not null\n");
	}

	printf("Writing y\n");	
	y->write_csr( "test.txt" );
	printf("Transposing y\n");	
	y->t();
	printf("Writing transposed y\n");	
	y->write_csr( "test_transpose.txt");
	printf("Done!\n");
	return(0);
}
