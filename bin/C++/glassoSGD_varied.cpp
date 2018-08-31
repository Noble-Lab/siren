#include "glassoSGD.h"

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

	const char* y2file = argv[1];
        const char* bfile = argv[2];
	const char* x1file = argv[3];
	double o2 = atof( argv[4] );
	double bw2 = atof( argv[5] );

        const char* outdir = argv[6];
	double sparsity = atof( argv[7] );
	bool normb = atoi( argv[8] );

        dbl_matrix *y2 = new dbl_matrix( y2file );
        dbl_matrix *b = new dbl_matrix( bfile );
        dbl_matrix *x1 = new dbl_matrix( x1file );

        double pm = 741.237;
        double bw1 = 0.1;
        //double bw1 = 0.10005079;
        //double pm = 749;
        //double o2 = 0.4;
        //o2 = 0.0;
	//double bw2 = 0.10005079*10;
	
	//xy[0]->sparsify(xy[1]);	
	/*
	int t = 10;
	dbl_matrix * b = new dbl_matrix( xy[0]->ncol, t );
	int size = t*xy[0]->ncol;
	for( int i=0; i<size; i++ ){ (*b)[i] = 1; }
	dbl_matrix * y = new dbl_matrix( xy[0]->nrow, t );
	multiply( xy[0]->m, xy[0]->nrow, xy[0]->ncol, b->m, b->nrow, t, 1.0, y->m );
	y->print();
*/
	printf("Making glassoSGD.\n");

	y2->t();
	b->t();

	if( normb ){
		col_magnorm(b->m, b->nrow, b->ncol);
	}

	glassoSGD *sgd = new glassoSGD(y2,b,outdir,sparsity, x1, pm, bw1, o2, bw2 );
	printf("Done Making\n");
	sgd->maxit = 100;	
	printf("Learning!\n");
	sgd->learn();
	//sgd->outdir = "gltestdir2/";
	//sgd->lambda*=2;
	//sgd->learn();

	


//	delete sgd;
	//printf("Deleting.\n");
	return(0);
}
