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
	
	for( int c = 0; c<argc; c++ ){
		printf("Arg %i: %s\n", c, argv[c] );
	}

	const char* y2file = argv[1];
        const char* bfile = argv[2];
	const char* x1file = argv[3];
	double o2 = atof( argv[4] );
	double bw2 = atof( argv[5] );

        const char* outdir = argv[6];
	double sparsity = atof( argv[7] );
	int normb = atoi( argv[8] );
	int gpnorm = 0;
	if( argc >= 10 ){
		gpnorm = atoi( argv[9] );
	}
	int batchits = 0;
	if( argc >= 11 ){
		batchits = atoi( argv[10] );
		printf("Batchits: %i\n", batchits );
	}
	char * precursor_annfile = NULL;
	if( argc >= 14 ){
		precursor_annfile = argv[13];
	}

	printf("argc:%i\n", argc );
	printf("normb: %i, gpnorm: %i\n", normb, gpnorm);
	printf("Reading Y\n");
        dbl_matrix *y2 = new dbl_matrix( y2file );
	printf("Reading Profiles B\n");
        dbl_matrix *b = new dbl_matrix( bfile );
	dbl_matrix * x1 = NULL;
	if( precursor_annfile == NULL ){
		printf("Reading Isotope Distributions\n");
	        x1 = new dbl_matrix( x1file );
	}
	else{
		printf("isotope distributions not given!\n");
	}

	b->max_norm();
	y2->max_norm();

        double pm = 741.237;
	if( argc >= 15 ){
		pm = atof(argv[14]);
	}
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


	// Change the dimensions so that y-nrow is the same as x->nrow (where it's really Bt->nrow)	
	if( y2-> nrow == b-> nrow ){
          printf("X and Y dimensions match\n");
	}
	else if( y2-> nrow == b-> ncol ){
          printf("X and Y dimensions don't match, transposing X\n");
	  b->t(); 
	}
	else if( y2->ncol == b->nrow ){
          printf("X and Y dimensions don't match, transposing Y\n");
	  y2->t();
         printf("New Y: %ix%i\n", y2->nrow, y2->ncol );
	}
	else if( y2->ncol == b->ncol ){
          printf("X and Y dimensions match but in the wrong way, transposing X and Y\n");
          b->t();
	  y2->t();
	}
	else{
	  printf("X and Y's dimensions don't match and nothing can be done about it :(\n");
	}

	printf("normb: %i\n",normb);
	if( normb == 1 ){	
		col_magnorm(b->m, b->nrow, b->ncol);
	}
	else if( normb == 2 ){
		col_squaremagnorm(b->m, b->nrow, b->ncol);
	}

	glassoSGD *sgd = new glassoSGD(y2,b,outdir,sparsity, x1, pm, bw1, o2, bw2, NULL, precursor_annfile );

	if( argc >= 12 ){
		sgd->stepsize = atoi( argv[11] );
	}

	sgd->batchits = batchits;
	sgd->gpnorm = gpnorm;
	printf("Done Making\n");
	sgd->maxit = 100;
	if( batchits != 0 ){
		sgd->maxit=batchits;

	}
	sgd->reset_epsilon_period = 15;	
	if( argc >= 13 ){
		sgd->reset_epsilon_period = atoi( argv[12] );
	}
//	sgd->X.write_csr("testb.txt");
	sgd->learn();
//	sgd->fg->precursorgroups[2]->output_latent_group("latentgroup2.txt");
//	sgd->fg->precursorgroups[3]->output_latent_group("latentgroup3.txt");
	//sgd->outdir = "gltestdir2/";
	//sgd->lambda*=2;
	//sgd->learn();
	double * group_penalties = sgd->fg->group_penalty_pp();

	char * outfile = new char[500];
	sprintf( outfile, "%s/group_penalties.txt", outdir );
        FILE *  file = fopen(outfile,"w");
	for( int n=0; n<sgd->N; n++ ){
		fprintf(file,"GP\t%i\t%f\n",n,group_penalties[n]);
	}	
	fclose(file);

	delete [] outfile;
	delete [] group_penalties;
//	delete sgd;
	//printf("Deleting.\n");
	return(0);
}
