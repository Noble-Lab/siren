#include "glassoSGDMS1.h"
#include "deconvolution.h"

int main( int argc, char* argv[] ){
	
	for( int c = 0; c<argc; c++ ){
		printf("Arg %i: %s\n", c, argv[c] );
	}

	const char* y1file = argv[1];
	const char* x1file = argv[2];
	const char* groupfile = argv[3];
	const char* outdir = argv[4];
	double sparsity = atof(argv[5] );
	int numits = atoi( argv[6] );
	double stepsize = atof( argv[7] );
	bool compute_sse = atoi( argv[8] );


	printf("Reading Y\n");
        dbl_matrix *y1 = new dbl_matrix( y1file );
	printf("Reading Isotope Distributions\n");
        dbl_matrix *x1 = new dbl_matrix( x1file );

	printf("Making glassoSGD.\n");

	glassoSGDMS1 * g = new glassoSGDMS1(y1,x1,outdir, sparsity, groupfile );
	printf("Learning!\n");
	//g->compute_sse = true;
	g->maxit = numits;	
	g->stepsize = stepsize;
	g->compute_sse = compute_sse;
	g->learn();
	g->learn_batch();
	//g->learn();
	g->write_B();


	printf("Segmenting elution profiles\n");
	char filename[500];
	sprintf(filename,"%s/segment_annotations.txt",outdir);
	dbl_matrix *dB;
	dB = segment_elution_profiles( &g->B, filename );
		
	sprintf(filename,"%s/ms1_elutionB_segmented.txt",outdir);
	dB->write_csr(filename);	
	//col_magnorm( dB->m, dB->nrow, dB->ncol );

	//print_matrix( g->B.m, g->B.nrow, g->B.ncol, "B" );
	return(0);
}





