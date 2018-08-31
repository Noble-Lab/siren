#include "glassoSGDMS1.h"

int main( int argc, char* argv[] ){
	
	const char* y1file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-12-13/testgroupms1_y1.txt";
	const char* x1file = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-12-13/testgroupms1_x1.txt";
	//const char* groupfile = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-12-13/testgroupms1_annotation.txt";
	const char* groupfile = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-12-13/testgroupms1_annotation_ungrouped.txt";
        const char* outdir = "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2016-12-13/";
	double sparsity = -1;

	printf("Reading Y\n");
        dbl_matrix *y1 = new dbl_matrix( y1file );
	printf("Reading Isotope Distributions\n");
        dbl_matrix *x1 = new dbl_matrix( x1file );

	printf("Making glassoSGD.\n");

	glassoSGDMS1 * g = new glassoSGDMS1(y1,x1,outdir, sparsity, groupfile );
	g->stepsize = 1.0;
	g->maxit = 100;	
	printf("Learning!\n");
	g->compute_sse = true;
	g->learn();
	printf("WHAT?\n");	
	//print_matrix( g->B.m, g->B.nrow, g->B.ncol, "B" );
	g->write_B();
	return(0);
}
