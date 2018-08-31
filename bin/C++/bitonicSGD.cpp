/**
 * \file bitonicSGD.cpp
 * AUTHOR: Alex Hu
 * CREATE DATE: March 5, 2015
 * PROJECT: DIA Regression
 * \brief Returns peptide abundances in DIA spectra that optimize
 * the bitonic regression objective function
 *****************************************************************************/
#include "bitonicSGD.h"
using namespace std;

//The arguments are: Yfilename Xfilename Wfilename omega max number of scans to use
int main( int argc, char* argv[] ){
	
//        for( int i=1; i<argc; i++ ){ printf("%i: %s\n", i, argv[i] );}

	vector<dbl_matrix*> xy = get_xy_withms1("/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2015-12-10/SGS_ms1_0.1/991.27_ms1matrix.txt","/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2015-12-10/SGS_PeptideArt123_0.50025395/database_1004.txt", "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/ms2/FromSonia/napedro_L120224_010_SW_ms1_0.1/991.27_1017.27.txt","/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/ms2/FromSonia/napedro_L120224_010_SW/1004.27_0.50025395.txt" );


	printf("Making bitonicSGD.\n");
	bitonicSGD *sgd = new bitonicSGD(xy[0],xy[1],"moo",90 );
	
	printf("Learning!\n");
	sgd->learn( 1 );	
	
	printf("Done! Cleaning up\n");
	//delete sgd;
	return(0);
};


