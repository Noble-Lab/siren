#include "lassoSGD.h"

int main( int argc, char* argv[] ){
	
	//vector<dbl_matrix*> xy = get_xy_withms1("/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2015-12-10/SGS_ms1_0.1/991.27_ms1matrix.txt","/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2015-12-10/SGS_PeptideArt123_0.50025395/database_1004.txt", "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/ms2/FromSonia/napedro_L120224_010_SW_ms1_0.1/991.27_1017.27.txt","/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/ms2/FromSonia/napedro_L120224_010_SW/1004.27_0.50025395.txt" );
	vector<dbl_matrix*> xy = get_xy("/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2015-12-10/SGS_PeptideArt123_0.50025395/database_1004.txt","/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/ms2/FromSonia/napedro_L120224_010_SW/1004.27_0.50025395.txt" );
	//vector<dbl_matrix*> xy = get_y_rowmajorx("/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2015-12-10/SGS_PeptideArt123_0.50025395/database_1004.txt","/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/ms2/FromSonia/napedro_L120224_010_SW/1004.27_0.50025395.txt" );
//	vector<dbl_matrix*> xy = get_xy("/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2015-12-10/SGS_ms1_0.1/991.27_ms1matrix.txt","/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/ms2/FromSonia/napedro_L120224_010_SW_ms1_0.1/991.27_1017.27.txt");

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
	printf("Making lassoSGD.\n");


	lassoSGD *sgd = new lassoSGD(xy[1],xy[0],"moo",0.0 );
	sgd->maxit = 10;	
	printf("Learning!\n");
	sgd->learn();	


	delete sgd;
	//printf("Deleting.\n");
	return(0);
}
