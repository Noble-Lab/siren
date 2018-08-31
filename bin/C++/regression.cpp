#include "regression.h"

int main( int argc, char* argv[] ){
	vector<dbl_matrix*> xy = get_xy_withms1("/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2015-12-10/SGS_ms1_0.1/991.27_ms1matrix.txt","/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2015-12-10/SGS_PeptideArt123_0.50025395/database_1004.txt", "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/ms2/FromSonia/napedro_L120224_010_SW_ms1_0.1/991.27_1017.27.txt","/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/ms2/FromSonia/napedro_L120224_010_SW/1004.27_0.50025395.txt" );

	xy[0]->write_csr( "test_x.txt");
	xy[1]->write_csr( "test_y.txt");
	regression *reg = new regression(*xy[0],*xy[1],"moo",90 );
	return(0);
}
