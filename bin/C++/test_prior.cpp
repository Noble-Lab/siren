#include "lassoBGD.h"

int main( int argc, char* argv[] ){

	for( int i=1; i<argc; i++ ){
		printf("%i: %s\n",i,argv[i]);
	}

	char * yfile = argv[1];
	char * xfile = argv[2];
	char * efile = argv[3];
	char * outdir = argv[4];
        double psparsity = atof( argv[5] );

	//vector<dbl_matrix*> xy = get_xy_withms1("/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2015-12-10/SGS_ms1_0.1/991.27_ms1matrix.txt","/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2015-12-10/SGS_PeptideArt123_0.50025395/database_1004.txt", "/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/ms2/FromSonia/napedro_L120224_010_SW_ms1_0.1/991.27_1017.27.txt","/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/ms2/FromSonia/napedro_L120224_010_SW/1004.27_0.50025395.txt" );
	//vector<dbl_matrix*> xy = get_xy("/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/2015-12-10/SGS_ms1_0.1/991.27_ms1matrix.txt","/net/noble/vol3/user/alexhu/proj/OpenSWATH/data/ms2/FromSonia/napedro_L120224_010_SW_ms1_0.1/991.27_1017.27.txt");
	vector<dbl_matrix*> xy = get_xy(xfile,yfile);
	int n;
	int t;
	bool * e = read_rowmajor_bool( efile, &n, &t, (xy[1])->ncol );

	printf("Making lassoBGD.\n");
	lassoBGD *bgd = new lassoBGD(xy[0],xy[1],outdir,psparsity, e );
	bgd->maxit = 100;
	printf("iterating!\n");
	bgd->iterate();

	return(0);
}
