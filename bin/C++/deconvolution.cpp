#include "deconvolution.h"

int main( int argc, char* argv[] ){

	dbl_matrix * B = new dbl_matrix( "test_B.txt" );
	print_matrix( B->m, B->nrow, B->ncol, "Abundances" );

	dbl_matrix * profiles = segment_elution_profiles( B );
	print_matrix( profiles->m, profiles->nrow, profiles->ncol, "Profiles" );

}
