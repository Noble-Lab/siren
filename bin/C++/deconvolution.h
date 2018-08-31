#include "dbl_matrix.h"
#include <vector>

// For a single row of values, find the indices
// of the minima
// minima needs to be a pre-initialized int of length len
// -1 in minima signals the end of the list of minima
// The minima can be 0
void find_minima( double * m, int len, int * minima ){
	int j = 0;
	if( m[0] < m[1] ){
		minima[j] = 0;
		j++;
	}
	double v;
	for( int t=1; t<len-1; t++ ){
		v = m[t];
		if( m[t-1] >= v && m[t+1] >= v && ( v != m[t-1] || v != m[t+1] )  ){
			minima[j] = t;
			j++;
		}
	} 
	if( m[len-1] < m[len-2] ){
		minima[j] = len-1;
		j++;
	}
	minima[j] = -1;
}

// For a single row of values, find the indices
// of the bounds of elution profiles
// minima needs to be a pre-initialized int array of length len
// -1 in minima signals the end of the list of minima
// bounds needs to be a pre-initialized int array of length len
// The bounds of non-zero elution profiles come in pairs in a sequence
// here
int find_bounds( double * b, int len, int * minima, int * bounds, int minlen = 7 ){

	if( minlen < 2 ){
		printf("Bound minlen cannot be less than 2\n");
		return(0);
	}
	int numprofiles=0;
	int t1;
	int t2;
	bool nonzero;
	int j=0;
	find_minima( b, len, minima );
	for( int i=0; i<len-1; i++ ){
		t1 = minima[i];
		t2 = minima[i+1];
		if( t2 == -1 ){
			break;
		}
		// Check that the area between the minima is non-zero
		// If it is, add the two minima to bounds and consider it an
		// elution profile.
		nonzero = true;
		for( int t=(t1+1); t<t2; t++ ){
			if( b[t] == 0 ){
				nonzero = false;
				break;
			}
		}
		if( nonzero && (t1+minlen) <= t2 ){// && abs(1015-t1) < 10 ){
			numprofiles++;	
			bounds[j]=t1;
			bounds[j+1]=t2;
			j += 2;
		}
	}	
	bounds[j] = -1; // Signifies the end of the bounds
	return(numprofiles);
}

// Create a new matrix that contains an elution profile per row from B
dbl_matrix * segment_elution_profiles( dbl_matrix * B, char const * annotationfilename ){

	int T = B->ncol;
	int N = B->nrow;
	printf("N,T = %i,%i\n", N, T );
	int size = N*T;
	int * bounds = new int[size];
	int * minima = new int[T];
	int numprofiles = 0;
	// For each N, find the time bounds for each elution profile.
	for( int n=0; n<N; n++ ){
		numprofiles += find_bounds( &(B->m[n*T]), T, minima, &bounds[n*T] );
//		printf("Finding bounds for %i\n", n );
//		print_matrix( &bounds[n*T], 1, T, "Bounds" );
//		printf("%i: %i\n", n, numprofiles );
	}
	printf("Numprofiles: %i\n", numprofiles );

	//exit(0);

	// This should be TxN!	
	dbl_matrix * profiles = new dbl_matrix( T, numprofiles );
	int pnum = 0;
	int i;

	int start;
	int end;
	int pix;

	double total_abundance;
	FILE * file;
	file = fopen(annotationfilename,"w");
	fprintf(file,"candidate index\tstart\tend\n");
	// For each peptide n
	for( int n=0; n<N; n++ ){
		// For each elution profile
		i = n*T;
		while( bounds[i] != -1 ){
			start = n*T+bounds[i];
			end = n*T+bounds[i+1];
			fprintf(file,"%i\t%i\t%i\n",n,start-n*T,end-n*T);
			total_abundance = 0;
			// Find the total abundance of the elution profile
//			for( int j=start; j<=end; j++ ){
//				total_abundance += B->m[j]; 			
//			}	
			total_abundance = 1;
			// Fill in the new profile matrix
	//		pix = bounds[i] + pnum*T;
			pix = pnum + bounds[i]*numprofiles;
			for( int j=start; j<=end; j++ ){
				profiles->m[pix] = B->m[j]/total_abundance;
				pix = pix+numprofiles; 
			}	
			i += 2;
			pnum++;
		}	
	}

	fclose(file);
	delete [] bounds;
	delete [] minima;	
	return(profiles);
}












