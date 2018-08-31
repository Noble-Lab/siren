#ifndef DBL_MATRIX_ROWMAJOR_H
#define DBL_MATRIX_ROWMAJOR_H

#include <vector>
#include "dbl_matrix.h"

class dbl_matrix_rowmajor: public dbl_matrix {
public:
	double * X;
	int * nzx;
	int maxnz;

	int * nz_per_column;


	// Initialization
	dbl_matrix_rowmajor();
	dbl_matrix_rowmajor( char const * mfile, char const * type );
	dbl_matrix_rowmajor( int NROW, int NCOL );
	dbl_matrix_rowmajor( int NROW, int NCOL, int MAXNZ );

	// Initialization


	// Operations
	dbl_matrix_rowmajor * copy_dbl_matrix_rowmajor();
	void sqrt();
	//void sparsify(dlb_matrix * y);

	void colnorm();

	double norm2();
	void rbind( dbl_matrix_rowmajor * b );
	int nnz();
	double max();
//	dbl_matrix_rowmajor * multiply( dbl_matrix_rowmajor * B, double s );
//	dbl_matrix_rowmajor * Tmultiply( dbl_matrix_rowmajor * B, double s );

	// "Sparsify X and Y
	// Remove the rows that have no signal in either X or Y
	void sparsify( dbl_matrix_rowmajor * b );
	void t(); // Tranpose it
	void col_magnorm();
	void sum_norm();
	void max_norm();
	bool * match_isotopes( dbl_matrix_rowmajor * X, dbl_matrix_rowmajor * Y, int minmatch );
	int * oldMIndices; // Maps the current m indices to the old m indices
	// before the uninformative rows were removed

	// Output
	void print();
	int write_csr( char const * outfile );
	~dbl_matrix_rowmajor();

};

dbl_matrix_rowmajor::dbl_matrix_rowmajor( ){
	nrow=0;
	ncol=0;
	X=NULL;
	nzx=NULL;
	m=NULL;
	maxnz=0;
	nz_per_column = NULL;	
} 

dbl_matrix_rowmajor::dbl_matrix_rowmajor(int NROW, int NCOL ){
	nrow=0;
	ncol=0;
	m=NULL;
	X=NULL;
	nzx=NULL;
	maxnz=0;
	nz_per_column = NULL;	
} 

dbl_matrix_rowmajor::dbl_matrix_rowmajor(int NROW, int NCOL, int MAXNZ ){
	nrow=NROW;
	ncol=NCOL;
	m=NULL;

	X=zeros(NROW,MAXNZ);
	nzx=zeros_int(NROW,MAXNZ,-1);

	maxnz=MAXNZ;
	nz_per_column = zeros_int(NCOL,1);	
} 

// This file contains the matrix in sparse form: i j value
// The first line contains the dimensions of the matrix: nrow ncol
dbl_matrix_rowmajor::dbl_matrix_rowmajor( char const * mfile, char const * type ="rowmajor" ){ 

	int nnz=0;
	maxnz = nnz;
	m=NULL;
	ifstream file( mfile );
	string line;
	getline( file, line);
	istringstream iss(line);
	iss >> nrow; iss >> ncol;
	oldMIndices = NULL;

	// Read it through one time just to see how many nonzero values
	// are in each row
	// rowmajor format
	int i;
	int j;
	char meh[50];
	while(file.good()){
		getline( file, line);
		istringstream s(line);
		if( line[0] == '>' ){
			if( nnz > maxnz ){
				maxnz = nnz;
			}
			nnz = 0;
			continue;	
		}
		s >> j;
		if( line.length() > 2 ){
			nnz += 1;	
		}
	}
	if( nnz > maxnz ){
		maxnz = nnz;
	}
	file.close();
	file.clear();
	file.open( mfile );

	X = zeros(nrow,maxnz);
	int size = nrow*maxnz;
	nzx = new int[size];
	for( i=0; i<size; i++ ){
		nzx[i] = -1;
	}

	getline( file, line);
	double v;

	// Read it to actually store the values now
	// rowmajor format
	i=0;
	j=0;
	int ix;
	printf("Reading rowmajor\n");
	nnz = 0;
	nz_per_column = zeros_int( ncol, 1 );	
	while(file.good()){
		getline( file, line);
		istringstream s(line);
		if( line.length() < 1 ) {
			continue;
		}
		if( line[0] == '>' ){
			s >> meh;
			s >> meh;
			i = atoi(meh);	
			ix = i*maxnz;
			continue;	
		}
		s >> j;
		if( i < nrow && j < ncol ){
			s >> v;
			X[ix] = v;
			nzx[ix] = j;
			nz_per_column[j] += 1;
			ix += 1;
			nnz += 1;
		}
	}
	file.close();
	printf("%s type: %s, (%ix%i)\n", mfile, type, nrow, ncol );
	printf("nnz: %i\n", nnz );
}

//int write_csr_matrix( double *M, int nrow, int ncol, char const * outfile ){
int dbl_matrix_rowmajor::write_csr(  char const * outfile ){
	FILE * file;
	printf("WHAT\n");
	file = fopen(outfile,"w");
	printf("IS THIS\n");
	int ix;
	int col;
	int nnz = 0;
	int end;
	fprintf(file,"%i\t%i\n",nrow,ncol);
	for( int r=0; r<nrow; r++ ){
		fprintf(file,"> %i\n", r );	
		ix = r*maxnz;
		end = ix+maxnz;
		col = nzx[ix];
		while( col != -1 && ix < end ){
			fprintf(file,"%i\t%f\n",col,X[ix]);
			ix += 1;
			col = nzx[ix];
			nnz += 1;
		}
	}
	fclose(file);
	return(nnz);
}

void dbl_matrix_rowmajor::t(){
	/*
	int nrow;
	int ncol;
	double * X;
	int * nzx;
	int maxnz;
	*/
	int * nzs = zeros_int( ncol, 1 );	
	// Find the new maxnz, which is the maximum number of rows each column
	// has a nonzero value in
	int end;
	int ix;
	for( int r=0; r<nrow; r++ ){
		ix = r*maxnz;
		end = ix+maxnz;	
		while( ix < end && nzx[ix] != -1 ){
			nzs[nzx[ix]] += 1;
			ix += 1;
		}	
	}
	int new_maxnz = maxint( nzs, ncol );
	printf("New maxnz: %i\n", new_maxnz );
	int new_nrow = ncol;
	int new_ncol = nrow;

	int new_size = new_maxnz*new_nrow;
	int * new_nzx = new int[new_size];

	for( int i=0; i<new_size; i++ ){ new_nzx[i] = -1; }
	double * new_X = zeros( new_maxnz, new_nrow );

	// Now fill the new matrices up!
	// nzs tells you where to put the new column numbers in the new_nzix
	// indexed by the new row number/old column number, the element is
	// the position within the new matrices of where you should put the 
	// new column number/old row number
	nzs[0] = 0;
	for( int i=1; i<ncol; i++ ){ nzs[i] = nzs[i-1] + new_maxnz; } 

	int new_row;
	int oldix;
	int newix;

	if( nz_per_column != NULL){
		delete [] nz_per_column;
		nz_per_column = zeros_int( nrow, 1 );
	}

	for( int r=0; r<nrow; r++ ){
		oldix = r*maxnz; // This index the old matrices
		end = oldix+maxnz;	
		while( oldix < end && nzx[oldix] != -1 ){
			new_row = nzx[oldix];
			newix = nzs[new_row];
			new_nzx[newix] = r;
			nz_per_column[r]+= 1;
			new_X[newix] = X[oldix];
			nzs[new_row] += 1;	
			oldix += 1;
		}	
	}

	delete [] X;
	delete [] nzx;
	delete [] nzs;
	X = new_X;
	nzx = new_nzx;
	maxnz = new_maxnz;
	nrow = new_nrow;
	ncol = new_ncol;	
}

/*
		x->sqrt();
		y->sqrt();
	}
	else{ printf("\n\n\nNot square rooting.\n"); }

	x->max_norm();
	y->max_norm();

	if ( rescale ){
		printf("Normalizing columns of x.\n");
		col_magnorm( x->m, x->nrow, x->ncol );
*/

void dbl_matrix_rowmajor::sqrt(){
	sqrt_mat( X, nrow, maxnz );
}

double dbl_matrix_rowmajor::max(){
	int size = nrow*maxnz;
	int max = 0; 
	for( int i=0; i<size; i++ ){
		if( nzx[i] != -1 && max < X[i]){
				max = X[i];
		}				
	}	
	return(max);
}

void dbl_matrix_rowmajor::max_norm(){
	int maxval = max(); 
	if( maxval == 0 ){
		return;
	}
	int size = nrow*maxnz;
	for( int i=0; i<size; i++ ){
		if( nzx[i] != -1){
			X[i] = X[i]/maxval;
		}				
	}	
}

void dbl_matrix_rowmajor::col_magnorm(){
	double * mags = zeros(ncol,1);
	int size = nrow*maxnz;
	int colnum;
	// Compute the sum of each column
	for( int i=0; i<size; i++ ){
		colnum = nzx[i];
		if( colnum != -1){
			mags[colnum] += X[i]*X[i];
		}				
	}
	for( int i=0; i<ncol; i++ ){
		mags[i] = pow(mags[i],0.5);
	}

	// Normalize each value by the sum of its column
	for( int i=0; i<size; i++ ){
		colnum = nzx[i];
		if( colnum != -1){
			X[i] /= mags[colnum];
		}				
	}
	delete [] mags;	
}

void dbl_matrix_rowmajor::colnorm(){
	double * sums = zeros(ncol,1);
	int size = nrow*maxnz;
	int colnum;
	// Compute the sum of each column
	for( int i=0; i<size; i++ ){
		colnum = nzx[i];
		if( colnum != -1){
			sums[colnum] += X[i];
		}				
	}
	// Normalize each value by the sum of its column
	for( int i=0; i<size; i++ ){
		colnum = nzx[i];
		if( colnum != -1){
			X[i] /= sums[colnum];
		}				
	}
	delete [] sums;	
}

double dbl_matrix_rowmajor::norm2(){
	double sums = 0;
	int size = nrow*maxnz;
	int colnum;
	// Compute the sum of each column
	for( int i=0; i<size; i++ ){
		colnum = nzx[i];
		if( colnum != -1){
			sums += X[i]*X[i];
		}				
	}
	return(sums);
}

/*
void dbl_matrix_rowmajor::sparsify(dlb_matrix * y){

	clock_t begin = clock();
	// Identify the m's to keep
	//bool occupied[M] = {false};
	int M = nrow;
	bool * occupied = new bool[M];
	for( int m=0; m<M; m++ ){ occupied[m] = false; };
	int nz;
	new_maxnz = 0;
	xnnz = 0;

	for( int m=0; m<M; m++ ){
		nz=0;
		for( int n=0; n<N; n++ ){
			if( X[m*N+n] > 0 ){
				if( ! occupied[m] ){
					newM++;
				}
				occupied[m] = true;
				peaksperpeptide[n] += 1;
				nz+=1;		
			}
		}
		xnnz += nz;
		if( nz > new_maxnz ){ maxnz = nz; }
	}
	printf(" M=%i\n", newM);
	
	// Create the new X and Y
	dbl_matrix *newX = new dbl_matrix_rowmajor(newM,N,maxnz);
	dbl_matrix *newY = new dbl_matrix(newM,T);
	nzx = zeros_int(newM,maxnz); for(int i=0; i<newM*maxnz; i++ ){ nzx[i] = -1; }
	oldMIndices = zeros_int(newM,1);
	Xcsr = zeros(newM,maxnz); 
	// Copy into them only the m's to keep
	int newm = 0;
	int xnoffset;
	int nzix;

	for( int m=0; m<M; m++ ){
		nzix = newm*maxnz;
		if( !occupied[m] ){ continue; }
		for( int n=0; n<N; n++ ){
			if( X[m*N+n] > 0.0 ){
				newX->m[newm*N+n] = X[m*N+n];
				nzx[nzix] = n;
				Xcsr[nzix] = X[m*N+n];
				nzix++;
			}


};*/

dbl_matrix_rowmajor::~dbl_matrix_rowmajor(){
	printf("Deleting dbl_matrix_rowmajor (%ix%i)\n",nrow,ncol);
	if( X != NULL ){
		delete [] X;
	}
	if( nzx != NULL ){
		delete [] nzx;
	}
}




#endif
