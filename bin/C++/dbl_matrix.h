#include <vector>

#ifndef DBL_MATRIX_H
#define DBL_MATRIX_H

#include "matrix.h"

class dbl_matrix {
public:
	int nrow;
	int ncol;
	double * m;
	double scale;

	// Initialization
	dbl_matrix();
	dbl_matrix( int NROW, int NCOL, double scale );
	dbl_matrix( char const * mfile, char const * type, double scale );

	// Destruction
	virtual ~dbl_matrix();

	// Operations
	dbl_matrix * copy_dbl_matrix();
	void sqrt();
	double norm2();
	void rbind( dbl_matrix * b );
	double &operator[](int i){
		return m[i];
	}
	int nnz();
	dbl_matrix * multiply( dbl_matrix * B, double s );
	dbl_matrix * Tmultiply( dbl_matrix * B, double s );

	// "Sparsify X and Y
	// Remove the rows that have no signal in either X or Y
	void sparsify( dbl_matrix * b );
	void t(); // Tranpose it
	void col_magnorm();
	void sum_norm();
	void max_norm();
	bool * match_isotopes( dbl_matrix * X, dbl_matrix * Y, int minmatch );
	int * oldMIndices; // Maps the current m indices to the old m indices
	// before the uninformative rows were removed

	// Output
	void print();
	int write_csr( char const * outfile );

};

bool * dbl_matrix::match_isotopes( dbl_matrix * X, dbl_matrix * Y, int minmatch = 2 ){

	int N = X->ncol;
	int M = X->nrow;
	int T = Y->ncol;
	int NT = N*T;
	bool * matched = new bool[NT]; for( int i=0; i<NT; i++ ) { matched[i] = false; }

	int count; // Counts how many isotope peaks have been considered for each peptide
	int nummatched;
	bool * Xt = transpose_bool( X->m, M, N );
	bool * Xtr; // A partiular row of Xt
	for( int n=0; n<N; n++ ){
		Xtr = &(Xt[n*M]);
		for( int t=0; t<T; t++ ){
			count = 0;
			nummatched = 0;	
			for( int m=0; m<M; m++ ){
				if( Xtr[m] ){
					count++;
					if( Y->m[m*T+t] != 0  ){
						nummatched++;
					}
					if( count == minmatch ){
						break;
					}
				}	
			}
			if( nummatched == minmatch ){
				matched[n*T+t] = true;	
			}	
		}
	} 
	return(matched);
}

void dbl_matrix::t(){

	double * newm = transpose( m, nrow, ncol );
	int temp = nrow;
	nrow = ncol;
	ncol = temp;
	delete [] m;
	m = newm;

}

void dbl_matrix::print(){
	printf("dbl_matrix dimensions: %ix%i with %i nnz\n",nrow,ncol, nnz());
}


// Looks for m/z bins that are unpopulated by X and removes them from
// both X and Y
void dbl_matrix::sparsify( dbl_matrix * y ){
	
	// Identify the m's to keep
	int M = nrow;
	int N = ncol;
	int T = y->ncol;
	bool * occupied = new bool[M]; for( int k=0; k<M; k++ ){ occupied[k] = false; }
	int newM = 0;
	int offset;
	int nz;
	int maxnz = 0;
	int xnnz = 0;
	printf("Sparsifying: from M=%i to", M );
	for( int r=0; r<M; r++ ){
		nz=0;
		for( int n=0; n<N; n++ ){
			if( m[r*N+n] > 0 ){
				if( ! occupied[r] ){
					newM++;
				}
				occupied[r] = true;
				nz+=1;		
			}
		}
		xnnz += nz;
		if( nz > maxnz ){ maxnz = nz; }
	}
	printf(" M=%i\n", newM);
	
	// Create the new X and Y
	double *newX = zeros(newM,N);
	double *newY = zeros(newM,T);
	printf("Setting oldMIndices to something!\n");
	oldMIndices = zeros_int( newM,1 );
	// Copy into them only the m's to keep
	int newm = 0;
	int xnoffset;
	for( int r=0; r<M; r++ ){
		if( !occupied[r] ){ continue; }
		for( int n=0; n<N; n++ ){
			if( m[r*N+n] != 0.0 ){
				newX[newm*N+n] = m[r*N+n];
			}
		}
		// Copy Ys
		for( int t=0; t<T; t++ ){
			if( (*y)[r*T+t] == 0 ){ continue; }
			newY[newm*T+t] = (*y)[r*T+t];
		}
		oldMIndices[newm] = r;
		newm++;
	}

	delete [] m;
	delete [] y->m;
	
	m = newX;
	nrow = newM;
	y->m = newY;
	y->nrow = newM;
	y->oldMIndices = copy( oldMIndices, newM, 1 );


	printf("maxnz=%i\n",maxnz);
	

};

void dbl_matrix::sum_norm(){

	double sum = 0;
	int size=nrow*ncol;
	for( int i=0; i<size; i++ ){
		if( m[i] != 0 ){
			sum += fabs(m[i]);			
		}
	}
	printf("sumnorm: %f\n",  sum );
	for( int i=0; i<size; i++ ){
		if( m[i] != 0 ){
			m[i] = m[i]/sum;	
		}
	}

}

void dbl_matrix::max_norm(){

	double max=0;
	int size=nrow*ncol;
	for( int i=0; i<size; i++ ){
		if( m[i] != 0 && fabs(m[i]) > max ){
			max = fabs(m[i]);		
		}
	}
	
	printf("Maxnorm max: %f\n", max );
	for( int i=0; i<size; i++ ){
		if( m[i] != 0 ){
			m[i] = m[i]/max;	
		}
	}
	scale = scale*max;
	printf("Scale = %f\n",scale);
}

int dbl_matrix::nnz(){

	int nnz = 0;
	int size=nrow*ncol;
	for( int i=0; i<size; i++ ){
		if( m[i] != 0 ){
			nnz += 1;
		}
	}
	return nnz;
}

dbl_matrix::dbl_matrix(){
	nrow = 1;
	ncol = 1;
	scale = 1;
	m = zeros(nrow,ncol);
	oldMIndices = NULL;
}

void dbl_matrix::sqrt(){
	sqrt_mat( m, nrow, ncol );
}

// Replaces its elements with a concatenated matrix
// of the original and b.
void dbl_matrix::rbind(dbl_matrix *b){

	int sizeA = nrow*ncol;
	int sizeB = b->nrow*b->ncol;
	int sizeC = sizeA+sizeB;

	double * C = zeros(nrow+b->nrow,b->ncol);
	for( int i=0; i<sizeA; i++ ){
		C[i] = m[i];
	}
	int cix = sizeA;
	for( int i=0; i<sizeB; i++ ){
		C[cix] = b->m[i];
		cix++;
	}

        delete [] m;
        m = C;
        nrow = nrow+b->nrow;
}


dbl_matrix::dbl_matrix( int NROW, int NCOL, double sc = 1.0 ){
	nrow = NROW;
	ncol = NCOL;
	scale = sc;
	m = zeros( nrow, ncol );
	oldMIndices = NULL;
};

double dbl_matrix::norm2(){
	double n2 = 0;
	int size = nrow*ncol;
	for( int i=0; i<size; i++ ){
		n2 += m[i]*m[i];
	}
	return(n2);
};

dbl_matrix * dbl_matrix::copy_dbl_matrix(){

	dbl_matrix * n = new dbl_matrix( nrow, ncol );
	int size = nrow*ncol;
	for( int i=0; i<size; i++ ){
		if( m[i] != 0 ){
			(*n)[i] = m[i];
		}
	}
	return n;

}

// This file contains the matrix in sparse form: i j value
// The first line contains the dimensions of the matrix: nrow ncol
//
dbl_matrix::dbl_matrix( char const * mfile, char const * type ="rowmajor", double sc = 1.0 ){ 

	ifstream file( mfile );
	string line;
	getline( file, line);
	istringstream iss(line);

	char d;

	iss >> nrow; iss >> ncol;
	iss >> d;

	char const * dense = "dense";
	char const * sparse = "sparse";
	if( d == 'D' || d == 'd' ){
		type = dense;
	}
	else if( d == 's' || d == 'S' ){
		type = sparse;
	}

	m = zeros((nrow),(ncol));
	oldMIndices = NULL;
	double v;
	scale = sc;
	// rowmajor format
	if( strcmp( type, "rowmajor") == 0 ){
		int i;
		int j;
		char meh[50];
		printf("Reading rowmajor\n");
		while(file.good()){
			getline( file, line);
			if( line.length() == 0 ){
				continue;
			}
			istringstream s(line);
			if( line[0] == '>' ){
				s >> meh;
				s >> meh;
				i = atoi(meh);	
				continue;	
			}
			s >> j;
			if( i < nrow && j < ncol ){
				s >> v;
				m[i*(ncol)+j] += v;
			}
		}
	}
	// sparse format
	else if( strcmp( type, "sparse" ) == 0 ){
		int i;
		int j;

		// I have a format of matrix that is just a column vector with ncol == 1,
		// so in vectors like that, only the row index and intensity is given.

		if( ncol == 1 ){
			while(file.good()){
				getline( file, line);
				istringstream s(line);
				s >> i;
				s >> v;
				if( i < nrow  ){
					s >> v; 
					m[i] += v;
				}
			}
		}

		else{
			while(file.good()){
				getline( file, line);
				istringstream s(line);
				s >> i;
				s >> j;
				if( i < nrow && j < ncol ){
					s >> v; 
					m[i*(ncol)+j] += v;
				}
			}
		}
	}
	
	// dense format
	else if( strcmp( type, "dense")  == 0){
		int size = (nrow)*(ncol);	
		int i=0;
		while(file.good()){
			getline( file, line);
			istringstream s(line);
			for(int c=0; c<ncol; c++ ){ 
				if( i == size ){ break; }	
				s >> m[i];
				i++;
			}
		}
	}
	file.close();


	printf("%s type: %s, (%ix%i)\n", mfile, type, nrow, ncol );
	printf("nnz: %i\n", nnz() );
}

//int write_csr_matrix( double *M, int nrow, int ncol, char const * outfile ){
int dbl_matrix::write_csr(  char const * outfile ){
	FILE * file;
	file = fopen(outfile,"w");
	int offset=0;
	fprintf(file,"%i\t%i\n",nrow,ncol);
	int nnz = 0;
	int end;
	for( int r=0; r<nrow; r++ ){
		end = offset+ncol;
		fprintf(file,">\t%i\n", r );	
		for( int c=offset; c<end; c++ ){
			if( m[c] > 0.0 ){
				fprintf(file,"%i\t%f\n", c-offset, m[c] );	
				nnz += 1;
			}
		}
		offset = end;
	}
	fclose(file);
	return(nnz);
}

dbl_matrix * dbl_matrix::multiply( dbl_matrix * B,  double s = 1.0 ){

        double sum=0;
	int bcols = B->ncol;
        int size = nrow*bcols;
	double * bm = B->m;

        dbl_matrix * p = new dbl_matrix(nrow,bcols);

        for( int i=0; i<nrow; i++ ){
                for( int j=0; j<bcols; j++ ){
                        sum=0;
                        for( int k=0; k<ncol; k++ ){
				sum += m[i*ncol+k]*bm[k*bcols+j];
                        }
                        (*p)[i*bcols+j] = sum;
                }
        }

        if( s != 1.0 ){
                for( int i=0; i<size; i++ ){
                        (*p)[i] *= s;
                }
        }
        return(p);
}

// void col_magnorm( double * M, int nrow, int ncol ){
/*
void dbl_matrix::col_magnorm(){
	col_magnorm( m, nrow, ncol );
}*/

// Here, arows == brows
// Returns At x B
dbl_matrix * dbl_matrix::Tmultiply( dbl_matrix * b, double s = 1.0 ){

	this->print(); b->print();
	int bcols = b->ncol;
	dbl_matrix * p = new dbl_matrix(ncol,bcols);	
	double * bm = b->m;
	for( int r=0; r<nrow; r++ ){
		for( int t=0; t<bcols; t++ ){	
			for( int n=0; n<ncol; n++ ){
				if( m[r*ncol+n] == 0 || bm[r*bcols+t] == 0 ){ continue; }
				(*p)[n*bcols+t] += m[r*ncol+n]*bm[r*bcols+t];	
			}
		}
	}
	int size = ncol*bcols;
        if( s != 1.0 ){
                for( int i=0; i<size; i++ ){
                        (*p)[i] *= s;
                }
        }

	return(p);
}

dbl_matrix::~dbl_matrix(){
	printf("Deleting dbl_matrix (%ix%i)\n",nrow,ncol);
	delete [] m;
	if( oldMIndices != NULL ){
		delete [] oldMIndices;
	}
}


#endif
