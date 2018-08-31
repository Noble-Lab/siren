/**
 * \file matrix.h
 * AUTHOR: Alex Hu
 * CREATE DATE: March 5, 2015
 * PROJECT: DIA Regression
 * \brief Returns peptide abundances in DIA spectra that optimize
 * the bitonic regression objective function
 *****************************************************************************/
#ifndef MATRIX_H
#define MATRIX_H

#include "memory.h"
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

// Writes out a sparse form of a matrix
// Given a sparse form of a matrix
// nrow and ncol are the real values of the matrix
// maxcol is the number of columns in inxes and vals
void write_csr_matrix( int * ixes, double * vals, int nrow, int ncol, int maxcol, char const * outfile, double scalar = 1.0 ){
	FILE * file;
	file = fopen(outfile,"w");	
	fprintf(file,"%i\t%i\n",nrow,ncol);
	int ix;
	int cnum;
	int end;
	double val;
	for( int r=0; r<nrow; r++ ){
		ix = r*maxcol;
		end = ix+maxcol;
		cnum = ixes[ix];
		if(  cnum != -1 ){
			fprintf(file,">\t%i\n",r);
		}
		else{
			continue;
		}
		while( cnum != -1 && ix < end ){
			if( vals[ix] > 0 ){
				fprintf(file,"%i\t%f\n", cnum, vals[ix]*scalar );
			}
			ix++;
			cnum = ixes[ix];
		}
	}
	fclose(file);
}

// Writes out a sparse form of a matrix
// Given a sparse form of a matrix
// nrow and ncol are the real values of the matrix
// maxcol is the number of columns in inxes and vals
void write_csr_ixes( int * ixes, int nrow, int ncol, int maxcol, char const * outfile ){
	FILE * file;
	printf("Writing (%ix%i) maxcol: %i, %s\n", nrow, ncol, maxcol, outfile );
	file = fopen(outfile,"w");	
	fprintf(file,"%i\t%i\n",nrow,ncol);
	int ix;
	int cnum;
	int end;
	for( int r=0; r<nrow; r++ ){
		ix = r*maxcol;
		end = ix+maxcol;
		cnum = ixes[ix];
		if(  cnum != -1 ){
			fprintf(file,">\t%i\n",r);
		}
		else{
			continue;
		}
		while( cnum != -1 && ix < end ){
			fprintf(file,"%i\t1\n", cnum );
			ix++;
			cnum = ixes[ix];			
		}
	}
	fclose(file);
	printf("Wrote %s\n", outfile );
}

// Performs XB where Xcsr is a row-major sparse matrix
// version of X
// where Xcsr contains the values in X and nzx
// containes the indices of the columns of the values in Xcsr
// and maxnz is the maximum number of columns associated with X
// M is the rows of X, T is the columns of T, N is the actually columns 
// of X, and maxnz is the number of columns of Xcsr and nzx
// B is NxT
// product is M by T
double sparseMultiply( double * Xcsr, int * nzx, double * B, int M, int T, int maxnz, double * product ){
	int nzix;
	int n;
	int nzend;	
	int Bix;
	int pix;

	int size = M*T;
	for( int i=0; i<size; i++ ){
		product[i] = 0;
	}
	for( int m=0; m<M; m++ ){
		nzix = m*maxnz;
		n=nzx[nzix];
		nzend = nzix+maxnz;	
		while( n != -1 ){
			Bix = n*T;
			pix = m*T;
			for( int t=0; t<T; t++ ){
				if( B[Bix] != 0 ){
					product[pix] += B[Bix]*Xcsr[nzix];
				}
				Bix++;
				pix++;
			}
			nzix++;
			if( nzix == nzend ){ break; }
			n=nzx[nzix];
		}
	}
};

void print_stats( double * m, int nrow, int ncol, const char *s  ){

	double min = m[0];
	double max = m[0];
	int nnz = 0;

	int size = nrow*ncol;
	for( int i=0; i<size; i++ ){

		if (m[i] < min){
			min = m[i];
		}
		if (m[i] > max){
			max = m[i];
		}
		if(m[i] != 0){
			nnz ++;
		}
	}

	printf("%s: min: %f, max: %f, nnz: %i\n", s, min, max, nnz );
		
}

double * zeros( int nrow, int ncol, double val = 0 ){
	int size=nrow*ncol;
	double * zero = new double[size];
	if( val == 0 ){
		for( int i=0; i<size; i++ ){ zero[i] = 0.0; }
	}
	else{
		for( int i=0; i<size; i++ ){ zero[i] = val; }
	}
	return( zero );
}

int * zeros_int( int nrow, int ncol, int val=0 ){
	int size=nrow*ncol;
	int * zero = new int[size];
	if( val == 0 ){
		for( int i=0; i<size; i++ ){ zero[i] = 0; }
	}
	else{
		for( int i=0; i<size; i++ ){ zero[i] = val; }
	}
	return( zero );
}

int nnz( double *M, int nrow, int ncol ){
	int nonzero = 0;
	int size = nrow*ncol;
	for(int i=0; i<size; i++ ){
		if( M[i] != 0 ){ nonzero++; }	
	}
	return( nonzero );
}

int nnz( bool *M, int nrow, int ncol ){
	int nonzero = 0;
	int size = nrow*ncol;
	for(int i=0; i<size; i++ ){
		if( M[i] ){ nonzero++; }	
	}
	return( nonzero );
}

// square roots the magnitudes of the negative numbers too
void sqrt_mat( double*M, int nrow, int ncol ){
	int size = nrow*ncol;
	for( int i=0; i<size; i++ ){
		if( M[i] >= 0 ){
			M[i] = sqrt(M[i]);
		}
		else{
			M[i] = -1*sqrt(-1*M[i]);
		}
	}
}


bool DELETE = true;

double * copy( double * m, int nrow, int ncol ){
	int size = nrow*ncol;
	double * n = new double[size];
	for( int i=0; i<size; i++ ){ n[i] = m[i]; }
	return(n); 	
}

int * copy( int * m, int nrow, int ncol ){
	int size = nrow*ncol;
	int * n = new int[size];
	for( int i=0; i<size; i++ ){ n[i] = m[i]; }
	return(n); 	
}

bool * falses( int nrow, int ncol ){
	int size=nrow*ncol;
	bool * falses = new bool[size];
	for( int i=0; i<size; i++ ){ falses[i] = false; }
	return( falses );
}

void print_matrix( int *M, int nrow, int ncol, const char * s ){
	//return;	
	int offset = 0;
	printf("Matrix %s %ix%i\n", s, nrow, ncol );
	for( int r=0; r<nrow; r++ ){
		offset = r*ncol;
		for( int c=0; c<ncol; c++ ){
			printf("%i\t",M[offset+c]);
		}
		printf("\n");
	}
}

void print_matrix( double *M, int nrow, int ncol, const char * s ){
	//return;	
	int offset = 0;
	printf("Matrix %s %ix%i\n", s, nrow, ncol );
	for( int r=0; r<nrow; r++ ){
		offset = r*ncol;
		for( int c=0; c<ncol; c++ ){
			printf("%f\t",M[offset+c]);
		}
		printf("\n");
	}
}

// Returns the index of array m that contains the minimum value
int argmin( double * m, int nrow, int ncol ){
	int size = nrow*ncol;
	int mini=0;
	for( int i=0; i<size; i++ ){
		if( m[i] < m[mini] ){
			mini = i;
		}
	}		
	return mini;
}

void deldoublestarstar( double ** s, int size ){
	for( int i=0; i<size; i++ ){
		delete [] s[i];
	}
	delete [] s;
}

void write_matrix( double *M, int nrow, int ncol, char const * outfile ){
	FILE * file;
	file = fopen(outfile,"w");
	int offset = 0;
	for( int r=0; r<nrow; r++ ){
		offset = r*ncol;
		for( int c=0; c<ncol; c++ ){
			if( M[offset+c] == 0.0 ){
				fprintf(file,"%i\t", 0 );	
			}
			else{
				fprintf(file,"%f\t",M[offset+c]);
			}
		}
		fprintf(file,"\n");
	}
	fclose(file);
}

// The format of this is:
// Returns the percent sparsity of the matrix
// nrow\tncol\n
// rownum\tcolnum\tvalue\n
// ...


double write_sparse_matrix( double *M, int nrow, int ncol, char const * outfile ){
	FILE * file;
	file = fopen(outfile,"w");
	int offset=0;
	fprintf(file,"%i\t%i\n",nrow,ncol);
	int nnz = 0;
	int end;
	for( int r=0; r<nrow; r++ ){
		end = offset+ncol;
		for( int c=offset; c<end; c++ ){
			if( M[c] > 0.0 ){
				fprintf(file,"%i\t%i\t%f\n", r, c-offset, M[c] );	
				nnz += 1;
			}
		}
		offset = end;
	}
	fclose(file);
	return(100.0*nnz/nrow/ncol);
}

int write_csr_matrix( double *M, int nrow, int ncol, char const * outfile, double scalar = 1.0 ){
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
			if( M[c] != 0.0 ){
				fprintf(file,"%i\t%e\n", c-offset, M[c]*scalar );	
				nnz += 1;
			}
		}
		offset = end;
	}
	fclose(file);
	return(nnz);
}

int write_csr_matrix( int *M, int nrow, int ncol, char const * outfile, double scalar = 1.0 ){
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
			if( M[c] != 0 ){
				fprintf(file,"%i\t%i\n", c-offset, M[c]*scalar );	
				nnz += 1;
			}
		}
		offset = end;
	}
	fclose(file);
	return(nnz);
}

int write_csr_matrix( bool *M, int nrow, int ncol, char const * outfile, double scalar = 1.0 ){
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
			if( M[c] != 0){
				fprintf(file,"%i\t%i\n", c-offset, 1 );	
				nnz += 1;
			}
		}
		offset = end;
	}
	fclose(file);
	return(nnz);
}

int write_csr_matrix( double *M, int nrow, int ncol, char const * outfile, int * oldMIndices, double scalar = 1.0 ){
	FILE * file;
	file = fopen(outfile,"w");
	int offset=0;
	fprintf(file,"%i\t%i\n",nrow,ncol);
	int nnz = 0;
	int end;

	if( scalar != 1.0 ){
		for( int r=0; r<nrow; r++ ){
			end = offset+ncol;
			fprintf(file,">\t%i\n", oldMIndices[r] );	
			for( int c=offset; c<end; c++ ){
				if( M[c] != 0.0 ){
					fprintf(file,"%i\t%e\n", c-offset, M[c]*scalar );	
					nnz += 1;
				}
			}
			offset = end;
		}
	}
	else{
		for( int r=0; r<nrow; r++ ){
			end = offset+ncol;
			fprintf(file,">\t%i\n", oldMIndices[r] );	
			for( int c=offset; c<end; c++ ){
				if( M[c] != 0.0 ){
					fprintf(file,"%i\t%e\n", c-offset, M[c] );	
					nnz += 1;
				}
			}
			offset = end;
		}
	}

	fclose(file);
	return(nnz);
}


int * row_nnz( double *M, int nrow, int ncol ){

	int * nonzero = new int[ncol];
	int size = nrow*ncol;
	int c;
	for(int i=0; i<size; i++ ){
		if( M[i] != 0 ){ 
			c = i%ncol;
			nonzero[c]++; 
		}	
	}
	return( nonzero );
}

int * col_nnz( double *M, int nrow, int ncol ){
	int * nonzero = new int[nrow];
	int size = nrow*ncol;
	int r;
	for(int i=0; i<size; i++ ){
		if( M[i] != 0 ){ 
			r = i/ncol;
			nonzero[r]++; 
		}	
	}
	return( nonzero );
}

double * transpose( double *M, int nrow, int ncol ){
	double * T = zeros( nrow, ncol );
	for( int r=0; r<nrow; r++ ){
		for( int c=0; c<ncol; c++ ){
			T[c*nrow+r] = M[r*ncol+c];
		}
	}
	return(T);
}

bool * transpose_bool( double *M, int nrow, int ncol ){
	bool * T = new bool[nrow*ncol];
	for( int r=0; r<nrow; r++ ){
		for( int c=0; c<ncol; c++ ){
			if( M[r*ncol+c] != 0 ){
				T[c*nrow+r] = true;
			}
			else{
				T[c*nrow+r] = false;
			}
		}
	}
	return(T);
}

bool * transpose( bool * M, int nrow, int ncol ){
	bool * T = new bool[nrow*ncol];
	int Mix;
	for( int r=0; r<nrow; r++ ){
		Mix = r*ncol;
		for( int c=0; c<ncol; c++ ){
			T[c*nrow+r] = M[Mix+c];
		}
	}
	return(T);
}

// Sets all the values in matrix M to zero
void zero( int * M, int nrow, int ncol, int val=0 ){
	int size = nrow*ncol;
	for( int i=0; i<size; i++ ){ M[i] = val; }
}

// Sets all the values in matrix M to zero
void zero( double * M, int nrow, int ncol, double val=0.0 ){
	int size = nrow*ncol;
	for( int i=0; i<size; i++ ){ M[i] = val; }
}

void transpose( double *M, int nrow, int ncol, double * T ){
	for( int r=0; r<nrow; r++ ){
		for( int c=0; c<ncol; c++ ){
			T[c*nrow+r] = M[r*ncol+c];
		}
	}
}

// Here, acols == brows
double * multiply( double*A, int arows, int acols, double*B, int brows, int bcols, double s = 1.0 ){

	double sum=0;
	int size = arows*bcols;
	double * P = new double[size];
	for( int i=0; i<arows; i++ ){
		for( int j=0; j<bcols; j++ ){
			sum=0;
			for( int k=0; k<acols; k++ ){
					sum += A[i*acols+k]*B[k*bcols+j];
			}
			P[i*bcols+j] = sum;
		}
	}
	if( s != 1.0 ){
		for( int i=0; i<size; i++ ){
			P[i] *= s;
		}
	}
	return(P);
}

// Here, acols == brows
double * multiply_csr( double*Acsr, int * nza,  int arows, int acols, double*B, int brows, int bcols, double * P ){
	printf("multiply_csr...\n");
	printf("arows:%i, acols:%i, brows:%i, bcols:%i\n", arows, acols, brows, bcols);
	double * sum;
	int nzix;
	int aix;
	int br;
//	print_matrix( nza, arows, acols, "nza");
	// For each row of a, ra
	for( int ra=0; ra<arows; ra++ ){	
//		printf("ra:%i\n",ra);
		// For each column of b, cb
		for( int cb=0; cb<bcols; cb++ ){
//			printf("-cb:%i\n",cb);
			sum = &P[ra*bcols+cb];
			*sum=0;
			aix=ra*acols; // Index of first (ra,0) in Acsr
			br=nza[aix]; // First non-zero column of A (dense), 
			// which is the row of b

			// For each row of a and column of b
			while( br != -1 ){
//				printf("--aix:%i, br:%i\n", aix, br);
				*sum += Acsr[aix]*B[br*bcols+cb];
				aix++;
				br=nza[aix];
			}	
		}
	}
	printf("done.\n");
}

void scale( double* M, int nrows, int ncols, double scalar ){
	int size= nrows*ncols;
	for( int i=0; i<size; i++ ){
		if( M[i] != 0 ){
			M[i] = M[i]*scalar;
		}
	}
}


void scale_col( double* M, int nrows, int ncols, int colnum, double scalar ){

	int size = nrows*ncols;
	while( colnum < size ){
		if( M[colnum] != 0 ){
			M[colnum] = M[colnum]*scalar;
		}
		colnum += ncols;
	}
}

// Here, acols == brows
void multiply( double*A, int arows, int acols, double*B, int brows, int bcols, double s, double * P ){
	double sum;
	for( int i=0; i<arows; i++ ){
		for( int j=0; j<bcols; j++ ){
			sum=0;
			for( int k=0; k<acols; k++ ){
				sum += A[i*acols+k]*B[k*bcols+j];
			}
			P[i*bcols+j] = s*sum;
		}
	}
}

// Here, arows == brows
// Returns At x B
//double sparseMultiply( double * Xcsr, int * nzx, double * B, int M, int T, int maxnz, double * product ){

void sparseTmultiply( double * Xcsr, int * nzx, double * Y , int M, int T, int N, int maxnz, double scalar, double * P ){
        int nzix;
        int n;
        int nzend;
        int Yix;
        int pix;

	printf("This is happening\n");

	int NT = N*T;
	for( int i=0; i<NT; i++ ){
		P[i] = 0;
	}

	for( int m=0; m<M; m++ ){
                nzix = m*maxnz;
                n=nzx[nzix];
                nzend = nzix+maxnz;
		while( n != -1 ){
			pix = n*T;
			Yix = m*T;
			for( int t=0; t<T; t++ ){	
				if( Y[Yix] != 0 ){
					P[pix] += scalar*Y[Yix]*Xcsr[nzix];	
				}
				pix++;
				Yix++;
			}
                        nzix++;
                        if( nzix == nzend ){ break; }
                        n=nzx[nzix];
		}
	}
	printf("SPARSETMULTIPLYDONE!\n");
}

// Here, arows == brows
// Returns At x B
//double sparseMultiply( double * Xcsr, int * nzx, double * B, int M, int T, int maxnz, double * product ){
/*
void sparseTmultiplyX( double * Xcsr, int * nzx, int M, int N, int maxnz, double scalar, double * P ){
        int nzix;
        int n;
        int nzend;
	int nzix2;
	int n2;
	int row;

	printf("This is happening\n");

	int NN = N*N;
	for( int i=0; i<NN; i++ ){
		P[i] = 0;
	}

	for( int m=0; m<M; m++ ){
                nzix = m*maxnz;
                n=nzx[nzix];
                nzend = nzix+maxnz;
		row = n*N;
		while( n != -1 ){
			nzix2 = nzix + 1;
			n2 = nzx[nzix2]
			while( n2 != -1 ){
				row = 					
				nzix2++;
				n2 = nzx[nzix];
				P[row+n2] += scalar*Xcsr[nzix]*Xcsr[nzix2];
			}	

			for( int t=0; t<T; t++ ){	
				if( Y[Yix] != 0 ){
					P[pix] += scalar*Y[Yix]*Xcsr[nzix];	
				}
				pix++;
				Yix++;
			}
                        nzix++;
                        if( nzix == nzend ){ break; }
                        n=nzx[nzix];
		}
	}
	printf("SPARSETMULTIPLYDONE!\n");
}
*/
// Here, arows == brows
// Returns At x B
//double sparseMultiply( double * Xcsr, int * nzx, double * B, int M, int T, int maxnz, double * product ){
double * sparseTmultiply( double * Xcsr, int * nzx, double * Y , int M, int T, int N, int maxnz, double scalar ){
	double * P = zeros( N, T );
	sparseTmultiply( Xcsr, nzx, Y , M,  T, N, maxnz, scalar, P );
	return(P);
}

// Here, arows == brows
// Returns At x B
double * Tmultiply( double*A, int arows, int acols, double*B, int brows, int bcols, double scalar ){

	double * P = zeros( acols, bcols );
	for( int m=0; m<arows; m++ ){
		for( int t=0; t<bcols; t++ ){	
			for( int n=0; n<acols; n++ ){
				if( A[m*acols+n] == 0 || B[m*bcols+t] == 0 ){ continue; }
				P[n*bcols+t] += scalar*A[m*acols+n]*B[m*bcols+t];	
			}
		}
	}
	return(P);
}

// Here, arows == brows
// Returns At x B
void * Tmultiply( double*A, int arows, int acols, double*B, int brows, int bcols, double scalar, double * P ){
	
	int size = acols*bcols;
	for( int i=0; i<size; i++ ){ P[i] = 0; }
	for( int m=0; m<arows; m++ ){
		for( int t=0; t<bcols; t++ ){	
			for( int n=0; n<acols; n++ ){
				P[n*bcols+t] += scalar*A[m*acols+n]*B[m*bcols+t];	
				if( isnan( P[n*bcols+t] ) ){
					printf("NAN!: (%i,%i, %i): %f, %f\n", m, n, t, A[m*acols+n]*B[m*bcols+t] );
				}
			}
		}
	}
}

double * add( double*A, double*B, int nrows, int ncols ){

	double * P = zeros( nrows, ncols );
	int size = nrows*ncols;
	for( int i=0; i<size; i++ ){
		P[i] = A[i] + B[i];
	}
	return(P);
}

double add( double*A, double*B, int nrows, int ncols, double * P ){
	int size = nrows*ncols;
	for( int i=0; i<size; i++ ){
		P[i] = A[i] + B[i];
	}
}

// Same as add except allows weights to scale each operand
// It also returns the unscaled l1 norm of the diff between A and B!
// Also do not allow any value of P to be less than 0!
double linearcombo( double*A, double*B, int nrows, int ncols, double * P, double cA, double cB ){
	int size = nrows*ncols;
	double l1norm=0;
	double diff;
	for( int i=0; i<size; i++ ){
		diff = A[i] - B[i];
		if( diff > 0 ){ l1norm += diff; }
		else{ l1norm -= diff; }
		P[i] = cA*A[i] + cB*B[i];
		if( P[i] < 0 ){
			P[i] = 0;
		}
	}
	return(l1norm);
}


void  subtract( double*A, double*B, int nrows, int ncols, double * P ){
	int size = nrows*ncols;
	for( int i=0; i<size; i++ ){
		P[i] = A[i]-B[i];
	}
}


void softthresh( double*A, int nrows, int ncols, double val ){
	int size = nrows*ncols;
	for( int i=0; i<size; i++ ){
		if( A[i] < val ){
			A[i] = 0;
		}
		else{
			A[i] = A[i]-val;
		}
	}
}


// Here, acols == bcols
// Returns A x Bt
double * multiplyT( double*A, int arows, int acols, double*B, int brows, int bcols, double scalar ){

	double * P = zeros( arows, brows );
	for( int m=0; m<arows; m++ ){
		for( int t=0; t<brows; t++ ){	
			for( int n=0; n<acols; n++ ){
				if( A[m*acols+n] != 0.0 && B[t*bcols+n] != 0.0  ){
					P[m*brows+t] += scalar*A[m*acols+n]*B[t*bcols+n];	
				}
			}
		}
	}
	return(P);
}

double * multiplyT( double*A, int arows, int acols, double*B, int brows, int bcols, double scalar, double * P ){

	for( int m=0; m<arows; m++ ){
		for( int t=0; t<brows; t++ ){
			P[m*brows+t] = 0;	
			for( int n=0; n<acols; n++ ){
				if( A[m*acols+n] != 0.0 && B[t*bcols+n] != 0.0  ){
					P[m*brows+t] += scalar*A[m*acols+n]*B[t*bcols+n];	
				}
				//P[m*brows+t] += scalar*A[m*acols+n]*B[t*bcols+n];	
			}
		}
	}
}

int maxint( int * intarray, int size ){
	if( size < 1 ){
		return 0;
	}
	int mint = intarray[0];
	for( int i=1; i<size; i++ ){
		if( intarray[i] > mint ){
			mint = intarray[i];
		}
	}	
	return( mint );
}

// This file contains the matrix in sparse form: i j value
// The first line contains the dimensions of the matrix: nrow ncol
double* read_dblmatrix( char const * mfile, int * nrow, int * ncol, int maxcols ){ 

	ifstream file( mfile );
	string line;
	getline( file, line);
	istringstream iss(line);
	iss >> *nrow; iss >> *ncol;
	if( maxcols > 0 && *ncol > maxcols ){
		*ncol = maxcols;
	} 	

	double* m = zeros((*nrow),(*ncol));
	int i;
	int j;
	double v;
	while(file.good()){
		getline( file, line);
		istringstream s(line);
		s >> i;
		s >> j;
		if( i < *nrow && j < *ncol ){
			s >> v;
			m[i*(*ncol)+j] += v;
		}
	}
	return(m);
}


// This file contains the matrix in sparse form: i j value
// The first line contains the dimensions of the matrix: nrow ncol
// The line and column numbers start at 1, not 0
double* read_dblmatrix_indexat1( char const * mfile, int * nrow, int * ncol, int maxcols ){ 

	ifstream file( mfile );
	string line;
	getline( file, line);
	istringstream iss(line);
	iss >> *nrow; iss >> *ncol;
	if( maxcols > 0 && *ncol > maxcols ){
		*ncol = maxcols;
	} 	

	double* m = zeros((*nrow),(*ncol));
	int i;
	int j;
	double v;
	while(file.good()){
		getline( file, line);
		istringstream s(line);
		s >> i;
		s >> j;
		i--; j--;
		if( i < *nrow && j < *ncol ){
			s >> v;
			m[i*(*ncol)+j] += v;
		}
	}
	return(m);
}

// This file contains the matrix in sparse form: i j value
// The first line contains the dimensions of the matrix: nrow ncol
double* read_rowmajor_dblmatrix( char const * mfile, int * nrow, int * ncol, int maxcols ){ 

	ifstream file( mfile );
	string line;
	getline( file, line);
	istringstream iss(line);
	iss >> *nrow; iss >> *ncol;
	if( maxcols > 0 && *ncol > maxcols ){
		*ncol = maxcols;
	} 	
	double* m = zeros((*nrow),(*ncol));
	int i;
	int j;
	double v;
	char meh[50];
	while(file.good()){
		getline( file, line);
		istringstream s(line);
		if( line[0] == '>' ){
			s >> meh;
			s >> meh;
			i = atoi(meh);	
			continue;	
		}
		s >> j;
		if( i < *nrow && j < *ncol ){
			s >> v;
			m[i*(*ncol)+j] += v;
		}
	}
	return(m);
}

// This file contains the matrix in sparse form: i j value
// The first line contains the dimensions of the matrix: nrow ncol
bool * read_rowmajor_bool( char const * mfile, int * nrow, int * ncol, int maxcols ){ 

	ifstream file( mfile );
	string line;
	getline( file, line);
	istringstream iss(line);
	iss >> *nrow; iss >> *ncol;
	if( maxcols > 0 && *ncol > maxcols ){
		*ncol = maxcols;
	} 	
	bool * m = falses( *nrow, *ncol );
	int i;
	int j;
	char meh[50];
	while(file.good()){
		getline( file, line);
		istringstream s(line);
		if( line[0] == '>' ){
			s >> meh;
			s >> meh;
			i = atoi(meh);	
			continue;	
		}
		s >> j;
		if( i < *nrow && j < *ncol ){
			m[i*(*ncol)+j]=true;
		}
	}
	return(m);
}

// This file contains the matrix in sparse form: i j value
// The first line contains the dimensions of the matrix: nrow ncol
double* read_rowmajor_dblmatrix( char const * mfile, int * nrow, int * ncol, int maxrows, int maxcols, bool transpose=false ){ 

	ifstream file( mfile );
	string line;
	getline( file, line);
	istringstream iss(line);

	if( transpose ){
		iss >> *ncol; iss >> *nrow;
	}	
	else{
		iss >> *nrow; iss >> *ncol;
	}

	if( maxcols > 0 && *ncol > maxcols ){
		*ncol = maxcols;
	} 
	if( maxrows > 0 && *nrow > maxrows ){
		*nrow = maxrows;
	} 
	double* m = zeros((*nrow),(*ncol));
	double v;
	int i;
	int j;
	char meh[50];

	if( transpose ){
		while(file.good()){
			getline( file, line);
			istringstream s(line);
			if( line[0] == '>' ){
				s >> meh;
				s >> meh;
				i = atoi(meh);	
				continue;	
			}
			s >> j;
			if( j < *nrow && i < *ncol ){
				s >> v;
				m[j*(*ncol)+i] += v;
			}
		}
	}
	else{
		while(file.good()){
			getline( file, line);
			istringstream s(line);
			if( line[0] == '>' ){
				s >> meh;
				s >> meh;
				i = atoi(meh);	
				continue;	
			}
			s >> j;
			if( i < *nrow && j < *ncol ){
				s >> v;
				m[i*(*ncol)+j] += v;
			}
		}
	}
	
	file.close();
	return(m);
}

// This file contains the matrix in sparse form: i j value
// The first line contains the dimensions of the matrix: nrow ncol
// acol must equal bcol

// This file contains the matrix in dense form
double* read_dblmatrix_dense( char const * mfile, int * nrow, int * ncol ){ 

	ifstream file( mfile );
	string line;
	
	getline( file, line);
        istringstream iss(line);
        iss >> *nrow; iss >> *ncol;
	
	int size = (*nrow)*(*ncol);	
	double* m = new double[size];
	int i=0;
	while(file.good()){
		getline( file, line);
		istringstream s(line);
		for(int c=0; c<*ncol; c++ ){ 
			if( i == size ){ break; }	
			s >> m[i];
			i++;
		}
	}
	return(m);
}

double norm1( double*M, int rows, int cols ){
	int size = rows*cols;
	float norm1=0;
	for(int i=0; i<size; i++){
		if( M[i] > 0 ){
			norm1 += M[i];
		}
		else if( M[i] < 0 ){
			norm1 -= M[i];
		}
	}
	return( norm1 );
}

double *compareM;

int compare( const void *a, const void *b ){
	int x = *(int*)a;
	int y = *(int*)b;
	return( compareM[x] > compareM[y] );
}

double max( double* M, int nrow, int ncol ){
//	printf("computing max...\n");
	int size=nrow*ncol;
	double max = M[0];
	for( int i=0; i<size; i++ ){
		if( max < M[i] ){ max=M[i]; }
	}
//	printf("Max = %f\n", max);
	return max;
}

double mean( double*M, int nrow, int ncol ){

	int size = nrow*ncol;
	double sum = 0.0;
	for( int i=0; i<size; i++ ){
		sum = sum+M[i];
	}		
	return( sum/size );
}

// Returns the value at the percentile highest
double percentile( double *M, int nrow, int ncol, double percentile, bool use_mean=false ){

	if( use_mean ){
		return( mean( M, nrow, ncol ) );
	}

	int size = nrow*ncol;
	//print_matrix( M, nrow, ncol, "Percentile M");
	int * indices = new int[size];
	compareM = M;
	percentile = percentile/100.0;
	for(int i=0; i<size; i++ ){ indices[i] = i; }
	// Sort indices...
	printf("Sorting...\n");
	qsort( &indices[0] , size, sizeof(int), compare);
	int index = size*percentile;
	printf( "%i index M: %f\n", indices[0], M[indices[0]] );
//	printf("Percentile: %f\tIndex: %i\tTotal: %i\n", percentile, index, size );
	if( index == size ){ --index; }
	index = indices[index];
	delete [] indices;
	return( M[index] );	
}

void subtract( double *M, int nrow, int ncol, double val ){
	int size = nrow*ncol;
	for( int i=0; i<size; i++ ){
		M[i] -= val;
	}
}

double * subtract( double *A, double* B, int nrow, int ncol ){
	int size = nrow*ncol;
	double * diff = new double[size];
	for( int i=0; i<size; i++ ){
		diff[i] -= A[i]-B[i];
	}
	return( diff );
}

void empty( double * A, int nrow, int ncol ){
	int size = nrow*ncol;
	for( int i=0; i<size; i++ ){ A[i] = 0.0; }
}

bool MAGNORM = false;
/*
 * This is a slow version of the method! Below is a faster version maybe
void col_magnorm( double * M, int nrow, int ncol ){
	// Each column is normalized by a scalar such that their magnitudes all equal 1
	double scalar;
	int i;
	for( int c=0; c<ncol; c++ ){
		scalar = 0;
		for( int r=0; r<nrow; r++ ){
			i = r*ncol + c;
			scalar += M[i]*M[i];
		}
		scalar = 1.0/scalar;
		scalar = sqrt(scalar);
		for( int r=0; r<nrow; r++ ){
			i = r*ncol + c;
			M[i] = scalar*M[i];
		}
	}

}*/

double * column_magnitude( double * M, int nrow, int ncol ){
	double * magnitude = zeros( ncol, 1 );
	int i=0;
	// Compute the magnitudes for every column
	for( int r=0; r<nrow; r++ ){
		for( int c=0; c<ncol; c++ ){
			magnitude[c] += M[i]*M[i];
			i++;
		}
	}
	for( int c=0; c<ncol; c++ ){
		magnitude[c] = sqrt(magnitude[c]);
	}
	return(magnitude);
}

// Each column is normalized by a scalar such that their magnitudes all equal 1
void col_magnorm( double * M, int nrow, int ncol ){
	double * scalar = column_magnitude(M,nrow,ncol);
	for( int c=0; c<ncol; c++ ){
		scalar[c] = 1.0/scalar[c];
	}
	// Scale every value based on the scalar!
	int i=0;
	for( int r=0; r<nrow; r++ ){
		for( int c=0; c<ncol; c++ ){
			M[i] *= scalar[c];
			i++;
		}
	}
	delete [] scalar;
}

// Each column is normalized by a scalar such that their magnitudes all equal 1
void col_squaremagnorm( double * M, int nrow, int ncol ){
	double * scalar = zeros( ncol, 1 );
	int i=0;
	// Compute the magnitudes for every column
	for( int r=0; r<nrow; r++ ){
		for( int c=0; c<ncol; c++ ){
			scalar[c] += M[i]*M[i];
			i++;
		}
	}
	// Compute the scalar for every column with which to multiply the 
	// values in the columns
	for( int c=0; c<ncol; c++ ){
		scalar[c] = 1.0/scalar[c];
	}
	// Scale every value based on the scalar!
	i=0;
	for( int r=0; r<nrow; r++ ){
		for( int c=0; c<ncol; c++ ){
			M[i] *= scalar[c];
			i++;
		}
	}
	delete [] scalar;
}

void hardthresh( double * M, int nrow, int ncol, double thresh ){
	int size = nrow*ncol;
	for( int i=0; i<size; i++ ){
		if( M[i] != 0 && M[i] < thresh ){
			M[i] = 0;
		}
	}
}


void checknan( double * M, int nrow, int ncol ){
	int size = nrow*ncol;
	int r;
	int c;
	for( int i=0; i<size; i++ ){
		if( isnan(M[i]) ){
			r = i/ncol;
			c = i-(ncol*r);
			printf("(%i,%i) = nan\n",r,c);		
		}
	}
}

int * mincolperrow( double *M, int nrow, int ncol ){

	int * mcpr = new int[nrow];
        int begin;
        int end;
	double min;
	for( int r=0; r<nrow; r++ ){
		begin = r*ncol;
		end = begin+ncol;
		min = M[begin];
		mcpr[r] = 0;
		for( int c=(begin+1); c<end; c++ ){
			if( M[c] < min ){
				min = M[c];
				mcpr[r] = c-begin;
			}
		}		
	}
	return(mcpr);

}

// acols must equal bcols
double * rbind( double * A, int arows, int acols, double * B, int brows, int bcols ){

	int sizeA = arows*acols;
	int sizeB = brows*bcols;
	int sizeC = sizeA+sizeB;

	double * C = zeros(arows+brows,bcols);
	for( int i=0; i<sizeA; i++ ){
		C[i] = A[i];
	}
	int cix = sizeA;
	for( int i=0; i<sizeB; i++ ){
		C[cix] = B[i];
		cix++;
	}
	return(C);
}
// Scales each column so that it sums to 1
void colnorm( double * M, int nrows, int ncols ){

	int size = nrows*ncols;
	double * sums = zeros(1,ncols);	
	int offset;

	for( int r=0; r<nrows; r++ ){
		offset = r*ncols;
		for( int c=0; c<ncols; c++ ){
			if( M[offset] != 0 ){
				sums[c] += M[offset];
			}
			offset++;	
		}
	}
	
	for( int r=0; r<nrows; r++ ){
		offset = r*ncols;
		for( int c=0; c<ncols; c++ ){
			if( M[offset] != 0 ){
				M[offset] = M[offset]/sums[c];
			}
			offset++;	
		}
	}
}

// Scales each column so that it sums to 1
void rownorm( double * M, int nrows, int ncols ){

	int size = nrows*ncols;
	double sum;	
	int offset;

	for( int r=0; r<nrows; r++ ){
		offset = r*ncols;
		sum = 0.0;
		for( int c=0; c<ncols; c++ ){
			if( M[offset] != 0 ){
				sum += M[offset];
			}
			offset++;	
		}
		if( sum == 0.0 ){
			continue;
		}
		offset = r*ncols;
		for( int c=0; c<ncols; c++ ){
			if( M[offset] != 0 ){
				M[offset] /= sum;
			}
			offset++;	
		}
	}
	
}


double norm2( double * m, int nrow, int ncol ){

	double n2=0;
	int size = nrow*ncol;
	for( int i=0; i<size; i++ ){
		if( m[i] != 0 ){ n2 += m[i]*m[i]; }
	}
	return n2;
}

// Normalizes the columns in m to be the same magnitude as the columns in standard
void magnorm( double * standard, int srows, int scols, double * M, int mrows, int mcols, double a=1.0 ){

	printf("Computing standard magnitudes...\n");
	// Compute the squared magnitude of each column in standard
	double * sumsquares = zeros(1,scols);
	int ix;
	for( int m=0; m<srows; m++ ){
		ix = m*scols;
		for( int t=0; t<scols; t++ ){
			sumsquares[t] += standard[ix]*standard[ix];
			ix++;
		}
	}
	printf("Computing m magnitudes...\n");
	// Compute the squared magnitude of each column in m
	double * msumsquares = zeros(1,mcols);
	for( int m=0; m<mrows; m++ ){
		ix = m*mcols;
		for( int t=0; t<mcols; t++ ){
			msumsquares[t] += M[ix]*M[ix];
			ix++;
		}
	}
	// Compute the scalar for each column
	printf("Computing scalars...\n");
	for( int t=0; t<mcols; t++ ){
		if( msumsquares[t] != 0 ){
			msumsquares[t] = sqrt(sumsquares[t]/msumsquares[t]);	
		}
	}
	// Multiply each element by the appropriate scalar
	printf("Scaling...\n");
	for( int m=0; m<mrows; m++ ){
		ix = m*mcols;
		for( int t=0; t<mcols; t++ ){
			M[ix] = M[ix]*msumsquares[t]*a;
			ix++; 
		}
	}
}

class csr {
	public:
		double * Mcsr;
		int * nzm;
		int nrow;
		int maxnz;
	
	csr( double * M, int nrowm, int ncol ){
	
		// Delete whatever's here if it's not null
	//	if( Mcsr != NULL){ delete [] Mcsr; }
	//	if( nzm != NULL){ delete [] nzm; }
	
		maxnz=0;
		int nz;
		int ix;
		nrow = nrowm;
		// Find the maximum number of non-zero columns per row
		for( int r=0; r<nrow; r++ ){
			nz = 0;
			for( int ix=r*ncol; ix<(r*(ncol+1)); ix++ ){
				if( M[ix] ){
					nz++;
				}
			}
			if( nz > maxnz ){
				maxnz = nz;
			}
		}
		// Create the matrix the stores the values
		printf("nrow=%i, ncol= %i, maxnz=%i\n", nrow, ncol, maxnz);
		nzm = zeros_int(nrow,maxnz); // For each row, stores the column number that was nonzero
		for( int i=0; i<(nrow*(maxnz)); i++ ){ nzm[i] = -1; } // Fill it with negative numbers
		Mcsr = zeros(nrow,maxnz); // For each row, stores the value that was nonzero

		// Fill these in!
		int nzix; //
		for( int r=0; r<nrow; r++ ){
			nzix = r*maxnz;
			ix = r*ncol;
			for( int c=0; c<ncol; c++ ){
				if( M[ix] ){
					nzm[nzix] = c;
					Mcsr[nzix] = M[ix];	
					nzix++;	
				}					
				ix++;
			}	
		}
	}

};

#endif
