/**
 * \file scorePeptides.cpp
 * AUTHOR: Alex Hu
 * CREATE DATE: March 7, 2015
 * PROJECT: DIA Regression
 * \brief Returns peptide abundances in DIA spectra that optimize
 * the bitonic regression objective function
 *****************************************************************************/
#include "memory.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

float theoretical_score = 0;

/* This came from http://www.johndcook.com/blog/cpp_phi/
 * It is NormalCDF(x,mean=0,var=1) 
*/ 
double phi(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
 
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
 
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
 
    return 0.5*(1.0 + sign*y);
}

// Computes the pearson correlation R
// between two vectors A and B of size N
double pearsonR( double * A, double * B, int N ){

	// Compute Means
	double muA = 0;
	double muB = 0;
	for( int i=0; i<N; i++ ){
		muA += A[i];
		muB += B[i];
	}	
	// Compute other values for r	
	double numerator=0;
	double denomA=0;
	double denomB=0;
	double diffA;
	double diffB;
	for( int i=0; i<N; i++ ){
		diffA = A[i]-muA;
		diffB = B[i]-muB;	
		numerator += diffA*diffB;
		denomA += diffA*diffA;
		denomB += diffB*diffB;
	}
	double r = numerator/sqrt(denomA)/sqrt(denomB); 
	return(r);
}

// Computes the p-value for the correlation between the top K peaks
// in the theoretical and observed spectra.
double theoretical_corrscore( double * Y, double * X, int N, int T, int * indices, int K, int start, int end, int n ){
	// The observed signal will be thes AUC for each of the
	// top K fragment ion chromatograms
	double * yvector = zeros(1,K);
	for( int t=start; t<=end; t++ ){
		for( int k=0; k<K; k++ ){
			yvector[k] += Y[indices[k]*T+t];
		}
	}
	// The theoretical signal is the intensities of the 
	// top K peaks.
	double * xvector = zeros(1,K);
	for( int k=0; k<K; k++ ){
		xvector[k] = X[indices[k]*N+n];
	}	
	double r = pearsonR( yvector, xvector, K );
	delete [] yvector;
	delete [] xvector;
	return(r);
}

void incorporate_theoretical_corrscore( double * npoints, double * meanr, double theoreticalr, int K ){
	// Weight the theoretical correlation as strongly
	printf("meanr: %f, theoreticalr: %f, theoretical_score: %f\n", meanr, theoreticalr, theoretical_score );
	*meanr = (*meanr+(theoreticalr*theoretical_score))/(1.0+theoretical_score);
	*npoints *= (1.0+theoretical_score);
	// Weight the theoretical correlation as strongly
	//meanr = (*meanr+theoreticalr)/2;
	//npoints *= 2;
	// Weight the theoretical correlation proportionately
	//*meanr = ( (*meanr)*(*npoints)+theoreticalr*K )/(*npoints+K);
	//*npoints += K;
}

/*
 * Copy and paste this thing to do theoretical scoring into a corr_score function
if( theoretical_score > 0 ){
	double tr = theoretical_corrscore( Y, X, N, T, indices, K, start, end, n );
	incorporate_theoretical_corrscore( &npoints, &meanr, tr, K );
}
*/

// This came from http://stackoverflow.com/questions/2328258/cumulative-normal-distribution-function-in-c-c
double phi2(double x)
{
  static const double RT2PI = sqrt(4.0*acos(0.0));

  static const double SPLIT = 7.07106781186547;

  static const double N0 = 220.206867912376;
  static const double N1 = 221.213596169931;
  static const double N2 = 112.079291497871;
  static const double N3 = 33.912866078383;
  static const double N4 = 6.37396220353165;
  static const double N5 = 0.700383064443688;
  static const double N6 = 3.52624965998911e-02;
  static const double M0 = 440.413735824752;
  static const double M1 = 793.826512519948;
  static const double M2 = 637.333633378831;
  static const double M3 = 296.564248779674;
  static const double M4 = 86.7807322029461;
  static const double M5 = 16.064177579207;
  static const double M6 = 1.75566716318264;
  static const double M7 = 8.83883476483184e-02;

  const double z = fabs(x);
  double c = 0.0;

  if(z<=37.0)
  {
    const double e = exp(-z*z/2.0);
    if(z<SPLIT)
    {
      const double n = (((((N6*z + N5)*z + N4)*z + N3)*z + N2)*z + N1)*z + N0;
      const double d = ((((((M7*z + M6)*z + M5)*z + M4)*z + M3)*z + M2)*z + M1)*z + M0;
      c = e*n/d;
    }
    else
    {
      const double f = z + 1.0/(z + 2.0/(z + 3.0/(z + 4.0/(z + 13.0/20.0))));
      c = e/(RT2PI*f);
    }
  }
  return x<=0.0 ? c : 1-c;
}

double ** auc_score_W( double* B, double *W, int N, int T, double thresh=0 ){ //This B should be NxT, not TxN!

	double ** scores_starts_ends_peaks = new double*[4];
	double * scores = new double[N];
	double * starts = new double[N];
	double * ends = new double[N];
	double * peaks = new double[N];
	scores_starts_ends_peaks[0] = scores;
	scores_starts_ends_peaks[1] = starts;
	scores_starts_ends_peaks[2] = ends;
	scores_starts_ends_peaks[3] = peaks;

	double score;
	double start;
	double max_score;
	double max_start;
	double max_end;
	int Boffset;
	int Bend;
	for( int n=0; n<N; n++ ){
		Boffset = n*T;
		Bend = Boffset+T;
		
		int maxix=Boffset;
		for( int Wix = Boffset; Wix < Bend; Wix ++ ){
			if( W[Wix] > W[maxix] ){
				maxix = Wix;
			}
		}
		peaks[n] = maxix-Boffset;
		
		// Start at the peak, then spread out from there.
		max_score = B[maxix];
		int Bix = maxix-1;
		while( Bix >= Boffset && B[Bix] > thresh ){
			max_score += B[Bix];
			Bix--;
		}
		starts[n] = Bix+1-Boffset;
		
		Bix = maxix+1;
		while( Bix < Bend && B[Bix] > thresh ){
			max_score += B[Bix];
			Bix++;
		}
		ends[n] = Bix-1-Boffset;	
		scores[n] = max_score;
	}
	//The first array is scores, the second array is 
	//elution peak start.
	return(scores_starts_ends_peaks);
}

double ** auc_score( double* B, int N, int T ){ //This B should be NxT, not TxN!

	double ** scores_starts_ends = new double*[3];
	double * scores = new double[N];
	double * starts = new double[N];
	double * ends = new double[N];
	scores_starts_ends[0] = scores;
	scores_starts_ends[1] = starts;
	scores_starts_ends[2] = ends;

	double score;
	double start;
	double max_score;
	double max_start;
	double max_end;
	int Boffset;
	int Bend;
	for( int n=0; n<N; n++ ){
		score = 0; start = 0; max_score = 0; max_start = 0; max_end = 0;
		Boffset = n*T;
		Bend = Boffset+T;
		for( int Bix=Boffset; Bix<Bend; Bix++ ){
			if( B[Bix] == 0 ){
				if( score > max_score ){
					max_score = score;
					max_start = start;
					max_end = Bix-Boffset-1; // the non-zero sequence is [start,end]
				}
				score = 0;
			}
			else{
				if( score == 0 ){
					start = Bix-Boffset;
				}
				score = score+B[Bix];
			}
		}
		if( score > max_score ){
			max_score = score;
			max_start = start;
			max_end = Bend;
		}
		scores[n] = max_score;
		starts[n] = max_start;	
		ends[n] = max_end;	
	}
	//The first array is scores, the second array is 
	//elution peak start.
	return(scores_starts_ends);
}

// Fills the array "indices" with the m/z bins of the
// top K peaks in peptide n
// spectrum is just an array of the same size as X to store intermediate data
void topKPeaks( double*X, int M, int N, int n, int K, int * indices, double * spectrum ){

	// Copy X[:,n] into spectrum
	// And find the first maximum value
	int m=0;
	int MN = N*M;
	double max=-1;
	for( int xix=n; xix<MN; xix+=N ){
		spectrum[m] = X[xix];
		if( spectrum[m] > max ){
			max = spectrum[m];
			indices[0] = m;
		}

		m++;
	}
	spectrum[indices[0]] = -1;
	
	for( int k=1; k<K; k++ ){
		max =-1;
		for(int m=0; m<M; m++){
			if( spectrum[m] > max || ( spectrum[m] == max && m > indices[k-1]) ){
				max = spectrum[m];
				indices[k] = m;
			}
		}
		spectrum[indices[k]] = -1;
	}
}

// Fills the array "indices" with the m/z bins of the
// top K peaks in peptide n
void topKPeaks_ms1( double*X, int M1, int M2, int N, int n, int K1, int K2, int * indices, double * spectrum ){
	// The first M1 rows are MS1 data, and
	// the next M2 rows are the MS2 data.
	int M = M1 + M2;

	double * X1 = X;
	double * X2 = &X[M1*N];

	int * indices1 = indices;
	int * indices2 = &indices[K1];	

	double * spectrum1 = spectrum;
	double * spectrum2 = &spectrum[M1];

	topKPeaks( X1, M1, N, n, K1, indices1, spectrum1 );
	topKPeaks( X2, M2, N, n, K2, indices2, spectrum2 );

	for( int k=K1; k<(K1+K2); k++ ){
		indices[k] = indices[k] + M1;	
	}

}


// This is a fisher transformation into a z-test  
// for pearson correlation r
// in observed matrix Y[(m1,m2),start,end]
// taken from wikipedia
// returned is r
double pairwise_corrscore( double*Y, int T, int m1, int m2, int start, int end ){
	//printf( "pairwise_corrscore (%i,%i)\n", start, end );
	int ix1 = m1*T+start;
	int ix2 = m2*T+start;
	int n = end-start+1;
	
	if( n < 4 ){ // Need at least 4 points to do the test
		return(0);
	}
	// Compute Means
	double mu1=0;
	double mu2=0;
	while( end >= start ){
		mu1 += Y[ix1];
		mu2 += Y[ix2];
		end--;
		ix1++;	
		ix2++;
	}
	mu1 = mu1/n;
	mu2 = mu2/n;

	// Compute other values for r	
	end = n+start-1;
	ix1 = m1*T+start;
	ix2 = m2*T+start;
	double numerator=0;
	double denom1=0;
	double denom2=0;
	double diff1;
	double diff2;
	while( end >= start ){
		diff1 = Y[ix1]-mu1;
		diff2 = Y[ix2]-mu2;	
		numerator += diff1*diff2;
		denom1 += diff1*diff1;
		denom2 += diff2*diff2;
		end--;
		ix1++;	
		ix2++;
	}
	double r = numerator/sqrt(denom1)/sqrt(denom2); 
	//printf("This is R: %f\n",r);
	if( isnan(r) ){
		end = n+start-1;
		printf("Infinite Correlation!\nn:%i, r: %f\n", n, r); 
		printf("m1: %i, m2: %i, start: %i, end: %i\n", m1, m2, start, end );	
		return(0); 	
	}
	if( r == 1.0 ){
		return(0.999);
	}
	return(r);
}

// This is a fisher transformation into a z-test  
// for pearson correlation r
// in observed matrix Y[(m1,m2),start,end]
// taken from wikipedia
// returned is r
//double pairwise_resid_corrscore( double*Y, int T, int m1, int m2, int start, int end ){
double pairwise_resid_corrscore( double* resid, double *X, double *Bt, int T, int M, int N, int m1, int m2, int pepnum, int start, int end ){
	//printf( "pairwise_corrscore (%i,%i)\n", start, end );
	int ix1 = m1*T+start;
	int ix2 = m2*T+start;
	int n = end-start+1;
	int Bix = pepnum*T+start;
	double X1 = X[m1*N+pepnum];
	double X2 = X[m2*N+pepnum];
	if( n < 4 ){ // Need at least 4 points to do the test
		return(0);
	}
	// Compute Means
	double mu1=0;
	double mu2=0;
	while( end >= start ){
		mu1 += resid[ix1]+Bt[Bix]*X1;
		mu2 += resid[ix2]+Bt[Bix]*X2;
//		mu1 += Bt[Bix]*X1;
//		mu2 += Bt[Bix]*X2;
//		if( resid[ix1] != 0 or resid[ix2] != 0 ){
//			printf("Something's up with resid %e, %e\n", resid[ix1], resid[ix2] );
//		}
		end--;
		Bix++;
		ix1++;	
		ix2++;
	}
	mu1 = mu1/n;
	mu2 = mu2/n;

	// Compute other values for r	
	end = n+start-1;
	ix1 = m1*T+start;
	ix2 = m2*T+start;
	double numerator=0;
	double denom1=0;
	double denom2=0;
	double diff1;
	double diff2;
	Bix = pepnum*T+start;
	while( end >= start ){
		diff1 = resid[ix1]+Bt[Bix]*X1-mu1;
		diff2 = resid[ix2]+Bt[Bix]*X2-mu2;	
		//diff1 = Bt[Bix]*X1-mu1;
		//diff2 = Bt[Bix]*X2-mu2;	
		numerator += diff1*diff2;
		denom1 += diff1*diff1;
		denom2 += diff2*diff2;
		end--;
		ix1++;	
		ix2++;
	}
	double r = numerator/sqrt(denom1)/sqrt(denom2); 
	
	//printf("This is R: %f\n",r);
	if( isnan(r) ){
		end = n+start-1;
		printf("Infinite Correlation!\nn:%i, r: %f\n", n, r); 
		printf("m1: %i, m2: %i, start: %i, end: %i\n", m1, m2, start, end );	
		return(0); 	
	}
	if( r == 1.0 ){

		return(0.999);

	}
	return(r);
}

double * get_resid( double*Y, double*X, double*Bt, int M, int T, int N ){
	double * resid = multiply( X, M, N, Bt, N, T );
	int end = M*T;
	for( int i=0; i<end; i++ ){
		resid[i] = Y[i] - resid[i];
		if( resid[i] < 0 ){
			resid[i] = 0;
		}
	}
	return( resid );	
}

//double ** corrscore_W_variable( double *Y, double *X, double *B, double* W, int M, int T, int N, int win, char* wintype, int K, int ms1size, int ms2size ){

// This B should also be NxT, not TxN!
//double ** corrscore( double *Y, double *X, double *B, int M, int T, int N, residual ){
double ** corrscore_W_backgroundsubtracted( double *Y, double *X, double *B, double* W, int M, int T, int N ){
	int K=min(10,M);
	printf("\n\nBACKGROUND SUBTRACTED!\n\n");
	// Get and save the AUC scores and starts	
	double ** auc = auc_score_W( B, W, N, T);		
	double ** corrscores = new double*[6];
	corrscores[0] = auc[0]; corrscores[1] = auc[1]; corrscores[2] = auc[2]; corrscores[3] = auc[3]; corrscores[4] = new double[N];
	corrscores[5] = new double[N];
	delete [] auc;

	int * indices = new int[K]; // These will contain the top K peaks
	double * spectrum = new double[M]; // This will contain a copy of a column of X
	int start;
	int end;
	int peak;
	int Boffset;
	double score;
	double z;
	double z_resid;
	double npoints;
	double meanr;
	double meanr_resid;
	double * resid = get_resid( Y, X, B, M, T, N );

	// for each X
	for( int n=0; n<N; n++ ){
		// Find the bounds of the elution profile
		start = corrscores[1][n];
		end = corrscores[2][n];
		peak = corrscores[3][n];
		if( start+7 < peak ){ start = peak-7; }
		if( end-7 > peak ){ end = peak+7; }
		Boffset = n*T;
		meanr = 0;
		meanr_resid= 0;

		topKPeaks( X, M, N, n, K, indices, spectrum ); // Find the top K peaks in spectrum X[:,n]

		for(int kix1=0; kix1<(K-1); kix1++ ){
			for( int kix2=kix1+1; kix2<K; kix2++ ){
				meanr += pairwise_corrscore( Y, T, indices[kix1], indices[kix2], start, end );
				meanr_resid += pairwise_resid_corrscore( resid, X, B, T, M, N, indices[kix1], indices[kix2], n, start, end );
			}
		}
		meanr = meanr/(K-1)/K*2;
		meanr_resid = meanr_resid/(K-1)/K*2;	
		npoints = (end-start+1)*(K-1)*K/2;
		if( theoretical_score > 0 ){
			double tr = theoretical_corrscore( Y, X, N, T, indices, K, start, end, n );
			incorporate_theoretical_corrscore( &npoints, &meanr, tr, K );
		}
		z = 0.5*log((1+meanr)/(1-meanr))*sqrt(npoints-3);
		z_resid = 0.5*log((1+meanr_resid)/(1-meanr_resid))*sqrt(npoints-3);
		if( isnan(z) ){ printf("Nan! r: %f, npoints: %i\n", meanr, npoints ); } 
		if( isnan(z_resid) ){ printf("Nan! r: %f, npoints: %i\n", meanr_resid, npoints ); } 
		corrscores[4][n] = z;
		corrscores[5][n] = z_resid;
	}
	delete [] resid;
	delete [] indices;
	delete [] spectrum;
	return(corrscores);
}

// This B should also be NxT, not TxN!
//double ** corrscore( double *Y, double *X, double *B, int M, int T, int N, residual ){
// Sets a minimum radius around the elution peak as minwin.
double ** corrscore_W_minwin( double *Y, double *X, double *B, double* W, int M, int T, int N, int minwin ){
	int K=min(10,M);

	// Get and save the AUC scores and starts	
	double ** auc = auc_score_W( B, W, N, T);		
	double ** corrscores = new double*[5];
	corrscores[0] = auc[0]; corrscores[1] = auc[1]; corrscores[2] = auc[2]; corrscores[3] = auc[3]; corrscores[4] = new double[N];
	delete [] auc;

	int * indices = new int[K]; // These will contain the top K peaks
	double * spectrum = new double[M]; // This will contain a copy of a column of X
	int start;
	int end;
	int peak;
	int Boffset;
	double score;
	double z;
	double npoints;
	double meanr;
	// for each X
	for( int n=0; n<N; n++ ){
		// Find the bounds of the elution profile
		start = corrscores[1][n];
		end = corrscores[2][n];
		peak = corrscores[3][n];
		if( start+7 < peak ){ start = peak-7; }
		if( end-7 > peak ){ end = peak+7; }
		if( peak-start < minwin ){ start = peak-minwin; }
		if( end-peak < minwin ){ end = peak+minwin; }
		Boffset = n*T;
		meanr = 0;
		topKPeaks( X, M, N, n, K, indices, spectrum ); // Find the top K peaks in spectrum X[:,n]
		for(int kix1=0; kix1<(K-1); kix1++ ){
			for( int kix2=kix1+1; kix2<K; kix2++ ){
				meanr += pairwise_corrscore( Y, T, indices[kix1], indices[kix2], start, end );
			}
		}
		meanr = meanr/(K-1)/K*2;
		npoints = (end-start+1)*(K-1)*K/2;
		if( theoretical_score > 0 ){
			double tr = theoretical_corrscore( Y, X, N, T, indices, K, start, end, n );
			incorporate_theoretical_corrscore( &npoints, &meanr, tr, K );
		}
		z = 0.5*log((1+meanr)/(1-meanr))*sqrt(npoints-3);
		if( isnan(z) ){ printf("Nan! r: %f, npoints: %i\n", meanr, npoints ); } 
		score = z;
		corrscores[4][n] = score;
	}
	delete [] indices;
	delete [] spectrum;
	return(corrscores);
}

// This B should also be NxT, not TxN!
//double ** corrscore( double *Y, double *X, double *B, int M, int T, int N, residual ){
// Sets a minimum radius around the elution peak as minwin.
double ** corrscore_W_fixedwin( double *Y, double *X, double *B, double* W, int M, int T, int N, int win ){
	int K=min(10,M);

	// Get and save the AUC scores and starts	
	double ** auc = auc_score_W( B, W, N, T);		
	double ** corrscores = new double*[5];
	corrscores[0] = auc[0]; corrscores[1] = auc[1]; corrscores[2] = auc[2]; corrscores[3] = auc[3]; corrscores[4] = new double[N];
	delete [] auc;

	int * indices = new int[K]; // These will contain the top K peaks
	double * spectrum = new double[M]; // This will contain a copy of a column of X
	int start;
	int end;
	int peak;
	int Boffset;
	double score;
	double z;
	double npoints;
	double meanr;
	// for each X
	for( int n=0; n<N; n++ ){
		// Find the bounds of the elution profile
		peak = corrscores[3][n];
		start = max(0,peak-win);
		end = min(T,peak+win);
		Boffset = n*T;
		meanr = 0;
		topKPeaks( X, M, N, n, K, indices, spectrum ); // Find the top K peaks in spectrum X[:,n]
		for(int kix1=0; kix1<(K-1); kix1++ ){
			for( int kix2=kix1+1; kix2<K; kix2++ ){
				meanr += pairwise_corrscore( Y, T, indices[kix1], indices[kix2], start, end );
			}
		}
		meanr = meanr/(K-1)/K*2;
		npoints = (end-start+1)*(K-1)*K/2;
		if( theoretical_score > 0 ){
			double tr = theoretical_corrscore( Y, X, N, T, indices, K, start, end, n );
			incorporate_theoretical_corrscore( &npoints, &meanr, tr, K );
		}
		z = 0.5*log((1+meanr)/(1-meanr))*sqrt(npoints-3);
		if( isnan(z) ){ printf("Nan! r: %f, npoints: %i\n", meanr, npoints ); } 
		score = z;
		corrscores[4][n] = score;
	}
	delete [] indices;
	delete [] spectrum;
	return(corrscores);
}

// This B should also be NxT, not TxN!
//double ** corrscore( double *Y, double *X, double *B, int M, int T, int N, residual ){
double ** corrscore_W( double *Y, double *X, double *B, double* W, int M, int T, int N ){
	int K=min(10,M);

	// Get and save the AUC scores and starts	
	double ** auc = auc_score_W( B, W, N, T);		
	double ** corrscores = new double*[5];
	corrscores[0] = auc[0]; corrscores[1] = auc[1]; corrscores[2] = auc[2]; corrscores[3] = auc[3]; corrscores[4] = new double[N];
	delete [] auc;

	int * indices = new int[K]; // These will contain the top K peaks
	double * spectrum = new double[M]; // This will contain a copy of a column of X
	int start;
	int end;
	int peak;
	int Boffset;
	double score;
	double z;
	double npoints;
	double meanr;
	// for each X
	for( int n=0; n<N; n++ ){
		// Find the bounds of the elution profile
		start = corrscores[1][n];
		end = corrscores[2][n];
		peak = corrscores[3][n];
		if( start+7 < peak ){ start = peak-7; } //Maximum window size is 7
		if( end-7 > peak ){ end = peak+7; } // Maximum window size is 7
		Boffset = n*T;
		meanr = 0;
		topKPeaks( X, M, N, n, K, indices, spectrum ); // Find the top K peaks in spectrum X[:,n]
		for(int kix1=0; kix1<(K-1); kix1++ ){
			for( int kix2=kix1+1; kix2<K; kix2++ ){
				meanr += pairwise_corrscore( Y, T, indices[kix1], indices[kix2], start, end );
			}
		}
		meanr = meanr/(K-1)/K*2;
		npoints = (end-start+1)*(K-1)*K/2;
		if( theoretical_score > 0 ){
			double tr = theoretical_corrscore( Y, X, N, T, indices, K, start, end, n );
			incorporate_theoretical_corrscore( &npoints, &meanr, tr, K );
		}
		z = 0.5*log((1+meanr)/(1-meanr))*sqrt(npoints-3);
		if( isnan(z) ){ printf("Nan! r: %f, npoints: %i\n", meanr, npoints ); } 
		score = z;
		corrscores[4][n] = score;
	}
	delete [] indices;
	delete [] spectrum;
	return(corrscores);
}


void write_scores( double**score, int numscores, int numpeps, const char * outfile ){
	FILE * file;
        file = fopen(outfile,"w");
        int s;
	if( numscores == 3 ){
		fprintf(file,"AUC\tstart\tend\n");
	}
	else if( numscores == 4 ){
		fprintf(file,"AUC\tstart\tend\tcorrelation\n");
	}
	else if( numscores == 5 ){
		fprintf(file,"AUC\tstart\tend\tpeak\tcorrelation\n");
	}
	else{
		fprintf(file,"AUC\tstart\tend\tpeak\tcorrelation\tbackground corrected correlation");
		for( int i=6; i<numscores; i++ ){
			fprintf(file,"\tfield%i", i );
		}
		fprintf(file,"\n");
	}
	for( int p=0; p<numpeps; p++ ){
		fprintf(file,"%f", score[0][p]);
		s=1;
		while( s < numscores ){
			if( s == 1 || s == 2 ){
				fprintf(file,"\t%i",int(score[s][p]));	
			}
			else if( numscores >= 5 && s == 3 ){
				fprintf(file,"\t%i",int(score[s][p]));	
			}
			else{
				fprintf(file,"\t%f",score[s][p]);	
			}
			s++;
		}
		fprintf(file,"\n");
	}
        fclose(file);
}



// This B should also be NxT, not TxN!
// Sets a fixed radius around the elution peak as fixedwin.
double ** corrscore_W_variable( double *Y, double *X, double *B, double* W, int M, int T, int N, int win, char* wintype, int K, int ms1size, int ms2size, bool background_subtracted=false ){
	K=min(K,M);

	// Get and save the AUC scores and starts	
	double ** auc = auc_score_W( B, W, N, T);		
	double ** corrscores = new double*[5];
	corrscores[0] = auc[0]; corrscores[1] = auc[1]; corrscores[2] = auc[2]; corrscores[3] = auc[3]; corrscores[4] = new double[N];
	delete [] auc;

	int K1 = 0;
	if( ms1size > 0 ){
		K1 = 3;
	}
	int K2 = K;
	K = K1 + K2;

	int * indices = new int[K]; // These will contain the top K peaks
	double * spectrum = new double[M]; // This will contain a copy of a column of X
	int start;
	int end;
	int peak;
	int Boffset;
	double score;
	double z;
	double z_resid;
	double npoints;
	double meanr;
	double meanr_resid ;

	double * resid;
	if( background_subtracted ){
		 resid = get_resid( Y, X, B, M, T, N );
	}
	// for each X
	for( int n=0; n<N; n++ ){
		// Find the bounds of the elution profile
		peak = corrscores[3][n];

		if( strcmp( "fixed", wintype ) == 0 ){
			start = max(0,peak-win);
			end = min(T,peak+win);
		}	
		else if( strcmp( "minimum", wintype ) == 0 ){
			start = corrscores[1][n];
			end = corrscores[2][n];
			if( start+7 < peak ){ start = peak-7; }
			if( end-7 > peak ){ end = peak+7; }
			if( peak-start < win ){ start = peak-win; }
			if( end-peak < win ){ end = peak+win; }
		}	
		else{
			start = corrscores[1][n];
			end = corrscores[2][n];
		}

		Boffset = n*T;
		meanr = 0;
		if( ms1size > 0 ){
			topKPeaks_ms1( X, ms1size, ms2size, N, n, K1, K2, indices, spectrum );
		}
		else{
			topKPeaks( X, M, N, n, K, indices, spectrum ); // Find the top K peaks in spectrum X[:,n]
		}
		for(int kix1=0; kix1<(K-1); kix1++ ){
			for( int kix2=kix1+1; kix2<K; kix2++ ){
				meanr += pairwise_corrscore( Y, T, indices[kix1], indices[kix2], start, end );
			}
		}
		meanr = meanr/(K-1)/K*2;
		npoints = (end-start+1)*(K-1)*K/2;
		if( theoretical_score > 0 ){
			double tr = theoretical_corrscore( Y, X, N, T, indices, K, start, end, n );
			incorporate_theoretical_corrscore( &npoints, &meanr, tr, K );
		}
		z = 0.5*log((1+meanr)/(1-meanr))*sqrt(npoints-3);
		if( isnan(z) ){ printf("Nan! r: %f, npoints: %i\n", meanr, npoints ); } 
		score = z;
		corrscores[4][n] = score;
	}
	if( resid != NULL ){ delete [] resid; }
	delete [] indices;
	delete [] spectrum;
	return(corrscores);
}
