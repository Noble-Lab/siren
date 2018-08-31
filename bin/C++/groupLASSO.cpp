#include <stdio.h>
#include "matrix.h"
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

double maxmass = 2400;
double res = 0.50025395;
double offset = 0.4;


int discretize(double mz,double bw,double offset){
        return( (int)((mz/bw) + 1.0 - offset) );
}

int M = discretize( maxmass, res, offset );
int ix( int row, int column ){
	return( row*M+column );
}


int main(){
	
	// Read in the amino acid masses
	std::vector<double> masses;
	std::vector<int> imasses;
	std::vector<char> alphabet;

	ifstream file( "/net/noble/vol3/user/alexhu/proj/OpenSWATH/bin/2016-08-03/aamasses.txt" );
	string line;
	getline( file, line );
	
	char letter[50];
	while(file.good()){

		getline( file, line);
		istringstream s(line);

		s >> letter;
		alphabet.push_back(letter[0]);	
		s >> letter; s >> letter; s >> letter;
		masses.push_back(atof(letter));
		imasses.push_back( discretize(atof(letter),res,offset) );
		//printf("%c\t%f\n", alphabet.back(), masses.back() );
	}

	std::sort(masses.begin(), masses.end());
	std::sort(imasses.begin(), imasses.end());

	//printf("%i %i %i\n", imasses[0], imasses[1], imasses[2] );
	//return(0);
	// Build the matrix of fragment edges
	bool * edges = new bool[M*M];
	int i = discretize( 0.0, res, offset );
	edges[ix(i,i)] = true;
	int j;
	int k;
	int m;
	while( i < M ){
		for( j=0; j<=i; j++ ){
			if( edges[ix(i,j)] ){
				for( k=0; k<imasses.size(); k++ ){
					m = imasses[k];
					if( (m+i) >= M ){
						break;
					}
					edges[ ix(i,i+m) ] = true;
					edges[ ix(i+m,i) ] = true;
				}
				break;
			}
		}
		i++;	
	}
	printf("Nonzero in edges: %i\n", nnz( edges, M, M ) );
	
	// Find all b-ions for precursor at charge 1 mass 1000
	m = 2000;
	vector<int> peaks;
	int lm = M+1;
	int * to_visit = new int[lm]; for( i = 0; i<lm; i++ ){ to_visit[i] = -1; } // Initialized with values of -1

	int vi = 0; //This is the index for the next peak to visit
	to_visit[vi] = m;
	bool * added = new bool[lm]; for( i=0; i<lm; i++ ){ added[i] = false; }
	added[m] = true;
	int ei = 1;
	int new_m;
	while( to_visit[vi] != -1 && vi < lm ){
		m = to_visit[vi];
		//printf("%i,",m);
		peaks.push_back(m);		
		for(new_m=0; new_m<m; new_m++ ){
			if( edges[ ix(m,new_m) ] ){
				if( !added[new_m] ){
					to_visit[ei] = new_m;
					added[ new_m ] = true;
					printf("%i %i %i\n", vi, ei, m );
					ei++;
				}
			}
		}
		vi++;
	}
	printf("Numpeaks: %i\n", peaks.size() );
	std::sort(peaks.begin(), peaks.end());

	for( i=0; i<peaks.size(); i++ ){
		printf("%i,",peaks[i] );
	}

	printf("\n");
	return(0);
}








