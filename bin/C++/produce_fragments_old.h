#include "regression.h" 

using namespace std;

/*
int discretize(double mz,double bw,double offset){
    //    int moo = (int)((mz/bw) + 1.0 - offset);
//	printf("moo: %i\n",moo);
	return( (mz-offset)/bw );
	//return( moo );
}*/

int discretize(double mz,double bw,double offset){
        return( (int)((mz/bw) + 1.0 - offset) );
}

double undiscretize(int mz,double bw,double offset){
	return( (0.5+mz+offset-1.0)*bw );
}


double massNH3Mono = 17.02655;     // mass of NH3 (monoisotopic)
double massCOMono =  27.9949;      // mass of CO (monoisotopic) // lost only by bmasses
double massH2OMono = 18.010564684; // mass of water (monoisotopic)
double massHMono = 1.0078246; // Mass of hydrogen (monoisotopic)

class fragmentgroup{
	public:
	// Data parameters
	int N; // Number of candidate peptides
 	int M1; // # Precursor m/z bins
	int M2; // # Fragment m/z bins
	double maxpmass;
	int maxpmassi;

	// Precursor m/z binning info
	double offset1;
	double binwidth1;
	// Fragment m/z binning minfo
	double offset2;
	double binwidth2;
	// Amino acid masses
	std::vector<double> amasses;
	std::vector<int> aimasses;
	std::vector<char> alphabet;

	// Integer masses of neutral losses
	int nh3;
	int co; 
	int h2o; 
	int h;

	// Input
	dbl_matrix * X1; // Theoretical MS1 spectra, MxN

	// Initialization
	fragmentgroup();
	fragmentgroup(dbl_matrix * x1, regression * r, double start_mz );
	fragmentgroup(dbl_matrix * x1, double m2, double o1, double bw1, double o2, double bw2,  char const * od ); // If you have the theoretical matrix and you want to infer what their m/z's are just from that

	// If instead you already have the m/zs of the precursors in a file called "annfile"
	fragmentgroup(char const * annfile, int N, double m2, double o1, double bw1, double o2, double bw2,  char const * od );


	// Variables to store analysis
	double * masses; // precursor masses
	int * charges; // precursor charges
	bool * edges; // edges
	bool * bmasses; // represent the masses of b-ions of the precursors
	FILE * logfile;
	FILE * groupfile;
	char const * outdir;

	// Functions
	void infer_mass_charge();
	void buildbmasses();
	void addpeaks( int m1, int m2 ); // Adds the vector of b-ions from m1 to the 
	double group_penalty();
	double  * group_penalty_pp();
	//vector of b-ions from m2

	void buildgroupmatrices();
	class groupcsrmatrix {
		public:
		fragmentgroup * fg;
		int n;

		// Data parameters
		double mass; // precursor mass
		int P; // number of groups
		int charge;

		// Data
		int maxp; // Maximum number of nonzero groups per precursor mass
		int * nz;
		int * nzperm; // Number of nonzero precursor groups per m
		double * csr;	
		double * csr_gradients;
		int maxm; // maximum m/z bin with at least one group
		double * penalties;	

		// Output stuff
		void output_latent_group( char const * outfile );

		// Constructor	
		groupcsrmatrix( fragmentgroup* x, int n );
		//groupcsrmatrix( int n, double pm ); // A vector of booleans that says
		//whether a b-ion exists at a particular m/z 
		~groupcsrmatrix();
	};

	vector<fragmentgroup::groupcsrmatrix*> precursorgroups;

	~fragmentgroup();
};

fragmentgroup::~fragmentgroup(){
	delete [] masses;
	delete [] charges;
	delete [] edges;
	delete [] bmasses;
	fclose(logfile);
//	fclose(groupfile);
}

fragmentgroup::fragmentgroup(char const * annfile, int n, double m2, double o1, double bw1, double o2, double bw2,  char const * od ){
	
	printf("X1 not given. N: %i\n", n);

	M2 = m2;
	offset1 = o1;
	binwidth1 = bw1;//0.10005079;
	offset2 = o2;//0.4;
	binwidth2 = bw2;//0.50025395;
	N = n;

	logfile = fopen("fglog.txt","w");
//	groupfile = fopen("groupannotations.txt","w");
//	fprintf(groupfile,"n\tgn\tb-nh3\tb-h2o\tb-co\tb\ty-nh3\ty-h2o\ty\n");
	outdir = od;	

	// Integer masses off neutral losses
	nh3 = (massNH3Mono)/binwidth2;
	co = (massCOMono)/binwidth2;  // lost only only by b ions 
	h2o = (massH2OMono)/binwidth2; 
	h = (massHMono)/binwidth2;
	
	printf("nh3: %i, co: %i, h2o: %i, h: %i \n", nh3, co, h2o, h );

	// Read in the masses and charges of the peptides
	charges = new int[N];
	masses = new double[N];
	char letter[50];

	ifstream afile( annfile ); // To read in the peptide annotations
	string line;
	getline( afile, line );
	maxpmass = 0;
	// These require a header whose values are:
	//  m/z | charge  
	n=0;
	while(afile.good()){
		if( n >= N ){
			break;
		}
		getline( afile, line);
		istringstream s(line);
		s >> letter;
		masses[n] = atof(letter);
		s >> letter; 
		charges[n] = atoi(letter);
		masses[n] = (masses[n]-massHMono)*charges[n];
		if( masses[n] > maxpmass ){
			maxpmass = masses[n];
		}
		n++;
	}
	maxpmassi = discretize( maxpmass, binwidth2, offset2 );
	if( M2 > maxpmassi ){ maxpmassi = M2; }

	// Read in the masses of the amino acids (using whatever alphabet is available
	ifstream file( "/net/noble/vol3/user/alexhu/proj/OpenSWATH/bin/2016-08-03/aamasses.txt" );
	getline( file, line );
	
	while(file.good()){
		getline( file, line);
		if( line.size() == 0 ){
			break;
		}	
		istringstream s(line);

		s >> letter;
		alphabet.push_back(letter[0]);	
		s >> letter; s >> letter; s >> letter;
		amasses.push_back(atof(letter));
		aimasses.push_back( (int)(atof(letter)/binwidth2) );
	}

	std::sort(amasses.begin(), amasses.end());
	std::sort(aimasses.begin(), aimasses.end());
	//exit(0);
	buildbmasses();
	buildgroupmatrices();
}

// x1 is the theoretical MS1 matrix
// r is the SGD object on which deconvolution should be performed
// Here, x1 is M by N.
fragmentgroup::fragmentgroup(dbl_matrix * x1, double m2, double o1, double bw1, double o2, double bw2,  char const * od ){
	
	X1 = x1;
	M1 = X1->nrow;
	M2 = m2;
	//M2 = r->M;
	N = x1->ncol;

	printf("x1 is (%ix%i)\n",M1,N);

	offset1 = o1;
	binwidth1 = bw1;//0.10005079;
	offset2 = o2;//0.4;
	binwidth2 = bw2;//0.50025395;

	logfile = fopen("fglog.txt","w");
//	groupfile = fopen("groupannotations.txt","w");
//	fprintf(groupfile,"n\tgn\tb-nh3\tb-h2o\tb-co\tb\ty-nh3\ty-h2o\ty\n");
	outdir = od;	

	// Integer masses off neutral losses
	nh3 = (massNH3Mono)/binwidth2;
	co = (massCOMono)/binwidth2;  // lost only only by b ions 
	h2o = (massH2OMono)/binwidth2; 
	h = (massHMono)/binwidth2;
	
	charges = new int[N];
	masses = new double[N];

	infer_mass_charge();

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
		amasses.push_back(atof(letter));
		aimasses.push_back( (int)(atof(letter)/binwidth2) );
	}

	std::sort(amasses.begin(), amasses.end());
	std::sort(aimasses.begin(), aimasses.end());
	//exit(0);
	buildbmasses();
	buildgroupmatrices();
}

void fragmentgroup::infer_mass_charge(){
	printf("infer mass charge()\n");
	bool * x1t = transpose_bool( X1->m, X1->nrow, X1->ncol );
	bool * x1trow;
	double diffs;
	int prevm;
	int npeaks;
	int startbin;
	double mz;
	maxpmass = 0;

	FILE * f = fopen( "masscharge.txt", "w" );
	fprintf(f,"mass\tcharge\n");
	for( int n=0; n<N; n++ ){
		x1trow = &(x1t[n*M1]);
		prevm = -1;
		npeaks = 0;
		diffs = 0;
		for( int m=0; m<M1; m++ ){
			if( x1trow[m] ){
				npeaks++;
				if( npeaks > 1 ){
					diffs += m-prevm;
				}
				else{
					startbin = m;
				}
				prevm = m;
			}
		}
		diffs = diffs/(npeaks-1);

		if( npeaks < 2 ){
			charges[n] = 1;
		}
		else{
			charges[n] = round(1.0/diffs/binwidth1);
		}
		mz = offset1 + ( 0.5+startbin )*binwidth1;
		masses[n] = (mz-massHMono)*charges[n];
		if( masses[n] > maxpmass ){
			maxpmass = masses[n];
		}
		fprintf(f,"%f\t%i\n",masses[n],charges[n]);
	}		
	//maxpmassi = (maxpmass-offset2)/binwidth2;
	maxpmassi = discretize( maxpmass, binwidth2, offset2 );
	if( M2 > maxpmassi ){ maxpmassi = M2; }
	delete [] x1t;
	fclose(f);
}

// bmasses represent the masses of the b-ions!
void fragmentgroup::buildbmasses(){
	// Build the matrix of fragment edges
	edges = falses(maxpmassi,M2);
	bool * erow;
	// Build matrix of fragment positions
	bmasses = falses(maxpmassi,M2);
	printf("Buildbmasses M2: %i, maxpmassi: %i, maxpmass: %f \n", M2, maxpmassi, maxpmass );
	printf("binwidth2: %f, offset2: %f\n", binwidth2, offset2);
	int i = discretize( 0.0, binwidth2, offset2 );
	printf(" -The zeroth index:%i\n",i);
	edges[i*M2+i] = true;
	int j;
	int k;
	int m;
	while( i < maxpmassi ){
		erow = &(edges[i*M2]);
		for( j=0; j<=min(i,M2-1); j++ ){
			if( erow[j] ){
				for( k=0; k<aimasses.size(); k++ ){
					m = aimasses[k];
					if( (m+i) < M2 ){
						erow[i+m] = true;
					}
					if( (m+i) < maxpmassi ){
						addpeaks( i, i+m );			
						if( i < M2 ){
							edges[ (i+m)*M2+i ] = true;
						}
					}
				}
			}
		}
		i++;
	}
	write_csr_matrix( bmasses, maxpmassi, M2, "bmasses.txt" );
	printf("Nonzero in edges: %i\n", nnz( edges, M2, M2 ) );
}

void fragmentgroup::buildgroupmatrices(){
	printf("Building group matrices...\n");
	groupcsrmatrix * gm;
	for( int n=0; n<N; n++ ){
		if( n % 250 == 0 ){
			printf("Creating group %i out of %i\n", n, N );
		}	
		gm = new fragmentgroup::groupcsrmatrix( this, n );
		precursorgroups.push_back( gm );
	}
}

// The assumption is that bmasses[m1,:] has already been correctly
// filled. This adds m1 to m2.
void fragmentgroup::addpeaks(int m1, int m2){
	bool * fp1 = &(bmasses[m1*M2]);
	bool * fp2 = &(bmasses[m2*M2]);

	if( m1 < M2 ){
		fp2[m1] = true;
	}
	for( int m=0; m<m1; m++ ){
		if( fp1[m] ){
			fp2[m] = true;
		}
	}
}

fragmentgroup::groupcsrmatrix::groupcsrmatrix( fragmentgroup * FG, int n_index  ){
	fg = FG;
	n = n_index; // index of precursor
	mass = fg->masses[n]; // precursor mass
	charge = fg->charges[n];

	double offset2 = fg->offset2;
	double binwidth2 = fg->binwidth2;
	bool * bmasses = fg->bmasses;
	int M2 = fg->M2;
	int massix = discretize( mass, binwidth2, offset2 );//(mass-offset2)/binwidth2;
	double pmz1 = mass + massHMono; // charge 1 precursor m/z
	double pmz = mass/(charge)+massHMono;

	//printf("mass: %f, charge: %i offset2: %f binwidth2: %f M2: %i\n", mass, fg->charges[n], offset2, binwidth2, M2 );

	// Count how many groups there are
	P=1; // 1 for the precursor!
	bmasses = &(bmasses[massix*M2]);
	for( int m=1; m<M2; m++ ){
		if( bmasses[m] ){
			P++;
		}	
	}
	//fprintf(fg->logfile,"%i ngroups: %i mass: %f charge %i pmz %f \n", n, P, mass, fg->charges[n], pmz );
	//printf("%i ngroups: %i mass: %f charge %i pmz %f \n", n, P, mass, fg->charges[n], pmz );
	

	//Output the annotations
	char outfile[200];
	sprintf(outfile,"fragment_annotations.txt");
	FILE * file;
	file = fopen(outfile,"w");
	fprintf(file,"precursor mass: %f, binwidth2: %f, offset2: %f\n", mass, binwidth2, offset2 );
	// Fill in a boolean matrix of bins vs groups 
	bool * groups = falses(P,M2); // It begins as P by M and then it will be transposed
	bool * row;
	int y1ix; // y-ion bin
	int b1ix;
	// Assume that bix=0 corresponds to a m/z of 0, and we should ignore it.
	// So we start at bix=1
	int pix = 0; // Index of the group 
	double fmz;
	bool print = true;
	for( int bix= discretize( 0.0, binwidth2, offset2 ) + 1; bix<(M2-fg->h); bix++ ){
	//for( int bix=1; bix<(M2-fg->h); bix++ ){
		row = &(groups[pix*M2]);
		if( bmasses[bix] ){
			fprintf(file,"Group %i\n",pix);
			//printf("bix: %i, pix: %i \n",bix, pix);
			b1ix = bix+fg->h; // charge 1 b ion (b1 ion)
			if( b1ix-fg->nh3 < M2 ){
				row[b1ix-fg->nh3] = true; // b1 ion ammonia loss
				fmz = undiscretize(b1ix-fg->nh3,binwidth2,offset2);	
				fprintf(file,"b ammonia loss %f %i\n", fmz, b1ix-fg->nh3);
			}
			if( b1ix-fg->h2o < M2 ){
				row[b1ix-fg->h2o] = true; // b1 ion water loss	
				fmz = undiscretize(b1ix-fg->h2o,binwidth2,offset2);	
				fprintf(file,"b water loss %f %i\n", fmz, b1ix-fg->h2o);
			}
			if( b1ix-fg->co < M2 ){
				row[b1ix-fg->co] = true; // b1 ion co loss	
				fmz = undiscretize(b1ix-fg->co,binwidth2,offset2);	
				fprintf(file,"b co loss %f %i\n", fmz, b1ix-fg->co);
			}
			if( b1ix < M2 ){
				row[b1ix] = true; // b1 ion	
				fmz = undiscretize(b1ix,binwidth2,offset2);	
				fprintf(file,"b %f %i\n", fmz, b1ix);
			}
			//y1ix = (pmz1-offset2)/binwidth2-bix;

			y1ix = discretize( pmz1-undiscretize(bix,binwidth2,offset2), binwidth2, offset2 );

			if( y1ix-fg->nh3 < M2 ){
				row[y1ix-fg->nh3] = true; // y1 ion ammonia loss	
				fmz = undiscretize(y1ix-fg->nh3,binwidth2,offset2);	
				fprintf(file,"y ammonia loss %f %i\n", fmz, y1ix-fg->nh3);
			}
			if( y1ix-fg->h2o < M2 ){
				row[y1ix-fg->h2o] = true; // y1 ion water loss	
				fmz = undiscretize(y1ix-fg->h2o,binwidth2,offset2);	
				fprintf(file,"y water loss %f %i\n", fmz, y1ix-fg->h2o);
			}
			if( y1ix < M2 ){
				row[y1ix] = true; // y1 ion 
				fmz = undiscretize(y1ix,binwidth2,offset2);	
				fprintf(file,"y %f %i\n", fmz, y1ix);
			}
//			fprintf(fg->groupfile,"%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n",N,pix,b1ix-fg->nh3,b1ix-fg->h2o,b1ix-fg->co,b1ix,y1ix-fg->nh3,y1ix-fg->h2o,y1ix);
			pix++;
		}
	}	
	fclose(file);
	// Add a group for the precursor and its losses
	//printf("Adding precursor to groups. %ith group\n", pix );
	row = &(groups[pix*M2]);
/* Old way of discretization
	row[(int)((pmz-offset2)/binwidth2)] = true;	
	row[(int)((pmz-offset2-massNH3Mono/charge)/binwidth2)] = true;	
	row[(int)((pmz-offset2-massH2OMono/charge)/binwidth2)] = true;	
	printf("Adding precursor to groups. %ith group, ix: %i\n", pix,(int)((pmz-offset2)/binwidth2) );
	fprintf(fg->groupfile,"%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n",N,pix,(int)((pmz-offset2-massNH3Mono/charge)/binwidth2),(int)((pmz-offset2-massH2OMono/charge)/binwidth2),-1,(int)((pmz-offset2)/binwidth2),-1,-1,-1);
*/		

	row[discretize(pmz,binwidth2,offset2)] = true;	
	row[discretize(pmz-(massNH3Mono/charge),binwidth2,offset2)] = true;
	row[discretize(pmz-(massH2OMono/charge),binwidth2,offset2)] = true;
	//printf("Adding precursor to groups. %ith group, ix: %i\n", pix,discretize(pmz,binwidth2,offset2) );

//	fprintf(fg->groupfile,"%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n",N,pix,discretize(pmz-(massNH3Mono/charge),binwidth2,offset2),discretize(pmz-(massH2OMono/charge),binwidth2,offset2),-1,discretize(pmz,binwidth2,offset2),-1,-1,-1);
	
	// Transpose it! make groups MxP
	bool * temp = groups;
	groups = transpose( temp, P, M2 ); 
	delete [] temp;	
	// Count ihe maximum number of groups that each mass bin is in
	maxp = 0;
	int mp;
	maxm = 0;
	for( int m=0; m<M2; m++ ){
		row = &(groups[m*P]);
		mp = 0;
		for( int p=0; p<P; p++ ){
			if( row[p] ){
				mp++;
			}
		}
		if( mp > 0 ){
			maxm = m;
		}
		if(mp > maxp){
			maxp = mp;
		}
	}
	maxm++;
	//printf("Maxp: %i, Maxm: %i\n",maxp);
	// Now make sparse versions
	csr = new double[maxm*maxp]; zero(csr,maxm,maxp,-1);
	csr_gradients = new double[maxm*maxp]; zero(csr,maxm,maxp,-1);
	nz = new int[maxm*maxp]; zero(nz,maxm,maxp,-1);
	nzperm = zeros_int(maxm,1);
	int nzix;
	for( int m=0; m<maxm; m++ ){
		nzix = m*maxp;
		row = &(groups[m*P]);
		for( int p=0; p<P; p++ ){
			if( row[p] ){
				nz[nzix] = p;
				csr[nzix] = 0.0;	
				csr_gradients[nzix] = 0.0;	
				nzperm[m] += 1;
				nzix++;
			}
		}
	}
	penalties = zeros(1,P);	
	//sprintf( filename, "groups_bool_%i.txt", fg->outdir, n );
	//sprintf( filename, "%s/groups_bool_%i.txt", fg->outdir, n );
	//write_csr_matrix( groups, M2, P, filename );	
	write_csr_matrix( groups, M2, P,"fragmentgroup.txt" );
	delete [] groups;
};

double fragmentgroup::group_penalty(){
	double penalty = 0;
	groupcsrmatrix * gm;
	for( int n=0; n<N; n++ ){
		gm = precursorgroups[n];
		for( int p=0; p<gm->P; p++ ){
			penalty += gm->penalties[p];
		} 
	}
	printf("-group penalty %f ", penalty);
	return( penalty );
}

double * fragmentgroup::group_penalty_pp(){
	double * penalty = zeros(1,N);
	groupcsrmatrix * gm;
	for( int n=0; n<N; n++ ){
		gm = precursorgroups[n];
		for( int p=0; p<gm->P; p++ ){
			penalty[n] += gm->penalties[p];
		} 
	}
	return( penalty );
}

void fragmentgroup::groupcsrmatrix::output_latent_group( char const * outfile ){
	//write_csr_matrix( nz, csr, maxm, P, maxp, outfile );
	write_csr_ixes( nz, maxm, P, maxp, outfile );
}

fragmentgroup::groupcsrmatrix::~groupcsrmatrix(){
	delete [] penalties;
	delete [] nz;
	delete [] csr;
	delete [] nzperm;
}


