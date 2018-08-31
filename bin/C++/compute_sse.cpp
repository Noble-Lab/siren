using namespace std;

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "alternateBW.h"
#include "deconvolution.h"
#include "bitonicSGD.h"
#include "lassoBGD.h" 
//#include "lassoADAM.h" 
#include "lassoSGD.h" 
#include "io.h" 

extern int opterr;
