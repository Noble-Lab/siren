# Siren pipeline

## Overview:

Siren can perform three functions:
- Identify elution peak scans, masses, and charges for peptide ions in MS1 data.
- Extract elution profiles for these peptide ions (abundance across multiple scans)
- Provide input for a modified version of DIA-Umpire (also included) to create pseudospectra
   for peptide ions from MS2 data accompanying the MS1 data.  

The inputs are:
- A .ms1 file (required), converted using msconvert.
- An mzXML (optional) associated with the .ms1 file. This can be used by modified DIA-Umpire in conjunction with Siren's extracted elution profile.

The outputs will appear in the same directory as the .ms1 file (or subdirectories created by Siren) and are:
1. Precursor ion annotations: tab-delimited text files containing the mass, charge, and peak scan for peptide ions annotated by Siren. The steps of the Siren pipeline will create multiple versions of this file with different pieces of information. 
2. Elution profiles in a row-major matrix: 
	- a text file containing a row-major sparse matrix containing elution profiles for each peptide ion contained in the preceding text files.
3. A text file containing the same elution profiles formatted for use by the modified form of DIA-Umpire, with a coefficient assigned to each elution peak that denotes its position in the L1 regularization path used to control FDR. 
4. Pseudospectra


# REQUIREMENTS

## System:

- Siren was designed and tested on Linux systems. It should probably work on a mac. The Siren pipeline is uses python and C++. The modified DIA-Umpire software is in a java .jar file.

## Dependencies:

Python 2.7.3 or later
Numpy 1.13.3 or later
Scipy 1.0.0rc1 or later
Scikit-learn 0.19.0 or later

## Installation

Siren relies on C++ programs for the computation of theoretical isotope distributions and 
a subset of the linear modeling. Please compile these programs by doing the following:

1. In bin/C++/Mercury8, type "make"
2. In bin/C++, type "g++ regressionApplication.cpp -o ratest_batch"

## Directions

To get outputs 1 and 2, run the following command within the directory ./bin/python/:
python siren.py FILENAME.ms1

To get outputs 1, 2, 3, and 4, run the following command within the directory ./bin/python/:
python siren.py MS1FILENAME.ms1 MZXMLFILENAME.mzXML

## Detailed description of pipeline

The script siren.py takes the input .ms1 file and runs a sequence of scripts that perform:

1. Noise reduction and data-driven generation of hypothetical ion species per spectrum
2. Generating decoy features for subsequent FDR estimation
3. Linking hypothetical ion species and decoys across multiple scans
4. Inference of abundances of hypothetical ion species via ordinary least squares regression
5. Extraction of hypothetical elution profiles from inferred abundances
6. Smoothing profiles and identifying elution peaks
7. Computing regularization paths to subsequently score elution peaks
8. Score elution peaks
9. Format elution peaks into modified DIA-Umpire pipeline. This file should appear as a tab-delimited "_precursorprofiles.txt" in the same directory as the .ms1 file. Each row is an elution profile for a unique ion species. The elution profile is in the last column, with the format "scan_number:inferred_abundance,scanum_number:inferred_abundance,scanum_number:inferred_abundance" for a contiguous sequence of scans.


