# Siren pipeline

Siren can perform three functions:
- Identify elution peak scans, masses, and charges for peptide ions in MS1 data.
- Extract elution profiles for these peptide ions (abundance across multiple scans)
- Provide input for a modified version of DIA-Umpire (also included) to create pseudospectra
   for peptide ions from MS2 data accompanying the MS1 data.  


The inputs are:
- A .ms1 file (required), converted using msconvert.
- An mzXML (optional) associated with the .ms1 file. This can be used by modified DIA-Umpire in conjunction with Siren's extracted elution profile.


The outputs are:
- Precursor ion annotations: a tab-delimited text file containing the mass, charge, and peak scan for peptide ions annotated by Siren
- Elution profiles in a row-major matrix: 
	- a text file containing a row-major sparse matrix containing elution profiles for each peptide ion contained in the preceding text file.
- a text file containing the same elution profiles formatted for use by the modified form of DIA-Umpire
- Pseudospectra


# REQUIREMENTS

## System:

- Siren was designed and tested on Linux systems. It should probably work on a mac.



The Siren pipeline is constructed using python and C++. The modified DIA-Umpire software is in a java .jar file.

To get outputs 1 and 2, run the following command:

To get output 3,



