import os
import sys 

ms1file = sys.argv[1]
ms2file = ''
if len(sys.argv) > 2:
	ms2file = sys.argv[2]

sys.stdout.write('Generating features from ' + ms1file )
os.system('python generate_theoretical_features_from_ms1_combine_individual.py ' + ms1file )
ms1dir = ms1file.replace('.ms1','_ms1_sparsebinned_singlescans_tannotated/')

sys.stdout.write('Generating decoy features (reversing sequences of top isotope intensities)' + ms1file )
os.system('python make_siren_decoys.py ' + ms1dir )

sys.stdout.write('Linking features over elution time')
os.system('python align_ms1_annotations_individual.py ' + ms1dir )
 
bdir = ms1dir.replace('_ms1_sparsebinned_singlescans_tannotated/','_peptide_abundances/')
os.system('python regression_justms1_sparse_job_sgd ' + ms1dir + ' ' + bdir )

sys.stdout.write('Generating un-regularized solutions')
os.system('python compute_paths_job.py '  + ms1dir + ' ' + bdir )

sys.stdout.write('Combining solutions into elution profiles')
alignmentfile = ms1dir + 'joint_annfile_minlen4.txt'
os.system('python combine_ms1_matrices_from_alignment_individual.py' + ms1dir + ' ' + alignmentfile + ' ' + bdir)

sys.stdout.write('Finding elution peaks')
bfile = bdir+'combined_B_sp0.txt'
os.system('python smooth_extract_withintensities_ms1matrices.py ' + bfile + ' 4 ' + alignmentfile ) 

sys.stdout.write('Generating regularization paths')
os.system('python compute_paths_job.py ' + ms1dir + ' ' + bdir ) 

sys.stdout.write('Estimating FDRs for elution peaks')
os.system('python compute_fdrs_from_paths.py ' + ms1dir + ' ' + bdir + ' ' + alignmentfile) 









