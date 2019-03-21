import os
import sys
import numpy as np
from scipy import stats
execfile('sg_filter.py')

# Arguments:
# 1. MS2 file in text file format
# 2. Same MS2 spectra but in original .ms2 format

def smooth_array(a,z,curve):
	halfa = int(len(curve)/2)
	halfb = len(curve)-halfa
	
	#np.convolve(a,curve,'same')
	
	for i in range(0,halfa):
		#print i, a[:(i+halfa)],curve[(halfb-i):]
		if a[i] == 0:
			continue
		z[i] = np.dot(a[:(i+halfa)],curve[(halfb-i):])
	for i in range(halfa,a.shape[0]-halfb):
		#print i,a[(i-halfa):(i+halfb)],curve
		if a[i] == 0:
			continue
		z[i] = np.dot(a[(i-halfa):(i+halfb)],curve)
	for i in range(a.shape[0]-halfb,a.shape[0]):
		#print i,a[(i-halfa):],curve[:len(a[(i-halfa):])]
		if a[i] == 0:
			continue
		z[i] = np.dot(a[(i-halfa):],curve[:len(a[(i-halfa):])])

# The dictionary has keys of tuple (discrete m/z, cycle number) and the elements are list of tuples (m/z, intensity)
# It replaces the values in the dictionary with the ones in an array indexed by cycle number that contains intensities
def update_dictionary( dictionary, a, fmz ):
	for i in range(0,a.size):
		if a[i] == 0:
			continue
		t = fmz,i+1
		s = sum(peak[1] for peak in dictionary[t])
		for j in range(0,len(dictionary[t])):
			dictionary[t][j] = dictionary[t][j][0], a[i]*dictionary[t][j][1]/s, dictionary[t][j][2]
			
def get_cycle( rt, cycle_times ):
	i = 0
	for i in range(0,len(cycle_times)):
		if rt < cycle_times[i]:
			break
	return i
			
# Read in the MS2 text file and populate the dictionary
def get_dictionaries_by_precursor(filename):

	spectra = 0
	dictionaries_by_precursor = {}
	fin = open( filename, 'r' )	
	cycle = 0
	cycle_times = [] # The element in cycle_times is the earliest retention time that falls in that index's cycle
	prev_pmz = float("inf")
	rt = 0
	sys.stderr.write('Constructing dictionary...\n')
	for line in fin:
		# Read in the line
		tokens = line.split()
		for i in range(0,len(tokens)):
			tokens[i] = float(tokens[i])
		[rt, pmz, fmz, intensity, charge] = tokens
		
		# Count spectra
		if pmz != prev_pmz:
			spectra = spectra+1
		
		# Update the cycle number
		if pmz < prev_pmz:
			cycle = cycle+1
			cycle_times.append( float(rt) )
			if cycle % 200 == 0:
				sys.stderr.write('Cycle ' + str(cycle)+'\n')
		prev_pmz = pmz	
	
		if not pmz in dictionaries_by_precursor:
			dictionaries_by_precursor[pmz] = {}
				
		# Populate dictionary
		t = int(fmz), cycle
		if t not in dictionaries_by_precursor[pmz]:
			dictionaries_by_precursor[pmz][t] = []
		dictionaries_by_precursor[pmz][t].append( (fmz,intensity,rt) )
		#print pmz,t,dictionaries_by_precursor[pmz][t]
	sys.stderr.write("Read %d spectra.\n" % spectra)
	fin.close()
	return [dictionaries_by_precursor,cycle_times,rt]


# iptimes are the points you want values for
# inputs is a sorted list of tuples (intensity, retention time) sorted by retention time in increasing order
# If a point to be interpolated is more than threshhold_time away from either boundary, then the boundary is set to 0.
def interpolate_points( inputs, iptimes, threshhold_time ): 
	peak_index = 0
	intensities = [ (0,iptime) for iptime in iptimes ]
	
	for i in range(0,len(iptimes)):
		rt = iptimes[i]
		left_rt = inputs[peak_index][1]
		if rt < left_rt:
			continue
		right_rt = inputs[peak_index+1][1]
		while not( rt >= left_rt and rt <= right_rt ): 
			peak_index = peak_index+1
			if peak_index+1 >= len(inputs):
				return( intensities )
			left_rt = inputs[peak_index][1]
			right_rt = inputs[peak_index+1][1]
		interval = right_rt - left_rt
		right_int = inputs[peak_index+1][0]
		left_int = inputs[peak_index][0]
		
		if (right_rt-rt) >= threshhold_time:
			right_int = 0
		if (rt-left_rt) >= threshhold_time:
			left_int = 0
		# Linear interpolation
		intensities[i] = ( ( right_int*(rt-left_rt)+left_int*(right_rt-rt) )/interval, iptimes[i] )

	
	return intensities 


#peak_pieces = [(0.0, 3.0), (1.0, 6.0), (2.0, 9.0), (3.0, 12.0), (4.0, 15.0), (5.0, 18.0), (6.0, 21.0)]
#iptimes =  [ float(x)+0.1 for x in range(2,22) ]
#x = interpolate_points( peak_pieces, iptimes, 3 )

def plot_pieces( pieces, col ):
	plt.plot( [piece[1] for piece in pieces],[piece[0] for piece in pieces], color = col )


# The values in the dictionary are triples consisting of fmz, intensity, and rt
# ep is in seconds
def sg_smooth_dictionary( dictionaries_by_precursor, cycle_times, rt, ep ):
	sys.stderr.write('Smoothing...\n')
	cycle = len(cycle_times)
	cycle_time = rt/cycle*60
	interval = int(ep/cycle_time)+1

	peak_num = 0
	
	N = interval/2
	pol_degree = 4
	curve = calc_coeff(N, pol_degree, 0)
	#print 'Curve', curve
	
	for pmz in dictionaries_by_precursor.keys():
		sys.stderr.write(str(pmz)+'\n')
		dictionary = dictionaries_by_precursor[pmz]

		prev_mz = -1		
		for mz in range(0,2000):
			exists = False
			peak_pieces = []
			for cycle_num in range(1,cycle+1):
				t = mz,cycle_num
	
				if t in dictionary:	
					# Put the peak into the chromatogram
					# Combines peaks in the same binned fmz window
					exists = True
					intensity = 0.0
					for peak_piece in dictionary[t]:
						intensity = intensity + peak_piece[1]
						peak_num = peak_num + 1
					peak_pieces.append( (intensity, dictionary[t][0][2]) )
			
			if not exists:
				continue
			# Prepare info for interpolation
			[start,end] = [peak_pieces[0][1],peak_pieces[len(peak_pieces)-1][1]]
			rt_interval = end-start
			if start == end:
				continue
			times = np.arange( start, end, rt_interval/cycle )
			#plt.clf()
			#plot_pieces( peak_pieces, "r" )
			intensities = interpolate_points(peak_pieces, times, rt_interval/cycle )
			interpolated_array = np.array( [ intensity[0] for intensity in intensities ] )
			#plot_pieces( intensities, "g" )
			#print 'Sum interpolated', sum( abs( interpolated_array ) )
			
			# Smooth the array
			smoothed_array = np.zeros(( len(interpolated_array) ))
			smooth_array( interpolated_array, smoothed_array, curve )
			
			# Re-interpolate
			## Put smoothed array back into tuple form into intensities
			#print 'Sum smoothed', sum( abs( smoothed_array) )
			intensities = [ (smoothed_array[i],times[i]) for i in range(0,len(intensities)) ] 
			## re-interpolate
			#print 'Sum intensity', sum( [ intensity[0] for intensity in intensities  ] )
			#print 'len(peak_pieces)', len(peak_pieces)
			new_values = interpolate_points( intensities, [ piece[1] for piece in peak_pieces ], rt_interval/cycle )
			#plot_pieces( new_values, "b" )
			#if mz == 1000:
			#	plt.savefig( 'smoothing_visualization.png' )
			#	return
			#print 'Sum', sum( [ new_value[0] for new_value in new_values  ] )
			
			
			# Put the new values into the dictionary
			index = 0
			for cycle_num in range(1,cycle+1):
				t = mz,cycle_num
				
				if t in dictionary:	
					intensity = new_values[index][0]
					#print 'intensity', intensity
					# Take the peak out of the smoothed chromatogram and put it into the dictionary
					# Combines peaks in the same binned fmz window
					s = sum( [ peak_piece[1] for peak_piece in dictionary[t]  ] )
					for peak_index in range(0,len(dictionary[t])):
						#print 'old', dictionary[t][peak_index]
						dictionary[t][peak_index] = ( dictionary[t][peak_index][0], intensity*dictionary[t][peak_index][1]/s, dictionary[t][peak_index][2] )
						#print 'new', dictionary[t][peak_index]
						
					index = index+1	
			
	sys.stderr.write("Smoothed %d peaks.\n" % peak_num)
			

# Smooth dictionary
def smooth_dictionary( dictionaries_by_precursor, cycle_times, rt, std, ep, smoothing_function = 'savitzky-golay' ):
	sys.stderr.write('Smoothing...\n')
	cycle = len(cycle_times)
	cycle_time = rt/cycle*60
	interval = int(ep/cycle_time)+1
	peak_num = 0
	
	if smoothing_function == 'savitzky-golay':
		N = interval/2
		pol_degree = 4
		curve = calc_coeff(N, pol_degree, 0)
	else: 
		curve = stats.norm(0,std)
		span = 3
		curve = curve.pdf(np.linspace(-span,span,interval))
		curve = curve/sum(curve)
	
	for pmz in dictionaries_by_precursor.keys():
		sys.stderr.write(str(pmz)+'\n')
		dictionary = dictionaries_by_precursor[pmz]
		
		prev_mz = -1
		a = np.zeros((cycle))
		z = np.zeros((cycle))
		
		for mz in range(0,2000):
			exists = False
			for cycle_num in range(1,cycle+1):
				t = mz,cycle_num
				if not t in dictionary:
					continue
					
				# Put the peak into the chromatogram
				for peak_piece in dictionary[t]:
					a[cycle_num-1] = a[cycle_num-1]+peak_piece[1]
					peak_num = peak_num+1
					exists = True
			
			if exists:
				# Smooth the array
				z = smooth_array(a,z,curve)
						
				# Put the chromatogram into the dictionary
				update_dictionary( dictionary, z, mz )
					
				# Re-initialize array
				a = np.zeros((cycle))
				z = np.zeros((cycle))
				
	sys.stderr.write("Smoothed %d peaks.\n" % peak_num)
			
# Output dictionary into new MS2 format
def dictionary_to_ms2( ms2filename, dictionaries_by_precursor, outfilename):
	ms2file = open(ms2filename, "r")
	fout = open(outfilename,'w')

	lineNum = 0
	numPeaks = 0
	numSkipped = 0
	rt = 0 #
	prev_precursor = 'boo'
	first_precursor = 'boo' #
	first_scan = 'boo'#
	second_stan = 'boo'#
	header = True
	max_fmz = 2000
	for line in ms2file:
		lineNum += 1
		line = line.rstrip()
	
		words = line.split()
		if ((words[0] == "H") or (words[0] == "Z") or (words[0] == "D")):
			fout.write( line + '\n' )
			header = True      
			continue
	
		elif (words[0] == "S"):
			header = True
			if (len(words) != 4):
				sys.stderr.write("Misformatted line %d.\n%s\n" % (lineNum, line))
				sys.exit(1)
			precursor = words[3]
			dictionary = dictionaries_by_precursor[float(precursor)]

			if first_precursor == 'boo':
				prev_precursor = precursor
				first_precursor = precursor
				first_scan = words[1]
				second_scan = words[2]
			elif float(precursor) < float(prev_precursor):
				rt = rt+1
				if rt % 100 == 0:	
					sys.stderr.write(str(rt)+'\n')
			prev_precursor = precursor		
			fout.write( line + '\n' )
		
		elif words[0] == "I":
			header = True
			fout.write( line + '\n' )
		
		elif len(words) == 2:
			tokens = line.split()
			if header == True:
				header = False
				pmz = float(precursor)
			
				for fmz in range(0,max_fmz):
					t = fmz,(rt+1)
					#print t
					if t in dictionary:
						for peak in dictionary[t]:
							fout.write( "%s\t%s\n" % ( str(peak[0]).rstrip('0'), str(peak[1]).rstrip('0') ) )
									
		else:
			sys.stderr.write("len(words) == %d\n" % len(words))
			sys.stderr.write("Misformatted line %d.\n%s\n" % (lineNum, line))
			sys.exit(1)
	fout.close()
	sys.stderr.write("Read %d lines.\n" % lineNum)

if False:
	for ms2txtfilename in sys.argv[1:]:
		sys.stderr.write(ms2txtfilename+'\n')
		[dictionaries_by_precursor,cycle_times,rt] = get_dictionaries_by_precursor(ms2txtfilename)
		smooth_dictionary( dictionaries_by_precursor, cycle_times, rt,0.5, 20 )
		ms2filename = ms2txtfilename[0:-4]+'.ms2'
		outfilename = 'ep20_'+ms2filename
		sys.stderr.write( 'Converting...\n' )
		dictionary_to_ms2( ms2filename, dictionaries_by_precursor, '/usr/tmp/'+outfilename)

if False:
	print 'Reading...'
	[dictionaries_by_precursor,cycle_times,rt] = get_dictionaries_by_precursor('/net/gs/vol1/home/alexhu/proj/crux-projects/2013MouseHeartProteome/data/charged/centroided/excerpt_525to650.txt')
	print 'Smoothing..'
	sg_smooth_dictionary( dictionaries_by_precursor, cycle_times, rt, 20 )
	print 'Converting'
	ms2filename = '/net/gs/vol1/home/alexhu/proj/crux-projects/2013MouseHeartProteome/data/charged/centroided/real_excerpt_525to650.ms2'
	outfilename = '/net/gs/vol1/home/alexhu/proj/crux-projects/2013MouseHeartProteome/data/charged/centroided/excerpt_525to650.ms2'
	sys.stderr.write( 'Converting...\n' )
	dictionary_to_ms2( ms2filename, dictionaries_by_precursor, outfilename)

#[dictionaries_by_precursor,cycle_times,rt] = get_dictionaries_by_precursor('/net/gs/vol1/home/alexhu/proj/crux-projects/2013MouseHeartProteome/data/charged/centroided/82593_lv_mcx_DIA_5mz_525to650.txt')

if False:
	for ms2txtfilename in sys.argv[1:]:
		sys.stderr.write(ms2txtfilename+'\n')
		[dictionaries_by_precursor,cycle_times,rt] = get_dictionaries_by_precursor(ms2txtfilename)
		sg_smooth_dictionary( dictionaries_by_precursor, cycle_times, rt, 20 )
		ms2filename = ms2txtfilename[0:-4]+'.ms2'
		outfilename = 'sg20_'+ms2filename
		sys.stderr.write( 'Converting...\n' )
		dictionary_to_ms2( ms2filename, dictionaries_by_precursor, outfilename)


