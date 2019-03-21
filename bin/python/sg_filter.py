# -*- coding: utf-8 -*-

"""
A module that implements the useful Savitzky–Golay filter.
See http://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_smoothing_filter .
"""

from math import *
from numpy import *

def calc_coeff(num_points, pol_degree, diff_order=0):
    """ Calculates filter coefficients for symmetric Savitzky–Golay filter.
    see: http://www.nrbook.com/a/bookcpdf/c14-8.pdf

    :param num_points: means that 2*num_points+1 values contribute to the smoother.
    :type num_points: int.
    :param pol_degree: is degree of fitting polynomial
    :type pol_degree: int.
    :param diff_order: is degree of implicit differentiation.
      0 means that filter results in smoothing of function.
      1 means that filter results in smoothing the first derivative of function.
    :type diff_order: int.
    :returns:  numpy.array -- The coefficients for the Savitzky–Golay filter.
    """

    # Setup interpolation matrix. You might use other interpolation
    # points and maybe other functions than monomials.

    x = arange(-num_points, num_points+1, dtype=int)
    monom = lambda x, deg : pow(x, deg)

    A = zeros((2*num_points+1, pol_degree+1), float)
    for i in range(2*num_points+1):
        for j in range(pol_degree+1):
            A[i,j] = monom(x[i], j)

    # calculate diff_order-th row of inv(A^T A)
    ATA = dot(A.transpose(), A)
    rhs = zeros((pol_degree+1,), float)
    rhs[diff_order] = (-1)**diff_order
    wvec = linalg.solve(ATA, rhs)

    # calculate filter-coefficients
    coeff = dot(A, wvec)

    return coeff

def smooth(signal, coeff):
    """ Applies coefficients calculated by :py:func:`calc_coeff` to signal. """
    N = (size(coeff)-1)/2
    res = convolve(signal, coeff)
    return res[N:-N]

def smooth_array(a,z,curve):
	halfa = int(len(curve)/2)
	halfb = len(curve)-halfa
	
	for n in xrange(0,a.shape[0]):
		#np.convolve(a,curve,'same')	
		for i in range(0,halfa):
			#print i, a[:(i+halfa)],curve[(halfb-i):]
			if a[n,i] == 0:
				continue
			z[i] = np.dot(a[n,:(i+halfa)],curve[(halfb-i):])
		for i in range(halfa,a.shape[0]-halfb):
			#print i,a[n,(i-halfa):(i+halfb)],curve
			if a[n,i] == 0:
				continue
			z[i] = np.dot(a[n,(i-halfa):(i+halfb)],curve)
		for i in range(a.shape[0]-halfb,a.shape[0]):
			#print i,a[n,(i-halfa):],curve[:len(a[n,(i-halfa):])]
			if a[n,i] == 0:
				continue
			z[i] = np.dot(a[n,(i-halfa):],curve[:len(a[n,(i-halfa):])])


def sparse_convolve( m, curve ):
	print '\nHello!\n'
	print m.shape
	nz = m.nonzero()[0]
	(minnz,maxnz) = ( nz.min(), nz.max() )
	return np.convolve( m.todense()[ minnz:(maxnz+1) ] )


def convolve_nonzero_columns( a, curve ):
	n = curve.shape[0]
	#return np.apply_along_axis(lambda m: np.convolve(m, curve, mode='same'), axis=0, arr=a)	 
	return np.apply_along_axis(lambda m: sparse_convolve(m, curve, mode='same'), axis=0, arr=a)	 

def zero_short_segments( a, seglen ):
	# For each row in a
	for r in xrange(0,a.shape[0]):	
		# Find the zero indices
		z = np.where( a[r,:] == 0 )[0]
		# For each zero index
		for i in xrange(0,z.shape[0]-1):
			# If the intervening sequence of zeros
			# is less than seglen long, zero it
			if (z[i+1]-z[i]-1) < seglen: 
				a[r,z[i]:z[i+1]] = 0

def smooth_array(a,z,curve):
	a_zeroed = np.copy(a)
	zero_short_segments( a_zeroed, 5 )
	for n in xrange(0,z.shape[0]):
		z[n,:] = np.convolve(a_zeroed[n,:],curve,'same')	
	z[ z <= 0] = 0
	z[ a == 0] = 0

#if __name__ == "__main__":

#    print "Sample run to demonstrate the SG Filter"
#    x=arange(0, 10., .1)
#    y=sin(x)
#    coeff = calc_coeff(5, 2, diff_order=1)
#    yPrime=smooth(y, coeff)

    # calibrate due to unknown correction factor:
#    N = len(x)
#    dy =cos(x)  # exakt
#    corrfac = dy[N/2]/yPrime[N/2]

#    import pylab
#    pylab.plot(x,y)
#    pylab.plot(x,corrfac*yPrime)
#    pylab.show()
