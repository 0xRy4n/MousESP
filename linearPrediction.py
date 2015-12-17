from numpy import *
from numpy import zeros as zeroes
from scipy.fftpack import fft, ifft
from scipy.linalg import *
import math


def frange(x, y, change):
	while x < y:
		yield x
		x += change


# data is in form of: [ [dx, dy, t], [dx, dy, t] ... ]
def interpolate(data):
	assert len(x) == len(y)
	xArray = []
	yArray = []
	for entry in data:
		xArray.append(entry[0])
		yArray.append(entry[1])

	interpolated = numpy.interp(frange(1, data[len(data) - 1][2], 0.02), xArray, yArray)

	retVal = [[]]

	for i, val in enumerate(interpolated):
		retVal[i][0] = i
		 
		retVal[i][1] = val
		retVal                                                 

	return(interpolated)

#
def resample(data):
	resampled_x = scipy.signal.resample(x,)


# Returns the function of a polynomial using the given coefficients
def getPolynomial(coeffs):

	# Defines the template of a polynomial
	def polynomial_template(x):
		retVal = 0

		# Algorithmic version of '(a_1 * x^n) + (a_2 * x^n-1) + ... (a_n * x^0)'
		for i in reversed(range(0, len(coeffs) - 1)):
			retVal  += coeffs[i] * (x**i)
		return(retVal)

	return(function_template)


# Gets the coefficients of a polynomial that passes through all points given by arrays x and y
def getCoefficients(x, y):
	assert(len(x) == len(y))

	deg = len(x) - 1
	coeffs = numpy.polyfit(x, y, deg) # Finds coefficients that best describes the polynomnial of the points within the x and y arrays with a degree of deg.

	
	return(coeffs)

# Autocorrelation is defined as R(i) = E[f(n)f(n-i)] where E is the expected value.
def autocorrelate(f, n, i):
	retVal = scipy.stats.rv_continous.expect(f(n)f(n - i)) 
	return(retVal)


def nextpow2(n):
	retVal = 2 * (ceil(log2(n)))
	return(retVal)

def lpc(x, p = None):
	x = array(x)
	m = 0
	n = 0

	try:
		m, n = x.shape
	except:
		m = 1
		n = x.shape[0]
	
	if n > 1 and m == 1:
		x = x.reshape(n, m)
		m, n = x.shape

	if p == None:
		p = m - 1
	elif p < 0:
		print("Order must be positive.")
		raise ValueError
	elif p > m:
		print("Order is too large, must be smaller than length of x.")
		raise ValueError

	X = fft(x, 2 ** nextpow2(2* m - 1))
	R = ifft(power(absolute(X), 2)) / m

	a, b, c = levinson(R, p)

	return(a)


def levinson(r, p = None):
	r_0 = real(r[0])
	row_1 = r[1:]
	
	if p == None:
		p = len(row_1) - 1

	reals = isrealobj(r)
	if reals:
		z = zeroes(p, dtype=float)
		z_orig = z
	else:
		z = zeroes(p, dtype=complex)
		z_orig = zeroes(p, dtype=complex)

	for i in range(0, p):
		base = row_1[i]
		if i == 0:
			branch = -base / r_0
		else:
			for k in range(0, i):
				base += z[k] * row[i - k - 1]
			branch = -base / r_0
		if reals:
			r_0 *= 1.0 - branch ** 2
		else:
			r_0 *= (1.0 - (branch.real ** 2 + branch.imag ** 2))

		z[i] = branch
		z_orig[i] = branch

		i_half = (i+1)/2
		if reals:
			for k in range(0, i_half):
				ik = i - k - 1
				base = z[k]
				z[k] += branch * z[ik]
				if k != ik:
					z[ik] += branch*base
		else:
			for k in range(0, i_half):
				ik = i - k - 1
				base = z[k]
				z[k] += branch * z[ik].conjugate()
				if k != ik:
					z[ik] += branch * base.conjugate()

		return(z, z_0, z_orig)






# Test
data = [1,2,3,4,5,6,7,8,200]
thing = lpc(data, 4)
print(thing)






