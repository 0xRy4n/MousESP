from numpy import *
from numpy import zeros as zeroes
from scipy.fftpack import fft, ifft
from scipy.linalg import *
import math


# For thing in range x to y with a jump of change
def frange(x, y, change):
        while x < y:
                yield x
                x += change


""" Polynomial Things """

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



""" LPC Things """

# Autocorrelation is defined as R(i) = E[f(n)f(n-i)] where E is the expected value.

# Gets next power of 2 for n
def nextpow2(n):
        retVal = 2 * (ceil(log2(n)))
        return(retVal)

# Gets the coefficients from the toeplitz matrix 
# x is 1 dimensional array, p is the order
def lpc(x, p, axis = -1):
        x = array(x)
        m = 0
        n = 0

        try:
                m, n = x.shape
        except:
                m = 1
                n = x.shape[0]
        
        if n > 1 and m == 1:
                x = x.reshape(n, m) # Transponds the x array
                m, n = x.shape

        if p == None:
                p = m - 1
        elif p < 0:
                print("Order must be positive.")
                raise ValueError
        elif p > m:
                print("Order is too large, must be smaller than length of x.")
                raise ValueError


        X = fft(x, 2 ** nextpow2(2* m - 1)) # Gets x's fourier transform, with a bin-width of 2^nextpow2(2 * m - 1)
        r = ifft(power(absolute(X), 2)) / m
        r = r[:,0] # reduces the r matrix to it's first column

        # Swaps axes around
        if axis != -1:
                r = swapaxes(r, axis, -1)

        # Performs levinson recursion r with an order of p
        a, e, k = levinson_1d(r, p)

        # Swaps axes around
        if axis != -1:
                a = swapaxes(a, axis, -1)
                e = swapaxes(e, axis, -1)
                k = swapaxes(k, axis, -1)
        
        return(a, e, k)


def levinson_1d(r, p):
        r = atleast_1d(r) # Makes sure r is only one dimension
        if r.ndim < 1:
                raise ValueError("Thing has to little dimensions.")

        n = r.size
        if n < 1:
                raise ValueError("The thing is empty.")
        elif p > n - 1:
                raise ValueError("P needs to be less than (r.size - 1)")

        if not isreal(all(r[0])):
                raise ValueError("First thing in r must be real.")
        elif not isfinite(1/all(r[0])):
                raise ValueError("First thing in r can't be 0.")

        a = empty(p + 1, r.dtype)
        t = empty(p + 1, r.dtype)
        k = empty(p, r.dtype)

        a[0] = 1 # Sets first thing in a in 
        e = r[0]

        for i in xrange(1, p+1):
                acc = r[i]
                for j in range(1, i):
                    acc += a[j] * r[i-j]
                k[i-1] = -acc / e
                a[i] = k[i-1]

                for j in range(p):
                    t[j] = a[j]

                for j in range(1, i):
                        a[j] += k[i-1] * conj(t[i-j])

                e *= 1 - k[i-1] * conj(k[i-1])

        return(a, e, k)



# Test
data = [1, 20, 2, 30, 3, 40, 4, 50, 5, 60, 7, 70, 8, 600]
thing, thing2, thing3 = lpc(data, 4)
print(thing)
print(thing2)
print(thing3)







