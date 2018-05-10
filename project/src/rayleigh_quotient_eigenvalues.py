'''
Power Iteration method of finding dominant eigenvector.
	source : https://www.wikiwand.com/en/Power_iteration
'''
import numpy as np
from numpy.linalg import norm


def rayleigh_quotient( A, b ):
	bt = conjugate_transpose( b )
	return ( np.dot( np.dot( bt, A ), b ) )/( np.dot( bt, b ) )

def conjugate_transpose( b ):
	return np.conj( b.T )

def __test():
	# same A as in power_iteration_eigenvector.py
	A = np.array( [ [ 0.5, 0.5 ], [ 0.2, 0.8 ] ] )

	# value estimated by test case in power_iteration_eigenvector.py
	b = np.array([0.70710679, 0.70710677])
	return rayleigh_quotient( A, b )#, np.linalg.eig( A )

if __name__ == '__main__':
	print __test()
	