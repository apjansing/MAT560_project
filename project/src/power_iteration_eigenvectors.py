'''
Power Iteration method of finding dominant eigenvector.
	source : https://www.wikiwand.com/en/Power_iteration
'''
import numpy as np
from numpy.linalg import norm


def power_iteration( A, b0, tol = 0.0000001 ):
	n = len( b0 )
	b1 = getBi( A, b0, n, tol )
	while norm( np.subtract( b1, b0 ), ord = n ) > tol:
	# for i in range(10000):
		b0 = b1
		b1 = getBi( A, b0, n )
	return b1

def getBi( A, b, n, tol = 0.0000001 ):
	b = np.dot( A, b )
	b_norm = norm( b, ord = n )
	if b_norm < tol:
		b /= 0.0000001
	else:
		b /= b_norm
	return b

def __test():
	A = np.array( [ [ 0.5, 0.5 ], [ 0.2, 0.8 ] ] )
	b = np.random.rand( 2 )
	return power_iteration( A, b )#, np.linalg.eig( A )[1]

def __random_test( n ):
	A = np.random.rand( n, n )
	b = np.random.rand( n )
	return power_iteration( A, b )#, np.linalg.eig( A )[1]

if __name__ == '__main__':
	print __test()
	print __test()
	for i in range(10, 151):
		__random_test( i )
	print 'Random square matrices of sizes 10-150 completed.'