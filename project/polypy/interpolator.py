from math import cos, pi, ceil, floor
import numpy as np
from sympy import symbols
from sympy.utilities.lambdify import lambdify


class interpolator(object):
    """docstring for interpolator"""

    def __init__(self, x=[0, 1, 2, 3, 4, -1, -2, -3, -4]):
        self.X = x

    '''
    Given the knows points of a function, F,
        Lagrange numerators, N, and demoninators, D,
        returns the sum all F[i] * [ N[i]/D[i] ]
    '''

    def getLagrangeFunction(self, F, N=None, D=None):
        x = symbols('x')
        if D is None:
            D = self.lagrangeDenominator(self.X)
        if N is None:
            N = self.lagrangeNumerator(self.getPolynomialParts(self.X))
        return lambdify(x, np.dot(F(self.X), np.divide(N, D)))

    '''
    Create Chebyshev n nodes from [a, b].
    Use this twice if the nodes cross the y-axis.
    '''

    def getChebyshevPoints(self, a, b, n=10):
        if a == -b:
            self.X = self.getChebyshevPoints(a, 0, int(floor(
                n / 2.))) + self.getChebyshevPoints(0, b, int(ceil(n / 2.)))
            return self.X
        else:
            x = []
            amb = (a - b) / 2.
            apb = (a + b) / 2.

            for i in range(1, n + 1):
                x += [apb + amb * (cos(pi * (2. * i - 1) / (2. * n)))]
            self.X = x
            return self.X

    '''
    Given a list of polynomal parts, returns the numerators for a
        Lagrange polynomial function.
        i.e. given [ x-1, x-2, x-3 ],
            returns [ (x-2)*(x-3), (x-1)*(x-3), (x-1)*(x-2) ]
    '''

    def lagrangeNumerator(self, N):
        numerators = []
        for i in range(len(N)):
            try:
                numerators += [np.prod(N[0:i]) * np.prod(N[i + 1:])]
            except:
                try:
                    numerators += [float(np.prod(N[0:i])) * np.prod(N[i + 1:])]
                except:
                    numerators += [np.prod(N[0:i]) * float(np.prod(N[i + 1:]))]
        return numerators

    '''
    Given a list of known values, returns the denominators for a
        Lagrange polynomial function.
        i.e. given [1, 2, 3], 
            returns [ (1-2)*(1-3) , (2-1)*(2-3), (3-1)*(3-1) ]
    '''

    def lagrangeDenominator(self, X=None):
        if X is None:
            X = self.X
        denominators = []
        for i in range(len(X)):
            deltaX = list(X[0:i]) + list(X[i + 1:])
            if len(X) > 0:
                denominator = []
                for dX in deltaX:
                    denominator += [X[i] - dX]
                denominators += [np.prod(denominator)]
        return denominators

    '''
    Given the numerator and denominator, N and D respectively, 
        returns an array of the [ n_i/d_i ]
    '''

    def lagrandreCoefficients(self, N, D):
        return np.divide(N, D)

    '''
    Using the lambda function p = lambda s, i : s - i,
        where i is a known value, and s is unknown,
        returns a portion of the polynomials needed for interpolation.
        i.e. return ( x - x_i )
    '''

    def getPolynomialParts(self, X=None):
        if X is None:
            X = self.X
        x = symbols('x')
        # polynomial portion
        p = lambda s, i: s - i
        return np.vectorize(p)(x, X)