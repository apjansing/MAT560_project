from dolfin import *
from dolfin.cpp.mesh import *
from mshr import *
from math import sin, cos
import numpy as np


class Simulator():
        """docstring for Simulator"""
        def __init__(self, mesh, U, eigenVectors, eigenValues,
         pExpr = Expression('sin(x[0])', degree=1), 
         vExpr = Expression('cos(x[0])', degree=1), C = 1):
                # Test for PETSc and SLEPc
                if not has_linear_algebra_backend("PETSc"):
                        print("DOLFIN has not been configured with PETSc. Exiting.")
                        exit()

                if not has_slepc():
                        print("DOLFIN has not been configured with SLEPc. Exiting.")
                        exit()
                
                parameters['reorder_dofs_serial'] = False
                self.timeScale = 1.0 / C
                self.fSpace = FunctionSpace(mesh, U)
                self.eVecs = eigenVectors
                self.eVals = eigenValues
                self.pVec = fem.interpolation.interpolate(pExpr, self.fSpace).vector()
                self.vVec = fem.interpolation.interpolate(vExpr, self.fSpace).vector()
                self.pEig = self.inEigenBasis(self.pVec)
                self.vEig = self.inEigenBasis(self.vVec)
                for i in range(len(self.eVals)):
                        self.vEig[i] /= self.eVals[i]

        def inEigenBasis(self, vec):
                return np.array([np.dot(vec, e) / np.dot(e, e) for e in self.eVecs])
        
        def vecAtTime(self, t):
                S = np.zeros(self.eVecs[0].shape)
                for i in range(len(self.eVecs)):
                        # position component
                        S += cos(self.eVals[i]*t)*self.pEig[i]*self.eVecs[i]
                        # velocity component
                        S += sin(self.eVals[i]*t)*self.vEig[i]*self.eVecs[i]
                return S

        def upperBound(self):
                m = np.zeros(self.eVecs[0].shape)
                for i in range(len(self.eVecs)):
                        m += np.abs(self.pEig[i]*self.eVecs[i])
                        m += np.abs(self.vEig[i]*self.eVecs[i])
                return np.max(m)
        
        def evaluate(self, t, x = None):
                f = Function(self.fSpace)
                f.vector().set_local(self.vecAtTime(t * self.timeScale))
                f.update()
                if x is None:
                        return f
                else:
                        return f(Point(x[0], x[1]))

