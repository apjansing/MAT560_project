'''
Reads an xml representation of 
code inspired from https://fenicsproject.org/qa/9414/fenics-mesh-generation-mark-inner-region/ 

Authors: Alexander Jansing, Aaron Gregory
Course: MAT 560 - Numerical Differential Equations
Date: April 29th, 2018
Project: Finding Eigenvalues and Eigenvectors of Plates
'''
from __future__ import print_function
from dolfin import *
from dolfin.cpp.mesh import *
from mshr import *
import pylab
from os import mkdir
from re import sub
from sys import argv
import numpy as np

def findEigenValuesVectorsAndPlot(filename):
    m = 3
    mesh = Mesh(filename)
    U = FiniteElement('CG', triangle, 1)
    V = FiniteElement('CG', triangle, 1)
    W = FunctionSpace(mesh, V * U)

    class DirichletBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary

    # Define boundary condition
    u0 = Constant(0.0)
    boundaryCondition1 = DirichletBC(W.sub(0), u0, DirichletBoundary())
    boundaryCondition2 = DirichletBC(W.sub(1), u0, DirichletBoundary())

    # Define the bilinear form
    (u, v) = TrialFunction(W)
    (f1, f2) = TestFunction(W)
    a = -(dot(grad(u), grad(f2)) + dot(grad(v), grad(f1)))*dx
    L = (u*f1 + v*f2)*dx

    # Create the matrices
    A = PETScMatrix()
    b = PETScMatrix()
    assemble(a, tensor = A)
    assemble(L, tensor = b)
    boundaryCondition1.apply(A)
    boundaryCondition2.apply(A)

    # Create eigensolver
    eigensolver = SLEPcEigenSolver(A, b)

    # Compute all eigenvalues of A x = \lambda x
    print("Computing eigenvalues. This can take a minute.")
    eigensolver.solve()

    print("found: ", eigensolver.get_number_converged())

    # get style name from filename
    baseName = sub('\.xml', '', filename)
    baseName = baseName[baseName.rfind('/')+1:]

    mkdirs(baseName)

    eigenVectors = []
    eigenValues = []

    u = Function(W)
    for i in reversed(range(eigensolver.get_number_converged())):
        r, c, rx, cx = eigensolver.get_eigenpair(i)
        print("eigenvalue %d = " % (i), r)
        eigenVectors += [ rx.vec().getArray() ]
        eigenValues += [ r ]
    eigenVectors = np.array(eigenVectors)
    eigenValues = np.array(eigenValues)

    print("eigenVectors", eigenVectors)
    print("eigenValues", eigenValues)

    np.save('eigen/%s/%sEigenVectors.npy' % (baseName, baseName), eigenVectors)
    np.save('eigen/%s/%sEigenValues.npy' % (baseName, baseName), eigenValues)

def mkdirs(s):
    try:
        mkdir('eigen')
    except Exception as e:
        pass
    try:
        mkdir('eigen/%s' % (s))
    except Exception as e:
        pass
    try:
        mkdir('eigen/%s/images' % (s))
    except Exception as e:
        pass
    try:
        mkdir('eigen/%s/vectors' % (s))
    except Exception as e:
        pass


if __name__ == "__main__":
    # Test for PETSc and SLEPc
    if not has_linear_algebra_backend("PETSc"):
        print("DOLFIN has not been configured with PETSc. Exiting.")
        exit()

    if not has_slepc():
        print("DOLFIN has not been configured with SLEPc. Exiting.")
        exit()

    findEigenValuesVectorsAndPlot(argv[1])